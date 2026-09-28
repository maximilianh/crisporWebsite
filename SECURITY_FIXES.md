# Security fixes in crispor.py - 28/09/26

Two fixes: (1) command injection through the CGI parameters, (2) a limit on the number
of queued jobs per IP address. Neither is committed yet.

## 1. Command injection through the CGI parameters

### Problem

- **Broken input filter.** Commit `2b1b54d` changed the whitelist regex to
  `[^+[]a-zA-Z0-9/:\n\r_. -]`. The unescaped `]` closes the character set, so the regex
  only matched the literal text `a-zA-Z0-9/:...` and let every character through.
  The same commit commented out the list of free-text parameters that skip the filter.
- **Remote command execution.** `getFreeEnergy()` and `showSecondaryStructure()` ran
  `echo <guideSeq> | RNAfold -T <temperature> ...` with `shell=True`. `guideSeq` and
  `temperature` come straight from the URL, so one GET request with `pamId`, `guideSeq`
  and `temperature` could run any command as the web server user.
- **Not only since `2b1b54d`.** The original regex also allowed newlines, and a newline
  starts a new shell command, so `temperature=37%0A<command>` worked before too.
- **Smaller issues.**
  - `extendAndGetSeq()` built a shell command with `chrom` inside quotes.
  - `org` allowed `..` and `/`, and is used in file paths and shell commands.

### Fix

- **Regex (`notOkChars`, `crispor.py:1381`).** The brackets are escaped:
  `[^+\[\]a-zA-Z0-9/:\n\r_. -]`.
- **Free-text parameters (`freeTextParams`, `crispor.py:1386`).** The exemption list
  is restored, plus `multiseq`, `pegPams` and `primers`, which were added since.
  - These parameters need characters the whitelist forbids, e.g. `>` in fasta headers
    or quotes in json.
  - They are no longer skipped entirely: they must not contain `` ` $ ; | & \ < `` or
    a null byte.
- **Line breaks (`multiLineParams`, `crispor.py:1412`).** Rejected in every other
  parameter, except the textareas `geneIds`, `addSeq`, `startSeq`, `endSeq` and
  `replaceInsertSeq`.
- **Per-parameter checks in `cgiGetParams()` (`crispor.py:1480`):**
  - `temperature`: a whole number from 0 to 100
  - `guideSeq`: only `ACGTUN`, upper or lower case
  - `org`: only `[A-Za-z0-9_.-]`, and no `..`
- **No shell.** `getFreeEnergy()`, `showSecondaryStructure()` and `extendAndGetSeq()`
  now start RNAfold, RNAplot and twoBitToFa with a list of arguments. The sequence is
  passed on stdin instead of through `echo`.
- **Not changed.** `libName` was already checked against the list of valid libraries.
  All SQL queries already use placeholders.

### Tests

- **Rejected:** `temperature=37;...`, `temperature=37<newline>...`,
  `guideSeq=ACGT$(...)` and `org=../../etc`. The marker file the test requests tried
  to create was never created.
- **Still working:**
  - the secondary structure page (SVG and free energy)
  - `getFreeEnergy()` values
  - genome sequence retrieval in `extendAndGetSeq()` (ce11)
  - fasta `seq`, json `geneModel` / `pegPams`, `pamId=s12+` and multi-line `geneIds`

## 2. Limit on the number of jobs per IP address

### Problem

Nothing limited the job queue. Every new sequence, genome, PAM or job name creates a
new batchId, so a bot could fill the queue, the disk and the RAM (`bwa bwasw` runs for
every new submission). The only protection was a few hard-coded blocked IPs.

### Fix

- **`MAXJOBSPERIP = 10` (`crispor.py:472`).** The maximum number of jobs, waiting or
  running, that one client can have in the queue.
- **`JobQueue.addJob(..., ip=None)` (`crispor.py:1111`).** When an IP is given,
  `_addJobLimited()` (`crispor.py:1164`) does two things in one `BEGIN IMMEDIATE`
  transaction, so parallel requests cannot slip past the limit:
  - A batchId that is already in the queue is always accepted and never counted.
  - Otherwise the job is refused if the client already has 10 jobs. The user is
    asked to reload the page later, and the job is queued on that reload.
- **Crashed jobs are not counted.** They keep `stepName='crash'` and are never removed
  from the queue.
- **No limit without an IP.** Command-line runs (`noIp`) and calls without `ip` behave
  as before.
- **All three places that add jobs pass the IP:** multiseq/multipam, the off-target
  search and the mutPeg (silent bystander) jobs.
- **IPv6 grouping (`clientKey()`, `crispor.py:1000`, `ipJobCount()`,
  `crispor.py:1152`).** Jobs are counted per client key:
  - IPv4: the address itself
  - IPv6: the /64 network, since one household, lab or phone usually has a whole /64
    and can use any address in it
  - IPv4-mapped IPv6 (`::ffff:1.2.3.4`): counted as the IPv4 address
  - The full address is still written to the queue and to `doneJobs.tsv`.

### Effect on users

- **Reopening a batch.** Opening an existing batchId is never refused, whether the job
  is done, waiting or running, and whoever opens it.
- **Shared networks.**
  - Behind one IPv4 address (NAT), or on one IPv6 /64, everyone shares 10 slots for
    new jobs.
  - This is stricter for sites that give each computer a public IPv4 address but put a
    whole building on one IPv6 /64.
- **IPv4 and IPv6 users.** They have one limit per protocol, so up to 20 jobs.
- **If a real institute hits the limit:** raise `MAXJOBSPERIP` or exempt its network.

### Tests

- **Refused:**
  - the 11th job of an IPv4 address, including in the `::ffff:` form
  - the 11th job from 10 different addresses in the same /64
- **Accepted:**
  - the neighbouring IPv4 address, another /64, and `1.2.3.44` (not counted as
    `1.2.3.4`)
  - a batchId already in the queue
  - a new job after one finished or crashed
  - `noIp`
- **`clientKey()` doctests pass.**

### Before deploying: check the client IP

In this Docker container, `REMOTE_ADDR` is always `172.17.0.1` (the Docker gateway):
every job in `doneJobs.tsv` and every request in the Apache access log has that IP.
As it is, the limit would apply to all users together. On the production server, make
sure `REMOTE_ADDR` is the real client address:

- **Behind a reverse proxy:** enable `mod_remoteip` (installed, not enabled) with
  `RemoteIPHeader X-Forwarded-For` and `RemoteIPInternalProxy <proxy IP>`.
- **Docker port mapping:** set `"userland-proxy": false` in the Docker daemon
  configuration, or use host networking.

## Still to do (not implemented)

- A global cap on the queue length, using `JobQueue.waitCount()`.
- Limit how many `bwa bwasw` processes can run at once in the CGI, or move this step to
  the job queue.
- Rate limiting in Apache (`mod_evasive` / `mod_qos`).
- A captcha on job submission, shown only once a client goes over the limit, if bots
  are still a problem after the limits above.
- Replace `pipes.quote` with `shlex.quote`: the `pipes` module was removed in
  Python 3.13.
- `html.escape()` on all CGI values printed in pages and in `errAbort()` messages.
