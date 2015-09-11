#!/usr/bin/env python
# the tefor crispr tool
# can be run as a CGI or from the command line

# todo: check symlink to /var/www

# python std library
import subprocess, tempfile, optparse, logging, atexit, glob, shutil
import Cookie, time, math, sys, cgi, re, array, random, platform, os
import hashlib, base64, string, logging, operator, urllib, sqlite3, time
import traceback, json, pwd, types

from datetime import datetime
from collections import defaultdict, namedtuple
from sys import stdout
from os.path import join, isfile, basename, dirname, getmtime, isdir
from StringIO import StringIO

# out own eff scoring library
import crisporEffScores

# don't report print as an error
# pylint: disable=E1601

# optional module for Excel export as native .xls files
# install with 'apt-get install python-xlwt' or 'pip install xlwt'
xlwtLoaded = True
try:
    import xlwt
except:
    sys.stderr.write("crispor.cgi - warning - the python xlwt module is not available\n")
    xlwtLoaded = False

# write debug output to stdout
DEBUG = False
#DEBUG = True

# system-wide temporary directory
TEMPDIR = os.environ.get("TMPDIR", "/tmp")
#TEMPDIR = "/dev/shm/"

# prefix in html statements before the directories "image/", "style/" and "js/" 
HTMLPREFIX =  ""
# alternative directory on local disk where image/, style/ and js/ are located
HTMLDIR = "/usr/local/apache/htdocs/crispor/"

# directory of crispor.cgi
baseDir = dirname(__file__)

# the segments.bed files use abbreviated genomic region names
segTypeConv = {"ex":"exon", "in":"intron", "ig":"intergenic"}

# directory for processed batches of offtargets ("cache" of bwa results)
batchDir = join(baseDir,"temp")

# the file where the job queue is stored
JOBQUEUEDB = join(baseDir, "jobs.db")

# directory for platform-independent scripts (e.g. Heng Li's perl SAM parser)
scriptDir = join(baseDir, "bin")

# directory for helper binaries (e.g. BWA)
binDir = join(baseDir, "bin", platform.system())

# directory for genomes
genomesDir = join(baseDir, "genomes")

DEFAULTORG = 'hg19'
DEFAULTSEQ = 'cttcctttgtccccaatctgggcgcgcgccggcgccccctggcggcctaaggactcggcgcgccggaagtggccagggcgggggcgacctcggctcacagcgcgcccggctattctcgcagctcaccatgGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCC'

pamDesc = [ ('NGG','NGG - Streptococcus Pyogenes'),
         ('NGA','NGA - S. Pyogenes mutant VQR'),
         ('NGCG','NGCG - S. Pyogenes mutant VRER'),
         ('NNAGAA','NNAGAA - Streptococcus Thermophilus'),
         ('NNGRRT','NNG(A/G)(A/G)T - Streptococcus Aureus'),
         ('NNNNGMTT','NNNNG(A/C)TT - Neisseiria Meningitidis'),
         ('NNNNACA','NNNNACA - Campylobacter jejuni'),
       ]

DEFAULTPAM = 'NGG'

# maximum size of an input sequence 
MAXSEQLEN = 1000

# BWA: allow up to X mismatches
maxMMs=4
# for guides that have a gap, allow fewer mismatches
maxGapMMs=2

# BWA: allow a single gap in the alignment?
allowGap = False

# maximum number of occurences in the genome to get flagged as repeats. 
# This is used in bwa samse, when converting the same file
# and for warnings in the table output.
MAXOCC = 40000

# the BWA queue size is 2M by default. We derive the queue size from MAXOCC
MFAC = 2000000/MAXOCC

# Highly-sensitive mode: MAXOCC is increased in processSubmission() and in the html UI if only one guide seq
# is run
HIGH_MAXOCC=600000

# minimum off-target score of standard off-targets (those that end with NGG)
#MINSCORE = 1.0
# minimum off-target score for alternative PAM off-targets
ALTPAMMINSCORE = 1.0

# for some PAMs, we allow other alternative motifs when searching for offtargets
# MIT and eCrisp do that, they use the motif NGG + NAG, we add one more, based on the
# on the guideSeq results in Tsai et al, Nat Biot 2014
# The NGA -> NGG rule was described by Kleinstiver...Young 2015 "Improved Cas9 Specificity..."
offtargetPams = {"NGG" : "NAG,NGA", "NGA" : "NGG" }

# global flag to indicate if we're run from command line or as a CGI
commandLineMode = False

# names/order of efficacy scores to show in UI
scoreNames = ["chariRank", "ssc", "doench", "wang", "crisprScan", "oof"]
# how many digits shall we show for each score? default is 0
scoreDigits = {
    "ssc" : 2,
}

# labels and descriptions of eff. scores
scoreDescs = {
    "doench" : ("Doench", "Doench et al score. Ranges from 0-100, higher values are better"),
    "ssc" : ("SSC", "SSC, Xu et al. score. Ranges mostly -2 to +2. Should be positive."),
    "crisprScan" : ("CrisprScan", "CrisprScan Zebrafish score by Moreno-Mateos et al"),
    "wang" : ("Wang", "SVM score by Wang/Wei/Sabatini/Lander, range 0-100"),
    "chariRank" : ("Chari", "Chari et al. SVM Rank-Percent, range: 0-100"),
    "oof" : ("Out-of-Frame", "Bae et al. Out-of-Frame score, range: 0-100")
}

# the headers for the guide and offtarget output files
guideHeaders = ["guideId", "guideSeq", "specScore", "effScore", "offtargetCount", "guideGenomeMatchGeneLocus"]
offtargetHeaders = ["guideId", "guideSeq", "offtargetSeq", "mismatchCount", "offtargetScore", "chrom", "start", "end", "locusDesc"]
# ====== END GLOBALS ============


# ==== CLASSES =====
class JobQueue:
    """
    simple job queue, using a db table as a backend
    jobs have different types and status. status can be updated while they run
    job running times are kept and old job info is kept in a separate table
    
    >>> q = JobQueue("/tmp/tempCrisporTest.db")
    >>> q.clearJobs()
    >>> q.waitCount()
    0
    >>> q.addJob("search", "abc123", "myParams")
    True

    only one job per jobId
    >>> q.addJob("search", "abc123", "myParams")
    False
    >>> q.waitCount()
    1
    >>> q.getStatus("abc123")
    u'Waiting'
    >>> q.startStep("abc123", "bwa", "Alignment with BWA")
    >>> q.getStatus("abc123")
    u'Alignment with BWA'
    >>> jobType, jobId, paramStr = q.popJob()

    >>> q.waitCount()
    0
    >>> q.jobDone("abc123")
    >>> q.waitCount()
    0

    can't pop from an empty queue
    >>> q.popJob()
    (None, None, None)
    >>> os.system("rm /tmp/tempCrisporTest.db")
    0
    """

    _queueDef = (
    'CREATE TABLE IF NOT EXISTS %s '
    '('
    '  jobType text,' # either "index" or "search"
    '  jobId text %s,' # unique identifier
    '  paramStr text,' # parameters for jobs, like db, options, etc. 
    '  isRunning int DEFAULT 0,' # indicates steps have started, done jobs are moved to doneJobs table
    '  stepName text,' # currently step, internal step name for timings
    '  stepLabel text,' # current step, human-readable status of job, for UI
    '  lastUpdate float,' # time of last update
    '  stepTimes text,' # comma-sep list of whole msecs, one per step
    '  startTime text' # date+time when job was put into queue
    ')')

    def __init__(self, dbName):
        self.dbName = dbName
        self.conn = sqlite3.Connection(dbName)
        self.conn.execute(self._queueDef % ("queue", "PRIMARY KEY"))
        self.conn.execute(self._queueDef % ("doneJobs", ""))
        self.conn.commit()
 
    def addJob(self, jobType, jobId, paramStr):
        " create a new job, returns False if not successful  "
        sql = 'INSERT INTO queue (jobType, jobId, isRunning, lastUpdate, ' \
            'stepTimes, paramStr, stepName, stepLabel, startTime) VALUES (?, ?, ?, ?, ?, ?, ?, ?, datetime("now"))'
        now = "%.3f" % time.time()
        try:
            self.conn.execute(sql, (jobType, jobId, 0, now, "", paramStr, "wait", "Waiting"))
            self.conn.commit()
            return True
        except sqlite3.IntegrityError:
            return False
        except sqlite3.OperationalError:
            errAbort("Cannot open DB file %s. Please contact services@tefor.net" % self.dbName)

    def getStatus(self, jobId):
        " return current job status label or None if job is not in queue"
        sql = 'SELECT stepLabel FROM queue WHERE jobId=?'
        try:
            status = self.conn.execute(sql, (jobId,)).next()[0]
        except StopIteration:
            status = None
        return status

    def dump(self):
        " for debugging, write the whole queue table to stdout "
        sql = 'SELECT * FROM queue'
        for row in self.conn.execute(sql):
            print "\t".join([str(x) for x in row])

    def jobInfo(self, jobId, isDone=False):
        " for debugging, return all job info as a tuple "
        if isDone:
            sql = 'SELECT * FROM doneJobs WHERE jobId=?'
        else:
            sql = 'SELECT * FROM queue WHERE jobId=?'
        try:
            row = self.conn.execute(sql, (jobId,)).next()
        except StopIteration:
            return []
        return row

    def startStep(self, jobId, newName, newLabel):
        " start a new step. Update lastUpdate, status and stepTime "
        self.conn.execute('BEGIN IMMEDIATE') # lock db
        sql = 'SELECT lastUpdate, stepTimes, stepName FROM queue WHERE jobId=?'
        lastTime, timeStr, lastStep = self.conn.execute(sql, (jobId,)).next()
        lastTime = float(lastTime)

        # append a string in format "stepName:milliSecs" to the timeStr
        now = time.time()
        timeDiff = "%d" % int((1000.0*(now - lastTime)))
        newTimeStr = timeStr+"%s=%s" % (lastStep, timeDiff)+","

        sql = 'UPDATE queue SET lastUpdate=?, stepName=?, stepLabel=?, stepTimes=?, isRunning=? WHERE jobId=?'
        self.conn.execute(sql, (now, newName, newLabel, newTimeStr, 1, jobId))
        self.conn.commit()

    def jobDone(self, jobId):
        " remove the job from the queue and add it to the queue log"
        self.conn.execute('BEGIN IMMEDIATE') # lock db
        sql = 'SELECT * FROM queue WHERE jobId=?'
        try:
            row = self.conn.execute(sql, (jobId,)).next()
        except StopIteration:
            # return if the job has already been removed
            logging.warn("jobDone - jobs %s has been removed already" % jobId)
            return

        sql = 'DELETE FROM queue WHERE jobId=?'
        self.conn.execute(sql, (jobId,))

        sql = 'INSERT INTO doneJobs VALUES (?,?,?,?,?,?,?,?,?)'
        self.conn.execute(sql, row)
        self.conn.commit()

    def waitCount(self):
        " return number of waiting jobs "
        sql = 'SELECT count(*) FROM queue WHERE isRunning=0'
        count = self.conn.execute(sql).next()[0]
        return count

    def popJob(self):
        " return (jobType, jobId, params) of first waiting job and set it to running state "
        self.conn.execute('BEGIN IMMEDIATE') # lock db
        sql = 'SELECT jobType, jobId, paramStr FROM queue WHERE isRunning=0 ORDER BY lastUpdate LIMIT 1'
        try:
            jobType, jobId, paramStr = self.conn.execute(sql).next()
        except StopIteration:
            self.conn.commit() # unlock db
            return None, None, None

        sql = 'UPDATE queue SET isRunning=1 where jobId=?'
        self.conn.execute(sql, (jobId,))
        self.conn.commit() # unlock db
        return jobType, jobId, paramStr

    def clearJobs(self):
        " clear the job table, removing running jobs, too "
        self.conn.execute("DELETE from queue")
        #self.conn.execute("DROP TABLE queue")
        self.conn.commit()

# ====== FUNCTIONS =====

def getParams():
    " get CGI parameters and return as dict "
    form = cgi.FieldStorage()
    params = {}

    #for key in ["pamId", "batchId", "pam", "seq", "org", "showAll", "download", "sortBy", "format", "ajax]:
    for key in form.keys():
        val = form.getfirst(key)
	if val!=None:
            params[key] = val

    if "pam" in params:
        if len(set(params["pam"])-set("ACTGNMKRY"))!=0:
            errAbort("Illegal character in PAM-sequence. Only ACTGMKRY and N allowed.")
    return params

transTab = string.maketrans("-=/+_", "abcde")

def makeTempBase(seq, org, pam):
    "create the base name of temp files using a hash function and some prettyfication "
    hasher = hashlib.sha1(seq+org+pam)
    batchId = base64.urlsafe_b64encode(hasher.digest()[0:20]).translate(transTab)[:20]
    return batchId

def saveSeqOrgPamToCookies(seq, org, pam):
    " create a cookie with seq, org and pam and print it"
    cookies=Cookie.SimpleCookie()
    expires = 365 * 24 * 60 * 60
    cookies['lastseq'] = seq
    cookies['lastseq']['expires'] = expires
    cookies['lastorg'] = org
    cookies['lastorg']['expires'] = expires
    cookies['lastpam'] = pam
    cookies['lastpam']['expires'] = expires
    print cookies

def debug(msg):
    if commandLineMode:
        logging.debug(msg)
    elif DEBUG:
        print msg
        print "<br>"

def errAbort(msg):
    " print err msg and exit "
    if commandLineMode:
        raise Exception(msg)

    print('<div style="float:left; text-align:left; width: 800px">')
    print(msg+"<p>")
    print('</div>')
    sys.exit(0)  # cgi must not exit with 1

def matchNuc(pat, nuc):
    " returns true if pat (single char) matches nuc (single char) "
    if pat in ["A", "C", "T", "G"] and pat==nuc:
        return True
    elif pat=="M" and nuc in ["A", "C"]:
        return True
    elif pat=="K" and nuc in ["T", "G"]:
        return True
    elif pat=="R" and nuc in ["A", "G"]:
        return True
    elif pat=="Y" and nuc in ["C", "T"]:
        return True
    else:
        return False

def gcContent(seq):
    " return GC content as a float "
    c = 0
    for x in seq:
        if x in ["G","C"]:
            c+= 1
    return (float(c)/len(seq))
            
def findPat(seq, pat):
    """ yield positions where pat matches seq, stupid brute force search 
    """
    for i in range(0, len(seq)-len(pat)+1):
        #print "new pos", i, seq[i:i+len(pat)],"<br>"
        found = True
        for x in range(0, len(pat)):
            #print "new step", x, "<br>"
            if pat[x]=="N":
                #print "N","<br>"
                continue
            seqPos = i+x
            if seqPos == len(seq):
                found = False
                break
            if not matchNuc(pat[x], seq[seqPos]):
                #print i, x, pat[x], seq[seqPos], "no match<br>"
                found = False
                break
            #print "match", i, x, found, "<br>"
        if found:
            #print "yielding", i, "<br>"
            yield i

def rndSeq(seqLen):
    " return random seq "
    seq = []
    alf = "ACTG"
    for i in range(0, seqLen):
        seq.append(alf[random.randint(0,3)])
    return "".join(seq)

def cleanSeq(seq):
    """ remove fasta header, check seq for illegal chars and return (filtered seq, user message) 
    special value "random" returns a random sequence.
    """
    #print repr(seq)
    if seq.startswith("random"):
        seq = rndSeq(800)
    lines = seq.strip().splitlines()
    #print "<br>"
    #print "before fasta cleaning", "|".join(lines)
    if len(lines)>0 and lines[0].startswith(">"):
        line1 = lines.pop(0)
    #print "<br>"
    #print "after fasta cleaning", "|".join(lines)
    #print "<br>"

    newSeq = []
    nCount = 0
    for l in lines:
        if len(l)==0:
            continue
        for c in l:
            if c not in "actgACTG":
                nCount +=1 
            else:
                newSeq.append(c)
    seq = "".join(newSeq)

    msgs = []
    if len(seq)>MAXSEQLEN:
        msgs.append("<strong>Sorry, this tool cannot handle sequences longer than 1kbp</strong><br>Below you find the results for the first %d bp of your input sequence.<br>" % MAXSEQLEN)
        seq = seq[:MAXSEQLEN]

    if nCount!=0:
        msgs.append("Sequence contained %d non-ACTG letters. They were removed." % nCount)

    return seq, "<br>".join(msgs)

def revComp(seq):
    " rev-comp a dna sequence with UIPAC characters "
    revTbl = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N' , 'M' : 'K', 'K' : 'M', 
        "R" : "Y" , "Y":"R" }
    newSeq = []
    for c in reversed(seq):
        newSeq.append(revTbl[c])
    return "".join(newSeq)

def findPams(seq, pam, strand, startDict, endSet):
    """ return two values: dict with pos -> strand of PAM and set of end positions of PAMs
    Makes sure to return only values with at least 20 bp left (if strand "+") or to the 
    right of the match (if strand "-")
    >>> findPams("GGGGGGGGGGGGGGGGGGGGGGG", "NGG", "+", {}, set())
    ({20: '+'}, set([23]))
    >>> findPams("CCAGCCCCCCCCCCCCCCCCCCC", "CCA", "-", {}, set())
    ({0: '-'}, set([3]))

    """

    # -------------------
    #          OKOKOKOKOK
    minPosPlus  = 20
    # -------------------
    # OKOKOKOKOK
    maxPosMinus = len(seq)-(20+len(pam))

    #print "new search", seq, pam, "<br>"
    for start in findPat(seq, pam):
        # need enough flanking seq on one side
        #print "found", start,"<br>"
        if strand=="+" and start<minPosPlus:
            #print "no, out of bounds +", "<br>"
            continue
        if strand=="-" and start>maxPosMinus:
            #print "no, out of bounds, -<br>"
            continue

        #print "match", strand, start, end, "<br>"
        startDict[start] = strand
        end = start+len(pam)
        endSet.add(end)
    return startDict, endSet

def rulerString(maxLen):
    " return line with positions every 10 chars "
    texts = []
    for i in range(0, maxLen, 10):
        numStr = str(i)
        texts.append(numStr)
        spacer = "".join([" "]*(10-len(numStr)))
        texts.append(spacer)
    return "".join(texts)

def showSeqAndPams(seq, startDict, pam, guideScores):
    " show the sequence and the PAM sites underneath in a sequence viewer "
    lines, maxY = distrOnLines(seq.upper(), startDict, len(pam))

    print "<div class='substep'>"
    print '<a id="seqStart"></a>'
    print "Found %d possible guide sequences in input (%d bp). Click on a PAM %s match to show its guide sequence.<br>" % (len(guideScores), len(seq), pam)
    print "Shown below are the PAM site and the nucleotide at position -3 5' of it.<br>"
    print '''Colors <span style="color:#32cd32; text-shadow: 1px 1px 1px #bbb">green</span>, <span style="color:#ffff00; text-shadow: 1px 1px 1px #888">yellow</span> and <span style="text-shadow: 1px 1px 1px #f01; color:#aa0014">red</span> indicate high, medium and low specificity of the PAM's guide sequence in the genome.'''
    print "</div>"
    print '''<div style="text-align: left; overflow-x:scroll; width:100%; background:#DDDDDD; border-style: solid; border-width: 1px">'''

    print '<pre style="font-size: 80%; display:inline; line-height: 0.95em; text-align:left">'+rulerString(len(seq))
    print seq

    for y in range(0, maxY+1):
        #print "y", y, "<br>"
        texts = []
        lastEnd = 0
        for start, end, name, strand, pamId  in lines[y]:
            spacer = "".join([" "]*((start-lastEnd)))
            lastEnd = end
            texts.append(spacer)
            score = guideScores[pamId]
            color = scoreToColor(score)

            #print score, opacity
            #texts.append('''<a style="text-shadow: 1px 1px 1px #bbb; color: %s" id="list%s" href="#%s" onmouseover="$('.hiddenExonMatch').show('fast');$('#show-more').hide();$('#show-less').show()" onfocus="window.location.href = '#seqStart'" >''' % (color, pamId,pamId))
            texts.append('''<a style="text-shadow: 1px 1px 1px #bbb; color: %s" id="list%s" href="#%s">''' % (color, pamId,pamId))
            texts.append(name)
            texts.append("</a>")
        print "".join(texts)
    print("</pre><br>")

    print '''</div>'''
    
def iterOneDelSeqs(seq):
    """ given a seq, create versions with each bp removed. Avoid duplicates 
    yields (delPos, seq)
    >>> list(iterOneDelSeqs("AATGG"))
    [(0, 'ATGG'), (2, 'AAGG'), (3, 'AATG')]
    """
    doneSeqs = set()
    for i in range(0, len(seq)):
        delSeq = seq[:i]+seq[i+1:]
        if delSeq not in doneSeqs:
            yield i, delSeq
        doneSeqs.add(delSeq)

def flankSeqIter(seq, startDict, pamLen, gapped=False):
    """ given a seq and dictionary of pos -> strand and the length of the pamSite
    yield tuples of (name, startPos, strand, flankSeq, pamSeq)
    """
    startList = sorted(startDict.keys())
    for startPos in startList:
        strand = startDict[startPos]

        if strand=="+":
            flankSeq = seq[startPos-20:startPos]
            pamSeq = seq[startPos:startPos+pamLen]
        else: # strand is minus
            flankSeq = revComp(seq[startPos+pamLen:startPos+pamLen+20])
            pamSeq = revComp(seq[startPos:startPos+pamLen])

        if "N" in flankSeq:
            continue

        if not gapped:
            yield "s%d%s" % (startPos, strand), startPos, strand, flankSeq, pamSeq
        else:
            # construct 20 sequences, each one with one nucleotide deleted
            # but always adding one more nucleotide 5' of the whole guide
            # if the sequence is not long enough, adding an "A", not elegant,
            # but this keeps the length at 20 bp which makes everything easier
            #if startPos>20:
                #firstNucl = seq[startPos-21]
            #else:
                #firstNucl = "A"
            for i, gapSeq in iterOneDelSeqs(flankSeq):
                pamId = "s%d%sg%d" % (startPos, strand, i)
                #yield pamId, startPos, strand, firstNucl+flankSeq[:i]+flankSeq[i+1:], pamSeq
                yield pamId, startPos, strand, gapSeq, pamSeq

def makeBrowserLink(dbInfo, pos, text, title, cssClasses=[]):
    " return link to genome browser (ucsc or ensembl) at pos, with given text "
    if dbInfo.server.startswith("Ensembl"):
        baseUrl = "www.ensembl.org"
        if dbInfo.server=="EnsemblPlants":
            baseUrl = "plants.ensembl.org"
        elif dbInfo.server=="EnsemblMetazoa":
            baseUrl = "metazoa.ensembl.org"
        elif dbInfo.server=="EnsemblProtists":
            baseUrl = "protists.ensembl.org"
        org = dbInfo.scientificName.replace(" ", "_")
        url = "http://%s/%s/Location/View?r=%s" % (baseUrl, org, pos)
    elif dbInfo.server=="ucsc":
        if pos[0].isdigit():
            pos = "chr"+pos
        url = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s" % (dbInfo.name, pos)
    else:
        #return "unknown genome browser server %s, please email services@tefor.net" % dbInfo.server
        url = "javascript:void(0)"

    classStr = ""
    if len(cssClasses)!=0:
        classStr = ' class="%s"' % (" ".join(cssClasses))

    return '''<a title="%s"%s target="_blank" href="%s">%s</a>''' % (title, classStr, url, text)

def makeAlnStr(seq1, seq2, pam, score, posStr):
    " given two strings of equal length, return a html-formatted string that highlights the differences "
    lines = [ [], [], [] ]
    last12MmCount = 0
    for i in range(0, len(seq1)-len(pam)):
        if seq1[i]==seq2[i]:
            lines[0].append(seq1[i])
            lines[1].append(seq2[i])
            lines[2].append(" ")
        else:
            lines[0].append("<b>%s</b>"%seq1[i])
            lines[1].append("<b>%s</b>"%seq2[i])
            lines[2].append("*")
            if i>7:
                last12MmCount += 1
    lines[0].append(" <i>"+seq1[-len(pam):]+"</i>")
    lines[1].append(" <i>"+seq2[-len(pam):]+"</i>")
    lines = ["".join(l) for l in lines]

    if len(posStr)>1 and posStr[0].isdigit():
        posStr = "chr"+posStr

    htmlText = "<small><pre>guide:      %s<br>off-target: %s<br>            %s</pre>Off-target score: %.2f<br>Position: %s</small>" % (lines[0], lines[1], lines[2], score, posStr)
    hasLast12Mm = last12MmCount>0
    return htmlText, hasLast12Mm
        
def parsePos(text):
    " parse a string of format chr:start-end:strand and return a 4-tuple "
    if text!=None and len(text)!=0 and text!="?":
        fields = text.split(":")
        if len(fields)==2:
            chrom, posRange = fields
            strand = "+"
        else:
            chrom, posRange, strand = fields
        start, end = posRange.split("-")
        start, end = int(start), int(end)
    else:
        chrom, start, end, strand = "", 0, 0, "+"
    return chrom, start, end, strand

def makePosList(countDict, guideSeq, pam, inputPos):
    """ for a given guide sequence, return a list of tuples that 
    describes the offtargets sorted by score and a string to describe the offtargets in the
    format x/y/z/w of mismatch counts
    inputPos has format "chrom:start-end:strand". All 0MM matches in this range
    are ignored from scoring ("ontargets")
    Also return the same description for just the last 12 bp and the score 
    of the guide sequence (calculated using all offtargets).
    """
    inChrom, inStart, inEnd, inStrand = parsePos(inputPos)
    count = 0
    otCounts = []
    posList = []
    scores = []
    last12MmCounts = []
    ontargetDesc = ""
    subOptMatchCount = 0

    # for each edit distance, get the off targets and iterate over them
    foundOneOntarget = False
    for editDist in range(0, maxMMs+1):
        #print countDict,"<p>"
        matches = countDict.get(editDist, [])

        #print otCounts,"<p>"
        last12MmOtCount = 0

        # create html and score for every offtarget
        otCount = 0
        for chrom, start, end, otSeq, strand, segType, geneNameStr, x1Count, gapPos in matches:
            # skip on-targets
            if segType!="":
                segTypeDesc = segTypeConv[segType]
                geneDesc = segTypeDesc+":"+geneNameStr
                geneDesc = geneDesc.replace("|", "-")
            else:
                geneDesc = geneNameStr

            # is this the on-target?
            # if we got a genome position, use it. Otherwise use a random off-target with 0MMs
            # as the on-target
            if editDist==0 and \
                ((chrom==inChrom and start >= inStart and end <= inEnd and x1Count < MAXOCC) \
                or (inChrom=='' and foundOneOntarget==False and x1Count < MAXOCC)):
                foundOneOntarget = True
                ontargetDesc = geneDesc
                continue

            otCount += 1
            guideNoPam = guideSeq[:len(guideSeq)-len(pam)]
            otSeqNoPam = otSeq[:len(otSeq)-len(pam)]
            if len(otSeqNoPam)==19:
                otSeqNoPam = "A"+otSeqNoPam # should not change the score a lot, weight0 is very low
            score = calcHitScore(guideNoPam, otSeqNoPam)
            scores.append(score)

            posStr = "%s:%d-%s" % (chrom, int(start)+1,end)
            alnHtml, hasLast12Mm = makeAlnStr(guideSeq, otSeq, pam, score, posStr)
            if not hasLast12Mm:
                last12MmOtCount+=1
            posList.append( (otSeq, score, editDist, posStr, geneDesc, alnHtml, gapPos) )
            # taking the maximum is probably not necessary, 
            # there should be only one offtarget for X1-exceeding matches
            subOptMatchCount = max(int(x1Count), subOptMatchCount)

        last12MmCounts.append(str(last12MmOtCount))
        # create a list of number of offtargets for this edit dist
        otCounts.append( str(otCount) )

    if subOptMatchCount > MAXOCC:
        guideScore = 0
        posList = []
        ontargetDesc = ""
        last12DescStr = ""
        otDescStr = ""
    else:
        guideScore = calcMitGuideScore(sum(scores))
        otDescStr = "&thinsp;-&thinsp;".join(otCounts)
        last12DescStr = "&thinsp;-&thinsp;".join(last12MmCounts)

    posList.sort(reverse=True, key=operator.itemgetter(1)) # sort by offtarget score

    return posList, otDescStr, guideScore, last12DescStr, ontargetDesc, subOptMatchCount

# --- START OF SCORING ROUTINES 

# MIT offtarget scoring

# aka Matrix "M"
hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]

def calcHitScore(string1,string2):
    " see 'Scores of single hits' on http://crispr.mit.edu/about "
    # The Patrick Hsu weighting scheme
    assert(len(string1)==len(string2)==20)

    dists = [] # distances between mismatches, for part 2
    mmCount = 0 # number of mismatches, for part 3
    lastMmPos = None # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(0, len(string1)):
        if string1[pos]!=string2[pos]:
            mmCount+=1
            if lastMmPos!=None:
                dists.append(pos-lastMmPos)
            score1 *= 1-hitScoreM[pos]
            lastMmPos = pos
    # 2nd part of the score
    if mmCount<2: # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists)/len(dists)
        score2 = 1.0 / (((19-avgDist)/19.0) * 4 + 1)
    # 3rd part of the score
    if mmCount==0: # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount**2)

    score = score1 * score2 * score3 * 100
    return score

def calcMitGuideScore(hitSum):
    """ Sguide defined on http://crispr.mit.edu/about 
    Input is the sum of all off-target hit scores. Returns the specificity of the guide.
    """
    score = 100 / (100+hitSum)
    score = int(round(score*100))
    return score

# --- END OF SCORING ROUTINES 

def parseChromSizes(fname):
    " return chrom sizes as dict chrom -> size "
    ret = {}
    for line in open(fname).read().splitlines():
        fields = line.split()
        chrom, size = fields[:2]
        ret[chrom] = int(size)
    return ret

def extendAndGetSeq(db, chrom, start, end, strand, flank=100):
    """ extend (start, end) by flank and get sequence for it using twoBitTwoFa.
    Return None if not possible to extend.
    >>> extendAndGetSeq("hg19", "chr21", 10000000, 10000005, "+", flank=3)
    'AAGGAATGTAG'
    """
    genomeDir = genomesDir
    sizeFname = "%(genomeDir)s/%(db)s/%(db)s.sizes" % locals()
    chromSizes = parseChromSizes(sizeFname)
    maxEnd = chromSizes[chrom]+1

    start -= flank
    end += flank
    if start < 0 or end > maxEnd:
        return None

    twoBitFname = "%(genomeDir)s/%(db)s/%(db)s.2bit" % locals()
    progDir = binDir
    genome = db
    cmd = "%(progDir)s/twoBitToFa %(genomeDir)s/%(genome)s/%(genome)s.2bit stdout -seq=%(chrom)s -start=%(start)s -end=%(end)s" % locals()
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    seqStr = proc.stdout.read()
    faFile = StringIO(seqStr)
    seqs = parseFasta(faFile)
    assert(len(seqs)==1)
    seq = seqs.values()[0].upper()
    if strand=="-":
        seq = revComp(seq)
    return seq

def getExtSeq(seq, start, end, strand, extUpstream, extDownstream, extSeq=None, extFlank=100):
    """ extend (start,end) by extUpstream and extDownstream and return the subsequence
    at this position in seq.
    Return None if there is not enough space to extend (start, end).
    extSeq is a sequence with extFlank additional flanking bases on each side. It can be provided
    optionally and is used if needed to return a subseq.
    Careful: returned sequence might contain lowercase letters.
    >>> getExtSeq("AACCTTGG", 2, 4, "+", 2, 4)
    'AACCTTGG'
    >>> getExtSeq("CCAACCTTGGCC", 4, 6, "-", 2, 3)
    'AAGGTTG'
    >>> getExtSeq("AA", 0, 2, "+", 2, 3)
    >>> getExtSeq("AA", 0, 2, "+", 2, 3, extSeq="CAGAATGA", extFlank=3)
    'AGAATGA'
    >>> getExtSeq("AA", 0, 2, "-", 2, 3, extSeq="CAGAATGA", extFlank=3)
    'CATTCTG'
    """
    assert(start>=0)
    assert(end<=len(seq))
    # check if the extended sequence really contains the whole input seq 
    # e.g. when user has added nucleotides to a otherwise matching sequence
    if extSeq!=None and (seq not in extSeq):
        debug("seq is not in extSeq")
        extSeq = None

    # extend
    if strand=="+":
        extStart, extEnd = start-extUpstream, end+extDownstream
    else:
        extStart, extEnd = start-extDownstream, end+extUpstream

    # check for out of bounds and get seq
    if extStart >= 0 and extEnd <= len(seq):
        subSeq = seq[extStart:extEnd]
    else:
        if extSeq==None:
            return None
        # lift to extSeq coords and get seq
        extStart += extFlank
        extEnd += extFlank
        assert(extStart >= 0)
        assert(extEnd <= len(extSeq))
        subSeq = extSeq[extStart:extEnd]

    if strand=="-":
        subSeq = revComp(subSeq)

    return subSeq

def pamStartToGuideRange(startPos, strand, pamLen):
    """ given a PAM start position and its strand, return the (start,end) of the guide.
    Coords can be negative or exceed the length of the input sequence.
    """
    if strand=="+":
        return (startPos-20, startPos)
    else: # strand is minus
        return (startPos+pamLen, startPos+pamLen+20)

def htmlHelp(text):
    " show help text with tooltip or modal dialog "
    print '''<img style="height:1.1em; width:1.0em" src="%simage/info-small.png" class="help tooltipster" title="%s" />''' % (HTMLPREFIX, text)

def htmlWarn(text):
    " show help text with tooltip "
    print '''<img style="height:1.1em; width:1.0em" src="%simage/warning-32.png" class="help tooltipster" title="%s" />''' % (HTMLPREFIX, text)

def readEnzymes():
    " parse restrSites.txt and return as dict length -> list of (name, seq) "
    fname = "restrSites.txt"
    enzList = {}
    for line in open(join(baseDir, fname)):
        name, seq1, seq2 = line.split()
        seq = seq1+seq2
        enzList.setdefault(len(seq), []).append( (name, seq) )
    return enzList
        
def patMatch(seq, pat, notDegPos=None):
    """ return true if pat matches seq, both have to be same length 
    do not match degenerate codes at position notDegPos (0-based)
    """
    assert(len(seq)==len(pat))
    for x in range(0, len(pat)):
        patChar = pat[x]
        nuc = seq[x]

        assert(patChar in "MKYRACTGN")
        assert(nuc in "MKYRACTGN")

        if notDegPos!=None and x==notDegPos and patChar!=nuc:
            #print x, seq, pat, notDegPos, patChar, nuc, "<br>"
            return False

        if patChar=="N":
            continue
        if patChar=="M" and nuc in ["A", "C"]:
            continue
        if patChar=="K" and nuc in ["T", "G"]:
            continue
        if patChar=="R" and nuc in ["A", "G"]:
            continue
        if patChar=="Y" and nuc in ["C", "T"]:
            continue
        if patChar!=nuc:
            return False
    return True

def findSite(seq, restrSite):
    """ return the positions where restrSite matches seq 
    seq can be longer than restrSite
    Do not allow degenerate characters to match at position len(restrSite) in seq
    """
    posList = []
    for i in range(0, len(seq)-len(restrSite)+1):
        subseq = seq[i:i+len(restrSite)]
        #print subseq==restrSite, subseq, restrSite,"<br>"

        # JP does not want any potential site to be suppressed
        #if i<len(restrSite):
            #isMatch = patMatch(subseq, restrSite, len(restrSite)-i-1)
        #else:
            #isMatch = patMatch(subseq, restrSite)
        isMatch = patMatch(subseq, restrSite)

        if isMatch:
            posList.append( (i, i+len(restrSite)) )
    return posList

def matchRestrEnz(allEnzymes, guideSeq, pamSeq):
    """ return list of enzymes that overlap the -3 position in guideSeq
    returns dict name -> list of matching positions
    """
    matches = {}
    #print guideSeq, pamSeq, "<br>"
    fullSeq = guideSeq+pamSeq
    for siteLen, sites in allEnzymes.iteritems():
        startSeq = len(fullSeq)-len(pamSeq)-3-(siteLen)+1
        seq = fullSeq[startSeq:]
        #print "restrEnz with len %d"% siteLen, seq, "<br>"
        for name, restrSite in sites:
            posList = findSite(seq, restrSite)
            if len(posList)!=0:
                liftOffset = startSeq+len(guideSeq)
                posList = [(liftOffset+x, liftOffset+y) for x,y in posList]
                matches.setdefault(name, []).extend(posList)
    return matches

def mergeGuideInfo(seq, startDict, pamPat, otMatches, inputPos, effScores, sortBy=None):
    """ 
    merges guide information from the sequence, the efficacy scores and the off-targets.
    creates rows with fields:


    for each pam in startDict, retrieve the guide sequence next to it and score it
    sortBy can be "effScore" or "mhScore" or "oofScore"
    """
    allEnzymes = readEnzymes()

    guideData = []
    guideScores = {}
    hasNotFound = False

    pamSeqs = list(flankSeqIter(seq, startDict, len(pamPat)))
    # make gapped and ungapped versions of the pam sequences
    #if allowGap:
        #pamSeqs.extend( flankSeqIter(seq, startDict, len(pamPat), gapped=True) )

    for pamId, startPos, strand, guideSeq, pamSeq in pamSeqs:
        # matches in genome
        # one desc in last column per OT seq
        if pamId in otMatches:
            pamMatches = otMatches[pamId]
            guideSeqFull = guideSeq + pamSeq
            mutEnzymes = matchRestrEnz(allEnzymes, guideSeq, pamSeq)
            posList, otDesc, guideScore, last12Desc, ontargetDesc, subOptMatchCount =\
                makePosList(pamMatches, guideSeqFull, pamPat, inputPos)

            #gStart, gEnd = pamStartToGuideRange(startPos, strand, len(pamPat))
            #seq34Mer = getExtSeq(seq, gStart, gEnd, strand, 4, 10, extSeq)
            #if seq34Mer==None:
                #effScore = "Too close to end"
            #else:
                #effScore = int(round(100*calcDoenchScore(seq34Mer[:30])))

            # get 60mer and calc oof-score
            #seq60Mer = getExtSeq(seq, gStart, gEnd, strand, 20, 20, extSeq)
            #if seq60Mer == None:
                #mhScore, oofScore = "Too close to end", ""
            #else:
                #assert(len(seq60Mer)==60)
                #mhScore, oofScore = calcMicroHomolScore(seq60Mer, 30)
                #oofScore = str(oofScore)+" %"

        # no off-targets found?
        else:
            posList, otDesc, guideScore = None, "Not found", None
            last12Desc = ""
            #effScore = 0
            hasNotFound = True
            mutEnzymes = []
            ontargetDesc = ""
            #mhScore, oofScore = 0, 0
            subOptMatchCount = False
            seq34Mer = None

        #gStart, gEnd = pamStartToGuideRange(startPos, strand, len(pamPat))
        #seq100Mer = getExtSeq(seq, gStart, gEnd, strand, 30, 50, extSeq)

        guideData.append( [guideScore, effScores.get(pamId, {}), startPos, strand, pamId, guideSeq, pamSeq, posList, otDesc, last12Desc, mutEnzymes, ontargetDesc, subOptMatchCount] )
        guideScores[pamId] = guideScore

    #if USESSC:
        #guideData = addSscScores(guideData, seqFieldIdx=-1, scoreFieldIdx=4)
    #else:
        # just remove the last field (=34mer)
        #guideData = [x[:-1] for x in guideData]

    if sortBy is not None and sortBy!="spec":
        sortFunc = (lambda row: row[1][sortBy])
    else:
        sortFunc = operator.itemgetter(0)

    guideData.sort(reverse=True, key=sortFunc)

    #guideData.sort(reverse=True, key=lambda row: 3*row[0]+row[1])
    return guideData, guideScores, hasNotFound

def printDownloadTableLinks(batchId):
    print '<div style="text-align:right">'
    print '<small>'
    print "Download tables: "
    print '<a href="crispor.cgi?batchId=%s&download=guides&format=xls">Guides</a>&nbsp;' % batchId
    print '<a href="crispor.cgi?batchId=%s&download=offtargets&format=xls">Off-targets</a>' % batchId
    print '</small>'
    print '</div>'

def hasGeneModels(org):
    " return true if this organism has gene model information "
    geneFname = join(genomesDir, org, org+".segments.bed")
    return isfile(geneFname)

def printTableHead(batchId, chrom, org):
    " print guide score table description and columns "
    # one row per guide sequence
    print '''<div class='substep'>Ranked by default from highest to lowest specificity score (<a target='_blank' href='http://dx.doi.org/10.1038/nbt.2647'>Hsu et al., Nat Biot 2013</a>) as on <a href="http://crispr.mit.org">http://crispr.mit.org</a>.'''
    print '''Also provided are efficacy scores (<a target="_blank" href="http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html">Doench et al., Nat Biot 2014</a>, <a href="http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design" target="_blank">tool</a> and <a target="_blank" href="http://genome.cshlp.org/content/early/2015/06/10/gr.191452.115">Xu et al, Gen Res 2014</a>, <a href="http://crispr.dfci.harvard.edu/SSC/" target="_blank">tool</a>), proximal GC content (<a target="_blank" href="http://www.cell.com/cell-reports/abstract/S2211-1247%2814%2900827-4">Ren et al., Cell Rep 2014</a>) and out-of-frame scores (<a target="_blank" href="http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html">Bae et al Nat Meth 2014</a>).<br>'''
    print '''Click on a column title to rank by its score.<br></div>'''

    printDownloadTableLinks(batchId)

    print """
    <script type="text/javascript">
    function onlyExons() {
        if ($("#onlyExonBox").prop("checked")) { 
            $(".otMore").show();
            $(".otMoreLink").hide();
            $(".otLessLink").hide();
            $(".notExon").hide();
            }
        else {
            if ($("#onlySameChromBox").prop("checked")) {
                $(".notExon:not(.diffChrom)").show();
            }
            else {
                $(".notExon").show();
                $(".otMoreLink").show();
                $(".otMore").hide();
            }
        }
    }
    function onlySameChrom() {
        if ($("#onlySameChromBox").prop("checked"))
            { 
            $(".otMore").show();
            $(".otMoreLink").hide();
            $(".otLessLink").hide();
            $(".diffChrom").hide();
            }
        else {
            if ($("#onlyExonBox").prop("checked")) {
                $(".diffChrom:not(.notExon)").show();
            }
            else {
                $(".diffChrom").show();
                $(".otMoreLink").show();
                $(".otMore").hide();
            }
        }
    }

    function showAllOts(classId) {
        $("#"+classId).show();
        $("#"+classId+"MoreLink").hide();
        $("#"+classId+"LessLink").show();
    }
    function showLessOts(classId) {
        $("#"+classId).hide();
        $("#"+classId+"MoreLink").show();
        $("#"+classId+"LessLink").hide();
    }
    </script>
    """

    print '<table id="otTable" style="background:white;table-layout:fixed; overflow:scroll; width:100%">'
    print '<colgroup span="1"><col></colgroup>'
    print '<colgroup span="1"><col></colgroup>'
    print '<colgroup span="1"><col></colgroup>'
    print '<colgroup span="3"><col><col><col></colgroup>'
    print '<colgroup span="1"><col></colgroup>'
    print '<colgroup span="1"><col></colgroup>'
    print '<colgroup span="1"><col></colgroup>'

    print '<thead>'
    print '<tr style="border-bottom:none; border-left:5px solid black; background-color:#F0F0F0">'
    
    print '<th style="width:80px; border-bottom:none">Position/<br>Strand'
    htmlHelp("You can click on the links in this column to highlight the <br>PAM site in the sequence viewer at the top of the page.")
    print '</th>'

    print '<th style="width:170px; border-bottom:none">Guide Sequence + <i>PAM</i><br>Restriction Enzymes'
    htmlHelp("Restriction enzymes potentially useful for screening mutations induced by the guide RNA.<br> These enzyme sites overlap cleavage site 3bp 5' to the PAM.<br>Digestion of the screening PCR product with this enzyme will not cut the product if the genome was mutated by Cas9.")

    print '<th style="width:70px; border-bottom:none"><a href="crispor.cgi?batchId=%s&sortBy=spec">Specificity Score</a>' % batchId
    htmlHelp("The higher the specificity score, the lower are off-target effects in the genome.&lt;br&gt;The specificity score ranges from 0-100 and measures the uniqueness of a guide in the genome. &lt;br&gt;See Hsu et al. Nat Biotech 2013")
    print "</th>"

    print '<th style="width:250px; border-bottom:none" colspan="%d">Predicted Efficacy' % (len(scoreNames)+1)
    htmlHelp("The higher the efficacy score, the more likely is cleavage at this position. &lt;br&gt;<br>The Chari, Doench, CrisprScan and Out-of-Frame scores range from 0-100. &lt;br&gt;The SSC score usually is in the range -2 to +2. 50-60% of inefficient guides have negative scores, see Xu et al.&lt;br&gt;Use the CrisprScan score for Zebrafish injections.<br>The last sub-column, prox. GC, indicates if the 6bp next to the PAM contain at least 4 Gs or Cs.<br>Guides like this lead to a heritable mutation rate of 60% in Drosophila (Ren et al., Cell Reports 2014)")

    #print '<th style="width:50">Prox. GC'
    #htmlHelp("At least four G or C nucleotides in the 6bp next to the PAM.<br>Ren, Zhihao, Jiang et al (Cell Reports 2014) showed that this feature is correlated with Cas9 activity (P=0.625). <br>When GC>=4, the guide RNA tested in Drosophila induced a heritable mutation rate in over 60% of cases.")
    print '</th>'

    #print '<th style="width:50px; border-bottom:none"><a href="crispor.cgi?batchId=%s&sortBy=oofScore">Out-of- Frame Score</a>' % batchId
    #htmlHelp("The Out-of-Frame Score predicts the percentage of clones that will carry out-of-frame deletions. For details see Bae et al, Nat Met 2014.")
    #print '</th>'

    print '<th style="width:120px; border-bottom:none">Off-targets for <br>0-1-2-3-4 mismatches<br><span style="color:grey">+ next to PAM </span>'
    htmlHelp("For each number of mismatches, the number of off-targets is indicated.<br>Example: 1-3-20-50-60 means 1 off-target with 0 mismatches, 3 off-targets with 1 mismatch, <br>20 off-targets with 3 mismatches, etc.<br>Off-targets are considered if they are flanked by one of the motifs NGG, NAG or NGA.<br>Shown in grey are the off-targets that have no mismatches in the 12 bp <br>adjacent to the PAM. These are the most likely off-targets.")
    #print "</th>"

    #print '<th style="width:120">Off-targets with no mismatches next to PAM</i>'
    print "</th>"
    print '<th style="width:*; border-bottom:none">Genome Browser links to matches sorted by off-target score'
    htmlHelp("For each off-target the number of mismatches is indicated and linked to a genome browser. <br>Matches are ranked by off-target score (see Hsu et al) from most to least likely.<br>Matches can be filtered to show only off-targets in exons or on the same chromosome as the input sequence.")

    print '<br><small>'
    #print '<form id="filter-form" method="get" action="crispor.cgi#otTable">'
    print '<input type="hidden" name="batchId" value="%s">' % batchId

    if hasGeneModels(org):
        print '''<input type="checkbox" id="onlyExonBox" onchange="onlyExons()">exons only'''
    else:
        print '<small style="color:grey">No exons.</small>'

    if chrom!="":
        if chrom[0].isdigit():
            chrom = "chrom "+chrom
        print '''<input type="checkbox" id="onlySameChromBox" onchange="onlySameChrom()">%s only''' % chrom
    else:
        print '<small style="color:grey">&nbsp;No match, no chrom filter</small>'
    # a hidden submit button
    # print '<input  type="submit" name="submit" value="submit">'
    #print '<input  type="submit" name="submit" value="1" style="position: absolute; height: 0px; width: 0px; border: none; padding: 0px;" hidefocus="true" tabindex="-1"/>'
    #print '</form></small>'
    print "</small>"
    print "</th>"
    print "</tr>"

    print '<tr style="border-top:none; border-left: solid black 5px; background-color:#F0F0F0">'
    print '<th style="border-top:none"></th>'
    print '<th style="border-top:none"></th>'
    print '<th style="border-top:none"></th>'

    for scoreName in scoreNames:
        if scoreNames=="oof":
            continue
        scoreLabel, scoreDesc = scoreDescs[scoreName]
        print '<th style="border: none; border-top:none; border-right: none" class="rotate"><div><span><a title="%s" href="crispor.cgi?batchId=%s&sortBy=%s">%s</a></span></div></th>' % (scoreDesc, batchId, scoreName, scoreLabel)

    # the ProxGC score comes next
    print '<th style="border: none; border-top:none; border-right: none; border-left:none" class="rotate"><div><span style="border-bottom:none"><a title="Ren et al 2014 obtained the highest cleavage when the final 6bp contained &gt;= 4 GCs"href="crispor.cgi?batchId=%s&sortBy=finalGc6">Prox GC</span></div></th>' % (batchId)

    # now the off score
    scoreName = "oof"
    scoreLabel, scoreDesc = scoreDescs[scoreName]
    print '<th style="border: none; border-top:none; border-right: none" class="rotate"><div><span><a title="%s" href="crispor.cgi?batchId=%s&sortBy=%s">%s</a></span></div></th>' % (scoreDesc, batchId, scoreName, scoreLabel)

    print '<th style="border-top:none"></th>'
    print '<th style="border-top:none"></th>'
    print "</tr>"
    print '</thead>'

def scoreToColor(guideScore):
    if guideScore > 50:
        color = "#32cd32"
    elif guideScore > 20:
        color = "#ffff00"
    else:
        color = "#aa0114"
    return color

def makeOtBrowserLinks(otData, chrom, dbInfo, pamId):
    " return a list with the html texts of the offtarget links "
    links = []

    i = 0
    for otSeq, score, editDist, pos, gene, alnHtml, gapPos in otData:
        cssClasses = ["tooltipster"]
        if not gene.startswith("exon:"):
            cssClasses.append("notExon")
        if pos.split(":")[0]!=chrom:
            cssClasses.append("diffChrom")

        classStr =  ""
        if len(cssClasses)!=0:
            classStr = ' class="%s"' % " ".join(cssClasses)

        link = makeBrowserLink(dbInfo, pos, gene, alnHtml, cssClasses)
        editDist = str(editDist)
        links.append( '''<div%(classStr)s>%(editDist)s:%(link)s</div>''' % locals() )

        #i+=1
        #if i>=3:
            #break

    #if not showAll and len(otData)>3:
         #print '''... <br>&nbsp;&nbsp;&nbsp;<a href="crispor.cgi?batchId=%s&showAll=1">- show all offtargets</a>''' % batchId
    return links

def filterOts(otDatas, minScore):
    " remove all offtargets with score < minScore "
    newList = []
    for otData in otDatas:
        score = otData[1]
        if score > minScore:
            newList.append(otData)
    return newList

def findOtCutoff(otData):
    " try cutoffs 0.5, 1.0, 2.0, 3.0 until not more than 20 offtargets left "
    for cutoff in [0.3, 0.5, 1.0, 2.0, 3.0, 10.0, 99.9]:
        otData = filterOts(otData, cutoff)
        if len(otData)<=30:
            return otData, cutoff

    if len(otData)>30:
        return otData[:30], None

    return otData, 1000

def printNoEffScoreFoundWarn(effScoresCount):
    if effScoresCount==0:
        print('<div style="text-align:left"><strong>Note:</strong> No guide could be scored for efficacy.')
        print("This happens when the input sequence is shorter than 100bp and there is no genome available to extend it. Please add flanking 50bp on both sides of the input sequence and submit this new, longer sequence.</div><br>")


def showGuideTable(guideData, pam, otMatches, dbInfo, batchId, org, showAll, chrom):
    " shows table of all PAM motif matches "
    print "<br><div class='title'>Predicted guide sequences for PAMs</div>" 

    showPamWarning(pam)
    showNoGenomeWarning(dbInfo)
    printTableHead(batchId, chrom, org)

    count = 0
    effScoresCount = 0
    for guideRow in guideData:
        guideScore, effScores, startPos, strand, pamId, guideSeq, \
            pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, subOptMatchCount = guideRow

        color = scoreToColor(guideScore)
        print '<tr id="%s" style="border-left: 5px solid %s">' % (pamId, color)

        # position and strand
        #print '<td id="%s">' % pamId
        print '<td>'
        print '<a href="#list%s">' % (pamId)
        print str(startPos+1)+" /"
        if strand=="+":
            print 'fw'
        else:
            print 'rev'
        print '</a>'
        print "</td>"

        # sequence
        print "<td>"
        print "<small>"
        print "<tt>"+guideSeq + " <i>" + pamSeq+"</i></tt>"
        print "<br>"

        scriptName = basename(__file__)
        if otData!=None and subOptMatchCount <= MAXOCC:
            print '<a href="%s?batchId=%s&pamId=%s&pam=%s" target="_blank">PCR primers</a>' % (scriptName, batchId, urllib.quote(str(pamId)), pam)

        if gcContent(guideSeq)>0.75:
            text = "This sequence has a GC content higher than 75%.<br>In the data of Tsai et al Nat Biotech 2015, the two guide sequences with a high GC content had more off-targets than all other sequences combined.<br> We do not recommend using guide sequences with such a high GC content."
            print "<br>"
            htmlWarn(text)
            print ' High GC content<br>'

        if gcContent(guideSeq)<0.25:
            text = "This sequence has a GC content lower than 25%.<br>In the data of Wang/Sabatini/Lander Science 2014, guides with a very low GC content had low cleavage efficacy."
            print "<br>"
            htmlWarn(text)
            print ' Low GC content<br>'

        print "<br>"
        if len(mutEnzymes)!=0:
            print "Restr. Enzymes:"
            print ",".join(mutEnzymes)
        print "</small>"
        print "</td>"

        # off-target score, aka specificity score aka MIT score
        print "<td>"
        if guideScore==None:
            print "No matches"
        else:
            print "%d" % guideScore
        print "</td>"

        # eff scores
        if effScores==None:
            print '<td colspan="%d">Too close to end</td>' % len(scoreNames)
            htmlHelp("The efficacy scores require some flanking sequence<br>This guide does not have enough flanking sequence in your input sequence and could not be extended as it was not found in the genome.<br>")
        else:
            for scoreName in scoreNames:
                score = effScores.get(scoreName, None)
                if score!=None:
                    effScoresCount += 1
                if score==None:
                    print '''<td>--</td>'''
                elif scoreDigits.get(scoreName, 0)==0:
                    print '''<td>%d</td>''' % int(score)
                else:
                    print '''<td>%0.2f</td>''' % (float(score))
            #print "<!-- %s -->" % seq30Mer
        # close GC > 4
        print "<td>"
        gcCount = guideSeq[-6:].count("G")+guideSeq[-6:].count("C")
        if gcCount >= 4:
            print "+"
        else:
            print "-"
        # main motif is "NGG" and last nucleotides are GGNGG
        if pam=="NGG" and patMatch(guideSeq[-2:], "GG"):
            print "<br>"
            print "<small>GG</small>"
            htmlHelp("Farboud/Meyer 2015 (Genetics 199(4)) obtained the highest cleavage in <i>C. elegans</i> with guides that end with <tt>GG</tt>. ")
        print "</td>"

        # mismatch description
        print "<td>"
        #otCount = sum([int(x) for x in otDesc.split("/")])
        if otData==None:
            # no genome match
            print otDesc
            htmlHelp("Sequence was not found in genome.<br>If you have pasted a cDNA sequence, note that sequences that overlap a splice site cannot be used as guide sequences<br>This warning also occurs if you have selected the wrong genome.")
        elif subOptMatchCount > MAXOCC:
            print ("Repeat")
            htmlHelp("At <= 4 mismatches, %d hits were found in the genome for this sequence. <br>This guide is a repeated region, it is too unspecific.<br>Usually, CRISPR cannot be used to target repeats." % subOptMatchCount)
        else:
            print otDesc
            print "<br>"

            # mismatch description, last 12 bp
            print '<small style="color:grey">'+last12Desc+"</small><br>"
            otCount = len(otData)
            print "<br><small>%d off-targets</small>" % otCount
        print "</td>"

        # links to offtargets
        print "<td><small>"
        if otData!=None:
            if len(otData)>500 and len(guideData)>1:
                otData, cutoff = findOtCutoff(otData)
                if cutoff!=None:
                    print "Too many off-targets. Showing those with score &gt;%0.1f " % cutoff
                else:
                    print "Too man off-targets. Showing "+str(len(otData))

                htmlHelp("This guide sequence has a high number of off-targets, its use is discouraged.<br>To show all off-targets, paste only the guide sequence into the input sequence box.")

            otLinks = makeOtBrowserLinks(otData, chrom, dbInfo, pamId)

            print "\n".join(otLinks[:3])
            if len(otLinks)>3:
                cssPamId = pamId.replace("-","minus").replace("+","plus") # +/-: not valid in css
                cssPamId = cssPamId+"More"
                print '<div id="%s" class="otMore" style="display:none">' % cssPamId
                print "\n".join(otLinks[3:])
                print "</div>"

                print '''<a id="%sMoreLink" class="otMoreLink" onclick="showAllOts('%s')">''' % (cssPamId, cssPamId)
                print 'show all...</a>'

                print '''<a id="%sLessLink" class="otLessLink" style="display:none" onclick="showLessOts('%s')">''' % (cssPamId, cssPamId)
                print 'show less...</a>'

        print "</small></td>"

        print "</tr>"
        count = count+1

    print "</table>"
    printDownloadTableLinks(batchId)

    printNoEffScoreFoundWarn(effScoresCount)

def linkLocalFiles(listFname):
    """ write a <link> statement for each filename in listFname. Version them via mtime
    (-> browser cache)
    """
    for fname in open(listFname).read().splitlines():
        fname = fname.strip()
        if not isfile(fname):
            fname = join(HTMLDIR, fname)
            if not isfile(fname):
                print "missing: %s<br>" % fname
                continue
        mTime = str(os.path.getmtime(fname)).split(".")[0] # seconds is enough
        if fname.endswith(".css"):
            url = fname.replace("/var/www/", "http://tefor.net/")
            print "<link rel='stylesheet' media='screen' type='text/css' href='%s?%s'/>" % (url, mTime)

def printHeader(batchId):
    " print the html header "

    print "<html><head>"

    printFile("header.inc")
    linkLocalFiles("includes.txt")
    printFile("/var/www/main/specific/googleanalytics/script.php")

    print '<link rel="stylesheet" type="text/css" href="%sstyle/tooltipster.css" />' % HTMLPREFIX
    print '<link rel="stylesheet" type="text/css" href="%sstyle/tooltipster-shadow.css" />' % HTMLPREFIX

    # the UFD combobox, https://code.google.com/p/ufd/wiki/Usage
    # patched to allow mouse wheel
    # https://code.google.com/p/ufd/issues/detail?id=86&q=mouse%20wheel
    print '<script type="text/javascript" src="%sjs/jquery.ui.ufd.js"></script>' % HTMLPREFIX
    print '<link rel="stylesheet" type="text/css" href="%sstyle/ufd-base.css" />' % HTMLPREFIX
    print '<link rel="stylesheet" type="text/css" href="%sstyle/plain.css" />' % HTMLPREFIX
    print '<link rel="stylesheet" type="text/css"  href="http://code.jquery.com/ui/1.11.1/themes/smoothness/jquery-ui.css" />'
    print '<script type="text/javascript" src="js/jquery.tooltipster.min.js"></script>'

    # activate tooltipster
   #theme: 'tooltipster-shadow',
    print ("""
    <script> 
    $(document).ready(function() { 
        $('.tooltipster').tooltipster({ 
            contentAsHTML: true,
            speed : 0
        }); });
    </script> """)

    # activate jqueryUI tooltips
    print ("""
    <script>
    $(function () {
       $(".tooltip").tooltip({
       relative : true,
       tooltipClass : "alignStyle",
       content: function () {
       return '<div style="width:300px">'+$(this).prop('title')+"</div>";
       }
      });
    });
    </script>""")


    # style of Jquery UI tooltips, default style is div.ui-tooltip
    print("""<style>
        .alignStyle {
            background-color: #FFFFFF;
            width: 350px;
            max-width: 400px;
            height: 110px;
            position : absolute;
            text-align: left;
            border:1px solid #cccccc;
        }
            </style>""")

    # style from https://css-tricks.com/rotated-table-column-headers/ to rotate table headers
    print("""<style>
       th.rotate {
         /* Something you can count on */
         /* height: 10px; */
         white-space: nowrap;
       }
       
       th.rotate > div {
         float:left;
         white-space: nowrap;
         position: relative;
         border-style: none;
         -webkit-transform: rotate(-90deg);
         -moz-transform: rotate(270deg);
         -ms-transform: rotate(270deg);
         -o-transform: rotate(270deg);

         transform: 
           /* Magic Numbers */
           /* translate(25px, 51px) */
           /* 45 is really 360 - 45 */
           rotate(270deg);
         width: 25px;
       }
       th.rotate > div > span {
         /* border-bottom: 1px solid #ccc; */
         padding: 0px 3px;
         white-space: nowrap;
       }
    </style>""")


    print("</head>")

    print'<body id="wrapper">'
    
    print "<div id='fb-root'></div>"
    print('<script src="facebooklike.js"></script>')

def firstFreeLine(lineMasks, y, start, end):
    " recursively search for first free line to place a feature (start, end) "
    #print "called with y", y
    if y>=len(lineMasks):
        return None
    lineMask = lineMasks[y]
    for x in range(start, end):
        if lineMask[x]!=0:
            return firstFreeLine(lineMasks, y+1, start, end)
        else:
            return y
    return None

def distrOnLines(seq, startDict, featLen):
    """ given a dict with start -> (start,end,name,strand) and a motif len, create lines of annotations such that
        the motifs don't overlap on the lines 
    """
    # max number of lines in y direction to draw
    MAXLINES = 18
    # amount of free space around each feature
    SLOP = 2

    # bitmask, one per line, 1 = we have a feature here, 0 = no feature here
    lineMasks = []
    for i in range(0, MAXLINES):
        lineMasks.append( [0]* (len(seq)+10) )

    # dict with lineCount (0...MAXLINES) -> list of (start, strand) tuples
    ftsByLine = defaultdict(list)
    maxY = 0
    for start in sorted(startDict):
        end = start+featLen
        strand = startDict[start]

        ftSeq = seq[start:end] 
        if strand=="+":
            label = '%s..%s'%(seq[start-3].lower(), ftSeq)
            startFt = start - 3
            endFt = end
        else:
            #print seq, end, "<br>"
            label = '%s..%s'%(ftSeq, seq[end+2].lower())
            startFt = start
            endFt = end + 3

        y = firstFreeLine(lineMasks, 0, startFt, endFt)
        if y==None:
            errAbort("not enough space to plot features")

        # fill the current mask
        mask = lineMasks[y]
        for i in range(max(startFt-SLOP, 0), min(endFt+SLOP, len(seq))):
            mask[i]=1

        maxY = max(y, maxY)

        pamId = "s%d%s" % (start, strand)
        ft = (startFt, endFt, label, strand, pamId) 
        ftsByLine[y].append(ft )
    return ftsByLine, maxY

def writePamFlank(seq, startDict, pam, faFname, gapped=False):
    " write pam flanking sequences to fasta file, optionally with versions where each nucl is removed "
    #print "writing pams to %s<br>" % faFname
    faFh = open(faFname, "w")
    for pamId, startPos, strand, flankSeq, pamSeq in flankSeqIter(seq, startDict, len(pam), gapped):
        faFh.write(">%s\n%s\n" % (pamId, flankSeq))
    faFh.close()

def runCmd(cmd):
    " run shell command, check ret code, replaces BIN and SCRIPTS special variables "
    cmd = cmd.replace("$BIN", binDir)
    cmd = cmd.replace("$SCRIPT", scriptDir)
    cmd = "set -o pipefail; " + cmd
    debug("Running %s" % cmd)
    ret = subprocess.call(cmd, shell=True, executable="/bin/bash")
    if ret!=0:
        if commandLineMode:
            logging.error("Error: could not run command %s." % cmd)
            sys.exit(1)
        else:
            print "Server error: could not run command %s.<p>" % cmd
            print "please send us an email, we will fix this error as quickly as possible. services@tefor.net "
            sys.exit(0)


def parseOfftargets(bedFname):
    """ parse a bed file with annotataed off target matches from overlapSelect,
    has two name fields, one with the pam position/strand and one with the
    overlapped segment 
    
    return as dict pamId -> editDist -> (chrom, start, end, seq, strand, segType, segName)
    segType is "ex" "int" or "ig" (=intergenic)
    if intergenic, geneNameStr is two genes, split by |
    """
    # example input:
    # chrIV 9864393 9864410 s41-|-|5|ACTTGACTG|0    chrIV   9864303 9864408 ex:K07F5.16
    # chrIV   9864393 9864410 s41-|-|5|ACTGTAGCTAGCT|9999    chrIV   9864408 9864470 in:K07F5.16
    debug("reading %s" % bedFname)

    # first sort into dict (pamId,chrom,start,end,editDist,strand) 
    # -> (segType, segName) 
    pamData = {}
    for line in open(bedFname):
        fields = line.rstrip("\n").split("\t")
        chrom, start, end, name, segment = fields
        nameFields = name.split("|")
        pamId, strand, editDist, seq = nameFields[:4]
        # strip the gap info from the pamId
        gapPos = None
        if "g" in pamId:
            pamId, gapPos = pamId.split('g')

        if len(nameFields)>4:
            x1Count = int(nameFields[4])
        else:
            x1Count = 0
        editDist = int(editDist)
        # some gene models include colons
        if ":" in segment:
            segType, segName = string.split(segment, ":", maxsplit=1)
        else:
            segType, segName = "", segment
        start, end = int(start), int(end)
        otKey = (pamId, chrom, start, end, editDist, seq, strand, x1Count, gapPos)

        # if a offtarget overlaps an intron/exon or ig/exon boundary it will
        # appear twice; in this case, we only keep the exon offtarget
        if otKey in pamData and segType!="ex":
            continue
        pamData[otKey] = (segType, segName)

    # index by pamId and edit distance
    indexedOts = defaultdict(dict)
    for otKey, otVal in pamData.iteritems():
        pamId, chrom, start, end, editDist, seq, strand, x1Score, gapPos = otKey
        segType, segName = otVal
        otTuple = (chrom, start, end, seq, strand, segType, segName, x1Score, gapPos)
        indexedOts[pamId].setdefault(editDist, []).append( otTuple )

    return indexedOts

class ConsQueue:
    """ a pseudo job queue that does nothing but report progress to the console """
    def startStep(self, batchId, desc, label):
        logging.info("Progress %s - %s - %s" % (batchId, desc, label))

def annotateBedWithPos(inBed, outBed):
    """
    given an input bed4 and an output bed filename, add an additional column 5 to the bed file
    that is a descriptive text of the chromosome pos (e.g. chr1:1.23 Mbp).
    """
    ofh = open(outBed, "w")
    for line in open(inBed):
        chrom, start = line.split("\t")[:2]
        start = int(start)

        if start>1000000:
            startStr = "%.2f Mbp" % (float(start)/1000000)
        else:
            startStr = "%.2f Kbp" % (float(start)/1000)
        desc = "%s %s" % (chrom, startStr)

        ofh.write(line.rstrip("\n"))
        ofh.write("\t")
        ofh.write(desc)
        ofh.write("\n")
    ofh.close()

def calcGuideEffScores(seq, extSeq, pam):
    """ given a sequence and an extended sequence, get all potential guides
    with pam, extend them to 100mers and score them with various eff. scores. Return a
    list of rows [headers, (guideSeq, 100mer, score1, score2, score3,...), ... ]
    """
    print seq, extSeq, pam
    seq = seq.upper()
    if extSeq:
        extSeq = extSeq.upper()
    startDict, endSet = findAllPams(seq, pam)
    pamInfo = list(flankSeqIter(seq, startDict, len(pam)))

    guideIds = []
    guides = []
    longSeqs = []
    for pamId, startPos, strand, guideSeq, pamSeq in pamInfo:
        guides.append(guideSeq+pamSeq)
        gStart, gEnd = pamStartToGuideRange(startPos, strand, len(pam))
        longSeq = getExtSeq(seq, gStart, gEnd, strand, 30, 50, extSeq)
        if longSeq!=None:
            longSeqs.append(longSeq)
            guideIds.append(pamId)

    if len(longSeqs)>0:
        effScores = crisporEffScores.calcAllScores(longSeqs)
    else:
        effScores = {}
    scoreNames = effScores.keys()

    # reformat to rows, write all scores to file
    rows = []
    for i, (guideId, guide, longSeq) in enumerate(zip(guideIds, guides, longSeqs)):
        row = [guideId, guide, longSeq]
        for scoreName in scoreNames:
            row.append(effScores[scoreName][i])
        rows.append(row)

    headerRow = ["guideId", "guide", "longSeq"]
    headerRow.extend(scoreNames)
    rows.insert(0, headerRow)
    return rows

def writeRow(ofh, row):
    " write list to file as tab-sep row "
    row = [str(x) for x in row]
    ofh.write("\t".join(row))
    ofh.write("\n")

def createBatchEffScoreTable(batchId):
    """ annotate all potential guides with efficacy scores and write to file.
    tab-sep file for easier debugging, no pickling
    """
    outFname = join(batchDir, batchId+".effScores.tab")
    seq, org, pam, position, extSeq = readBatchParams(batchId)
    seq = seq.upper()
    if extSeq:
        extSeq = extSeq.upper()

    guideRows = calcGuideEffScores(seq, extSeq, pam)
    guideFh = open(outFname, "w")
    for row in guideRows:
        writeRow(guideFh, row)
    guideFh.close()
    logging.info("Wrote eff scores to %s" % guideFh.name)

def readEffScores(batchId):
    " parse eff scores from tab sep file and return as dict pamId -> dict of scoreName -> value "
    effScoreFname = join(batchDir, batchId)+".effScores.tab"
    # old batches during transition time don't have this file yet, so make one now
    if not isfile(effScoreFname):
        createBatchEffScoreTable(batchId)

    seqToScores = {}
    for row in lineFileNext(open(effScoreFname)):
        scoreDict = {}
        rowDict = row._asdict()
        # the first three fields are the pamId, shortSeq, longSeq, they are not scores
        allScoreNames = row._fields[3:]
        for scoreName in allScoreNames:
            score = rowDict[scoreName]
            if "." in score:
                score = float(score)
            else:
                score = int(score)
            scoreDict[scoreName] = score
        seqToScores[row.guideId] = scoreDict
    return seqToScores

def processSubmission(faFnames, genome, pam, bedFname, batchBase, batchId, queue):
    """ search fasta file against genome, filter for pam matches and write to bedFName 
    optionally write status updates to work queue.
    """
    genomeDir = genomesDir # make var local, see below
    pamLen = len(pam)

    # increase MAXOCC if there is only a single query
    if len(parseFasta(open(faFnames[0])))==1:
        global MAXOCC
        MAXOCC=max(HIGH_MAXOCC, MAXOCC)

    queue.startStep(batchId, "effScores", "Calculating guide efficacy scores")
    createBatchEffScoreTable(batchId)

    if genome=="noGenome":
        # skip off-target search
        open(bedFname, "w") # create a 0-byte file to signal job completion
        queue.startStep(batchId, "done", "Job completed")
        return

    saFname = batchBase+".sa"

    matchesBedFname = batchBase+".matches.bed"
    open(matchesBedFname, "w") # truncate to 0 size
    assert(len(faFnames) <= 2)

    # 2nd fa filename is optional and if present, are sequences with gaps
    for i, faFname in enumerate(faFnames):
        # ALIGNMENT
        if i==0:
            maxDiff = maxMMs
            queue.startStep(batchId, "bwa", "Alignment of potential guides, mismatches <= %d" % maxDiff)
            convertMsg = "Converting alignments"
            seqLen = 20
        elif i==1:
            maxDiff = maxGapMMs
            queue.startStep(batchId, "bwa", "Alignment of potential gapped guides, mismatches <= %d" % maxDiff)
            convertMsg = "Converting gapped alignments"
            seqLen = 19

        bwaM = MFAC*MAXOCC # -m is queue size in bwa
        cmd = "$BIN/bwa aln -o 0 -m %(bwaM)s -n %(maxDiff)d -k %(maxDiff)d -N -l %(seqLen)d %(genomeDir)s/%(genome)s/%(genome)s.fa %(faFname)s > %(saFname)s" % locals()
        runCmd(cmd)

        queue.startStep(batchId, "saiToBed", convertMsg)
        maxOcc = MAXOCC # create local var from global
        # EXTRACTION OF POSITIONS + CONVERSION + SORT/CLIP
        # the sorting should improve the twoBitToFa runtime
        cmd = "$BIN/bwa samse -n %(maxOcc)d %(genomeDir)s/%(genome)s/%(genome)s.fa %(saFname)s %(faFname)s | $SCRIPT/xa2multi.pl | $SCRIPT/samToBed %(pamLen)s | sort -k1,1 -k2,2n | $BIN/bedClip stdin %(genomeDir)s/%(genome)s/%(genome)s.sizes stdout >> %(matchesBedFname)s " % locals()
        runCmd(cmd)

    # arguments: guideSeq, mainPat, altPats, altScore, passX1Score
    filtMatchesBedFname = batchBase+".filtMatches.bed"
    queue.startStep(batchId, "filter", "Removing matches without a PAM motif")
    altPats = offtargetPams.get(pam, "na")
    bedFnameTmp = bedFname+".tmp"
    altPamMinScore = str(ALTPAMMINSCORE)
    # EXTRACTION OF SEQUENCES + ANNOTATION
    faFnameStr = ",".join(faFnames)
    cmd = "$BIN/twoBitToFa %(genomeDir)s/%(genome)s/%(genome)s.2bit stdout -bed=%(matchesBedFname)s | $SCRIPT/filterFaToBed %(faFnameStr)s %(pam)s %(altPats)s %(altPamMinScore)s %(maxOcc)d > %(filtMatchesBedFname)s" % locals()
    runCmd(cmd)

    segFname = "%(genomeDir)s/%(genome)s/%(genome)s.segments.bed" % locals()

    # if we have gene model segments, annotate them, otherwise just use the chrom position
    if isfile(segFname):
        queue.startStep(batchId, "genes", "Annotating matches with genes")
        cmd = "cat %(filtMatchesBedFname)s | $BIN/overlapSelect %(segFname)s stdin stdout -mergeOutput -selectFmt=bed -inFmt=bed | cut -f1,2,3,4,8 2> %(batchBase)s.log > %(bedFnameTmp)s " % locals()
        runCmd(cmd)
    else:
        queue.startStep(batchId, "chromPos", "Annotating matches with chromosome position")
        annotateBedWithPos(filtMatchesBedFname, bedFnameTmp)

    # make sure the final bed file is never in a half-written state, it is our signal that the job is complete
    shutil.move(bedFnameTmp, bedFname)
    queue.startStep(batchId, "done", "Job completed")


def lineFileNext(fh):
    """ 
        parses tab-sep file with headers as field names 
        yields collection.namedtuples
        strips "#"-prefix from header line
    """
    line1 = fh.readline()
    line1 = line1.strip("\n").strip("#")
    headers = line1.split("\t")
    Record = namedtuple('tsvRec', headers)
   
    for line in fh:
        line = line.rstrip("\n")
        fields = line.split("\t")
        try:
            rec = Record(*fields)
        except Exception, msg:
            logging.error("Exception occured while parsing line, %s" % msg)
            logging.error("Filename %s" % fh.name)
            logging.error("Line was: %s" % repr(line))
            logging.error("Does number of fields match headers?")
            logging.error("Headers are: %s" % headers)
            #raise Exception("wrong field count in line %s" % line)
            continue
        # convert fields to correct data type
        yield rec

def readGenomes():
    " return list of genomes supported "
    genomes = {}

    myDir = dirname(__file__)
    genomesDir = join(myDir, "genomes")
    for subDir in os.listdir(genomesDir):
        infoFname = join(genomesDir, subDir, "genomeInfo.tab")
        if isfile(infoFname):
            row = lineFileNext(open(infoFname)).next()
            # add a note to identify UCSC genomes
            if row.server.startswith("ucsc"):
                addStr="UCSC "
            else:
                addStr = ""
            genomes[row.name] = row.scientificName+" - "+row.genome+" - "+addStr+row.description
            #genomes[row.name] = row.genome+" - "+row.scientificName+" - "+row.description

    genomes = genomes.items()
    genomes.sort(key=operator.itemgetter(1))
    return genomes

def printOrgDropDown(lastorg):
    " print the organism drop down box "
    genomes = readGenomes()
    print '<select id="genomeDropDown" class style="float:left; max-width:350px" name="org" tabindex="2">'
    print '<option '
    if lastorg == "noGenome":
        print 'selected '
    print 'value="noGenome">No Genome, no off-target ranking, show only cleavage efficacy scores</option>'

    for db, desc in genomes:
        print '<option '
        if db == lastorg :
            print 'selected '
        print 'value="%s">%s</option>' % (db, desc)

    print "</select>"
    #print ('''
      #<script type="text/javascript">
      #$("#genomeDropDown").ufd({maxWidth:350, listWidthFixed:false});
      #</script>''')
    print ('''<br>''')

def printPamDropDown(lastpam):
    
    print '<select style="float:left" name="pam" tabindex="3">'
    for key,value in pamDesc:
        print '<option '
        if key == lastpam :
            print 'selected '
        print 'value="%s">%s</option>' % (key, value)
    print "</select>"

def printForm(params):
    " print html input form "
    scriptName = basename(__file__)

    # The returned cookie is available in the os.environ dictionary
    cookies=Cookie.SimpleCookie(os.environ.get('HTTP_COOKIE'))
    if "lastorg" in cookies and "lastseq" in cookies and "lastpam" in cookies:
       lastorg   = cookies['lastorg'].value
       lastseq   = cookies['lastseq'].value
       lastpam   = cookies['lastpam'].value
    else:
       lastorg = DEFAULTORG
       lastseq = DEFAULTSEQ
       lastpam = DEFAULTPAM

    print """
<form id="main-form" method="post" action="%s">

<div class="introtext">
 CRISPOR v2.0 (CRISPr selectOR, http://crispor.tefor.net) is a program that helps design and evaluate target sites for use with the CRISPR/Cas9 system.
    <div onclick="$('.about-us').toggle('fast');" class="title" style="cursor:pointer;display:inline;font-size:large;font-style:normal">
        <img src="%simage/info.png" class="infopoint" style="vertical-align:text-top;">
    </div>
    <div class="about-us"><br>
    CRISPOR v2.0 uses the BWA algorithm to identify guide RNA sequences for CRISPR mediated genome editing.<br>
    It searches for off-target sites (with and without mismatches), shows them in a table and annotates them with flanking genes.<br>
    For more information on principles of CRISPR-mediated genome editing, check the <a href="https://www.addgene.org/CRISPR/guide/">Addgene CRISPR guide</a>.</div>
</div>

<div class="windowstep subpanel" style="width:50%%;">
    <div class="substep">
        <div class="title" style="cursor:pointer;" onclick="$('#helptext1').toggle('fast')">
            Step 1
            <img src="%simage/info.png" class="infopoint" >
        </div>
       Submit a single sequence for guide RNA identification and analysis
    </div>

    <textarea tabindex="1" style="width:100%%" name="seq" rows="10"
              placeholder="Enter the sequence of the gene you want to target - example: %s">
    %s
    </textarea>
    <div style="text-align:left"><small>Text case is preserved, e.g. you can mark ATGs with lowercase letters.<br>Instead of a sequence, you can also paste a chromosome range, e.g. "chr4:111539573-111539808"</small></div>
    <div id="helptext1" class="helptext">CRISPOR conserves the lowercase and uppercase format of your sequence (allowing to highlight sequence features of interest such as ATG or STOP codons)</div>
    
</div>
<div class="windowstep subpanel" style="width:40%%">
    <div class="substep">
        <div class="title" style="cursor:pointer;" onclick="$('#helpstep2').toggle('fast')">
            Step 2
            <img src="image/info.png" class="infopoint">
        </div>
        Choose a species genome

    </div>
    """% (scriptName, HTMLPREFIX, HTMLPREFIX, lastseq,lastseq)

    printOrgDropDown(lastorg)
    #print '<small style="float:left">Type a species name to search for it</small>'

    print """<div id="helpstep2" class="helptext">More information on these species can be found on the <a href="http://www.efor.fr">EFOR</a> website.
To add your genome of interest to the list, send us 
<a href="mailto:services@tefor.net">an email</a>.</div>
"""
    print """
</div>
<div class="windowstep subpanel" style="width:40%%">
    <div class="substep">
        <div class="title" style="cursor:pointer;" onclick="$('#helpstep3').toggle('fast')">
            Step 3
            <img src="%simage/info.png" class="infopoint">
        </div>
        Choose a Protospacer Adjacent Motif (PAM)
    </div>
    """ % HTMLPREFIX
    printPamDropDown(lastpam)
    print """
    <div id="helpstep3" class="helptext">The most common system uses the NGG PAM recognized by Cas9 from S. <i>pyogenes</i>. The VRER and VQR mutants were described by <a href="http://www.nature.com/nature/journal/vaop/ncurrent/abs/nature14592.html" target="_blank">Kleinstiver et al</a>, Nature 2015</div>
</div>


    <input type="submit" name="submit" value="SUBMIT" tabindex="4"/>
    """    
    printFile("sponsors.inc")
    print """



<style>
      
   div.sponsors:hover
   {
        -webkit-animation-play-state: paused;
        opacity:1;
   }
   div.sponsors
   {
       -webkit-animation:           mymove 2s infinite; /* Chrome, Safari, Opera */
       -webkit-animation-direction: alternate; /* Chrome, Safari, Opera */

       animation:           mymove 2s infinite;        
       animation-direction: alternate;
   }
   /* Chrome, Safari, Opera */
   @-webkit-keyframes mymove
   {
       from
       {                
           opacity: 0.7;
       }
       to   
       {             
           opacity: 1;
       }
   }
   /* Standard syntax */
   @keyframes mymove
   {
       from
       {                
           opacity: 0.7;
       }
       to   
       {             
           opacity: 1;
       }
   }

</style>
</form>
    """

def readBatchParams(batchId):
    """ given a batchId, return the genome, the pam, the input sequence and the
    chrom pos and extSeq, a 100bp-extended version of the input sequence.
    Returns None for pos if not found. """

    batchBase = join(batchDir, batchId)
    inputFaFname = batchBase+".input.fa"
    if not isfile(inputFaFname):
        errAbort('Could not find the batch %s. We cannot keep Crispor runs for more than '
                'a few months. Please resubmit your input sequence via'
            ' <a href="crispor.cgi">the query input form</a>' % batchId)

    ifh = open(inputFaFname)
    ifhFields = ifh.readline().replace(">","").strip().split()
    if len(ifhFields)==2:
        genome, pamSeq = ifhFields
        position = None
    else:
        genome, pamSeq, position = ifhFields

    inSeq = ifh.readline().strip()

    ifh.seek(0)
    seqs = parseFasta(ifh)
    ifh.close()

    extSeq = None
    if "extSeq" in seqs:
        extSeq = seqs["extSeq"]

    # some older batch files don't include a position
    if position==None:
        position = coordsToPosStr(*findBestMatch(genome, inSeq))

    return inSeq, genome, pamSeq, position, extSeq

def findAllPams(seq, pam):
    """ find all matches for PAM and return as dict startPos -> strand and a set
    of end positions
    """
    startDict, endSet = findPams(seq, pam, "+", {}, set())
    startDict, endSet = findPams(seq, revComp(pam), "-", startDict, endSet)
    return startDict, endSet

def newBatch(seq, org, pam, skipAlign=False):
    """ obtain a batch ID and write seq/org/pam to their files. 
    Return batchId, position string and a 100bp-extended seq, if possible.
    """
    batchId = makeTempBase(seq, org, pam)
    if skipAlign:
        chrom, start, end, strand = None, None, None, None
    else:
        chrom, start, end, strand = findBestMatch(org, seq)
    # define temp file names
    batchBase = join(batchDir, batchId)
    # save input seq, pamSeq, genome, position for primer design later
    inputFaFname = batchBase+".input.fa"
    posStr = coordsToPosStr(chrom, start, end, strand)
    ofh = open(inputFaFname, "w")
    ofh.write(">%s %s %s\n%s\n" % (org, pam, posStr, seq))

    # try to get a 100bp-extended version of the input seq
    extSeq = None
    if chrom!=None:
        extSeq = extendAndGetSeq(org, chrom, start, end, strand)
        if extSeq!=None:
            ofh.write(">extSeq\n%s\n" % (extSeq))
    ofh.close()
    return batchId, posStr, extSeq

def readDbInfo(org):
    " return a dbInfo object with the columsn in the genomeInfo.tab file "
    myDir = dirname(__file__)
    genomesDir = join(myDir, "genomes")
    infoFname = join(genomesDir, org, "genomeInfo.tab")
    if not isfile(infoFname):
        return None
    dbInfo = lineFileNext(open(infoFname)).next()
    return dbInfo

def printQueryNotFoundNote(dbInfo):
    print "<div class='title'>Query sequence, not present in the genome of %s</div>" % dbInfo.scientificName
    print "<div class='substep'>"
    print "<em><strong>Note:</strong> The query sequence was not found in the selected genome."
    print "This can be a valid query, e.g. a GFP sequence.<br>"
    print "If not, you might want to check if you selected the right genome for your query sequence.<br>"
    print "When reading the list of guide sequences and off-targets below, bear in mind that the software cannot distinguish off-targets from on-targets now, so some 0-mismatch targets are expected. In this case, the scores of guide sequences are too low.<p>"
    print "</em></div>"

def getOfftargets(seq, org, pam, batchId, startDict, queue):
    """ write guides to fasta and run bwa or use older cached results.
    Return name of the BED file with the matches.
    Write progress status updates to queue object.
    """
    batchBase = join(batchDir, batchId)
    otBedFname = batchBase+".bed"
    flagFile = batchBase+".running"

    if isfile(flagFile):
       errAbort("This sequence is still being processed. Please wait for ~20 seconds "
           "and try again, e.g. by reloading this page. If you see this message for "
           "more than 1-2 minutes, please email services@tefor.net")

    if not isfile(otBedFname) or commandLineMode:
        # write potential PAM sites to file 
        faFnames = [batchBase+".fa"]
        writePamFlank(seq, startDict, pam, faFnames[0])
        # if needed: add a 2nd fasta file with gapped sequences
        if allowGap:
            faFnames.append(batchBase+".gap.fa")
            writePamFlank(seq, startDict, pam, faFnames[1], gapped=True)

        if commandLineMode:
            processSubmission(faFnames, org, pam, otBedFname, batchBase, batchId, queue)
        else:
            q = JobQueue(JOBQUEUEDB)
            ip = os.environ["REMOTE_ADDR"]
            wasOk = q.addJob("search", batchId, "ip=%s,org=%s,pam=%s" % (ip, org, pam))
            if not wasOk:
                print "Job is running"
            return None

    return otBedFname

def startAjaxWait(batchId):
    """ print the ajax script to stdout """
    scriptName = basename(__file__)
    #$.ajax({ url: "server", success: function(data){
    #//$.getJSON( "%(scriptName)s?batchId=%(batchId)s&ajaxStatus=1", gotStatus);
    print """
    Status: <div id="statusEl">Waiting</div>
    <script>
    function gotStatus( data ) {
      status = data["status"];
      if (status==null || status=="null")
          { window.location.href=window.location.href; }
      else
          { $("#statusEl").text(status); }
      };
     
    (function poll(){
       setTimeout(function(){
            $.getJSON( "%(scriptName)s?batchId=%(batchId)s&ajaxStatus=1", gotStatus);
            poll();
          } , 2000)})();

    </script>
    """ % locals()

def showPamWarning(pam):
    if pam!="NGG":
        print '<div style="text-align:left">'
        print "<strong>Note:</strong> Your query involves a Cas9 that is not from S. Pyogenes. "
        print "Please bear in mind that specificity end efficacy scores were designed using data with S. Pyogenes Cas9 and might not be applicable to this particular Cas9.<br>"
        print '</div>'

def showNoGenomeWarning(dbInfo):
    if dbInfo==None:
        print('<div style="text-align:left;"><strong>Note:</strong> As there is no genome that can be used to get flanking sequence for your sequence, efficacy scores close to the start and the end of your sequence cannot be calculated and will be shown as "--".</div>')

def getSeq(db, posStr):
    """
    given a database name and a string with the position as chrom:start-end, return the sequence as
    a string.
    """
    chrom, start, end, strand =  parsePos(posStr)
    if end-start > MAXSEQLEN:
        errAbort("Input sequence range too long. Please retry with a sequence range shorter than %d bp." % MAXSEQLEN)
    genomeDir = genomesDir # pull in global var
    twoBitFname = "%(genomeDir)s/%(db)s/%(db)s.2bit" % locals()
    cmd = ["twoBitToFa", twoBitFname, "-seq="+chrom, "-start="+str(start), "-end="+str(end), "stdout"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    seqStr = proc.stdout.read()
    # remove fasta header line
    lines = seqStr.splitlines()
    lines.pop(0)
    seq = "".join(lines)
    return seq

def crisprSearch(params):
    " do crispr off target search and eff. scoring "
    # retrieve sequence if not provided
    if "pos" in params and not "seq" in params:
        params["seq"] = getSeq(params["org"], params["pos"])

    if "batchId" in params:
        # if we're getting only the batchId, extract the parameters from the batch
        # this allows a stable link to a batch that is done
        batchId = params["batchId"]
        seq, org, pam, position, extSeq = readBatchParams(batchId)
        seq, warnMsg = cleanSeq(seq)
    else:
        seq, org, pam = params["seq"], params["org"], params["pam"]
        # the "seq" parameter can contain a chrom:start-end position instead of the sequence.
        if re.match(" *[a-zA-Z0-9_-]+: *[0-9, ]+ *- *[0-9,]+ *", seq):
            seq = getSeq(params["org"], seq.replace(" ","").replace(",",""))

        seq, warnMsg = cleanSeq(seq)
        batchId, position, extSeq = newBatch(seq, org, pam)
        print ("<script>")
        print ('''history.replaceState('crispor.cgi', document.title, '?batchId=%s');''' % (batchId))
        print ("</script>")

    if len(warnMsg)!=0:
        print warnMsg+"<p>"

    batchBase = join(batchDir, batchId)

    # read genome info tab file into memory
    dbInfo = readDbInfo(org)

    # search PAMs
    uppSeq = seq.upper()
    startDict, endSet = findAllPams(uppSeq, pam)
    otBedFname = getOfftargets(uppSeq, org, pam, batchId, startDict, None)

    if otBedFname is None:
        # this can happen only in CGI mode. Job has been added to the queue or is not done yet. 
        startAjaxWait(batchId)
        return

    # be more sensitive if only a single guide seq is run
    if len(startDict)==1:
        global MAXOCC
        MAXOCC=max(HIGH_MAXOCC, MAXOCC)

    if dbInfo==None:
        print "<div class='title'>No Genome selected, specificity scoring is deactivated</div>"
        print('<div style="text-align:left;"><strong>Note:</strong> There are no predicted off-targets below and all specificity scores will be shown in red as their score is 0. <br></div>')

    elif position=='?':
        printQueryNotFoundNote(dbInfo)
    else:
        genomePosStr = ":".join(position.split(":")[:2])
        print "<div class='title'><em>%s</em> sequence found at " % (dbInfo.scientificName)
        print '<span style="text-decoration:underline">'
        mouseOver = "link to UCSC or Ensembl Genome Browser"
        if dbInfo.server=="manual":
            mouseOver = "no genome browser link available for this organism"
        print makeBrowserLink(dbInfo, genomePosStr, genomePosStr, mouseOver, ["tooltipster"])
        print "</span></div>"
        #print " (link to Genome Browser)</div>"

    otMatches = parseOfftargets(otBedFname)
    effScores = readEffScores(batchId)
    sortBy = (params.get("sortBy", None))
    guideData, guideScores, hasNotFound = mergeGuideInfo(uppSeq, startDict, pam, otMatches, position, effScores, sortBy)

    if hasNotFound and not position=="?":
        print('<div style="text-align:left"><strong>Note:</strong> At least one of the possible guide sequences was not found in the genome. ')
        print("If you pasted a cDNA sequence, note that sequences with score 0, e.g. splice junctions, are not in the genome, only in the cDNA and are not usable as CRISPR guides.</div><br>")

    showSeqAndPams(seq, startDict, pam, guideScores)

    showAll = (params.get("showAll", 0)=="1")

    chrom, start, end, strand = parsePos(position)

    showGuideTable(guideData, pam, otMatches, dbInfo, batchId, org, showAll, chrom)

    print '<br><a class="neutral" href="crispor.cgi">'
    print '<div class="button" style="margin-left:auto;margin-right:auto;width:150;">New Query</div></a>'

def runPhp(script):
    " run a file through php, write result to stdout. accepts a full or a relative path "
    if "/" in script:
        path = script
    else:
        myDir = dirname(__file__)
        path = "%s/%s" % (myDir, script)

    if not isfile(path):
        return
    proc = subprocess.Popen("php "+path, shell=True, stdout=subprocess.PIPE)
    script_response = proc.stdout.read()
    print script_response

def printFile(fname):
    if "/" in fname:
        path = fname
    else:
        myDir = dirname(__file__)
        path = "%s/%s" % (myDir, fname)
#
    if not isfile(path):
        print "install error: %s not found" % path
        return
    print open(path).read()

def printTeforBodyStart():
    print "<div class='logo'><a href='http://tefor.net/main/'><img src='http://tefor.net/main/images/logo/logo_tefor.png' alt='logo tefor infrastructure'></a></div>"
    printFile("menu.inc")

    print '<div id="bd">'
    print '<div class="centralpanel" style="margin-left:0px">'
    #runPhp("networking.php")
    print '<div class="subpanel" style="background:transparent;box-shadow:none;">'
    print '<div class="contentcentral" style="margin-left:0px; width:100%">'

def printTeforBodyEnd():
    print '</div>'
    print '</div>'
    print '</div>'
    runPhp("footer.php")
    print '</div>'

def iterGuideRows(guideData, addHeaders=False):
    "yield rows from guide data "
    headers = ["guideId", "guideSeq", "specScore", "offtargetCount", "guideGenomeMatchGeneLocus"]
    for scoreName in scoreNames:
        headers.append(scoreName+"EffScore")
    if addHeaders:
        yield guideHeaders

    #print "\t".join(headers)

    pamIdRe = re.compile(r's([0-9]+)([+-])g?([0-9]*)')
    for guideRow in guideData:
        guideScore, effScores, startPos, strand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, subOptMatchCount = guideRow

        otCount = 0
        if otData!=None:
            otCount = len(otData)

        # s20+ or s20+g0 if gapped
        #pamPos = int(pamId[1:-1])+1
        #strand = pamId[-1]
        pamPos, strand, rest = pamIdRe.match(pamId).groups()
        if strand=="+":
            strDesc = 'fw'
        else:
            strDesc = 'rev'
        guideDesc = str(pamPos)+strDesc

        row = [guideDesc, guideSeq+pamSeq, guideScore, otCount, ontargetDesc]
        for scoreName in scoreNames:
            row.append(effScores.get(scoreName, "NotEnoughFlankSeq"))
        row = [str(x) for x in row]
        yield row

def iterOfftargetRows(guideData, addHeaders=False):
    " yield bulk offtarget rows for the tab-sep download file "
    if addHeaders:
        headers = list(offtargetHeaders) # clone list
        if allowGap:
            headers.append("gapPos")
        yield headers

    for guideRow in guideData:
        guideScore, effScores, startPos, strand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, subOptMatchCount = guideRow

        otCount = 0

        if otData!=None:
            otCount = len(otData)
            for otSeq, score, editDist, pos, gene, alnHtml, gapPos in otData:
                gene = gene.replace(",", "_").replace(";","-")

                chrom, start, end, strand = parsePos(pos)
                row = [pamId, guideSeq+pamSeq, otSeq, editDist, score, chrom, start, end, gene]
                if allowGap:
                    row.append(gapPos)
                row = [str(x) for x in row]
                yield row

def xlsWrite(rows, title, outFile, colWidths, fileFormat):
    """ given rows, writes a XLS binary stream to outFile, if xlwt is available 
    Otherwise writes a tab-sep file.
    colWidths is a list of widths of columns, in Arial characters.
    """
    if xlwtLoaded and not fileFormat=="tsv":
        charSize = 269 # see http://reliablybroken.com/b/2011/10/widths-heights-with-xlwt-python/
        wb = xlwt.Workbook()
        ws = wb.add_sheet(title)

        for rowCount, row in enumerate(rows):
            for colCount, col in enumerate(row):
                if col.isdigit():
                    col = int(col)
                ws.write(rowCount, colCount, col)

        # set sizes in characters per column
        for colId, colWidth in enumerate(colWidths):
            ws.col(colId).width = charSize*colWidth

        wb.save(outFile)
    else:
        for row in rows:
            outFile.write("\t".join(row))
            outFile.write("\n")
    outFile.flush()

def downloadFile(params):
    " "
    batchId = params["batchId"]
    seq, org, pam, position, extSeq = readBatchParams(batchId)
    uppSeq = seq.upper()

    startDict, endSet = findAllPams(uppSeq, pam)

    otBedFname = join(batchDir, batchId+".bed")
    effScoreFname = join(batchDir, batchId+".effScores.tab")

    otMatches = parseOfftargets(otBedFname)
    effScores = readEffScores(batchId)
    guideData, guideScores, hasNotFound = mergeGuideInfo(uppSeq, startDict, pam, otMatches, position, effScores)

    if position=="?":
        queryDesc = org+"_unknownLoc"
    else:
        queryDesc = org+"_"+position.strip(":+-").replace(":","_")
        print org, position, queryDesc

    fileFormat = params.get("format", "tsv")
    if params["download"]=="guides":
        print "Content-Disposition: attachment; filename=\"guides_%s.%s\"" % (queryDesc, fileFormat)
        print "" # = end of http headers
        xlsWrite(iterGuideRows(guideData, addHeaders=True), "guides", sys.stdout, [6,28,10,10], fileFormat)

    elif params["download"]=="offtargets":
        print "Content-Disposition: attachment; filename=\"offtargets-%s.%s\"" % (queryDesc, fileFormat)
        print "" # = end of http headers
        otRows = list(iterOfftargetRows(guideData, addHeaders=True))
        otRows.sort(key=operator.itemgetter(4), reverse=True)
        xlsWrite(otRows, "offtargets", sys.stdout, [6,28,28,5], fileFormat)

def printBody(params):
    " main dispatcher function "
    if len(params)==0:
        printForm(params)
    elif "batchId" in params and "pamId" in params and "pam" in params:
        primerDetailsPage(params)
    elif (("seq" in params or "pos" in params) and "org" in params and "pam" in params) or \
         "batchId" in params:
        crisprSearch(params)
    else:
        errAbort("Unrecognized CGI parameters.")

def parseBoulder(tmpOutFname):
    " parse a boulder IO style file, as output by Primer3 "
    data = {}
    for line in open(tmpOutFname):
        key, val = line.rstrip("\n").split("=")
        data[key] = val
    return data

def runPrimer3(seq, tmpInFname, tmpOutFname, targetStart, targetLen, prodSizeRange):
        """ return primers from primer3 in format seq1, tm1, pos1, seq2, tm2, pos2"""
        conf = """SEQUENCE_TEMPLATE=%(seq)s
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_PRODUCT_SIZE_RANGE=%(prodSizeRange)s
SEQUENCE_TARGET=%(targetStart)s,%(targetLen)s
=""" % locals()
        open(tmpInFname, "w").write(conf)

        cmdLine = "primer3_core %s > %s" % (tmpInFname, tmpOutFname)
        runCmd(cmdLine)

        p3 = parseBoulder(tmpOutFname)
        seq1 = p3["PRIMER_LEFT_0_SEQUENCE"]
        seq2 = p3["PRIMER_RIGHT_0_SEQUENCE"]
        tm1 = p3["PRIMER_LEFT_0_TM"]
        tm2 = p3["PRIMER_RIGHT_0_TM"]
        pos1 = int(p3["PRIMER_LEFT_0"].split(",")[0])
        pos2 = int(p3["PRIMER_RIGHT_0"].split(",")[0])
        return seq1, tm1, pos1, seq2, tm2, pos2

def parseFastaAsList(fileObj):
    " parse a fasta file, return list (id, seq) "
    seqs = []
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if seqId!=None:
                seqs.append( (seqId, "".join(parts)) )
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts)!=0:
        seqs.append( (seqId, "".join(parts)) )
    return seqs

def parseFasta(fileObj):
    " parse a fasta file, where each seq is on a single line, return dict id -> seq "
    seqs = {}
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if seqId!=None:
                seqs[seqId]  = "".join(parts)
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts)!=0:
        seqs[seqId]  = "".join(parts)
    return seqs

def coordsToPosStr(chrom, start, end, strand):
    " convert coords to a string "
    if chrom==None:
        return "?"
    locStr = "%s:%d-%d:%s" % (chrom, start, end, strand)
    return locStr

def findBestMatch(genome, seq):
    """ find best match for input sequence from batchId in genome and return as
    a string "chrom:start-end:strand or None if not found "
    """
    if genome=="noGenome":
        return None, None, None, None

    # write seq to tmp file
    tmpFaFh = tempfile.NamedTemporaryFile(prefix="crisporBestMatch", suffix=".fa")
    tmpFaFh.write(">tmp\n%s" % seq)
    tmpFaFh.flush()
    logging.debug("seq: %s" % open(tmpFaFh.name).read())
    faFname = tmpFaFh.name

    # create temp SAM file
    tmpSamFh = tempfile.NamedTemporaryFile(prefix="crisporBestMatch", suffix=".sam")
    samFname = tmpSamFh.name

    genomeDir = genomesDir # make local var
    cmd = "$BIN/bwa bwasw -T 20 %(genomeDir)s/%(genome)s/%(genome)s.fa %(faFname)s > %(samFname)s" % locals()
    runCmd(cmd)

    chrom, start, end = None, None, None
    for l in open(samFname):
        if l.startswith("@"):
            continue
        logging.debug("%s" % l)
        l = l.rstrip("\n")
        fs = l.split("\t")
        qName, flag, rName, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = fs[:11]
        if (int(flag) and 2) == 2:
            strand = "-"
        else:
            strand = "+"
        if not re.compile("[0-9]*").match(cigar):
            continue
        if cigar=="*":
            logging.debug("No best match found")
            return None, None, None, None
            #errAbort("Sequence not found in genome. Are you sure you have pasted the correct sequence and also selected the right genome?")
        # XX why do we get soft clipped sequences from BWA? repeats?
        cleanCigar = cigar.replace("M","").replace("S", "")
        if not cleanCigar.isdigit():
            logging.debug("Best match found, but cigar string was %s" % cigar)
            return None, None, None, None
        matchLen = int(cleanCigar)
        chrom, start, end =  rName, int(pos)-1, int(pos)-1+matchLen # SAM is 1-based
        #print chrom, start, end, strand

    # delete the temp files
    tmpSamFh.close()
    tmpFaFh.close()
    logging.debug("Found best match at %s:%d-%d:%s" % (chrom, start, end, strand))
    return chrom, start, end, strand

def designPrimer(genome, chrom, start, end, strand, guideStart, batchId):
    " create primer for region around chrom:start-end, write output to batch "
    " returns (leftPrimerSeq, lTm, lPos, rightPrimerSeq, rTm, rPos, amplified sequence)" 
    flankStart = start - 1000
    flankEnd = end + 1000

    if flankStart<0:
        errAbort("Not enough space on genome sequence to design primer. Please design it manually")

    flankFname = join(batchDir, batchId+".inFlank.fa")
    cmd = "twoBitToFa genomes/%(genome)s/%(genome)s.2bit %(flankFname)s -seq=%(chrom)s -start=%(flankStart)d -end=%(flankEnd)d" % locals()
    runCmd(cmd)

    flankSeq = parseFasta(open(flankFname)).values()[0]
    tmpFname = join(batchDir, batchId+".primer3.in")
    tmpOutFname = join(batchDir, batchId+".primer3.out")
    # the guide seq has to be at least 150bp away from the left PCR primer for agarose gels
    lSeq, lTm, lPos, rSeq, rTm, rPos = \
        runPrimer3(flankSeq, tmpFname, tmpOutFname, 1000+guideStart-150, 330, "300-600")
    targetSeq = flankSeq[lPos:rPos+1]
    return lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq

def markupSeq(seq, start, end):
    " print seq with start-end in bold "
    return seq[:start]+"<u>"+seq[start:end]+"</u>"+seq[end:]

def makeHelperPrimers(guideName, guideSeq):
    " return dict with various names -> primer for primer page "
    primers = defaultdict(list)
    guideRnaFw = "guideRna%sT7sense" % guideName
    guideRnaRv = "guideRna%sT7antisense" % guideName

    # T7 plasmids
    if guideSeq.lower().startswith("gg"):
        primers["T7"].append((guideRnaFw, "TA<b>%s</b>" % guideSeq))
        primers["T7"].append((guideRnaRv, "AAAC<b>%s</b>" % revComp(guideSeq[2:])))
    else:
        primers["T7"].append((guideRnaFw, "TAGG<b>%s</b>" % guideSeq))
        primers["T7"].append((guideRnaRv, "AAAC<b>%s</b>" % revComp(guideSeq)))

    # T7 in-vitro
    prefix = ""
    if not guideSeq.lower().startswith("gg"):
        prefix = "GG"
    specPrimer = "TAATACGACTCACTATA%s<b>%s</b>GTTTTAGAGCTAGAAATAGCAAG" % (prefix, guideSeq)

    primers["T7iv"].append(("guideRNA%sT7PromSense" % guideName, specPrimer))
    primers["T7iv"].append(("guideRNAallT7PromAntisense (constant primer used for all guide RNAs)", "AAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC"))

    # mammalian cells
    fwName = "guideRNA%sU6sense" % guideName
    revName = "guideRNA%sU6antisense" % guideName
    if guideSeq.lower().startswith("g"):
        primers["mammCells"].append((fwName, "ACACC<b>%s</b>G" % guideSeq))
        primers["mammCells"].append((revName, "AAAAC<b>%s</b>G" % revComp(guideSeq)))
    else:
        primers["mammCells"].append((fwName, "ACACC<u>G</u><b>%s</b>G" % guideSeq))
        primers["mammCells"].append((revName, "AAAAC<b>%s</b><u>C</u>G" % revComp(guideSeq)))
        primers["mammCellsNote"] = True

    return primers

def printPrimerTableAll(primers):
    print '<table class="primerTable">'
    for key, primerList in primers.iteritems():
        if key.endswith("Note"):
            continue
        for name, seq in primerList:
            name = name.split()[0]
            print '<tr>'
            print "<td>%s</td>" % name
            print "<td><tt>%s</tt></td>" % seq
            print "</tr>"
    print "</table>"

def printPrimerTable(primerList, onRows=False):
    " given a list of (name, seq) tuples, print a table "
    print '<table class="primerTable">'
    for name, seq in primerList:
        if onRows:
            print "<tr><td>%s</td></tr>" % name
            print "<tr><td><tt>%s</tt></td></tr>" % seq
        else:
            print '<tr>'
            print "<td>%s</td>" % name
            print "<td><tt>%s</tt></td>" % seq
            print "</tr>"
    print "</table>"

def printDropboxLink():
    print """
    <a id="dropbox-link"></a>
    <script type="text/javascript" src="https://www.dropbox.com/static/api/2/dropins.js" id="dropboxjs" data-app-key="lp7njij87kjrcxv"></script>");
    <script type="text/javascript">
    options = {
        success: function(files) {
           $('textarea[name=hgct_customText]').val(files[0].link);
           alert('Here is the file link: ' + files[0].link);
        },
        cancel: function() {},
        linkType: 'direct',
        multiselect: true,
    };
    var button = Dropbox.createChooseButton(options);
    document.getElementById('dropbox-link').appendChild(button);
    </script>
    """

def primerDetailsPage(params):
    """ create primers with primer3 around site identified by pamId in batch
    with batchId. Output primers as html
    """
    batchId, pamId, pam = params["batchId"], params["pamId"], params["pam"]
    inSeq, genome, pamSeq, position, extSeq = readBatchParams(batchId)
    seqLen = len(inSeq)
    batchBase = join(batchDir, batchId)

    # find position of guide sequence in genome at MM0
    otBedFname = batchBase+".bed"
    otMatches = parseOfftargets(otBedFname)
    if pamId not in otMatches or 0 not in otMatches[pamId]:
        errAbort("No perfect match found for guide sequence in the genome. Cannot design primers for a non-matching guide sequence.<p>Are you sure you have selected the right genome? <p> If you have selected the right genome and entered a cDNA as the query sequence, please note that sequences that overlap a splice site are not part of the genome and cannot be used as guide sequences.")

    matchList = otMatches[pamId][0] # = get all matches with 0 mismatches
    if len(matchList)!=1:
        errAbort("Multiple perfect matches for this guide sequence. Cannot design primer. Please select another guide sequences or email services@tefor.net to discuss your strategy or modifications to this software.")
        # XX we could show a dialog: which match do you want to design primers for?
        # But who would want to use a guide sequence that is not unique?

    chrom, start, end, seq, strand, segType, segName, x1Count = matchList[0]
    start = int(start)
    end = int(end)

    # retrieve guideSeq + PAM sequence from input sequence
    # XX guide sequence must not appear twice in there
    pamFname = batchBase+".fa"
    pams = parseFasta(open(pamFname))
    guideSeq = pams[pamId]
    guideStrand = pamId[-1]
    guideSeqWPam = seq

    if strand=="+":
        guideStart = inSeq.find(guideSeq)
        highlightSeq = guideSeqWPam
    else:
        guideStart = inSeq.find(revComp(guideSeq))
        highlightSeq = revComp(guideSeqWPam)

    lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq = \
        designPrimer(genome, chrom, start, end, strand, 0, batchId)

    guideStart = targetSeq.upper().find(highlightSeq.upper())
    guideEnd = guideStart + len(highlightSeq)

    if not chrom.startswith("ch"):
        chromLong = "chr"+chrom
    else:
        chromLong = chrom

    seqParts = ["<i><u>%s</u></i>" % targetSeq[:len(lSeq)] ] # primer 1
    seqParts.append("&nbsp;")
    seqParts.append(targetSeq[len(lSeq):guideStart]) # sequence before guide

    seqParts.append("<strong>") 
    seqParts.append(targetSeq[guideStart:guideEnd]) # guide sequence including PAM
    seqParts.append("</strong>")

    seqParts.append(targetSeq[guideEnd:len(targetSeq)-len(rSeq)])# sequence after guide

    seqParts.append("&nbsp;")
    seqParts.append("<i><u>%s</u></i>" % targetSeq[-len(rSeq):]) # primer 2

    targetHtml = "".join(seqParts)

    # prettify guideSeqWPam to highlight the PAM
    guideSeqHtml = "%s %s" % (guideSeqWPam[:-len(pam)], guideSeqWPam[-len(pam):])

    print '''<div style='width: 80%; margin-left:10%; margin-right:10%; text-align:left;'>'''
    print "<h2>Guide sequence: %s</h2>" % (guideSeqHtml)

    print "<h3>Validation Primers</h3>"
    guidePos = int(pamId.strip("s+-"))+1
    guideStrand = pamId[-1]
    if guideStrand=="+":
        primerGuideName = str(guidePos)+"forw"
    else:
        primerGuideName = str(guidePos)+"rev"

    print '<table class="primerTable">'
    print '<tr>'
    print "<td>guideRna%sLeft</td>" % primerGuideName
    print "<td>%s</td>" % (lSeq)
    print "<td>Tm %s</td>" % (lTm)
    print "</tr><tr>"
    print "<td>guideRna%sRight</td>" % primerGuideName
    print "<td>%s</td>" % (rSeq)
    print "<td>Tm %s</td>" % (rTm)
    print '</tr></table>'

    print "<h3>Genome fragment with validation primers and guide sequence</h3>"
    if strand=="-":
        print("Your guide sequence is on the reverse strand relative to the genome sequence, so it is reverse complemented in the sequence below.<p>")

    print '''<div style='word-wrap: break-word; word-break: break-all;'>'''
    print "<strong>Genomic sequence %s:%d-%d including primers, forward strand:</strong><br> <tt>%s</tt><br>" % (chromLong, start, end, targetHtml)
    print '''</div>'''
    print "<strong>Sequence length:</strong> %d<p>" % (rPos-lPos)
    print '<small>Method: Primer3.2 with default settings, target length 400-600 bp</small>'

    # restriction enzymes
    allEnzymes = readEnzymes()
    pamSeq = seq[-len(pam):]
    mutEnzymes = matchRestrEnz(allEnzymes, guideSeq, pamSeq)
    if len(mutEnzymes)!=0:
        print "<h3>Restriction Enzyme Sites for PCR product validation</h3>"

        print "Cas9 induces mutations next to the PAM site."
        print "If a mutation is induced, then it is very likely that one of the followingenzymes no longer cuts your PCR product amplified from the mutant sequence."
        print "For each restriction enzyme, the guide sequence with the restriction site underlined is shown below.<p>"

        for enzName, posList in mutEnzymes.iteritems():
            print "<strong>%s</strong>:" % enzName
            for start, end in posList:
                print markupSeq(guideSeqWPam, start, end)
            print "<br>"

    # primer helper

    print """
    <style>
        table.primerTable {
            border-width: 1px;
            border-color: #DDDDDD;
            border-collapse: collapse;
        }
        table.primerTable td {
            border-width: 1px;
            border-color: #DDDDDD;
            border-collapse: collapse;
        }
    </style>
    """

    print "<hr>"
    print "<h2>Expression of guide RNA</h2>"

    print "<h4>Summary of all primers explained below</h4>"
    primers = makeHelperPrimers(primerGuideName, guideSeq)
    printPrimerTableAll(primers)

    print "<p>Depending on the biological system studied, different options are available for expression of Cas9 and guide RNAs. In zebrafish embryos, guide RNAs and Cas9 are usually made by in vitro transcription and co-injected. In mammalian cells, guide RNAs and Cas9 are usually expressed from transfected plasmid and typically driven by U6 and CMV promoters."


    # T7 from plasmids
    print "<h3>In vitro with T7 RNA polymerase from plasmid DNA</h3>"
    print 'To produce guide RNA by in vitro transcription with T7 RNA polymerase, the guide RNA sequence can be cloned in a variety of plasmids (see <a href="http://addgene.org/crispr/empty-grna-vectors/">AddGene website</a>).<br>'
    print "For the guide sequence %s, the following primers should be ordered for cloning into the BsaI-digested plasmid DR274 generated by the Joung lab<p>" % guideSeqWPam

    printPrimerTable(primers["T7"])

    # T7 from primers, in vitro
    print "<h3>In vitro with T7 polymerase using overlapping oligonucleotides</h3>"
    print "Template for in vitro synthesis of guide RNA with T7 RNA polymerase can be prepared by annealing and primer extension of the following primers:<p>"

    printPrimerTable(primers["T7iv"], onRows=True)

    print 'The protocol for template preparation from oligonucleotides and in-vitro transcription can be found in <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4038517/?report=classic">Gagnon et al. PLoS ONE 2014</a>.'

    # MAMMALIAN CELLS
    print "<h3>In mammalian cells from plasmid</h3>"
    if "tttt" in guideSeq.lower():
        print "The guide sequence %s contains the motif TTTT, which terminates RNA polymerase. This guide sequence cannot be transcribed in mammalian cells." % guideSeq
    else:
        print "The guide sequence %s does not contain the motif TTTT, which terminates RNA polymerase. This guide sequence can be transcribed in mammalian cells." % guideSeq

        print "<br>"
        print "To express guide RNA in mammalian cells, a variety of plasmids are available. For example, to clone the guide RNA sequence in the plasmid pM3636, where guide RNA expression is driven by a human U6 promoter, the following primers should be used :"
        print "<br>"
        if "mammCellsNote" in primers:
            print("<strong>Note:</strong> Efficient transcription from the U6 promoter requires a 5' G. This G has been added to the sequence below, where it is underlined.<br>")

        printPrimerTable(primers["mammCells"])
    print "<hr>"
    print '</div>'

def runTests():
    guideSeq = "CTCTTTACGCAGAGGGATGT"
    testRes = {"ATTTTTATGCAGAGTGATGT":     0.4, 
               "TTCTTTACCCGGAGGGATGA": 0.2, 
               "CTGTTTACACACAGGGATTT": 0.2, 
               "CTCTCTGTGCAGATGGATGT": 0.1, 
               "ATCTTAAAGCAGATGGATGT": 0.1, 
               "CTCTTTCCGCAGAGGCTTGT": 0.1, 
               "CTCGTAGCGCAGAGGGAGGT": 0.1, 
               "CTCTTTAAAGAGATGGATGT": 0.1, 
               "CACTTCACTCAGAGGCATGT": 0.1, 
               "CTTTTTTCTCAGAAGGATGT": 0.1, 
               "CTCTTTACACAGAGAGACGT": 0.1, 
               "CTCTTTTCTCAGAGAGATGG": 0.1, 
               "CTATTTACCCAAATGGATGT": 0.1, 
               "CTCTTTGCACAGGGGGAAGT": 0, 
               "CTCTTTGCACAGGGGGAAGT": 0, 
               "CTCTTCACACAGAGGAATGA": 0, 
               "CTCTTTCCACAGGGGAATGT": 0 }

    testRes2 = {
       "GAGTCTAAGCAGAAGAAGAA":     2.2,
       "GAGTCCTAGCAGGAGAAGAA": 1.8,
       "GAGAGCAAGCAGAAGAAGAA": 1.6,
       "GAGTACTAGAAGAAGAAGAA": 1.6,
       "ACGTCTGAGCAGAAGAAGAA": 1.5,
       "GCGACAGAGCAGAAGAAGAA": 1.5,
       "GAGTAGGAGGAGAAGAAGAA": 1.4,
       "GATGCCGTGAAGAAGAAGAA": 1.3,
       "GATTCCTACCAGAAGAAGAA": 1,
       "GAATCCAAGCAGAAGAAGAG": 1,
       "AAGTACTGGCAGAAGAAGAA": 0.9,
       "AGGTGCTAGCAGAAGAAGAA": 0.9,
       "GGGGCCAGGCAGAAGAAGAA": 0.9,
       "ATGTGCAAGCAGAAGAAGAA": 0.9,
       "ACCTCCCAGCAGAAGAAGAA": 0.9,
       "CCCTCCCAGCAGAAGAAGAA": 0.9,
       "TCATCCAAGCAGAAGAAGAA": 0.9,
       "TTCTCCAAGCAGAAGAAGAA": 0.9,
       "GGTGCCAAGCAGAAGAAGAA": 0.9,
       "GCACCCCAGCAGAAGAAGAA": 0.9,
       "CAGTCCAGGAAGAAGAAGAA": 0.9,
       "AAGCCCAAGGAGAAGAAGAA": 0.9,
       "CACTCCAAGTAGAAGAAGAA": 0.9,
       "GAGTCCGGGAAGGAGAAGAA": 0.9,
       "GGTTCCCAGGAGAAGAAGAA": 0.9,
       "AAGTCTGAGCACAAGAAGAA": 0.9,
       "GAGGACAAGAAGAAGAAGAA": 0.9,
       "GTCTGCGATCAGAAGAAGAA": 0.8,
       "GGTTCTGTGCAGAAGAAGAA": 0.8,
       "AGGTGGGAGCAGAAGAAGAA": 0.8,
       "AAGAGCGAGCGGAAGAAGAA": 0.8,
       "CAATTTGAGCAGAAGAAGAA": 0.8,
       "AATACAGAGCAGAAGAAGAA": 0.8,
       "CAAACGGAGCAGAAGAAGAA": 0.8,
       "AAGTGAGAGTAGAAGAAGAA": 0.8,
       "AAGTAGGAGAAGAAGAAGAA": 0.8,
       "AAGTTGGAGAAGAAGAAGAA": 0.8,
       "CAGGCTGAGAAGAAGAAGAA": 0.8,
       "TAGTCAGGGGAGAAGAAGAA": 0.8,
       "TAGTCAGGGGAGAAGAAGAA": 0.8,
       "AAGTGGGAGGAGAAGAAGAA": 0.8,
       "TAGTCAGGGGAGAAGAAGAA": 0.8,
       "TCTTCCGAGCTGAAGAAGAA": 0.8,
       "GCGGCCGATGAGAAGAAGAA": 0.8,
       "GCGTCCGCCAAGAAGAAGAA": 0.8,
       "GCTCCTGAGCAGAAGAAGAA": 0.8,
       "CACTCTGAGGAGAAGAAGAA": 0.8,
       "GTGTGGGAGGAGAAGAAGAA": 0.8,
       "GGGTAAGAGTAGAAGAAGAA": 0.8
    }
    #for seq, expScore in testRes.iteritems():
        #score = calcHitScore(guideSeq, seq)
        #print score, "%0.1f" % score, expScore

    guideSeq = "GAGTCCGAGCAGAAGAAGAA"
    for seq, expScore in testRes2.iteritems():
        score = calcHitScore(guideSeq, seq)
        #print score, "%0.1f" % score, expScore
    
def parseArgs():
    " parse command line options into args and options "
    parser = optparse.OptionParser("""usage: %prog [options] org fastaInFile guideOutFile 

Command line interface for the Crispor tool.

    org          = genome identifier, like hg19 or ensHumSap
    fastaInFile  = Fasta file
    guideOutFile = tab-sep file, one row per guide

    If many guides have to scored in batch: Add GGG to them to make them valid
    guides, separate these sequences by N characters and supply as a single
    fasta sequence, a few dozen to ~100 per file.
    """)

    parser.add_option("-d", "--debug", dest="debug", \
        action="store_true", help="show debug messages, do not delete temp directory")
    parser.add_option("-t", "--test", dest="test", \
        action="store_true", help="run internal tests")
    parser.add_option("-p", "--pam", dest="pam", \
        action="store", help="PAM-motif to use, default %default", default="NGG")
    parser.add_option("-o", "--offtargets", dest="offtargetFname", \
        action="store", help="write offtarget info to this filename")
    parser.add_option("-m", "--maxOcc", dest="maxOcc", \
        action="store", type="int", help="MAXOCC parameter, 20mers with more matches are excluded")
    parser.add_option("", "--mm", dest="mismatches", \
        action="store", type="int", help="maximum number of mismatches, default %default", default=4)
    parser.add_option("", "--gap", dest="allowGap", \
        action="store_true", help="allow a gap in the match")
    parser.add_option("", "--skipAlign", dest="skipAlign", \
        action="store_true", help="do not align the input sequence. The on-target will be a random match with 0 mismatches.")
    parser.add_option("", "--minAltPamScore", dest="minAltPamScore", \
        action="store", type="float", help="minimum off-target score for alternative PAMs, default %default", \
        default=ALTPAMMINSCORE)
    parser.add_option("", "--worker", dest="worker", \
        action="store_true", help="Run as worker process: watches job queue and runs jobs")
    parser.add_option("", "--user", dest="user", \
        action="store", help="for the --worker option: switch to this user at program start")
    parser.add_option("", "--clear", dest="clear", \
        action="store_true", help="clear the worker job table and exit")
    parser.add_option("-g", "--genomeDir", dest="genomeDir", \
        action="store", help="directory with genomes, default %default", default=genomesDir)

    #parser.add_option("-f", "--file", dest="file", action="store", help="run on file") 
    (options, args) = parser.parse_args()

    if len(args)==0 and not options.test and not options.worker and not options.clear:
        parser.print_help()
        sys.exit(0)

    if options.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    return args, options

def main():
    # detect if running under apache or not
    if 'REQUEST_METHOD' in os.environ and sys.argv[-1] != "--worker":
        mainCgi()
    else:
        mainCommandLine()
    
def delBatchDir():
    " called at program exit, in command line mode "
    if not isdir(batchDir):
        return
    logging.debug("Deleting dir %s" % batchDir)
    fnames = glob.glob(join(batchDir, "*"))
    if len(fnames)>50:
        raise Exception("cowardly refusing to remove many temp files")
    for fname in fnames:
        os.remove(fname)
    os.removedirs(batchDir)

def runQueueWorker(userName):
    " in an infinite loop, take jobs from the job queue in jobs.db and finish them "
    if userName!=None:
        uid =  pwd.getpwnam(userName)[2]
        os.setuid(uid)

    try:
       # Store the Fork PID
       pid = os.fork()

       if pid > 0:
         print 'PID: %d' % pid
         os._exit(0)

    except OSError, error:
      print 'Unable to fork. Error: %d (%s)' % (error.errno, error.strerror)
      os._exit(1)

    print("%s: Worker daemon started. Waiting for jobs." % datetime.ctime(datetime.now()))
    sys.stdout.flush()

    q = JobQueue(JOBQUEUEDB)
    while True:
        if q.waitCount()==0:
            #q.dump()
            time.sleep(1+random.random()/10)
            continue

        jobType, batchId, paramStr = q.popJob()
        if jobType=="search":
            print "found job"
            jobError = False
            try:
                seq, org, pam, position, extSeq = readBatchParams(batchId)
                uppSeq = seq.upper()
                startDict, endSet = findAllPams(uppSeq, pam)
                print "searching for offtargets:  ", seq, org, pam, position
                getOfftargets(uppSeq, org, pam, batchId, startDict, q)
            except:
                exStr = traceback.format_exc()
                print " - WORKER CRASHED WITH EXCEPTION -"
                print exStr
                q.startStep(batchId, "crash", exStr)
                print exStr
                jobError = True

            if not jobError:
                q.jobDone(batchId)
        elif jobType is None:
            logging.error("No job")
        else:
            #raise Exception()
            logging.error("Illegal jobtype: %s - %s. Marking as done." % (jobType, batchId))
            q.jobDone(batchId)

def clearQueue():
    " empty the job queue "
    q = JobQueue(JOBQUEUEDB)
    q.clearJobs()
    print("Worker queue now empty")

# this won't work, it's very hard to spawn from CGIs
#def spawnWorkers(minCount):
#    #" check if we're running minCount worker threads, if not, start them up "
#    workerCount = 0
#    for line in os.popen("ps aux"):
#        sys.stderr.write(line+"\n")
#        if line.endswith("--worker") and "python" in line and basename(__file__) in line: 
#            workerCount += 1
#
#    if workerCount > minCount:
#        sys.stderr.write("ERROR: too many workers\n")
#        return
#
#    if workerCount < minCount:
#        sys.stderr.write("Found %d running workers, need %d\n" % (workerCount, minCount))
#        #for i in range(workerCount+1, minCount+1):
#            ## fork processes. Need to use at because apache will kill all direct children
#            #os.system("at now <<< %s --worker" % __file__)
#            #logFname = join(batchDir, "worker%d.log" % i)
#            #os.system("echo '%s --worker >> %s 2>&1 ' | at now" % (__file__, logFname))
#            #sys.stderr.write("Started one worker, log file %s\n" % logFname)
#        # see http://stackoverflow.com/questions/6024472/start-background-process-daemon-from-cgi-script
#        sys.stdout.flush()
#        os.close(sys.stdout.fileno()) # Break web pipe
#        sys.stderr.flush()
#        os.close(sys.stderr.fileno()) # Break web pipe
#        if os.fork(): # Get out parent process
#           return
#
#runQueueWorker()

def mainCommandLine():
    " main entry if called from command line "
    global commandLineMode
    commandLineMode = True

    args, options = parseArgs()

    if options.test:
        runTests()
        import doctest
        doctest.testmod()
        sys.exit(0)

    #if options.ajax:
        #sendStatus(options.ajax)

    if options.worker:
        runQueueWorker(options.user)
        sys.exit(0)

    if options.clear:
        clearQueue()
        sys.exit(0)

    org, inSeqFname, outGuideFname = args

    # different genomes directory?
    if options.genomeDir != None:
        global genomesDir
        genomesDir = options.genomeDir

    # handle the alignment/filtering options
    if options.maxOcc != None:
        global MAXOCC
        MAXOCC = options.maxOcc
        HIGH_MAXOCC = options.maxOcc

    if options.minAltPamScore!=None:
        global ALTPAMMINSCORE
        ALTPAMMINSCORE = options.minAltPamScore

    if options.mismatches:
        global maxMMs
        maxMMs = options.mismatches

    if options.allowGap:
        global allowGap
        allowGap = True

    skipAlign = False
    if options.skipAlign:
        skipAlign = True

    # get sequence
    seqs = parseFasta(open(inSeqFname))

    # make a directory for the temp files
    # and put it into the global variable, so all functions will use it
    global batchDir
    batchDir = tempfile.mkdtemp(dir=TEMPDIR, prefix="crispor")
    logging.debug("Created directory %s" % batchDir)
    if options.debug:
        logging.info("debug-mode, temporary directory %s will not be deleted" % batchDir)
    else:
        atexit.register(delBatchDir)

    # prepare output files
    guideFh = open(join(batchDir, "guideInfo.tab"), "w")
    guideFh.write("\t".join(guideHeaders)+"\n")
    if options.offtargetFname:
        offtargetFh = open(join(batchDir, "offtargetInfo.tab"), "w")
        offtargetFh.write("\t".join(offtargetHeaders)+"\n")

    for seqId, seq in seqs.iteritems():
        seq = seq.upper()
        logging.info("running on sequence ID '%s'" % seqId)
        # get the other parameters and write to a new batch
        seq = seq.upper()
        pam = options.pam
        batchId, position, extSeq = newBatch(seq, org, pam, skipAlign)
        logging.debug("Temporary output directory: %s/%s" % (batchDir, batchId))

        if position=="?":
            logging.error("no match found for sequence %s in genome %s" % (inSeqFname, org))

        startDict, endSet = findAllPams(seq, pam)
        otBedFname = getOfftargets(seq, org, pam, batchId, startDict, ConsQueue())
        otMatches = parseOfftargets(otBedFname)

        effScores = readEffScores(batchId)

        guideData, guideScores, hasNotFound = mergeGuideInfo(seq, startDict, pam, otMatches, position, effScores)

        for row in iterGuideRows(guideData):
            guideFh.write("\t".join(row))
            guideFh.write("\n")

        if options.offtargetFname:
            for row in iterOfftargetRows(guideData):
                offtargetFh.write("\t".join(row))
                offtargetFh.write("\n")

    guideFh.close()
    shutil.move(guideFh.name, outGuideFname)
    logging.info("guide info written to %s" % outGuideFname)

    if options.offtargetFname:
        offtargetFh.close()
        shutil.move(offtargetFh.name, options.offtargetFname)
        logging.info("off-target info written to %s" % options.offtargetFname)

    if not options.debug:
       shutil.rmtree(batchDir)

def sendStatus(batchId):
    " send batch status as json "
    print "Content-type: application/json\n"
    q = JobQueue(JOBQUEUEDB)
    status = q.getStatus(batchId)
    if "Traceback" in status:
        d = {"status" : "An error occured. Please send an email to services@tefor.net and tell us that the failing batchId was %s. We can usually fix this quickly. Thanks!" % batchId}
    else:
        d = {"status":status}
    print json.dumps(d)

def mainCgi():
    " main entry if called from apache "
    # XX IS THE SCRIPT SYMLINKED ? XX
    if os.getcwd()!="/var/www/crispor":
        # only activate stackdumps if running on a development machine
        import cgitb
        cgitb.enable()

    # parse incoming parameters
    params = getParams()
    batchId = None

    if "batchId" in params and "download" in params:
        downloadFile(params)
        return

    if "ajaxStatus" in params and "batchId" in params:
        sendStatus(params["batchId"])
        return

    # save seq/org/pam into a cookie, if they were provided
    if "seq" in params and "org" in params and "pam" in params:
        seq, org, pam = params["seq"], params["org"], params["pam"]
        saveSeqOrgPamToCookies(seq, org, pam)
        #batchId = makeTempBase(seq, org, pam)

    # print headers
    print "Content-type: text/html\n"
    print "" # = end of http headers

    printHeader(batchId)
    printTeforBodyStart()

    printBody(params)     # main dispatcher, branches based on the params dictionary

    printTeforBodyEnd()
    print("</body></html>")

main()
