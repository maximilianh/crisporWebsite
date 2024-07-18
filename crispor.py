#!/data/www/crispor/venv/bin/python3
# if you do not want the hardcoded PATH above, delete this line and the one above to use the default Python3 interpreter
#!/usr/bin/env python3
# I know that this line looks unprofessional to you, but modifying the PATH on a shared Apache webserver is not obvious.

# the tefor crispr tool
# can be run as a CGI or from the command line

# python std library
import subprocess, tempfile, optparse, logging, atexit, glob, shutil
import http.cookies, time, sys, cgi, re, random, platform, os, pipes
import hashlib, base64, string, logging, operator, urllib.request, urllib.parse, urllib.error, time
import traceback, json, pwd, gzip, zlib

from io import StringIO
from collections import defaultdict, namedtuple
from datetime import datetime
from itertools import product
from os.path import abspath, basename, dirname, isdir, isfile, join, relpath

try:
    # prefer the pip package, it's more up-to-date than the native package
    import pysqlite3 as sqlite3
except:
    import sqlite3
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict # python2.6 users: run 'sudo pip install ordereddict'

# for matplotlib, improves "import" performance
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-cache"

# try to load external dependencies
# we're going into great lengths to create a readable error message
needModules = set(["pytabix", "twobitreader", "pandas", "matplotlib", "scipy"])
try:
    import tabix # if not found, install with 'pip install pytabix'
    needModules.remove("pytabix")
except:
    pass

try:
    import twobitreader # if not found, install with 'pip install twobitreader'
    needModules.remove("twobitreader")
except:
    pass

try:
    import pandas # required by doench2016 score. install with 'pip install pandas'
    needModules.remove("pandas")
    import scipy # required by doench2016 score. install with 'pip install scipy'
    needModules.remove("scipy")
    import matplotlib # required by doench2016 score. install with 'pip install matplotlib'
    needModules.remove("matplotlib")
    import numpy # required by doench2016 score. install with 'pip install numpy'
    needModules.remove("numpy")
except:
    pass

if len(needModules)!=0:
    print("Content-type: text/html\n")
    print(("Python interpreter path: %s<p>" % sys.executable))
    print(("These python modules were not found: %s<p>" % ",".join(needModules)))
    print("To install all requirements in one line, run: sudo pip install biopython numpy scikit-learn==0.16.1 pandas twobitreader<p>")
    sys.exit(0)

# our own eff scoring library
import crisporEffScores

# don't report print as an error
# pylint: disable=E1601

# optional module for Excel export as native .xls files
# install with 'apt-get install python-xlwt' or 'pip install xlwt'
xlwtLoaded = True
try:
    import xlwt
except:
    sys.stderr.write("crispor.py - warning - the python xlwt module is not available\n")
    xlwtLoaded = False

# optional module for mysql support
#try:
    #import MySQLdb
    #mysqldbLoaded = True
#except:
    #mysqldbLoaded = False

# version of crispor
versionStr = "5.2"

# contact email
contactEmail='crispor@tefor.net'

# url to this server
ctBaseUrl = "http://crispor-max.tefor.net/temp/customTracks"

# write debug output to stdout
DEBUG = False
#DEBUG = True

# use bowtie for off-target search?
useBowtie = False

# calculate the efficienc scores?
doEffScoring = True

# system-wide temporary directory
#TEMPDIR = os.environ.get("TMPDIR", "/var/tmp")
TEMPDIR = "/var/tmp"

# a hack for cluster jobs at UCSC:
# - default to ramdisk
if isdir("/scratch/tmp"):
    TEMPDIR = "/dev/shm/"

# skipAlign is useful if your input sequence is not in the genome at all
# - don't do bwasw
#    - this will trigger auto-ontarget: any perfect match is the on-target
# - do not calculate efficiency scores
skipAlign = False

# prefix in html statements before the directories "image/", "style/" and "js/"
HTMLPREFIX =  ""
# alternative directory on local disk where image/, style/ and js/ are located
HTMLDIR = "/usr/local/apache/htdocs/crispor/"

# directory of crispor.py
baseDir = dirname(__file__)

# filename of this script, usually crispor.py
myName = basename(__file__)

# the segments.bed files use abbreviated genomic region names
segTypeConv = {"ex":"exon", "in":"intron", "ig":"intergenic"}

# directory for processed batches of offtargets ("cache" of bwa results)
batchDir = join(baseDir,"temp")

# sqlite3 db with gzipped old json batch files, to avoid hitting the ext4 inode limits
batchArchive = "/data/crisporJobArchive.db"

# the file where the sqlite job queue is stored
#JOBQUEUEDB = join(TEMPDIR, "crisporJobs.db") # TEMPDIR is mapped away for security reasons under Redhat/Centos for CGIs
JOBQUEUEDB = "/data/www/temp/crisporJobs.db"

# alternatively: connection info for mysql
jobQueueMysqlConn = {"socket":None, "host":None, "user": None, "password" : None}

# directory for platform-independent scripts (e.g. Heng Li's perl SAM parser)
scriptDir = join(baseDir, "bin")

# directory for helper binaries (e.g. BWA)
binDir = abspath(join(baseDir, "bin", platform.system()))

# directory for genomes
genomesDir = join(baseDir, "genomes")

DEFAULTORG = 'hg19'
DEFAULTSEQ = 'cttcctttgtccccaatctgggcgcgcgccggcgccccctggcggcctaaggactcggcgcgccggaagtggccagggcgggggcgacctcggctcacagcgcgcccggctattctcgcagctcaccatgGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCC'

# used if hg19 is not available
ALTORG = 'sacCer3'
ALTSEQ = 'ATTCTACTTTTCAACAATAATACATAAACatattggcttgtggtagCAACACTATCATGGTATCACTAACGTAAAAGTTCCTCAATATTGCAATTTGCTTGAACGGATGCTATTTCAGAATATTTCGTACTTACACAGGCCATACATTAGAATAATATGTCACATCACTGTCGTAACACTCT'

pamDesc = [ ('NGG','20bp-NGG - Sp Cas9, SpCas9-HF1, eSpCas9 1.1'),
         ('NNG','20bp-NNG - Cas9 S. canis'),
         ('NGN','20bp-NGN - SpG'),
         ('NNGT','20bp-NNGT - Cas9 S. canis - high efficiency PAM, recommended'),
         ('NAA','20bp-NAA - iSpyMacCas9'),
         #('TTN','TTN-23bp - Cpf1 F. Novicida'), # Jean-Paul: various people have shown that it's not usable yet
         ('NNGRRT','21bp-NNG(A/G)(A/G)T - Cas9 S. Aureus'),
         ('NNGRRT-20','20bp-NNG(A/G)(A/G)T - Cas9 S. Aureus with 20bp-guides'),
         ('NGK','20bp-NG(G/T) - xCas9, recommended PAM, see notes'),
         #('NGN','20bp-NGN or GA(A/T) - xCas9 (low efficiency, not recommended)'),
         #('NGG-BE1','20bp-NGG - BaseEditor1, modifies C->T'),
         ('NNNRRT','21bp-NNN(A/G)(A/G)T - KKH SaCas9'),
         ('NNNRRT-20','20bp-NNN(A/G)(A/G)T - KKH SaCas9 with 20bp-guides'),
         ('NGA','20bp-NGA - Cas9 S. Pyogenes mutant VQR'),
         ('NNNNCC','24bp-NNNNCC - Nme2Cas9'),
         ('NGCG','20bp-NGCG - Cas9 S. Pyogenes mutant VRER'),
         ('NNAGAA','20bp-NNAGAA - Cas9 S. Thermophilus'),
         ('NGGNG','20bp-NGGNG - Cas9 S. Thermophilus'),
         ('NNNNGMTT','20bp-NNNNG(A/C)TT - Cas9 N. Meningitidis'),
         ('NNNNACA','20bp-NNNNACA - Cas9 Campylobacter jejuni, original PAM'),
         ('NNNNRYAC','22bp-NNNNRYAC - Cas9 Campylobacter jejuni, revised PAM'),
         ('NNNVRYAC','22bp-NNNVRYAC - Cas9 Campylobacter jejuni, opt. efficiency'),
         ('TTCN','TTCN-20bp - CasX'),
         ('TTTV','TTT(A/C/G)-23bp - Cas12a (Cpf1)  - recommended, 23bp guides'),
         ('TTTV-21','TTT(A/C/G)-21bp - Cas12a (Cpf1) - 21bp guides recommended by IDT'),
         ('TTTN','TTTN-23bp - Cas12a (Cpf1) - low efficiency'),
         ('ATTN','ATTN-23bp - BhCas12b v4'),
         ('NGTN','NGTN-23bp - ShCAST/AcCAST, Strecker et al, Science 2019'),
         ('TYCV','T(C/T)C(A/C/G)-23bp - TYCV As-Cpf1 K607R'),
         ('TATV','TAT(A/C/G)-23bp - TATV As-Cpf1 K548V'),
         ('TTTA','TTTA-23bp - TTTA LbCpf1'),
         ('TCTA','TCTA-23bp - TCTA LbCpf1'),
         ('TCCA','TCCA-23bp - TCCA LbCpf1'),
         ('CCCA','CCCA-23bp - CCCA LbCpf1'),
         ('GGTT','GGTT-23bp - CCCA LbCpf1'),
         ('YTTV','YTTV-20bp - MAD7 Nuclease, Lui, Schiel, Maksimova et al, CRISPR J 2020'),
         ('TTYN','TTYN- or VTTV- or TRTV-23bp - enCas12a E174R/S542R/K548R - Kleinstiver et al Nat Biot 2019'),
         ('NNNNCNAA','20bp-NNNNCNAA - Thermo Cas9 - Walker et al, Metab Eng Comm 2020'),
         ('NNN','20bp-NNN - SpRY, Walton et al Science 2020'), # https://science.sciencemag.org/content/368/6488/290.abstract
         ('NRN','20bp-NRN - SpRY (high efficiency PAM)'),
         ('NYN','20bp-NYN - SpRY (low efficiency PAM)'),
         #('VTTV','(A/C)TT(A/C)-23bp - enCas12a S542R - Kleinstiver et al Nat Biot 2019'),
         #('TRTV','T(A/G)T(A/C)-23bp - enCas12a K548R - Kleinstiver et al Nat Biot 2019'),
       ]

DEFAULTPAM = 'NGG'

# the default base editor modification window
DEFAULTBEWIN = "1-7"

# for some PAMs, there are alternative main PAMs. These are also shown on the main sequence panel
multiPams = {
    #"NGN" : ["GAW"],
    "TTYN" : ["VTTV", "TRTV"]
}

# these PAMs are not specific. Allow only short sequences for them.
slowPams = ["TTYN", "NNG"]

# allow only very short sequences for these
verySlowPams = ["NNN", "NRN", "NYN"]

# for some PAMs, we allow other alternative motifs when searching for offtargets
# MIT and eCrisp do that, they use the motif NGG + NAG, we add one more, based on the
# on the guideSeq results in Tsai et al, Nat Biot 2014
# The NGA -> NGG rule was described by Kleinstiver...Young 2015 "Improved Cas9 Specificity..."
# NNGTRRT rule for S. aureus is in the new protocol "SaCas9 User manual"
# ! the length of the alternate PAM has to be the same as the original PAM!
offtargetPams = {
    "NGG" : ["NAG","NGA"],
    #"NGN" : ["GAW"],
    "NGK" : ["GAW"],
    "NGA" : ["NGG"],
    "NNGRRT" : ["NNGRRN"],
    "TTTV" : ["TTTN"],
    'ATTN' : ["TTTN", "GTTN"],
    "TTYN" : ["VTTV", "TRTV"]
}

# maximum size of an input sequence
MAXSEQLEN = 2300
# maximum input size when specifying "no genome"
MAXSEQLEN_NOGENOME = 25000
# maximum input size when using xCas9 or sCanis
MAXSEQLEN2 = 600
# maximum input size for NNN SpRY or similar PAMs
MAXSEQLEN3 = 150

# BWA: allow up to X mismatches
maxMMs=4

# maximum number of occurences in the genome to get flagged as repeats.
# This is used in bwa samse, when converting the same file
# and for warnings in the table output.
MAXOCC = 60000

# the BWA queue size is 2M by default. We derive the queue size from MAXOCC
MFAC = 2000000/MAXOCC

# the length of the guide sequence, set by setupPamInfo
GUIDELEN=None
# length of the PAM sequence
PAMLEN=None

# the name of the base editor, if any. This is the flag to activate
# baseEditor mode in the UI
baseEditor = None

# input sequences are extended by X basepairs so we can calculate the efficiency scores
# and can better design primers
FLANKLEN=100

# the name of the currently processed batch, assigned only once
# in readBatchParams and only for json-type batches
batchName = ""

# are we doing a Cpf1 run?
# this variable changes almost all processing and
# has to be set on program start, as soon as we know
# the PAM we're running on
pamIsFirst=None
saCas9Mode=False

# Highly-sensitive mode (not for CLI mode):
# MAXOCC is increased in processSubmission() and in the html UI if only one
# guide seq is run
# Also, the number of allowed mismatches is increased to 5 instead of 4
#HIGH_MAXOCC=600000
#HIGH_maxMMs=5

# minimum off-target score of standard off-targets (those that end with NGG)
# This should probably be based on the CFD score these days
# But for now, I'll let the user do the filtering
MINSCORE = 0.0

# minimum off-target score for alternative PAM off-targets
# There is not a lot of data to support this cutoff, but it seems
# reasonable to have at least some cutoff, as otherwise we would show
# NAG and NGA like NGG and the data shows clearly that the alternative
# PAMs are not recognized as well as the main NGG PAM.
# so for now, I just filter out very degenerative ones. the best solution
# would be to have a special penalty on the CFD score, but CFS does not
# support non-NGG PAMs (is this actually true?)
ALTPAMMINSCORE = 1.0

# how much shall we extend the guide after the PAM to match restriction enzymes?
pamPlusLen = 5

# global flag to indicate if we're run from command line or as a CGI
commandLineMode = False

# names/order of efficiency scores to show in UI
cas9ScoreNames = ["fusi", "crisprScan", "rs3"]
allScoreNames = ["fusi", "fusiOld", "chariRank", "ssc", "wuCrispr", "doench", "wang", "crisprScan", "aziInVitro", "ccTop", "rs3"]

mutScoreNames = []
spCas9MutScoreNames = ["oof", 'lindel'] # lindel is only added for spCas9
otherMutScoreNames = ["oof"] # lindel is only added for spCas9

cpf1ScoreNames = ["seqDeepCpf1"]

saCas9ScoreNames = ["najm"]

# to make the CFD more comparable to the MIT score, Nicholas Parkinson suggests to multiply it with 100.
# can be switched on with the URL argument fixCfd=1
doCfdFix=False

# how many digits shall we show for each score? default is 0
scoreDigits = {
    "ssc" : 1,
}

# List of AddGene plasmids, their long and short names:
addGenePlasmids = [
("43860", ("MLM3636 (Joung lab)", "MLM3636")),
("49330", ("pAc-sgRNA-Cas9 (Liu lab)", "pAcsgRnaCas9")),
("42230", ("pX330-U6-Chimeric_BB-CBh-hSpCas9 (Zhang lab) + derivatives", "pX330")),
("52961", ("lentiCRISPR v2 (Zhang lab)", "lentiCrispr")),
("52963", ("lentiGuide-Puro (Zhang lab)", "lentiGuide-Puro")),
]

addGenePlasmidsAureus = [
("61591", ("pX601-AAV-CMV::NLS-SaCas9-NLS-3xHA-bGHpA;U6::BsaI-sgRNA (Zhang lab)", "pX601")),
("61592", ("pX600-AAV-CMV::NLS-SaCas9-NLS-3xHA-bGHpA (Zhang lab)", "pX600")),
("61593", ("pX602-AAV-TBG::NLS-SaCas9-NLS-HA-OLLAS-bGHpA;U6::BsaI-sgRNA (Zhang lab)", "pX602")),
("65779", ("VVT1 (Joung lab)", "VVT1"))
]

# list of AddGene primer 5' and 3' extensions, one for each AddGene plasmid
# format: prefixFw, prefixRw, u6-G-suffix, restriction enzyme, link to protocol
addGenePlasmidInfo = {
"43860" : ("ACACC", "AAAAC", "G", "BsmBI", "https://www.addgene.org/static/data/plasmids/43/43860/43860-attachment_T35tt6ebKxov.pdf"),
"49330" : ("TTC", "AAC", "", "Bsp QI", "http://bio.biologists.org/content/3/1/42#sec-9"),
"42230" : ("CACC", "AAAC", "", "Bbs1", "https://www.addgene.org/static/data/plasmids/52/52961/52961-attachment_B3xTwla0bkYD.pdf"),
"52961" : ("CACC", "AAAC", "", "BsmBI", "https://www.addgene.org/static/data/plasmids/52/52961/52961-attachment_B3xTwla0bkYD.pdf"),
"61591" : ("CACC", "AAAC", "", "BsaI", "https://www.addgene.org/static/data/plasmids/61/61591/61591-attachment_it03kn5x5O6E.pdf"),
"61592" : ("CACC", "AAAC", "", "BsaI", "https://www.addgene.org/static/data/plasmids/61/61592/61592-attachment_iAbvIKnbqNRO.pdf"),
"61593" : ("CACC", "AAAC", "", "BsaI", "https://www.addgene.org/static/data/plasmids/61/61592/61592-attachment_iAbvIKnbqNRO.pdf"),
"65779": ("CACC", "AAAC", "", "BsmBI (aka Esp3l)", "https://www.addgene.org/static/data/plasmids/65/65779/65779-attachment_G8oNyvV6pA78.pdf"),
"52963": ("CACC", "AAAC", "", "BsmBI (aka Esp3l)", "https://www.addgene.org/static/data/plasmids/52/52963/52963-attachment_IPB7ZL_hJcbm.pdf")
}

# the barcodes for subpool tagging for oligo pool tables
satMutBarcodes = [
   (0,  "No Subpool barcode"),
   (1,  "Subpool 1: CGGGTTCCGT/GCTTAGAATAGAA"),
   (2,  "Subpool 2: GTTTATCGGGC/ACTTACTGTACC"),
   (3,  "Subpool 3: ACCGATGTTGAC/CTCGTAATAGC"),
   (4,  "Subpool 4: GAGGTCTTTCATGC/CACAACATA"),
   (5,  "Subpool 5: TATCCCGTGAAGCT/TTCGGTTAA"),
   (6,  "Subpool 6: TAGTAGTTCAGACGC/ATGTACCC"),
   (7,  "Subpool 7: GGATGCATGATCTAG/CATCAAGC"),
   (8,  "Subpool 8: ATGAGGACGAATCT/CACCTAAAG"),
   (9,  "Subpool 9: GGTAGGCACG/TAAACTTAGAACC"),
   (10, "Subpool 10: AGTCATGATTCAG/GTTGCAAGTCTAG"),
]

# Restriction enzyme supplier codes
rebaseSuppliers = {
"B":"Life Technologies",
"C":"Minotech",
"E":"Agilent",
"I":"SibEnzyme",
"J":"Nippon Gene",
"K":"Takara",
"M":"Roche",
"N":"NEB",
"O":"Toyobo",
"Q":"Molecular Biology Resources",
"R":"Promega",
"S":"Sigma",
"V":"Vivantis",
"X":"EURx",
"Y":"SinaClon BioScience"
}

# labels and descriptions of eff. scores
scoreDescs = {
    "doench" : ("Doench '14", "Range: 0-100. Linear regression model trained on 880 guides transfected into human MOLM13/NB4/TF1 cells (three genes) and mouse cells (six genes). Delivery: lentivirus. The Fusi score can be considered an updated version this score, as their training data overlaps a lot. See <a target='_blank' href='http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html'>Doench et al.</a>"),
    "wuCrispr" : ("Wu-Crispr", "Range 0-100. Aka 'Wong score'. SVM model trained on previously published data. The aim is to identify only a subset of efficient guides, many guides will have a score of 0. Takes into account RNA structure. See <a target='_blank' href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0784-0'>Wong et al., Gen Biol 2015</a>"),
    "ssc" : ("Xu", "Range ~ -2 - +2. Aka 'SSC score'. Linear regression model trained on data from &gt;1000 genes in human KBM7/HL60 cells (Wang et al) and mouse (Koike-Yusa et al.). Delivery: lentivirus. Ranges mostly -2 to +2. See <a target='_blank' href='http://genome.cshlp.org/content/early/2015/06/10/gr.191452.115'>Xu et al.</a>"),
    "crisprScan" : ["Moreno-Mateos", "Also called 'CrisprScan'. Range: mostly 0-100. Linear regression model, trained on data from 1000 guides on &gt;100 genes, from zebrafish 1-cell stage embryos injected with mRNA. See <a target=_blank href='http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3543.html'>Moreno-Mateos et al.</a>. Recommended for guides transcribed <i>in-vitro</i> (T7 promoter). Click to sort by this score. Note that under 'Show all scores', you can find a Doench2016 model trained on Zebrafish scores, Azimuth in-vitro, which should be slightly better than this model for zebrafish."],
    "wang" : ("Wang", "Range: 0-100. SVM model trained on human cell culture data on guides from &gt;1000 genes. The Xu score can be considered an updated version of this score, as the training data overlaps a lot. Delivery: lentivirus. See <a target='_blank' href='http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3972032/'>Wang et al.</a>"),
    "chariRank" : ("Chari", "Range: 0-100. Support Vector Machine, converted to rank-percent, trained on data from 1235 guides targeting sequences that were also transfected with a lentivirus into human 293T cells. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v12/n9/abs/nmeth.3473.html'>Chari et al.</a>"),
    "fusi" : ("Doench '16", "Aka the 'Fusi-Score', since V4.4 using the version 'Azimuth', scores are slightly different than before April 2018 but very similar (click 'show all' to see the old scores). Range: 0-100. Boosted Regression Tree model, trained on data produced by Doench et al (881 guides, MOLM13/NB4/TF1 cells + unpublished additional data). Delivery: lentivirus. See <a target='_blank' href='http://biorxiv.org/content/early/2015/06/26/021568'>Fusi et al. 2015</a> and <a target='_blank' href='http://www.nature.com/nbt/journal/v34/n2/full/nbt.3437.html'>Doench et al. 2016</a> and <a target=_blank href='https://crispr.ml/'>crispr.ml</a>. Recommended for guides expressed in cells (U6 promoter). Click to sort the table by this score."),
    "fusiOld" : ("OldDoench '16", "The original implementation of the Doench 2016 score, as received from John Doench. The scores are similar, but not exactly identical to the 'Azimuth' version of the Doench 2016 model that is currently the default on this site, since Apr 2018."),
    "rs3" : ("Doench-RuleSet3", "The Doench Rule Set 3 (RS3) score (-200-+200). Similar to the Doench 2014 and Doench 2016/Fusi/Azimuth score, but updated and more accurate. See <a href='https://www.nature.com/articles/s41467-022-33024-2' target=_blank>. Scores shown are multiplied with 100 for easier display. RS3 is configured here to use the Hsu-TRACR sequence."),
    "najm" : ("Najm 2018", "A modified version of the Doench 2016 score ('Azimuth'), by Mudra Hegde for S. aureus Cas9. Range 0-100. See <a target=_blank href='https://www.nature.com/articles/nbt.4048'>Najm et al 2018</a>."),
    "ccTop" : ("CCTop", "The efficiency score used by CCTop, called 'crisprRank'."),
    "aziInVitro" : ("Azimuth in-vitro", "The Doench 2016 model trained on the Moreno-Mateos zebrafish data. Unpublished model, gratefully provided by J. Listgarden. This should be better than Moreno-Mateos, but we have not found the time to evaluate it yet."),
    "housden" : ("Housden", "Range: ~ 1-10. Weight matrix model trained on data from Drosophila mRNA injections. See <a target='_blank' href='http://stke.sciencemag.org/content/8/393/rs9.long'>Housden et al.</a>"),
    "proxGc" : ("ProxGCCount", "Number of GCs in the last 4pb before the PAM"),
    "seqDeepCpf1" : ("DeepCpf1", "Range: ~ 0-100. Convolutional Neural Network trained on ~20k Cpf1 lentiviral guide results. This is the score without DNAse information, 'Seq-DeepCpf1' in the paper. See <a target='_blank' href='https://www.nature.com/articles/nbt.4061'>Kim et al. 2018</a>"),
    "oof" : ("Out-of-Frame", "Range: 0-100. Out-of-Frame score, only for deletions. Predicts the percentage of clones that will carry out-of-frame deletions, based on the micro-homology in the sequence flanking the target site. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html'>Bae et al. 2014</a>. Click the score to show the predicted deletions."),
    "lindel": ("Lindel", "Wei Chen Frameshift ratio (0-100). Predicts probability of a frameshift caused by any type of insertion or deletion. See <a href='https://academic.oup.com/nar/article/47/15/7989/5511473'>Wei Chen et al, Bioinf 2018</a>. Click the score to see the most likely deletions and insertions."),
}

# the headers for the guide and offtarget output files
guideHeaders = ["guideId", "targetSeq", "mitSpecScore", "cfdSpecScore", "offtargetCount", "targetGenomeGeneLocus"]
offtargetHeaders = ["guideId", "guideSeq", "offtargetSeq", "mismatchPos", "mismatchCount", "mitOfftargetScore", "cfdOfftargetScore", "chrom", "start", "end", "strand", "locusDesc"]

# library descriptions
libLabels = [
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4486245/
    ("human_brunello" , "Human, Brunello, Doench Nat Bio 2016 (recommended)"),
    ("human_avana" , "Human, Avana, Doench Nat Bio 2016"),
    ("human_geckov2" , "Human, GeCKO V2, Sanjana Nat Meth 2014"),
    ("mouse_brie" , "Mouse, Brie, Doench Nat Bio 2016 (recommended)"),
    ("mouse_geckov2" , "Mouse, GeCKO V2, Sanjana Nat Meth 2014"),
    ("mouse_asiago" , "Mouse, Asiago, Doench Nat Bio 2016"),
]

# a file crispor.conf in the directory of the script allows to override any global variable
myDir = dirname(__file__)
confPath =join(myDir, "crispor.conf")
if isfile(confPath):
    exec(open(confPath).read())
    #execfile(confPath)

cgiParams = None

# ====== END GLOBALS ============

def setupPamInfo(pam):
    " modify a few globals based on the current pam "
    global GUIDELEN
    global pamIsFirst
    global addGenePlasmids
    global PAMLEN
    global scoreNames
    global baseEditor
    global saCas9Mode
    global mutScoreNames
    global isSpg

    pamOpt = None
    if "-" in pam:
        pam, pamOpt = pam.split("-")
        if pamOpt=="BE1":
            baseEditor = "BE1"
        elif pamOpt=="spg":
            isSpg = True

    PAMLEN = len(pam)

    pamIsFirst = False
    scoreNames = cas9ScoreNames

    if pamIsCasX(pam):
        logging.debug("switching on CasX mode, guide length is 20bp")
        GUIDELEN = 20
        pamIsFirst = True
        scoreNames = cpf1ScoreNames
    if pamIsCpf1(pam):
        logging.debug("switching on Cpf1 mode, guide length is 23bp")
        GUIDELEN = 23
        pamIsFirst = True
        scoreNames = cpf1ScoreNames
        if pamOpt:
            GUIDELEN=int(pamOpt)
    elif pam=="NGTN":
        logging.debug("switching on Cpf1 mode for ShCAST, guide length is 23bp")
        GUIDELEN = 23
        pamIsFirst = True
    elif pam=="NNNNRYAC" or pam=="NNNVRYAC":
        GUIDELEN = 22
    elif pam=="NNGRRT" or pam=="NNNRRT":
        logging.debug("switching on S. aureus mode, guide length is 21bp")
        addGenePlasmids = addGenePlasmidsAureus
        GUIDELEN = 21
        if pamOpt=="20":
            GUIDELEN=20
        saCas9Mode = True
        scoreNames = saCas9ScoreNames
    elif pam=="NNNNCC":
        GUIDELEN = 24
    else:
        GUIDELEN = 20

    if GUIDELEN==20 and pam=="NGG":
        mutScoreNames = spCas9MutScoreNames
    else:
        mutScoreNames = otherMutScoreNames

    logging.debug("Enzyme info: pam=%s, guideLen=%d, pamIsFirst=%s, saCas9Mode=%s" %
        (pam, GUIDELEN, pamIsFirst, saCas9Mode))

    return pam


# ==== CLASSES =====
class JobQueue:
    """
    simple job queue, using a db table as a backend
    jobs have different types and status. status can be updated while they run
    job running times are kept and old job info is kept in a separate table

    >>> q = JobQueue()
    >>> q.openSqlite()
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
    'Waiting'
    >>> q.startStep("abc123", "bwa", "Alignment with BWA")
    >>> q.getStatus("abc123")
    'Alignment with BWA'
    >>> jobType, jobId, paramStr = q.popJob()

    >>> q.waitCount()
    0
    >>> q.jobDone("abc123")
    >>> q.waitCount()
    0

    can't pop from an empty queue
    #>>> q.popJob()
    #(None, None, None)
    #>>> os.system("rm /tmp/tempCrisporTest.db")
    #0
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
    '  startTime text ' # date+time when job was put into queue
    ')')

    #def __init__(self):
        #" no inheritance needed here "
        #self.openSqlite(JOBQUEUEDB)

    def openSqlite(self, dbName=JOBQUEUEDB):
        self.dbName = dbName
        self.conn = sqlite3.connect(dbName, timeout=10)
        # trying to increase stability and locks?
        if sqlite3.version_info[0] >= 3:
            self.conn.execute("PRAGMA journal_mode=WAL;")
        #self.conn.set_trace_callback(print) # for debugging: print all sql statements

        try:
            self.conn.execute(self._queueDef % ("queue", "PRIMARY KEY"))
        except sqlite3.OperationalError:
            errAbort("cannot open the file %s" % JOBQUEUEDB)
        self.conn.commit()
        self._chmodJobDb()

    def _chmodJobDb(self):
        # umask is not respected by sqlite, bug http://www.mail-archive.com/sqlite-users@sqlite.org/msg59080.html
        try:
            os.chmod(JOBQUEUEDB, 0o666)
        except OSError:
            # if the file was created by other job, we can't chmod, as we're the CGI. Just silently ignore this
            pass

    def addJob(self, jobType, jobId, paramStr):
        " create a new job, returns False if not successful  "
        self._chmodJobDb()

        sql = 'INSERT INTO queue (jobType, jobId, isRunning, lastUpdate, ' \
            'stepTimes, paramStr, stepName, stepLabel, startTime) VALUES (:jobType, :jobId, :isRunning, :lastUpdate, ' \
                ':stepTimes, :paramStr, :stepName, :stepLabel, :startTime)'
        now = "%.3f" % time.time()
        values = {'jobType' : jobType, 'jobId' : jobId, 'isRunning' : 0, 'lastUpdate' : now, 'stepTimes':"", 'paramStr':paramStr, 'stepName':"wait",
                "stepLabel":"Waiting", "startTime":now}
        try:
            cur = self.conn.cursor()
            cur.execute(sql, values)
            self.commitRetry()
            return True
        except sqlite3.IntegrityError as ev:
            # if the job is already in the queue, on Python3 (not Python2!) an integrity error will be thrown.
            # so this actually means that all is fine.
            # instead of checking if the job exists and then add it, we just add it and tolerate the error, saves one query.
            # (if the pipeline crashed or server was restarted, jobs may be on disk but not in the queue)
            return True
        except sqlite3.OperationalError:
            errAbort("Cannot open DB file %s. Please contact %s" % (self.dbName, contactEmail))

    def getStatus(self, jobId):
        " return current job status label or None if job is not in queue"
        sql = 'SELECT stepLabel FROM queue WHERE jobId=?'
        try:
            status = self.conn.execute(sql, (jobId,)).fetchmany(1)[0][0]
        except StopIteration:
            logging.debug("getStatus got StopIteration")
            status = None
        self.commitRetry()
        return status

    def dump(self):
        " for debugging, write the whole queue table to stdout "
        sql = 'SELECT * FROM queue'
        for row in self.conn.execute(sql):
            print("\t".join([str(x) for x in row]))

    def jobInfo(self, jobId, isDone=False):
        " for debugging, return all job info as a tuple "
        print("job info<br>")
        if isDone:
            sql = 'SELECT * FROM doneJobs WHERE jobId=?'
        else:
            sql = 'SELECT * FROM queue WHERE jobId=?'
        try:
            row = next(self.conn.execute(sql, (jobId,)))
        except StopIteration:
            return []
        return row

    def commitRetry(self):
        " try to commit 10 times "
        tryCount = 0
        while tryCount < 10:
            try:
                self.conn.commit()
                break
            except sqlite3.OperationalError:
                time.sleep(3)
                tryCount += 1

        if tryCount >= 30:
            raise Exception("Database locked for a long time")

    def execute(self, cmd, *args):
        " try to execute command 10 times "
        tryCount = 0
        while tryCount < 10:
            try:
                res = self.conn.execute(cmd, *args)
                return res
            except sqlite3.OperationalError:
                time.sleep(3)
                tryCount += 1

        if tryCount >= 10:
            raise Exception("Cannot execute %s, Database locked for a long time" % cmd)

    def startStep(self, jobId, newName, newLabel):
        " start a new step. Update lastUpdate, status and stepTime "
        self.startTransaction()
        sql = 'SELECT lastUpdate, stepTimes, stepName FROM queue WHERE jobId=?'
        logging.debug(sql)
        result = self.conn.execute(sql, (jobId,)).fetchmany(1)[0]
        logging.debug(result)
        lastTime, timeStr, lastStep = result
        logging.debug("end")
        lastTime = float(lastTime)

        # append a string in format "stepName:milliSecs" to the timeStr
        now = time.time()
        timeDiff = "%d" % int((1000.0*(now - lastTime)))
        newTimeStr = timeStr+"%s=%s" % (lastStep, timeDiff)+","

        sql = 'UPDATE queue SET lastUpdate=?, stepName=?, stepLabel=?, stepTimes=?, isRunning=? WHERE jobId=?'
        self.conn.execute(sql, (now, newName, newLabel, newTimeStr, 1, jobId))

        self.commitRetry()

    def startTransaction(self):
        logging.debug("Starting transaction")
        tryCount = 0
        while tryCount < 10:
            try:
                self.conn.execute('BEGIN') # indicate that transaction should start now
                break
            except sqlite3.OperationalError:
                logging.debug("Waiting since transaction start failed")
                time.sleep(3)
                tryCount += 1

        if tryCount >= 10:
            raise Exception("Database locked for a long time")

    def jobDone(self, jobId):
        " remove the job from the queue and add it to the queue log"
        print("job done<br>")
        self.startTransaction()
        sql = 'SELECT * FROM queue WHERE jobId=?'
        try:
            row = next(self.conn.execute(sql, (jobId,)))
        except StopIteration:
            # return if the job has already been removed
            logging.warn("jobDone - jobs %s has been removed already" % jobId)
            self.commitRetry() # release lock
            return

        sql = 'DELETE FROM queue WHERE jobId=?'
        self.conn.execute(sql, (jobId,))
        self.commitRetry() # release lock

        # good to have a log file of the old jobs
        with open("doneJobs.tsv", "a") as ofh: # if this triggers an error: run 'touch doneJobs.tsv && chmod a+rw doneJobs.tsv' in the crispor dir.
            row = [str(x) for x in row]
            line = "\t".join(row)
            ofh.write(line)
            ofh.write("\n")

    def waitCount(self):
        " return number of waiting jobs, wait until database is ready "
        sql = 'SELECT count(*) FROM queue WHERE isRunning=0'
        count = None
        while count is None:
            try:
                count = self.conn.execute(sql).fetchall()[0][0]
            except sqlite3.OperationalError:
                logging.debug("OperationalError on waitCount()")
                time.sleep(1+random.random()/10)
        return count

    def popJob(self):
        " return (jobType, jobId, params) of first waiting job and set it to running state "
        print('pop job<br>')
        self.startTransaction()
        sql = 'SELECT jobType, jobId, paramStr FROM queue WHERE isRunning=0 ORDER BY lastUpdate LIMIT 1'
        try:
            jobType, jobId, paramStr = next(self.execute(sql))
        except StopIteration:
            logging.debug("No data for '%s'" % sql)
            self.commitRetry() # unlock db
            return None, None, None

        sql = 'UPDATE queue SET isRunning=1 where jobId=?'
        self.execute(sql, (jobId,))
        self.commitRetry() # unlock db

        return jobType, jobId, paramStr

    def clearJobs(self):
        " clear the job table, removing running jobs, too "
        self.conn.execute("DELETE from queue")
        self.commitRetry()

    def close(self):
        " "
        self.conn.close()

# ====== FUNCTIONS =====
contentLineDone = False

# the queue workers should be able to never abort
doAbort = True

def getTwoBitFname(db):
    " return the name of the twoBit file for a genome "
    # at UCSC, try to use local disk, if possible
    locPath = join("/scratch", "data", db, db+".2bit")
    if isfile(locPath):
        return locPath
    path = join(genomesDir, db, db+".2bit")
    return path

def errAbort(msg, isWarn=False):
    " print err msg and exit "
    if commandLineMode:
        raise Exception(msg)

    if not contentLineDone:
        print("Content-type: text/html\n")

    print('<div style="position: absolute; padding: 10px; left: 100; top: 100; border: 10px solid black; background-color: white; text-align:left; width: 800px; font-size: 18px">')

    if isWarn:
        print("<strong>Warning:</strong><p> ")
    else:
        print("<strong>Error:</strong><p> ")

    print((msg+"<p>"))
    print(("If you think this is a bug or you have any other suggestions, please do not hesitate to contact us %s<p>" % contactEmail))
    if isWarn:
        print("In the email, please also send us the full URL of the page.")
    else:
        print("Please also send us the full URL of the page where you see the error. Thanks!")
    print('</div>')

    if doAbort:
        sys.exit(0)  # cgi must not exit with 1

# allow only dashes, digits, characters, underscores and colons in the CGI parameters
# and +
notOkChars = re.compile(r'[^+a-zA-Z0-9/:\n\r_. -]')

def checkVal(key, inStr):
    """ remove special characters from input string, to protect against injection attacks """
    if key!="geneIds":
        if len(inStr) > 10000:
            errAbort("input parameter %s is too long" % key)
    else:
        if len(inStr) > 100000:
            errAbort("Pasting more than tens of thousands of gene IDs makes little sense. Copy/paste error?")

    matchObj =notOkChars.search(inStr)
    if matchObj!=None:
        errAbort("input parameter %s contains an invalid character %s (ASCII %d)" % (key, repr(matchObj.group()), ord(matchObj.group())))
    return inStr

def cgiGetParams():
    " get CGI parameters and return as dict "
    form = cgi.FieldStorage()
    global cgiParams
    cgiParams = {}

    # parameters are:
    #"pamId", "batchId", "pam", "seq", "org", "download", "sortBy", "format", "ajax
    for key in list(form.keys()):
        val = form.getfirst(key)
        if val!=None:
            # "seq" is cleaned by cleanSeq later
            val = urllib.parse.unquote(val)
            if key not in ["seq", "name"]:
                checkVal(key, val)
            cgiParams[key] = val

    if "pam" in cgiParams:
        legalChars = set("ACTGNMKRYVBE120345/-")
        illegalChars = set(cgiParams["pam"])-legalChars
        if len(illegalChars)!=0:
            errAbort("Illegal character in PAM-sequence. Only %s are allowed."+"".join(legalChars))

    if "batchId" in cgiParams:
        batchId = cgiParams["batchId"]
        if not batchId.isalnum() or len(batchId) > 30:
            errAbort("Invalid batchId")

    return cgiParams

def cgiGetStr(params, argName, default=None):
    val = params.get(argName, None)
    if val==None and default==None:
        errAbort("'%s' parameter must be specified" % argName)
    if val==None:
        return default
    return val

def cgiGetNum(params, argName, default):
    " get CGI parameter which must be a number "
    val = params.get(argName, None)
    if val==None:
        return default
    if not val.isdigit():
        errAbort("'%s' parameter must be a number" % argName)
    val = int(val)
    return val

transTab = str.maketrans("-=/+_", "abcde")

def makeTempBase(seq, org, pam, batchName):
    "create the base name of temp files using a hash function and some prettyfication "
    hasher = hashlib.sha1(seq.encode("latin1")+org.encode("latin1")+pam.encode("latin1")+batchName.encode("latin1"))
    shortHash = hasher.digest()[0:20]
    batchId = base64.urlsafe_b64encode(shortHash).decode('latin1').translate(transTab)[:20]
    return batchId

def makeTempFile(prefix, suffix):
    " return a temporary file that is deleted upon exit, unless DEBUG is set "
    if DEBUG:
        fname = join("/tmp", prefix+suffix)
        fh = open(fname, "wt")
    else:
        fh = tempfile.NamedTemporaryFile(mode="wt", dir=TEMPDIR, prefix="primer3In", suffix=".txt")
    return fh

def pamIsCpf1(pam):
    " if you change this, also change bin/filterFaToBed and bin/samToBed!!! "
    return (pam in ["TTN", "TTTN", "TYCV", "TATV", "TTTV", "TTTR", "ATTN", "TTTA", "TCTA", "TCCA", "CCCA", "YTTV", "TTYN"])

def pamIsCasX(pam):
    " if you change this, also change bin/filterFaToBed and bin/samToBed!!! "
    return (pam in ["TTCN"])

def pamIsSaCas9(pam):
    " only used for notes and efficiency scores, unlike its Cpf1 cousin function "
    return (pam.split("-")[0] in ["NNGRRT", "NNNRRT"])

def isSlowPam(pam):
    " do not allow input sequences > 500 bp "
    if pamIsXCas9(pam) or pam=="TTYN" or pam=="NNG":
        return True
    else:
        return False

def pamIsXCas9(pam):
    " "
    return (pam in ["NGK", "NGN"])

def pamIsSpCas9(pam):
    " only used for notes and efficiency scores, unlike its Cpf1 cousin function "
    return (pam in ["NGG", "NGA", "NGCG"])

def saveSeqOrgPamToCookies(seq, org, pam):
    " create a cookie with seq, org and pam and print it"
    cookies=http.cookies.SimpleCookie()
    expires = 365 * 24 * 60 * 60
    if len(seq)<3000:
        cookies['lastseq'] = seq
    else:
        cookies['lastseq'] = "(last sequence was too long, could not be saved in Internet Browser cookie)"

    cookies['lastseq']['expires'] = expires
    cookies['lastorg'] = org
    cookies['lastorg']['expires'] = expires
    cookies['lastpam'] = pam
    cookies['lastpam']['expires'] = expires
    print(cookies)

def debug(msg):
    if commandLineMode:
        logging.debug(msg)
    elif DEBUG:
        print(msg)
        print("<br>")

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
    seq = seq.upper()
    pat = pat.upper()
    patLen = len(pat)
    for i in range(0, len(seq)-patLen+1):
        subseq = seq[i:i+patLen]
        if patMatch(subseq, pat):
            yield i

def rndSeq(seqLen):
    " return random seq "
    seq = []
    alf = "ACTG"
    for i in range(0, seqLen):
        seq.append(alf[random.randint(0,3)])
    return "".join(seq)

def cleanSeq(seq, db):
    """ remove fasta header, check seq for illegal chars and return (filtered
    seq, user message) special value "random" returns a random sequence.
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
            if c not in "actgACTGNn":
                nCount +=1
            else:
                newSeq.append(c)
    seq = "".join(newSeq)

    msgs = []
    tooLongHint = """
    Please split your input sequence into shorter sequences or use
    the <a href='downloads/'>stand-alone version</a> on your own Linux or Mac server to process longer sequences in batch.<br>
    """

    if len(seq)>MAXSEQLEN and db!="noGenome":
        errMsg = "<strong>Sorry, this tool cannot handle sequences longer than %d bp</strong><br>" % (MAXSEQLEN)
        errAbort(errMsg+tooLongHint)
    if len(seq)>MAXSEQLEN_NOGENOME and db=="noGenome":
        errMsg = "<strong>Sorry, this tool cannot handle sequences longer than %d bp when using the 'No Genome' option.</strong><br>" % (MAXSEQLEN_NOGENOME)
        errAbort(errMsg+tooLongHint)

    if nCount!=0:
        msgs.append("Sequence contained %d non-ACTGN letters. They were removed." % nCount)

    return seq, "<br>".join(msgs)

revTbl = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N' , 'M' : 'K', 'K' : 'M',
    "R" : "Y" , "Y":"R" , "g":"c", "a":"t", "c":"g","t":"a", "n":"n", "V" : "B", "v":"b",
    "B" : "V", "b": "v", "W" : "W", "w" : "w"}

def revComp(seq):
    " rev-comp a dna sequence with UIPAC characters "
    newSeq = []
    for c in reversed(seq):
        newSeq.append(revTbl[c])
    return "".join(newSeq)

def docTestInit(isCpf1, guideLen):
    global pamIsFirst
    global GUIDELEN
    pamIsFirst=isCpf1
    GUIDELEN=guideLen

def findPams (seq, pam, strand, startDict, endSet):
    """ return two values: dict with pos -> strand of PAM and set of end positions of PAMs
    Makes sure to return only values with at least GUIDELEN bp left (if strand "+") or to the
    right of the match (if strand "-")
    If the PAM is cpf1, then this is inversed: pos-strand matches must have at least GUIDELEN
    basepairs to the right, neg-strand matches must have at least GUIDELEN bp on their left
    >>> docTestInit(False, 20)
    >>> findPams("GGGGGGGGGGGGGGGGGGGGGGG", "NGG", "+", {}, set())
    ({20: '+'}, {23})
    >>> findPams("CCAGCCCCCCCCCCCCCCCCCCC", "CCA", "-", {}, set())
    ({0: '-'}, {3})
    >>> docTestInit(True, 20)
    >>> findPams("TTTNCCCCCCCCCCCCCCCCCTTTN", "TTTN", "+", {}, set())
    ({0: '+'}, {4})
    >>> docTestInit(False, 20)
    >>> findPams("CCCCCCCCCCCCCCCCCCCCCAAAA", "NAA", "-", {}, set())
    ({}, set())
    >>> findPams("AAACCCCCCCCCCCCCCCCCCCCC", "NAA", "-", {}, set())
    ({0: '-'}, {3})
    >>> findPams("CCCCCCCCCCCCCCCCCCCCCCCCCAA", "NAA", "-", {}, set())
    ({}, set())
    >>> findPams("GTTGTGTTTTACAATGCAGAGAGTGGAGGATGCTTTTTATACATTGGTGAGAGAGATCCGACAGTACAGATTGAAAAAAATCAGCAAAGAAGAAAAGACTCCTGGCTGTGTGAAAATTAAAAAATGCGTTATAATGTAATCTGGTAAGTTGAGCATATTCATTCTGGTACAAAGCAGATGTCTTCAGAGGTAACA", "TATV", "-", {}, set())
    ({37: '-', 129: '-'}, {41, 133})
    >>> findPams("GTTGTGTTTTACAATGCAGAGAGTGGAGGATGCTTTTTATACATTGGTGAGAGAGATCCGACAGTACAGATTGAAAAAAATCAGCAAAGAAGAAAAGACTCCTGGCTGTGTGAAAATTAAAAAATGCGTTATAATGTAATCTGGTAAGTTGAGCATATTCATTCTGGTACAAAGCAGATGTCTTCAGAGGTAACA", "TATV", "+", {}, set())
    ({37: '+', 129: '+'}, {41, 133})

    """
    assert(pamIsFirst is not None)

    if pamIsFirst:
        maxPosPlus  = len(seq)-(GUIDELEN+len(pam))
        minPosMinus = GUIDELEN
    else:
        # -------------------
        #          OKOKOKOKOK
        minPosPlus  = GUIDELEN
        # -------------------
        # OKOKOKOKOK
        maxPosMinus = len(seq)-(GUIDELEN+len(pam))

    #print "new search", seq, pam, "minPosPlus=",minPosPlus, "guideLen=", GUIDELEN, "<br>"
    for start in findPat(seq, pam):
        if pamIsFirst:
            # need enough flanking seq on one side
            #return("Cpf1 mode found", start,"<br>")
            if strand == "+" and start > maxPosPlus:
                continue
            if strand == "-" and start < minPosMinus:
                continue
        else:
            # return("non-Cpf1 mode found", start,"<br>")
            if strand=="+" and start < minPosPlus:
                continue
            if strand=="-" and start > maxPosMinus:
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

def varDictToHtml(varDict, seq, varShortLabel):
    " make a list of one html string per position in the sequence "
    if varDict is None:
        return None

    varHtmls = []
    for i in range(0, len(seq)):
        if not i in varDict:
            varHtmls.append(".")
        else:
            varHooverLines = []
            showStar = False # show a star if change is non-simple SNP
            varInfos = varDict[i]
            for chrom, pos, refAll, altAll, infoDict in varInfos:
                varHooverLines.append("%s: %s &rarr; %s<br>" % (varShortLabel, refAll, altAll))
                if "freq" in infoDict:
                    varHooverLines.append("&nbsp;<b>Freq:</b> %s<br>" % infoDict["freq"])
                #if "dbg" in infoDict:
                    #varHooverLines.append("%s<br>" % infoDict["dbg"])
                if "varId" in infoDict:
                    varHooverLines.append("&nbsp;<b>ID:</b> %s<br>" % infoDict["varId"])
                if "ExAC" in varDict["label"]:
                        endPos = int(pos)+len(refAll)
                        varHooverLines.append('&nbsp;<a target=_blank href="http://exac.broadinstitute.org/region/%s-%s-%d">ExAC Browser</a><br>' % (chrom, pos, endPos))

                if len(refAll)!=1 or len(altAll)!=1:
                    showStar = True

            if len(varInfos)!=1:
                showStar = True

            varDesc = "".join(varHooverLines)

            if showStar:
                dispChar = "*"
            else:
                dispChar = altAll
            varHtmls.append("<u class='tooltipsterInteract' title='%s'>%s</u>" % (varDesc, dispChar))
    return varHtmls

def cssClassesFromSeq(guideSeq, suffix=""):
    " The CSS class of guide row and links in seq viewer depend on the first nucl of guide "
    classNames = ["guideRow"]
    if guideSeq[0].upper()!="G":
        classNames.append("guideRowNoPrefixG"+suffix)
    if not guideSeq.startswith("GG"):
        classNames.append("guideRowNoPrefixGG"+suffix)
    if guideSeq[0].upper()!="A":
        classNames.append("guideRowNoPrefixA"+suffix)
    classStr = " ".join(classNames)
    return classStr

def buildCodonTable():
    " from http://www.petercollingridge.co.uk/tutorials/bioinformatics/codon-table/ "
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(list(zip(codons, amino_acids)))
    return codon_table

def buildOneToThree():
    " return one-letter -> three-letter conversion table for amino acids "
    oneToThree = \
        {'C':'Cys', 'D':'Asp', 'S':'Ser', 'Q':'Gln', 'K':'Lys',
         'I':'Ile', 'P':'Pro', 'T':'Thr', 'F':'Phe', 'N':'Asn',
         'G':'Gly', 'H':'His', 'L':'Leu', 'R':'Arg', 'W':'Trp',
         'A':'Ala', 'V':'Val', 'E':'Glu', 'Y':'Tyr', 'M':'Met',
         'U':'Sec', '*':'Stop',
         'X':'Stop',  # is this really used like that?
         'Z':'Glx', # special case: asparagine or aspartic acid
         'B':'Asx'  # special case: glutamine or glutamic acid
         }
    return oneToThree

def makeExonLines(exonInfo, seq, selTransId):
    """ create text that draws exons, input is transId -> (exonNumber, exStart, exEnd, exFrame).
    returns a list of (transId (=label), symbol (=mouseover), ASCII-line) """
    lines = []
    #maxLabelLen = 0
    codonTable = buildCodonTable()
    #oneToThree = buildOneToThree()
    seqLen = len(seq)
    seq = seq.upper()

    for (transId, symbol), exRows in exonInfo.items():
        if selTransId!="allTrans" and transId!=selTransId:
            continue
        line = [" "]*seqLen
        mouseOvers = {} # position -> mouseOver-text or null, for end of mouse over
        for exIdx, (exNum, exStart, exEnd, exFrame, nextFrame, exStrand) in enumerate(exRows):
            if exFrame==-1:
                for i in range(exStart, exEnd):
                    line[i]="="
                exonLabel = "noncoding"
                if (exEnd-exStart)>len(exonLabel)+4:
                    # center the exon label on the exon
                    mid = exStart+int((exEnd-exStart)*0.5)
                    halfLen = int(len(exonLabel)*0.5)
                    labStart = mid-halfLen
                    for i in range(0, len(exonLabel)):
                        line[labStart+i] = exonLabel[i]
                    line[labStart-1] = " "
                    line[labStart+len(exonLabel)] = " "
            else:
                exonDesc = "gene %s<br>transcript %s<br>exon %d<br>start phase %s" % (symbol, transId, exNum+1, exFrame)
                if nextFrame is not None:
                    exonDesc += "<br>end phase %s" % nextFrame
                    if (exFrame+nextFrame) % 3 == 0:
                        exonDesc += "<br>Removing the exon retains the reading frame"
                    else:
                        exonDesc += "<br>Removing the exon will destroy the reading frame"

                mouseOvers[exStart] = exonDesc
                mouseOvers[exEnd] = None

                for i in range(exStart, exStart+exFrame):
                        line[i] = "-"
                for i in range(exStart+exFrame, exEnd, 3):
                    codon = seq[i:i+3]
                    if len(codon)==3:
                        shortAa = codonTable[codon]
                        if exStrand=="+":
                            longAa = shortAa+"]]"
                        else:
                            longAa = "[["+shortAa # highlighting rev. dir. more
                    else:
                        # codon is split by splice site
                        longAa = "-"
                    for j in range(0, len(longAa)):
                        line[i+j] = longAa[j]

        # now merge the mouse overs as span tags into the ASCII line
        newLine = []
        for pos, char in enumerate(line):
            if pos in mouseOvers:
                overString = mouseOvers[pos]
                if overString is not None:
                    newLine.append("<span class='tooltipsterInteract' title='%s'>" % overString)
                if char=="-":
                    newLine.append("<span class='tooltipsterInteract' title='This codon goes over a splice site. The nucleotides of the split codon are not translated to amino acids but shown as dashes.'>-</span>")
                else:
                    newLine.append(char)
                if overString is None:
                    newLine.append("</span>")
            else:
                newLine.append(char)

        # and fix up the < signs
        newLineStr = "".join(newLine)
        newLineStr = newLineStr.replace("[[", "<lc>&lt;&lt;</lc>")
        newLineStr = newLineStr.replace("]]", "<lc>&gt;&gt;</lc>")

        lines.append((symbol, transId, newLineStr))
        #maxLabelLen = max(maxLabelLen, len(transId))
    return lines

def getGeneModels(org):
    " read possible gene models for org and return as list (name, desc) or None if no gene models "
    mask = join(genomesDir, org, "*.bb")
    fnames = glob.glob(mask)

    descFname = join(genomesDir, org, "genes.tsv")
    if not isfile(descFname):
        return None

    geneDescs = {}
    for line in open(descFname):
        fname, desc = line.split(maxsplit=1)
        geneDescs[fname] = desc

    ret = []
    for fname in fnames:
        baseName = basename(fname)
        name = baseName.split('.')[0]
        desc = geneDescs.get(baseName, name)
        ret.append((name, desc))
    return ret

def getSelGeneModel(org):
    " return (list of (name, desc) of models, selected gene model name) "
    geneModels = getGeneModels(org)
    selGeneModel = None
    selTransId = None

    if geneModels:
        #selGeneModel = cgiParams.get("geneModel", geneModels[0][0])
        selGeneModel = cgiParams.get("geneModel", "noGenes")
        geneModels.insert(0, ("noGenes", "Do not show"))
        possNames = [x for x,y in geneModels]
        if not selGeneModel in possNames:
            errAbort("The gene model name specified with the argument geneModel is invalid")

        selTransId = cgiParams.get("transId", "allTrans")

    return geneModels, selGeneModel, selTransId

def printSeqForCopy(seq):
    " print a hidden text area so we can copy the sequence to the clipboard "
    print('<input id="seqAsText" type="text" style="display:none">')
    print(seq)
    print("</input>")

def calcKomorScore(guideSeq, pos):
    " return base editing score given the guide sequence and the position "
    return pos/7.0 # temporary hack

def makeEditLines(seq, pamSeqs, winStart, winEnd, guideScores):
    " create the lines that show the possible baseEditor edits "
    editInfos = []
    for i in range(0, len(seq)):
        editInfos.append(defaultdict(list))

    upSeq = seq.upper()
    for pamId, pamStart, guideStart, strand, guideSeq, pamSeq, pamPlusSeq in pamSeqs:
        specScore = guideScores[pamId]
        if strand=="+":
            fromPos = guideStart+winStart
            toPos = guideStart+winEnd
            fromNucl = "C"
            toNucl = "T"
        else:
            guideEnd = guideStart+GUIDELEN
            fromPos = guideEnd-winEnd
            toPos = guideEnd-winStart
            fromNucl = "G"
            toNucl = "A"

        for pos in range(fromPos, toPos):
            # position of mutated nucl on forw strand guide
            if strand=="+":
                mutPos = pos-guideStart
            else:
                mutPos = GUIDELEN - (pos - guideStart) - 1

            if upSeq[pos]==fromNucl:
                beScore = calcKomorScore(guideSeq, mutPos)
                editInfos[pos][toNucl].append((pamId, guideSeq, pamSeq, mutPos, beScore, specScore))

    altNucls = ["A", "T"]

    editLabels = []
    for an in altNucls:
        editLabels.append("Edits to "+an)

    editLines = []
    for i in range(0, len(altNucls)):
        editLines.append([" "]*len(seq))

    # rearrange into lines of text + JSON
    jsonData = defaultdict(list)
    for pos, eiDict in enumerate(editInfos):
        if not eiDict:
            continue
        jsonData[pos] = eiDict
        for nucl, guideData in eiDict.items():
            yPos = altNucls.index(nucl)
            editLines[yPos][pos] = "<d pos=%d>%s</d>" % (pos, nucl)

    ret = []
    for label, lineChars in zip(editLabels, editLines):
        ret.append( (label, None, "".join(lineChars)) )

    return ret, jsonData

def makePamLines(lines, maxY, pamIdToSeq, guideScores):
    for y in range(0, maxY+1):
        texts = []
        lastEnd = 0
        for start, end, name, strand, pamId  in lines[y]:
            guideSeq = pamIdToSeq.get(pamId)
            if guideSeq==None:
                # when there is an N in the guide, the PAM is valid, but the guide is not
                continue
            classStr = cssClassesFromSeq(guideSeq, suffix="Seq")

            spacer = "".join([" "]*((start-lastEnd)))
            lastEnd = end
            texts.append(spacer)

            score = guideScores[pamId]
            # XX How can this happen for non-Cpf1 enzymes? Can this ever happen?
            if score is None and not pamIsFirst:
                continue
            color = scoreToColor(score)[0]

            texts.append('''<a class='%s' style="text-shadow: 1px 1px 1px #bbb; color: %s" id="list%s" href="#%s">''' % (classStr, color, pamId,pamId))
            texts.append(name)
            texts.append("</a>")
        yield ("", None, ''.join(texts))

def getBeWin(winVal):
    " return (start, end) of base editor window given CGI variable "
    fs = winVal.split("-")
    if len(fs)!=2:
        errAbort("parameter beWin must contain only one dash")
    start = fs[0].strip()
    end = fs[1].strip()
    if not start.isdigit() or not end.isdigit():
        errAbort("parameter beWin must be two dash-separated numbers")
    start = int(start)
    end = int(end)
    return start, end

def printLines(lines, labelLen):
    " print list of (label, string) such that label is at least labelLen characters long "
    for label, mouseOver, line in lines:
        if mouseOver is not None:
            labelStr =('<span class="tooltipsterInteract" title="{:s}">{:'+str(labelLen)+'s} </span>').format(label, mouseOver)
        else:
            labelStr = ('{:'+str(labelLen)+'s} ').format(label)

        print((labelStr), end=' ')
        print(line)

def getMaxLen(lines):
    " given a list of tuples where first element is the label, return the longest label len "
    maxLen = 0
    for l in lines:
        label = l[0]
        maxLen = max(maxLen, len(label))
    return maxLen

def printJson(name, obj):
    print("<script>")
    print((name), end=' ')
    print(("="), end=' ')
    print((json.dumps(obj)))
    print("</script>")

def showSeqAndPams(org, seq, startDict, pam, guideScores, varHtmls, varDbs, varDb, minFreq, position, pamIdToSeq):
    " show the sequence and the PAM sites underneath in a sequence viewer "
    pamSeqs = list(flankSeqIter(seq, startDict, len(pam), True))

    lines, maxY = distrOnLines(seq.upper(), startDict, len(pam))

    posLabel = "Position"
    varLabel = "Variants"
    seqLabel = "Sequence"
    exonLabelLen = 0
    editLines = []
    exonLines = []

    geneModels, selGeneModel, selTransId = getSelGeneModel(org)
    #selGeneModel = None
    #geneModels = None

    if baseEditor:
        beWinStart, beWinEnd = getBeWin(cgiParams.get("beWin", DEFAULTBEWIN))
        editLines, jsonData = makeEditLines(seq, pamSeqs, beWinStart, beWinEnd, guideScores)
        printJson("editData", jsonData)

    pamLines = list(makePamLines(lines, maxY, pamIdToSeq, guideScores))

    labelLen = max(len(varLabel), len(seqLabel), len(posLabel), getMaxLen(pamLines))

    if selGeneModel!=None:
        exonInfo, maxTransIdLen = getExonInfo(org, selGeneModel, position)
        labelLen = max(labelLen, maxTransIdLen)

    if baseEditor:
        labelLen = max(labelLen, getMaxLen(editLines))
    if selGeneModel:
        labelLen = max(labelLen, exonLabelLen)

    print("<div class='substep'>")
    print('<a id="seqStart"></a>')
    print("Your input sequence is %d bp long. It contains %d possible guide sequences.<br>" % (len(seq), len(guideScores)))

    if not pamIsFirst:
        print("Shown below are their PAM sites and the expected cleavage position located -3bp 5' of the PAM site.<br>")
        print("Click on a match for the PAM %s below to show its %d bp-long guide sequence. " % (pam, GUIDELEN))
        print("(Need help? Look at the <a target=_blank href='manual/#annotseq'>CRISPOR manual</a>)<br>")
        print('''Colors <span style="color:#32cd32; text-shadow: 1px 1px 1px #bbb">green</span>, <span style="color:#ffff00; text-shadow: 1px 1px 1px #888">yellow</span> and <span style="text-shadow: 1px 1px 1px #f01; color:#aa0014">red</span> indicate high, medium and low specificity of the PAM's guide sequence in the genome.<p>''')
    else:
        print("Click on a match for the PAM %s below to show its %d bp-long guide sequence.<br>" % (pam, GUIDELEN))

    if baseEditor or varDb or selGeneModel:
        print(("""<form style="display:inline" id="paramForm" action="%s" method="GET">""" % basename(__file__)))

    if geneModels:
        print ("Gene Models:")
        printDropDown("geneModel", geneModels, selGeneModel, style="width:20em")
        if selGeneModel!="noGenes":
            print ("Transcript:")
            transIdInfo = [("allTrans", "All Transcripts")]
            for transId, sym in list(exonInfo.keys()):
                transIdInfo.append( (transId, sym+" / "+transId) )

            printDropDown("transId", transIdInfo, selTransId, style="width:20em")
            # XX XXXXXX
            exonLines = makeExonLines(exonInfo, seq, selTransId)
            #exonLines = []
        print("""<input style="height:18px;margin:0px;font-size:10px;line-height:normal" type="submit" name="submit" value="Update">""")
        print("""<br>""")

    if baseEditor:
        print ("Base Editor modification window:")
        print(("""<input type="text" name="beWin" size="10" value="%s">""" % DEFAULTBEWIN))
        print("""<input style="height:18px;margin:0px;font-size:10px;line-height:normal" type="submit" name="submit" value="Update">""")
        print("<br>")

    if varDb is not None:
        print ("Variant database:")
        varDbList = [(b,c) for a,b,c,d in varDbs] # only keep fname+label
        printDropDown("varDb", varDbList, varDb)

        if minFreq==0.0:
            minFreq="0.0"
        else:
            minFreq = str(minFreq)

        # pull out the hasAF field for this varDb
        varDbHasAF = False
        for shortLabel, fname, desc, hasAF in varDbs:
            if fname==varDb:
                varDbHasAF = hasAF
                break

        if varDbHasAF:
            print("""&nbsp; Min. frequency: """)
            print(("""<input type="text" name="minFreq" size="8" value="%s">""" % minFreq))
        print("""<input style="height:18px;margin:0px;font-size:10px;line-height:normal" type="submit" name="submit" value="Update">""")
        print(("<small style='margin-left:30px'><a href='mailto:%s'>Missing a variant database? We can add it.</a></small>" % contactEmail))

    if position=="?":
        print("<small style=''>Input sequence not in genome, cannot show genome variants.</small>")
    elif varDb is None:
        print(("<small style=''><a href='mailto:%s'>Suggest a genome variants database to show on this page</a></small>" % contactEmail))

    print("</div>")
    print('''<div class="blueHighlight" style="text-align: left; overflow-x:scroll; width:98vw; background:#DDDDDD; border-style: solid; border-width: 1px">''')

    print('''<pre style="font-family: Source Code Pro; font-size: 80%; display:inline; line-height: 0.95em; text-align:left">''')
    print(('{:'+str(labelLen)+'s} ').format(posLabel), end=' ')
    print(rulerString(len(seq)))

    if varHtmls is not None:
        print(('{:'+str(labelLen)+'s} ').format(varLabel), end=' ')
        print("".join(varHtmls))

    print(('{:'+str(labelLen)+'s} ').format(seqLabel), end=' ')
    print (seq)

    printLines(exonLines, labelLen)

    if baseEditor:
        printLines(editLines, labelLen)

    printLines(pamLines, labelLen)


    print("</pre><br>")

    print('''</div>''')

    #printSeqForCopy(seq)

    if pamIsCpf1(pam):
        print('<div style="line-height: 1.0; padding-top: 5px; font-size: 15px">Cpf1 has a staggered site: cleavage occurs usually - but not always - after the 18th base on the non-targeted strand which has the TTTV PAM motif (indicate by "\\" in the schema above). Cleavage mostly occurs after the 23rd base on the targeted strand which has the AAAN motif (indicated by "/" in the schema above). See <a target=_blank href="http://www.sciencedirect.com/science/article/pii/S0092867415012003">Zetsche et al 2015</a>, in particular <a target=_blank href="http://www.sciencedirect.com/science?_ob=MiamiCaptionURL&_method=retrieve&_eid=1-s2.0-S0092867415012003&_image=1-s2.0-S0092867415012003-gr3.jpg&_cid=272196&_explode=defaultEXP_LIST&_idxType=defaultREF_WORK_INDEX_TYPE&_alpha=defaultALPHA&_ba=&_rdoc=1&_fmt=FULL&_issn=00928674&_pii=S0092867415012003&md5=11771263f3e390e444320cacbcfae323">Fig 3</a>.</div>')
    elif pamIsCasX(pam):
        print('<div style="line-height: 1.0; padding-top: 5px; font-size: 15px">We have no description yet on how exactly the CasX cleavage looks like. Please contact crispor@tefor.net if you have an idea how to describe the cleavage site.</div>')

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

def flankSeqIter(seq, startDict, pamLen, doFilterNs):
    """ given a seq and dictionary of pamPos -> strand and the length of the pamSite
    yield tuples of (name, pamStart, guideStart, strand, flankSeq, pamSeq)

    flankSeq is the guide sequence (=flanking the PAM).

    if doFilterNs is set, will not return any sequences that contain an N character
    pamPlusSeq are the 5bp after the PAM. If not enough space, pamPlusSeq is None
    """

    startList = sorted(startDict.keys())
    for pamStart in startList:
        strand = startDict[pamStart]

        pamPlusSeq = None
        if pamIsFirst: # Cpf1: get the sequence to the right of the PAM
            if strand=="+":
                guideStart = pamStart+pamLen
                flankSeq = seq[guideStart:guideStart+GUIDELEN]
                pamSeq = seq[pamStart:pamStart+pamLen]
                if pamStart-pamPlusLen >= 0:
                    pamPlusSeq = seq[pamStart-pamPlusLen:pamStart]
            else: # strand is minus
                guideStart = pamStart-GUIDELEN
                flankSeq = revComp(seq[guideStart:pamStart])
                pamSeq = revComp(seq[pamStart:pamStart+pamLen])
                if pamStart+pamLen+pamPlusLen < len(seq):
                    pamPlusSeq = revComp(seq[pamStart+pamLen:pamStart+pamLen+pamPlusLen])
        else: # common case: get the sequence on the left side of the PAM
            if strand=="+":
                guideStart = pamStart-GUIDELEN
                flankSeq = seq[guideStart:pamStart]
                pamSeq = seq[pamStart:pamStart+pamLen]
                if pamStart+pamLen+pamPlusLen < len(seq):
                    pamPlusSeq = seq[pamStart+pamLen:pamStart+pamLen+pamPlusLen]
            else: # strand is minus
                guideStart = pamStart+pamLen
                flankSeq = revComp(seq[guideStart:guideStart+GUIDELEN])
                pamSeq = revComp(seq[pamStart:pamStart+pamLen])
                if pamStart-pamPlusLen >= 0:
                    pamPlusSeq = revComp(seq[pamStart-pamPlusLen:pamStart])

        if "N" in flankSeq and doFilterNs:
            continue

        yield "s%d%s" % (pamStart, strand), pamStart, guideStart, strand, flankSeq, pamSeq, pamPlusSeq

def makeBrowserLink(dbInfo, pos, text, title, cssClasses, ctUrl=None):
    " return link to genome browser (ucsc or ensembl) at pos, with given text "
    if dbInfo is None:
        errAbort("Your batchID relates to a genome that is not present anymore. You will have to change the version of the site. Or contact us and send us the full URL of this page.")

    if dbInfo.server.startswith("Ensembl"):
        baseUrl = "www.ensembl.org"
        urlLabel = "Ensembl"

        # link back to archive, if possible
        if dbInfo.description.startswith("Ensembl "):
            ensVersion = dbInfo.description.split()[1]
            if ensVersion.isdigit():
                baseUrl = "e%s.ensembl.org" % ensVersion

        elif dbInfo.server=="EnsemblPlants":
            baseUrl = "plants.ensembl.org"
        elif dbInfo.server=="EnsemblMetazoa":
            baseUrl = "metazoa.ensembl.org"
        elif dbInfo.server=="EnsemblProtists":
            baseUrl = "protists.ensembl.org"
        org = dbInfo.scientificName.replace(" ", "_")
        pos = pos.replace(":+","").replace(":-","") # remove the strand
        url = "http://%s/%s/Location/View?r=%s" % (baseUrl, org, pos)
    elif dbInfo.server=="ucsc" or dbInfo.name.startswith("GCA_") or dbInfo.name.startswith("GCF_"):
        urlLabel = "UCSC"
        if len(pos)>0 and pos[0].isdigit():
            pos = "chr"+pos
        # remove the strand
        pos = pos.replace(":+","").replace(":-","")
        url = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s" % (dbInfo.name, pos)
        if ctUrl is not None:
            url+= "&hgt.customText=%s" % ctUrl
    # some limited support for gbrowse
    elif dbInfo.server.startswith("http://"):
        urlLabel = "GBrowse"
        chrom, start, end, strand = parsePos(pos)
        start = start+1
        url = "%s/?name=%s:%d..%d" % (dbInfo.server, chrom, start, end)
    else:
        chrom, start, end, strand = parsePos(pos)
        if chrom is not None and chrom.startswith("NC_"):
            start = start+1
            url = "https://www.ncbi.nlm.nih.gov/nuccore/%s?report=graph&log$=seqview&v=%d-%d" % \
            (chrom, start, end)
            urlLabel = "NCBI "

        else:
            #return "unknown genome browser server %s, please email services@tefor.net" % dbInfo.server
            urlLabel = None
            url = "javascript:void(0)"

    classStr = ""
    if len(cssClasses)!=0:
        classStr = ' class="%s"' % (" ".join(cssClasses))

    if title is None:
        if urlLabel != None:
            title = "Link to %s Genome Browser" % urlLabel
        else:
            title = "No Genome Browser link available yet for this organism"
    return '''<a title="%s"%s target="_blank" href="%s">%s</a>''' % (title, classStr, url, text)

def highlightMismatches(guide, offTarget, pamLen):
    " return a string that marks mismatches between guide and offtarget with * "
    if pamLen!=0:
        if pamIsFirst:
            offTarget = offTarget[pamLen:]
        else:
            offTarget = offTarget[:-pamLen]
    assert(len(guide)==len(offTarget))

    s = []
    for x, y in zip(guide, offTarget):
        if x==y:
            s.append(".")
        else:
            s.append("*")
    return "".join(s)

def parseNewAlias(ifh):
    " part of parseAlias(): IGV-compatible format: first is UCSC, all other columns are aliases "
    toUcsc = {}
    for line in ifh:
        if line.startswith("#"):
            continue
        row = line.rstrip("\n").split("\t")
        for i in range(1, len(row)):
            toUcsc[row[i]] = row[0]
    return toUcsc

def parseAlias(fname):
    " parse tsv file with at least two columns, orig chrom name and new chrom name. copied from chromToUcsc script from the UCSC tools. "
    logging.debug("alias file is in IGV-format")
    toUcsc = {}
    if fname.startswith("http://") or fname.startswith("https://"):
        ifh = urlopen(fname)
        if fname.endswith(".gz"):
            data = gzip.GzipFile(fileobj=ifh).read().decode()
            ifh = data.splitlines()
    elif fname.endswith(".gz"):
        ifh = gzip.open(fname, "rt")
    else:
        ifh = open(fname)

    firstLine = True
    for line in ifh:
        if line.startswith("#") and firstLine:
            return parseNewAlias(ifh)
        if line.startswith("alias"):
            continue
        row = line.rstrip("\n").split("\t")
        toUcsc[row[0]] = row[1]
        firstLine = False
    return toUcsc

chromAlias = None

def applyChromAlias(db, chrom):
    " if chrom is in chromAlias, return the human-readable name "
    global chromAlias
    if chromAlias==-1: # == chromAlias file not present
        return chrom
    elif chromAlias is None:
        chromAliasFname = join("genomes", db, db+".chromAlias.txt")
        if not isfile(chromAliasFname):
            chromAlias = -1
            return chrom
        else:
            chromAlias = parseAlias(chromAliasFname)

    return chromAlias.get(chrom, chrom)

def makeAlnStr(org, seq1, seq2, pam, mitScore, cfdScore, posStr, chromDist):
    """ given two strings of equal length, return a html-formatted string of several lines
    that show the two sequences and a line that highlights where they differ
    """
    lines = [ [], [], [] ]
    last12MmCount = 0
    inLinkage = False
    hlSeed = False

    if pamIsSpCas9(pam):
        hlSeed = True

    if pamIsFirst:
        lines[0].append("<i>"+seq1[:len(pam)]+"</i> ")
        lines[1].append("<i>"+seq2[:len(pam)]+"</i> ")
        lines[2].append("".join([" "]*(len(pam)+1)))

    if pamIsFirst:
        guideStart = len(pam)
        guideEnd = len(seq1)
    else:
        guideStart = 0
        guideEnd = len(seq1)-len(pam)

    for i in range(guideStart, guideEnd):
        if hlSeed and i==10:
            lines[1].append("<u>")

        if seq1[i]==seq2[i]:
            lines[0].append(seq1[i])
            lines[1].append(seq2[i])
            lines[2].append(" ")
        else:
            lines[0].append("<b>%s</b>" % seq1[i])

            lines[1].append("<b>%s</b>" % seq2[i])

            lines[2].append("*")
            if i>7:
                last12MmCount += 1

        if hlSeed and i==guideEnd-1:
            lines[1].append("</u>")

    if not pamIsFirst:
        lines[0].append(" <i>"+seq1[-len(pam):]+"</i>")
        lines[1].append(" <i>"+seq2[-len(pam):]+"</i>")
    lines = ["".join(l) for l in lines]

    chrom, chromPos, strand = posStr.split(":")
    chrom = applyChromAlias(org, chrom)
    posStr = ":".join((chrom, chromPos, strand))

    if len(posStr)>1 and posStr[0].isdigit():
        posStr = "chr"+posStr

    htmlText1 = "<small><pre>guide:      %s<br>off-target: %s<br>            %s</pre>" \
        % (lines[0], lines[1], lines[2])

    if pamIsCpf1(pam) or pamIsCasX(pam):
        htmlText2 = "Cpf1/CasX: No off-target scores available</small>"
    elif saCas9Mode:
        htmlText2 = "SaCas9 Tycko Score: %s" % mitScore
    else:
        if cfdScore==None:
            cfdStr = "Cannot calculate CFD score on non-ACTG characters"
        else:
            cfdStr = "%f" % cfdScore

        htmlText2 = "CFD Off-target score: %s<br>MIT Off-target score: %.2f<br>Position: %s</small>" % (cfdStr, mitScore, posStr)
        if chromDist!=None and org!=None:
            htmlText2 += "<br><small>Distance from target: %.3f Mbp</small>" % (float(chromDist)/1000000.0)
            if org.startswith("mm") or org.startswith("hg") or org.startswith("rn"):
                if chromDist > 20000000:
                    htmlText2 += "<br><small>&gt;20Mbp = unlikely to be in linkage with target</small>"
                else:
                    htmlText2 += "<br><small>&lt;20Mbp= likely to be in linkage with "
                    "target! Even if no linkage: beware of chromosomal rearrangements "
                    "when using this guide!</small>"
                    inLinkage = True

    hasLast12Mm = last12MmCount>0
    return htmlText1+htmlText2, hasLast12Mm, inLinkage

def parsePos(text):
    """ parse a string of format chr:start-end:strand and return a 4-tuple
    Strand defaults to + and end defaults to start+23
    """
    if text!=None and len(text)!=0 and text!="?":
        fields = text.split(":")
        if len(fields)==2:
            chrom, posRange = fields
            strand = "+"
        else:
            chrom, posRange, strand = fields
        posRange = posRange.replace(",","")
        if "-" in posRange:
            start, end = posRange.split("-")
            start, end = int(start), int(end)
            if start > end:
                start, end = end, start
                strand = "-"
        else:
            # if the end position is not specified (as by default done by UCSC outlinks), use start+23
            start = int(posRange)
            end = start+23
    else:
        chrom, start, end, strand = "", 0, 0, "+"
    return chrom, start, end, strand

def annotateOfftargets(org, countDict, guideSeq, pam, inputPos):
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
    mitOtScores = []
    cfdScores = []
    last12MmCounts = []
    ontargetDesc = ""
    repCount = 0 # if repCount for a guide is !=0, then the guide should not be used. repCount is then the number
    # of matches for the guide in the genome (not looking at the PAM)

    # for each edit distance, get the off targets and iterate over them
    foundOneOntarget = False
    isSaCas9 = pamIsSaCas9(pam)
    isCpf1 = pamIsCpf1(pam)

    for editDist in range(0, maxMMs+1):
        #print countDict,"<p>"
        matches = countDict.get(editDist, [])

        #print otCounts,"<p>"
        last12MmOtCount = 0

        # create html and score for every offtarget
        otCount = 0
        for chrom, start, end, otSeq, strand, segType, geneNameStr, totalAlnCount, isRep in matches:
            # if repCount is > 0, then this means that the guide should not be used and we cannot
            # even get any off-targets
            if (totalAlnCount > MAXOCC) or (totalAlnCount > 1 and isRep):
                repCount = totalAlnCount # any off-target with this condition will trigger the whole guide to be suppressed

            # skip on-targets
            if segType!="":
                segTypeDesc = segTypeConv[segType]
                geneDesc = segTypeDesc+":"+geneNameStr
                geneDesc = geneDesc.replace("|", "-")
            else:
                geneDesc = geneNameStr

            # is this not an off-target but the on-target?
            # if we got a genome position, use it. Otherwise use a random off-target with 0MMs
            # as the on-target ("auto-ontarget" mode)
            if editDist==0 and \
                repCount==0 and \
                ((chrom==inChrom and start >= inStart and end <= inEnd) \
                or (inChrom=='' and foundOneOntarget==False)):
                foundOneOntarget = True
                ontargetDesc = geneDesc
                continue

            otCount += 1
            guideNoPam = guideSeq[:len(guideSeq)-len(pam)]
            otSeqNoPam = otSeq[:len(otSeq)-len(pam)]

            if len(otSeqNoPam)==19:
                otSeqNoPam = "A"+otSeqNoPam # should not change the score a lot, weight0 is very low
                guideNoPam = "A"+guideNoPam

            if isCpf1:
                # Cpf1 has no off-target scores yet
                mitScore=0.0
                cfdScore=0.0
            elif isSaCas9:
                mitScore = calcSaHitScore(guideNoPam, otSeqNoPam)
                cfdScore = -1

            else:
                # MIT score must not include the PAM
                mitScore = calcHitScore(guideNoPam, otSeqNoPam)
                # this is a heuristic based on the guideSeq data where alternative
                # PAMs represent only ~10% of all cleaveage events.
                # We divide the MIT score by 5 to make sure that these off-targets
                # are not ranked among the top but still appear in the list somewhat
                if pam=="NGG" and otSeq[-2:]!="GG":
                    mitScore = mitScore * 0.2

                # CFD score must include the PAM
                cfdScore = calcCfdScore(guideSeq, otSeq)

            mitOtScores.append(mitScore)
            if cfdScore != -1:
                cfdScores.append(cfdScore)

            posStr = "%s:%d-%s:%s" % (chrom, start+1,end, strand)
            if (chrom==inChrom):
                dist = abs(start-inStart)
            else:
                dist = None

            parNum = isInPar(org, chrom, start, end)
            if parNum is not None:
                posStr += " PAR%s" % parNum

            alnHtml, hasLast12Mm, inLinkage = makeAlnStr(org, guideSeq, otSeq,
                pam, mitScore, cfdScore, posStr, dist)
            if not hasLast12Mm:
                last12MmOtCount+=1
            posList.append( (otSeq, mitScore, cfdScore, editDist, posStr, geneDesc,
                alnHtml, inLinkage) )

        last12MmCounts.append(str(last12MmOtCount))
        # create a list of number of offtargets for this edit dist
        otCounts.append( str(otCount) )

    # calculate the guide scores
    if pamIsCpf1(pam):
        guideScore = -1
        guideCfdScore = -1
    else:
        if repCount>0:
            guideScore = 0
            guideCfdScore = 0
        else:
            guideScore = calcMitGuideScore(sum(mitOtScores))

            if doCfdFix:
                guideCfdScore = calcCfdGuideScore(sum(cfdScores))
            else:
                guideCfdScore = calcMitGuideScore(sum(cfdScores))

    # obtain the off-target info: coordinates, descriptions and off-target counts
    if repCount>0:
        posList = []
        ontargetDesc = ""
        last12DescStr = ""
        otDescStr = ""
    else:
        otDescStr = "&thinsp;-&thinsp;".join(otCounts)
        last12DescStr = "&thinsp;-&thinsp;".join(last12MmCounts)

    if pamIsCpf1(pam):
        # sort by edit dist if using Cfp1
        posList.sort(key=operator.itemgetter(3))
    else:
        # sort by CFD score if we have it
        posList.sort(reverse=True, key=operator.itemgetter(2))

    return posList, otDescStr, guideScore, guideCfdScore, last12DescStr, \
        ontargetDesc, repCount

# --- START OF SCORING ROUTINES

saGuide = None
saScorer = None
def calcSaHitScore(guideSeq, otSeq):
    """
    saCas9 offtarget scoring from Tycko et al, https://www.nature.com/articles/s41467-018-05391-2
    see bin/src/pairwise-library-screen/
    """
    global saScorer
    global saGuide
    if guideSeq!=saGuide:
        sys.path.append("bin/src/pairwise-library-screen")
        import predictSingle
        saGuide = guideSeq
        saScorer = predictSingle.SaCas9Scorer(len(guideSeq))

    # to be compatible with the MIT score, has to be in the range 0-100
    # for the MIT aggregate guide specificity score
    return 100.0*saScorer.calcScore(guideSeq, otSeq)

# MIT offtarget scoring, "Hsu score"

# aka Matrix "M"
hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]

def calcHitScore(string1,string2):
    " see 'Scores of single hits' on http://crispr.mit.edu/about "
    # The Patrick Hsu weighting scheme
    # S. aureus requires 21bp long guides. We fudge by using only last 20bp
    matrixStart = 0
    maxDist = 19

    assert(string1[0].isupper())
    assert(len(string1)==len(string2))
    #for nmCas9 and a few others with longer guides, we limit ourselves to 20bp
    if len(string1)>20:
        string1 = string1[-20:]
        string2 = string2[-20:]
    # for 19bp guides, we fudge a little, but first pos has no weight anyways
    elif len(string1)==19:
        string1 = "A"+string1
        string2 = "A"+string2
    # for shorter guides, I'm not sure if this score makes sense anymore, we force things
    elif len(string1)<19:
        matrixStart = 20-len(string1)
        maxDist = len(string1)-1

    assert(len(string1)==len(string2))

    dists = [] # distances between mismatches, for part 2
    mmCount = 0 # number of mismatches, for part 3
    lastMmPos = None # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(matrixStart, len(string1)):
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
        score2 = 1.0 / (((maxDist-avgDist)/float(maxDist)) * 4 + 1)
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

def calcCfdGuideScore(hitSum):
    " suggested by Nicholas Parkinson "
    norm_score = 100.* 100.0 / (hitSum)
    return norm_score

# === SOURCE CODE cfd-score-calculator.py provided by John Doench =====
# The CFD score is an improved specificity score

def get_mm_pam_scores():
    """
    """
    import pickle
    dataDir = join(dirname(__file__), 'CFD_Scoring')
    mm_scores = pickle.load(open(join(dataDir, 'mismatch_score.pkl'),'rb'))
    pam_scores = pickle.load(open(join(dataDir, 'pam_scores.pkl'),'rb'))
    return (mm_scores,pam_scores)

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

#Calculates CFD score
def calc_cfd(wt,sg,pam):
    #mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score)

mm_scores, pam_scores = None, None

def calcCfdScore(guideSeq, otSeq):
    """ based on source code provided by John Doench
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGAAAGGG")
    0.4635989007074176
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGGGGG")
    1.0
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "aaaaGaGaGGGGGGGGGGGGGGG")
    0.5140384614450001

    # mismatches:      *               !!
    >>> calcCfdScore("ATGGTCGGACTCCCTGCCAGAGG", "ATGGTGGGACTCCCTGCCAGAGG")
    0.5

    # mismatches:    *  ** *
    >>> calcCfdScore("ATGGTCGGACTCCCTGCCAGAGG", "ATGATCCAAATCCCTGCCAGAGG")
    0.53625000020625

    >>> calcCfdScore("ATGTGGAGATTGCCACCTACCGG", "ATCTGGAGATTGCCACCTACAGG")
    0.384615385

    """
    global mm_scores, pam_scores
    if mm_scores is None:
        mm_scores,pam_scores = get_mm_pam_scores()
    wt = guideSeq.upper()
    off = otSeq.upper()
    m_wt = re.search('[^ATCG]',wt)
    m_off = re.search('[^ATCG]',off)
    if (m_wt is None) and (m_off is None):
        pam = off[-2:]
        sg = off[:20]
        cfd_score = calc_cfd(wt,sg,pam)
        if doCfdFix:
            cfd_score = cfd_score*100.0
        return cfd_score

    return -1
# ==== END CFD score source provided by John Doench

# --- END OF SCORING ROUTINES

def getSizeFname(genome):
    " return name of chrom.sizes file "
    genomeDir = genomesDir # make local
    sizeFname = "%(genomeDir)s/%(genome)s/%(genome)s.sizes" % locals()
    return sizeFname

def parseChromSizes(genome):
    " return chrom sizes as dict chrom -> size "
    sizeFname = getSizeFname(genome)
    ret = {}
    for line in open(sizeFname).read().splitlines():
        fields = line.split()
        chrom, size = fields[:2]
        ret[chrom] = int(size)
    return ret

def extendAndGetSeq(db, chrom, start, end, strand, oldSeq, flank=FLANKLEN):
    """ extend (start, end) by flank and get sequence for it using twoBitTwoFa.
    Return None if not possible to extend.
    #>>> extendAndGetSeq("hg19", "chr21", 10000000, 10000005, "+", flank=3)
    #'AAGGAATGTAG'
    """
    assert("|" not in chrom) # we are using | to split info in BED files. | is not allowed in the fasta
    chromSizes = parseChromSizes(db)
    maxEnd = chromSizes[chrom]+1

    start -= flank
    end += flank
    if start < 0 or end > maxEnd:
        return None

    genomeDir = genomesDir
    twoBitFname = "%(genomeDir)s/%(db)s/%(db)s.2bit" % locals()
    progDir = binDir
    genome = db
    cmd = "%(progDir)s/twoBitToFa %(genomeDir)s/%(genome)s/%(genome)s.2bit stdout -seq='%(chrom)s' -start=%(start)s -end=%(end)s" % locals()
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, encoding="utf8")
    seqStr = proc.stdout.read()
    proc.wait()
    if proc.returncode!=0:
        errAbort("Could not run '%s'. Return code %s" % (cmd, str(proc.returncode)))
    faFile = StringIO(seqStr)
    seqs = parseFasta(faFile)
    assert(len(seqs)==1)
    seq = list(seqs.values())[0].upper()


    if strand=="-":
        seq = revComp(seq)

    genomeSeq = seq[FLANKLEN:(FLANKLEN+len(oldSeq))].upper()
    if oldSeq.upper() not in genomeSeq:
        logging.warn("Input sequence has SNPs compared to genome, not returning extended seq:")
        logging.warn("- Input sequence:  %s" % oldSeq)
        logging.warn("- Genome sequence: %s" % genomeSeq)
        logging.warn("- Diff String    : %s" % highlightMismatches(oldSeq, genomeSeq, 0))
        return None
    # ? make sure that user annotations, like added Ns, are retained in the long sequence
    #fixedSeq = seq[:100]+oldSeq+seq[-100:]
    #assert(len(fixedSeq)==len(seq))
    return seq

def getExtSeq(seq, start, end, strand, extUpstream, extDownstream, extSeq=None, extFlank=FLANKLEN):
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

    # extend
    if strand=="+":
        extStart, extEnd = start-extUpstream, end+extDownstream
    else:
        extStart, extEnd = start-extDownstream, end+extUpstream

    # check for out of bounds and get seq
    if extStart >= 0 and extEnd <= len(seq):
        logging.debug("using input seq, pos %d-%d" % (extStart, extEnd))
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
        logging.debug("using extended seq, pos %d-%d" % (extStart, extEnd))

    if strand=="-":
        logging.debug("revcomp'ing result")
        subSeq = revComp(subSeq)

    # check that the extended sequence really contains the whole input seq
    # e.g. when user has added nucleotides to a otherwise matching sequence
    #if seq.upper() not in subSeq.upper():
        #debug("seq is not in extSeq")
        #subSeq = None

    logging.debug("Got -%d/+%d-extended seq for (%d, %d, %s) = %s. Result: %s." %
        (extUpstream, extDownstream, start, end, strand, seq[start:end], subSeq))
    return subSeq

def pamStartToGuideRange(startPos, strand, pamLen):
    """ given a PAM start position and its strand, return the (start,end) of the guide.
    Coords can be negative or exceed the length of the input sequence.
    """
    if not pamIsFirst:
        if strand=="+":
            return (startPos-GUIDELEN, startPos)
        else: # strand is minus
            return (startPos+pamLen, startPos+pamLen+GUIDELEN)
    else:
        if strand=="+":
            return (startPos+pamLen, startPos+pamLen+GUIDELEN)
        else: # strand is minus
            return (startPos-GUIDELEN, startPos)

def htmlHelp(text):
    " show help text with tooltip or modal dialog "
    className = "tooltipster"
    if "href" in text:
        className = "tooltipsterInteract"

    print('''<img style="padding-bottom: 3px; height:1.1em; width:1.0em" src="%simage/info-small.png" class="help %s" title="%s" />''' % (HTMLPREFIX, className, text))

def htmlWarn(text):
    " show help text with tooltip "
    print('''<img style="height:0.9em; width:0.8em; padding-bottom: 2px" src="%simage/warning-32.png" class="help tooltipster" title="%s" />''' % (HTMLPREFIX, text))

def readRestrEnzymes():
    """ parse restrSites.txt and
    return as dict length -> list of (name, suppliers, seq) """
    fname = "restrSites.txt"
    enzList = {}
    for line in open(join(baseDir, fname)):
        if line.startswith("#"):
            continue
        seq, name, suppliers = line.rstrip("\n").rstrip("\r").split("\t")
        suppliers = tuple(suppliers.split(","))
        enzList.setdefault(len(seq), []).append( (name, suppliers, seq) )
    return enzList

def patMatch(seq, pat, notDegPos=None):
    """ return true if pat matches seq, both have to be same length
    do not match degenerate codes at position notDegPos (0-based)
    """
    assert(len(seq)==len(pat))
    for x in range(0, len(pat)):
        patChar = pat[x]
        nuc = seq[x]

        assert(patChar in "MKYRACTGNWSDVBH")
        assert(nuc in "MKYRACTGNWSDX")

        if notDegPos!=None and x==notDegPos and patChar!=nuc:
            return False

        if nuc=="X":
            return False

        if patChar=="N":
            continue

        if patChar=="H" and nuc in "ACT":
            continue
        if patChar=="D" and nuc in "AGT":
            continue
        if patChar=="B" and nuc in "CGT":
            continue
        if patChar=="V" and nuc in "ACG":
            continue

        if patChar=="W" and nuc in "AT":
            continue
        if patChar=="S" and nuc in "GC":
            continue
        if patChar=="M" and nuc in "AC":
            continue
        if patChar=="K" and nuc in "TG":
            continue
        if patChar=="R" and nuc in "AG":
            continue
        if patChar=="Y" and nuc in "CT":
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

        # JP does not want any potential site to be suppressed
        #if i<len(restrSite):
            #isMatch = patMatch(subseq, restrSite, len(restrSite)-i-1)
        #else:
            #isMatch = patMatch(subseq, restrSite)
        isMatch = patMatch(subseq, restrSite)

        if isMatch:
            posList.append( (i, i+len(restrSite)) )
    return posList

def matchRestrEnz(allEnzymes, guideSeq, pamSeq, pamPlusSeq, pamPat):
    """ return list of enzymes that overlap the -3 position in guideSeq
    returns dict (name, pattern, suppliers) -> list of matching positions
    """
    matches = defaultdict(set)

    if pamPlusSeq is None:
        pamPlusSeq = "XXXXX" # make sure that we never match a restriction site outside the seq boundaries

    fullSeq = concatGuideAndPam(guideSeq, pamSeq, pamPlusSeq)
    #print guideSeq, pamSeq, pamPlusSeq, fullSeq, "<br>"

    for siteLen, sites in allEnzymes.items():
        if pamIsCpf1(pamPat):
            # most modified position: 4nt from the end
            # see http://www.nature.com/nbt/journal/v34/n8/full/nbt.3620.html
            # Figure 1
            startSeq = len(fullSeq)-4-pamPlusLen-(siteLen)+1
        else:
            # most modified position for Cas9: 3bp from the end
            startSeq = len(fullSeq)-len(pamSeq)-3-pamPlusLen-(siteLen)+1

        seq = fullSeq[startSeq:].upper()
        for name, suppliers, restrSite in sites:
            posList = findSite(seq, restrSite)
            if len(posList)!=0:
                liftOffset = startSeq
                posList = [(liftOffset+x, liftOffset+y) for x,y in posList]
                matches.setdefault((name, restrSite, suppliers), set()).update(posList)
    return matches

def mergeGuideInfo(seq, startDict, pamPat, otMatches, inputPos, effScores, sortBy=None, org=None):
    """
    merges guide information from the sequence, the efficiency scores and the off-targets.
    creates rows with too many fields. needs refactoring.

    for each pam in startDict, retrieve the guide sequence next to it and score it
    sortBy can be "effScore", "mhScore", "oofScore" or "pos"
    """
    allEnzymes = readRestrEnzymes()

    guideData = []
    guideScores = {}
    hasNotFound = False
    pamIdToSeq = {}

    pamSeqs = list(flankSeqIter(seq.upper(), startDict, len(pamPat), True))

    for pamId, pamStart, guideStart, strand, guideSeq, pamSeq, pamPlusSeq in pamSeqs:
        # matches in genome
        # one desc in last column per OT seq
        if pamId in otMatches:
            pamMatches = otMatches[pamId]
            guideSeqFull = concatGuideAndPam(guideSeq, pamSeq)
            mutEnzymes = matchRestrEnz(allEnzymes, guideSeq, pamSeq, pamPlusSeq, pamPat)
            posList, otDesc, guideScore, guideCfdScore, last12Desc, ontargetDesc, \
               repCount = \
                   annotateOfftargets(org, pamMatches, guideSeqFull, pamPat, inputPos)
            if repCount!=0:
                guideScore = 0
                guideCfdScore = 0

        # no off-targets found?
        else:
            posList, otDesc, guideScore = None, "Not found", -1
            guideCfdScore = -1
            last12Desc = ""
            hasNotFound = True
            mutEnzymes = []
            ontargetDesc = ""
            repCount = 0
            seq34Mer = None

        guideRow = [guideScore, guideCfdScore, effScores.get(pamId, {}), pamStart, guideStart, strand, pamId, guideSeq, pamSeq, posList, otDesc, last12Desc, mutEnzymes, ontargetDesc, repCount]
        guideData.append( guideRow )
        guideScores[pamId] = guideScore
        pamIdToSeq[pamId] = guideSeq

    if sortBy == "pos":
        sortFunc = (lambda row: row[3])
        reverse = False
    elif sortBy == "offCount":
        sortFunc = (lambda row: len(row[9]))
        reverse = False
    elif sortBy == "cfdSpec":
        sortFunc = operator.itemgetter(1)
        reverse = True
    elif sortBy == "spec" or sortBy is None:
        sortFunc = (lambda row: row[0])
        reverse = True
    elif sortBy is not None and not sortBy.endswith("pec"):
        sortFunc = (lambda row: row[2].get(sortBy, 0))
        reverse = True
    else:
        errAbort("Unknown sortBy value. This is a bug. Please contact us.")

    guideData.sort(reverse=reverse, key=sortFunc)

    return guideData, guideScores, hasNotFound, pamIdToSeq

def printDownloadTableLinks(batchId, addTsv=False):
    print('<div id="downloads" style="text-align:left">')
    print("Download as Excel tables: ", end=' ')
    print('<a href="crispor.py?batchId=%s&download=guides&format=xls">Guides</a>&nbsp;/&nbsp;' % batchId, end=' ')
    if not pamIsFirst and not saCas9Mode:
        print('<a href="crispor.py?batchId=%s&showAllScores=1&download=guides&format=xls">Guides, all scores</a>&nbsp;/&nbsp;' % batchId, end=' ')
    print('<a href="crispor.py?batchId=%s&download=offtargets&format=xls">Off-targets</a>&nbsp;/&nbsp;' % batchId, end=' ')
    print(('<a href="crispor.py?batchId=%s&satMut=1">Saturating mutagenesis assistant</a><br>' % batchId))
    #print "<small>Plasmid Editor: ",
    #print '<a href="crispor.py?batchId=%s&download=genbank">Guides</a></small>' % batchId,

    if addTsv:
        print("<small>Tab-sep format: ", end=' ')
        print('<a href="crispor.py?batchId=%s&download=guides&format=tsv">Guides</a>&nbsp;/&nbsp;' % batchId, end=' ')
        print('<a href="crispor.py?batchId=%s&download=offtargets&format=tsv">Off-targets</a></small>' % batchId, end=' ')

    print('</div>')

def hasGeneModels(org):
    " return true if this organism has gene model information "
    geneFname = join(genomesDir, org, org+".segments.bed")
    return isfile(geneFname)

def printTableHead(pam, batchId, chrom, org, varHtmls, showColumns):
    " print guide score table description and columns "
    # one row per guide sequence
    if not pamIsCpf1(pam):
        print('''<div class='substep'>Ranked by default from highest to lowest specificity score (<a target='_blank' href='http://dx.doi.org/10.1038/nbt.2647'>Hsu et al., Nat Biot 2013</a>). Click on a column title to rank by a score.<br>''')
        #print("""<b>Our recommendation:</b> Use Fusi for in-vivo (U6) transcribed guides, Moreno-Mateos for in-vitro (T7) guides injected into Zebrafish/Mouse oocytes.<br>""")
        print('''If you use this website, please cite our <a href="https://academic.oup.com/nar/article/46/W1/W242/4995687">paper in NAR 2018</a>.''')
        print("Too much information? Look at the <a target=_blank href='manual/'>CRISPOR manual</a>.<p>")
        print('</div>')

    printDownloadTableLinks(batchId)

    print("""
    <script type="text/javascript">
    function allRows() {
        $("guideRow").show();
    }

    //function copySeq() {
        //var c = new ClipboardJS('#seqAsText');
        //var copyText = document.getElementById("seqAsText");
        //var selRes = copyText.select();
        //var val = copyText.value;
        //var res = document.execCommand("copy");
        //alert("The input sequence is now in your clipboard. You can paste it into other programs.");
    //}

    $(document).ready( function() {
        //#$('#copyLink').click( copySeq );
        var clipboard = new ClipboardJS('#copyLink');
        clipboard.on('success', function(e) {
            alert("The input sequence is now in your clipboard. You can paste it into other programs.");
            console.info('Action:', e.action);
            console.info('Text:', e.text);
            console.info('Trigger:', e.trigger);
            e.clearSelection();
        });

        $('d').mouseenter( onEditHover );
        $('d').mouseleave ( onEditOut );
    });

    function onEditOut() {
        $('#editHover').hide();
    }

    function colorChar(str, pos) {
    /* put a span-color tag around the char at pos in str and return result */
        var prefix = str.substring(0, pos);
        var hlChar = str[pos];
        var suffix = str.substring(pos+1);
        return prefix+"<mut>"+hlChar+"</mut>"+suffix;
    }

    function onEditHover(ev) {
    /* user hovers over an edit letter */
        ev.preventDefault();
        var oldEl = document.getElementById("editHover");
        if (oldEl)
            oldEl.remove();

        console.log(ev.target);
        const boundBox = ev.target.getBoundingClientRect();
        var x = boundBox.left;
        var y = boundBox.top;
        y += 14;

        var div = document.createElement('div');
        div.id = "editHover";
        div.style.width="400px";
        div.style.height="200px";
        div.style.border="1px solid black";
        div.style.padding="10px";
        div.style.position="fixed";
        div.style.backgroundColor="white";
        div.style.left=x+"px";
        div.style.top=y+"px";

        var pos = parseInt(this.getAttribute("pos"));
        var nucl = this.textContent;
        if (nucl.toUpperCase()==="T")
            origNucl = "C";
        else
            origNucl = "G";
        var htmls=[];
        htmls.push("The following guides can mutate "+origNucl+" to "+nucl+" at position "+pos+":<br>");
        htmls.push("<table class='editTable'>");
        htmls.push("<tr><th>Guide ID</th><th>Guide Sequence</th><th>Komor score</th><th>Spec. Score</th></tr>");

        var guides = editData[pos][nucl];
        guides.sort( function (a, b) { a[4] - b[4] } ); // sort by komor score
        for (var i=0; i<guides.length; i++) {
            guide = guides[i];
            pamId = guide[0];
            guideSeq = guide[1];
            pam = guide[2];
            mutPos = guide[3];
            beScore = guide[4];
            specScore = guide[5];

            htmls.push("<tr>");

            htmls.push("<td>"+pamId+"</td>");
            htmls.push("<td><tt>"+colorChar(guideSeq, mutPos)+" "+pam+"</tt></td>");
            htmls.push("<td>"+beScore.toFixed(2)+"</td>");
            htmls.push("<td>"+specScore+"</td>");

            htmls.push("</tr>");
        }

        htmls.push("</table>");
        $(div).append(htmls.join(""));
        document.body.appendChild(div);
    }

    function onlyWith(doPrefix) {
        /* show only guide rows and guide sequence viewer features that start with a prefix */

        if ($("#onlyWith"+doPrefix+"Box").prop("checked"))
            {
            $(".prefixBox").prop("checked", false);
            $("#onlyWith"+doPrefix+"Box").prop("checked", true);
            //$(".guideRow").show();
            $(".guideRow").css("visibility", "visible");
            $(".guideRowNoPrefix"+doPrefix).hide();
            // special handling for sequence viewer: hide() would destroy the layout there
            $(".guideRowNoPrefix"+doPrefix+"Seq").css("visibility", "hidden");
            }
        else
            {
            $(".prefixBox").prop("checked", false);
            $(".guideRow").show();
            $(".guideRowNoPrefix"+doPrefix+"Seq").css("visibility", "visible");
            }
    }

    function displayClass(className, dispVal) {
    /* hide in a loop, works around Safari stack size limits that crash jquery functions */
        var els = document.getElementsByClassName(className);
        for (var el of els) {
            el.style.display = dispVal;
        }
    }

    function onlyExons() {
    /* show only off-targets in exons */
        if ($("#onlyExonBox").prop("checked")) {
            $(".otMore").show();
            $(".otMoreLink").hide();
            $(".otLessLink").hide();
            displayClass("notExon", "none");
        }
        else {
            if ($("#onlySameChromBox").prop("checked")) {
                $(".notExon:not(.diffChrom)").show();
            }
            else {
                displayClass("notExon", "block");
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
    """)

    print('<table id="otTable" style="background:white;table-layout:fixed; overflow:scroll; width:100%">')

    print('<thead>')
    print('<tr style="border-bottom:none; border-left:5px solid black; background-color:#F0F0F0">')

    print('<th style="width:80px; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=pos" class="tooltipster" title="Click to sort the table by the position of the PAM site">Position/<br>Strand</a>' % batchId)
    htmlHelp("You can click on the links in this column to highlight the <br>PAM site in the sequence viewer at the top of the page.")
    print('</th>')

    print('<th style="width:235px; border-bottom:none">Guide Sequence + <i>PAM</i><br>')

    print ('+ Restriction Enzymes')
    htmlHelp("Restriction enzymes can be very useful for screening mutations induced by the guide RNA using PCR and Restrictrion frament length polymorphism (RFLP).<br>Enzyme sites shown here overlap the main cleavage site 3bp 5' to the PAM.<br>Digestion of the PCR product with these enzymes will not cut the product if the genome was mutated by Cas9. This is a lot easier than screening with the T7 assay, Surveyor or sequencing.")
    print('<br>')

    if varHtmls is not None:
        print(' + Variants')
        htmlHelp("Variants that overlap the guide sequence are shown. You can change the variant database with the drop-down box above the sequence viewer at the top of the page.")
        print('<br>')

    print('''<small>''')
    print('''<input type="checkbox" class="prefixBox" id="onlyWithGBox" onchange="onlyWith('G')">Only G-''')
    print('''<input type="checkbox" class="prefixBox" id="onlyWithGGBox" onchange="onlyWith('GG')">Only GG-''')
    print('''<input type="checkbox" class="prefixBox" id="onlyWithABox" onchange="onlyWith('A')">Only A-''')
    htmlHelp("The three checkboxes allow you to show only guides that start with GG-, G- or A-. While we recommend prefixing a 20bp guide with G for U6 expression with spCas9, some protocols recommend using only guides with a G- prefix for U6 and A- for U3.")
    print('''</small>''')

    if not pamIsCpf1(pam):
        print('<th style="width:80px; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=spec" class="tooltipster" title="Click to sort the table by specificity score. Hover over the (i) bubble on the right to get more information about the specificity score.">MIT Specificity Score</a>' % batchId)
        if pamIsSaCas9(pam):
            htmlHelp("The higher the specificity score, the lower are off-target effects in the genome.<br>This specificity score has been adapted for SaCas9 and based on the off-target scores shown on mouse-over. The algorithm was provided by Josh Tycko. Like the MIT score for spCas9, it is aggregated from all off-target scores and ranges 0-100. See <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6063963/'>Tycko et al. Nat Comm 2018</a> for details.")
        else:
            htmlHelp("The higher the specificity score, the lower are off-target effects in the genome.<br>The specificity score ranges from 0-100 and measures the uniqueness of a guide in the genome. See <a href='http://dx.doi.org/10.1038/nbt.2647'>Hsu et al. Nat Biotech 2013</a>. We recommend values &gt;50, where possible. See <a target=_blank href='manual/#offs'>the CRISPOR manual</a>")
        print("</th>")

    if "cfdGuideScore" in showColumns:
        print('<th style="width:60px; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=cfdSpec" class="tooltipster" title="Click to sort the table by CFD specificity score">CFD Spec. score</a>' % batchId)
        htmlHelp("The CFD specificity score, inspired like guidescan.com, behaves like the MIT specificity score, but it is based on the more accurate CFD off-target model, from <a href='http://www.nature.com/nbt/journal/v34/n2/full/nbt.3437.html'>Doench 2016</a>, which is also used by Crispor to rank the off-targets. The CFD specificity score correlates better than the MIT score with the total off-target cleavage fraction of a guide, see <a target=_blank href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6731277/'>Tycko et al, Nat Comm 2019</a> and also the <a target=_blank href='/manual/#faq'>CRISPOR manual</a>.")
        print("</th>")

    if len(scoreNames)==2 or pamIsCpf1(pam) or pamIsSaCas9(pam):
       print('<th style="width:150px; height:100px; border-bottom:none" colspan="%d">Predicted Efficiency' % (len(scoreNames)))
    else:
       print('<th style="width:270px; border-bottom:none" colspan="%d">Predicted Efficiency' % (len(scoreNames))) # -1 because proxGc is in scoreNames but has no column

    htmlHelp("The higher the efficiency score, the more likely is cleavage at this position. For details on the scores, mouseover their titles below.<br>Note that these predictions are not very accurate, they merely enrich for more efficient guides by a factor of 2-3 so you have to test a few guides to see the effect. <a target=_blank href='manual/#onEff'>Read the CRISPOR manual</a>")

    if not pamIsCpf1(pam) and not pamIsSaCas9(pam):
        if cgiParams.get("showAllScores", "0")=="0":
            print(("""<br><a style="font-size:12px" href="%s" class="tooltipsterInteract" title="By default, only the two most relevant scores are shown, based on our study <a href='http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2'>Haeussler et al. 2016</a>. Click this link to show all efficiency scores.">Show all scores</a>""" % cgiGetSelfUrl({"showAllScores":"1"}, anchor="otTable")))
            scoreDescs["crisprScan"][0] = "Mor.-Mateos"
        else:
            print(("""<br><a style="font-size:12px" href="%s" class="tooltipsterInteract" title="Show only the two main scores">Show main scores</a>""" % cgiGetSelfUrl({"showAllScores":None}, anchor="otTable")))

    print('</th>')

    mhColName="Outcome"
    if not baseEditor:
        if len(mutScoreNames)<=1:
            mhColName = ""
            oofWidth=45
            #oofDesc = "Click on score to show micro-homology"
            #oofDesc = ""
        else:
            oofWidth=67
            #oofDesc = "Click score for details"

        colSpan = len(mutScoreNames)
        print('<th colspan=%d style="width:%dpx; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=oof" class="tooltipster" title="Prediction of the DNA sequence after strand break repair. Click to sort the table by frameshift/out-of-frame scores. Hover over the score names to show information about a particular score. Click a score number to see the predicted indel pattern around the guide.">%s</a>' % (colSpan, oofWidth, batchId, mhColName))
        #htmlHelp(scoreDescs["oof"][1])
        #print "<small>%s</small>" % oofDesc
        print('</th>')

    print('<th style="width:117px; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=offCount" class="tooltipster" title="Click to sort the table by number of off-targets">Off-targets for <br>0-1-2-3-4 mismatches<br></a><span style="color:grey">+ next to PAM </span>' % (batchId))

    altPamsHelp = [pam]
    if pam in offtargetPams:
        altPamsHelp.extend(offtargetPams[pam])

    htmlHelp("For each number of mismatches, the number of off-targets is indicated.<br>Example: 1-3-20-50-60 means 1 off-target with 0 mismatches, 3 off-targets with 1 mismatch, <br>20 off-targets with 2 mismatches, etc.<br>The CRISPOR website only searches up to four mismatches (use the command line version for 5 or 6). Off-targets are considered if they are flanked by one of these motifs: %s .<br>Shown in grey are the off-targets that have no mismatches in the 12 bp adjacent to the PAM. These are the most likely off-targets." % (", ".join(altPamsHelp)))

    print("</th>")
    print('<th style="width:*; border-bottom:none">Genome Browser links to matches sorted by CFD off-target score')
    htmlHelp("For each off-target the number of mismatches is indicated and linked to a genome browser. <br>Matches are ranked by CFD off-target score (see Doench 2016 et al) from most to least likely.<br>Matches can be filtered to show only off-targets in exons or on the same chromosome as the input sequence.<br>On most organisms, you can click the links below to open a window with a genome browser at this position.")

    print('<br><small>')
    print('<input type="hidden" name="batchId" value="%s">' % batchId)

    if hasGeneModels(org):
        print('''<input type="checkbox" id="onlyExonBox" onchange="onlyExons()">exons only''')
    else:
        print('<small title="When this genome was loaded into CRISPOR, gene models were not available. Contact us if you want to filter for off-targets in exons and think that a gene models are now available for this genome." style="color:grey">No exons.</small>')

    if chrom!="":
        if chrom[0].isdigit():
            chrom = "chrom "+chrom
        print('''<input type="checkbox" id="onlySameChromBox" onchange="onlySameChrom()">%s only''' % chrom)
    else:
        print('<small style="color:grey">&nbsp;No match, no chrom filter</small>')

    print("</small>")
    print("</th>")
    print("</tr>")

    # subheaders
    print('<tr style="border-top:none; border-left: solid black 5px; background-color:#F0F0F0">')

    print('<th style="border-top:none"></th>')
    print('<th style="border-top:none"></th>')

    if "cfdGuideScore" in showColumns:
        print('<th style="border-top:none"></th>')

    if not pamIsCpf1(pam):
        print('<th style="border-top:none"></th>')

    for scoreName in scoreNames:
        if scoreName in ["oof", "proxGc"] or "oof" in scoreName:
            continue
        scoreLabel, scoreDesc = scoreDescs[scoreName]
        print('<th style="width: 10px; border: none; border-top:none; border-right: none" class="rotate"><div><span><a title="%s" class="tooltipsterInteract" href="crispor.py?batchId=%s&sortBy=%s">%s</a></span></div></th>' % (scoreDesc, batchId, scoreName, scoreLabel))

    if "proxGc" in scoreNames:
        # the ProxGC score comes next
        print('''<th style="border: none; border-top:none; border-right: none; border-left:none" class="rotate">''')
        print('''<div><span style="border-bottom:none">''')
        print('''<a title="This column shows two heuristics based on observations rather than computational models: <a href='http://www.cell.com/cell-reports/abstract/S2211-1247%2814%2900827-4'>Ren et al</a> 2014 obtained the highest cleavage in Drosophila when the final 6bp contained &gt;= 4 GCs, based on data from 39 guides. <a href='http://www.genetics.org/content/early/2015/02/18/genetics.115.175166.abstract'>Farboud et al.</a> obtained the highest cleavage in C. elegans for the 10 guides that ended with -GG, out of the 50 guides they tested.<br>The column contains + if the final GC count is &gt;= 4 and GG if the guide ends with GG." href="crispor.py?batchId=%s&sortBy=finalGc6" class="tooltipsterInteract">Prox GC</span></div></th>''' % (batchId))

    # these are empty cells to fill up the row and avoid white space
    for scoreName in mutScoreNames:
        scoreLabel, scoreDesc = scoreDescs[scoreName]
        print('<th style="width: 10px; border-top:none; border-right: none" class="rotate"><div><span><a title="%s" class="tooltipsterInteract" href="crispor.py?batchId=%s&sortBy=%s">%s</a></span></div></th>' % (scoreDesc, batchId, scoreName, scoreLabel))

    print('<th style="border-top:none"></th>')
    print('<th style="border-top:none"></th>')

    print("</tr>")
    print('</thead>')

def scoreToColor(guideScore):
    if guideScore is None:
        color = ("#000000", "black")
    elif guideScore > 50:
        color = ("#32cd32", "green")
    elif guideScore > 20:
        color = ("#ffff00", "yellow")
    elif guideScore==-1:
        color = ("#000000", "black")
    else:
        color = ("#aa0114", "red")
    return color

def hexToRgb(hexCode):
    " convert hex color to RGB in UCSC format, https://stackoverflow.com/questions/29643352/converting-hex-to-rgb-value-in-python "
    hexCode = hexCode.lstrip("#")
    return ",".join(tuple(str(int(hexCode[i:i+2], 16)) for i in (0, 2 ,4)))

def makeOtBrowserLinks(otData, chrom, dbInfo, pamId):
    " return a list with the html texts of the offtarget links "
    links = []

    i = 0
    for otSeq, score, cfdScore, editDist, pos, gene, alnHtml, inLinkage in otData:
        cssClasses = ["tooltipster"]
        if not gene.startswith("exon:"):
            cssClasses.append("notExon")
        if pos.split(":")[0]!=chrom:
            cssClasses.append("diffChrom")
        if inLinkage:
            cssClasses.append("inLinkage")

        classStr =  ""
        if len(cssClasses)!=0:
            classStr = ' class="%s"' % " ".join(cssClasses)

        link = makeBrowserLink(dbInfo, pos, gene, alnHtml, ["tooltipster"])
        editDist = str(editDist)
        links.append( '''<div%(classStr)s>%(editDist)s:%(link)s</div>''' % locals() )

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

def printNote(s):
    print('<div style="text-align:left; background-color: aliceblue; padding:5px; border: 1px solid black"><strong>Note:</strong>')
    print(s)
    print("</div>")

def printWarning(s):
    print('<div style="text-align:left; background-color: #FFDDDD; padding:5px; border: 1px solid black"><strong>Warning:</strong>')
    print(s)
    print('</div>')

def printNoEffScoreFoundWarn(effScoresCount, pam):
    if effScoresCount==0 and not pamIsCpf1(pam):
        note = "No guide could be scored for efficiency. This happens when the input sequence is shorter than 100bp and there is no genome available to extend it or if there is simply not guide socring method. In the first case, please add flanking 50bp on both sides of the input sequence and submit this new, longer sequence. For the second case, you can contact me and suggest an efficiency scoring method, send me the published paper in this case."
        printNote(note)

def showGuideTable(guideData, pam, otMatches, dbInfo, batchId, org, chrom, varHtmls):
    " shows table of all PAM motif matches "
    print("<br><div class='title'>Predicted guide sequences for PAMs</div>")

    global scoreNames
    if (cgiParams.get("showAllScores", "0")=="1"):
        scoreNames = allScoreNames

    showColumns = set()

    # show the CFD guide score?
    if pamIsSpCas9(pam):
        showColumns.add("cfdGuideScore")

    showPamWarning(pam)
    showNoGenomeWarning(dbInfo)
    printTableHead(pam, batchId, chrom, org, varHtmls, showColumns)

    count = 0
    effScoresCount = 0
    showProxGcCol = ("proxGc" in scoreNames)

    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, pamStart, guideStart, strand, pamId, guideSeq, \
            pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, repCount = guideRow

        color = scoreToColor(guideScore)[0]

        classStr = cssClassesFromSeq(guideSeq)
        print('<tr id="%s" class="%s" style="border-left: 5px solid %s">' % (pamId, classStr, color))

        # position and strand
        #print '<td id="%s">' % pamId
        print('<td>')
        print('<a href="#list%s">' % (pamId))
        print(str(pamStart+1)+" /")
        if strand=="+":
            print('fw')
        else:
            print('rev')
        print('</a>')
        print("</td>")

        # sequence with variants and PCR primer link
        print("<td>")
        print("<small>")

        # guide sequence + PAM sequence
        if pamIsFirst:
            fullGuideHtml = "<tt><i>"+pamSeq+"</i> " + guideSeq+"</tt>"
            spacePos = len(pamSeq)
        else:
            fullGuideHtml = "<tt>"+guideSeq + " <i>" + pamSeq+"</i></tt>"
            spacePos = len(guideSeq)
        print(fullGuideHtml)
        print("<br>")

        # variant-string
        if varHtmls is not None:
            varFound = False
            varStrs = []
            guideHtmlStart = min(guideStart, pamStart)
            guideHtmls = varHtmls[guideHtmlStart:guideHtmlStart+len(guideSeq)+len(pamSeq)]
            if strand=="-":
                guideHtmls = list(reversed(guideHtmls))

            for i in range(len(guideHtmls)):
                html = guideHtmls[i]
                if html!=".":
                    varFound = True
                if i==spacePos:
                    varStrs.append("&nbsp;")
                varStrs.append(html)
            print(("<tt style='color:#888888'>%s</tt><br>" % ("".join(varStrs))))

        if "TTTT" in guideSeq.upper():
            text = "This guide contains the sequence TTTT. It cannot be transcribed with a U6 or U3 promoter, as TTTT terminates the transcription."
            htmlWarn(text)
            print(' Not with U6/U3')
            print("<br>")

        if pam=="NGG":
            grafType = crisporEffScores.getGrafType(guideSeq)
            if grafType:
                if grafType=="tt":
                    grafText = "The guide ends with TTC or TTT or contains only T and C in the last four nucleotides and more than 2 Ts or at least one TT and one T or C ('TT-motif'). These guides should be avoided in polymerase III (Pol III)-based gene editing experiments requiring high sgRNA expression levels."
                elif grafType=="gcc":
                    grafText = "The guide ends with [AGT]GCC or GCCT ('GCC motif'). These sgRNAs appear to be inefficient irrespective of the delivery method and should thus be generally avoided."

                text = "This guide contains one of the motifs described by <a target=_blank href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6352712/'>Graf et al, Cell Reports 2019</a>. %s " % grafText
                htmlWarn(text)
                print(' Inefficient')
                print("<br>")

        if gcContent(guideSeq)>0.75:
            text = "This sequence has a GC content higher than 75%.<br>In the data of Tsai et al Nat Biotech 2015, the two guide sequences with a high GC content had almost as many off-targets as all other sequences combined. We do not recommend using guide sequences with such a high GC content."
            htmlWarn(text)
            print(' High GC content')
            print("<br>")

        if gcContent(guideSeq)<0.25:
            text = "This sequence has a GC content lower than 25%.<br>In the data of Wang/Sabatini/Lander Science 2014, guides with a very low GC content had low cleavage efficiency."
            htmlWarn(text)
            print(' Low GC content<br>')
            print("<br>")

        if len(mutEnzymes)!=0:
            print("<div style='margin-top: 3px'>Enzymes: <i>", end=' ')
            print(", ".join([x.split("/")[0] for x,y,z in list(mutEnzymes.keys())]))
            print("</i></div>")

        scriptName = basename(__file__)
        if otData!=None and repCount == 0:
            print(('&nbsp;<a href="%s?batchId=%s&pamId=%s&pam=%s" target="_blank"><strong>Cloning / PCR primers</strong></a>' % (scriptName, batchId, urllib.parse.quote(str(pamId)), pam) ))

        print("</small>")
        print("</td>")

        # off-target score, aka specificity score aka MIT score
        if not pamIsCpf1(pam):
            print("<td>")
            if guideScore==None:
                print("No matches")
            else:
                print("%d" % guideScore)
            print("</td>")

        # guide score based on CFD scores, aka guidescan score
        if "cfdGuideScore" in showColumns:
            print("<td>")
            if guideCfdScore==None:
                print("No matches")
            else:
                print("%d" % guideCfdScore)
            print("</td>")

        # eff scores
        if effScores==None:
            print('<td colspan="%d">Too close to end</td>' % len(scoreNames))
            htmlHelp("The efficiency scores require some flanking sequence<br>This guide does not have enough flanking sequence in your input sequence and could not be extended as it was not found in the genome.<br>")
        else:
            for scoreName in scoreNames:
                # out-of-frame and prox. gc need special treatment
                if scoreName in ["oof", "proxGc"]:
                    continue
                score = effScores.get(scoreName, None)
                if score!=None:
                    effScoresCount += 1
                if score==None:
                    print('''<td>--</td>''')
                elif scoreName=="ssc":
                    # save some space
                    numStr = '%.1f' % (float(score))
                    print('''<td style="font-size:small">%s</td>'''  % numStr)
                elif scoreDigits.get(scoreName, 0)==0:
                    print('''<td>%d</td>''' % int(score))
                else:
                    print('''<td>%0.1f</td>''' % (float(score)))
            #print "<!-- %s -->" % seq30Mer

        if showProxGcCol:
            print("<td>")
            # close GC > 4
            finalGc = int(effScores.get("finalGc6", -1))
            if finalGc==1:
                print("+")
            elif finalGc==0:
                print("-")
            else:
                print("--")

            # main motif is "NGG" and last nucleotides are GGNGG
            if int(effScores.get("finalGg", 0))==1:
                print("<br>")
                print("<small>-GG</small>")
            print("</td>")

        if not baseEditor:
            for mutScoreName in mutScoreNames:
                print("<td>")
                oofScore = str(effScores.get(mutScoreName, None))

                if mutScoreName=="oof":
                    scoreDesc = "out-of-frame deletions"
                else:
                    scoreDesc = "frameshift mutations"

                if oofScore==None or oofScore=="None":
                    print("--")
                else:
                    print("""<a href="%s?batchId=%s&pamId=%s&showMh=%s" target=_blank class="tooltipster" title="This score indicates how likely %s are. Click to show the induced deletions based on the micro-homology around the cleavage site.">%s</a>""" % (myName, batchId, urllib.parse.quote(pamId), mutScoreName, scoreDesc, oofScore))
                    #print """<br><br><small><a href="%s?batchId=%s&pamId=%s&showMh=1" target=_blank class="tooltipster">Micro-homology</a></small>""" % (myName, batchId, pamId)
                print("</td>")

        # mismatch description
        print("<td>")
        #otCount = sum([int(x) for x in otDesc.split("/")])
        if otData==None:
            # no genome match
            print(otDesc)
            htmlHelp("This exact sequence was not found in the genome.<br>If you have pasted a cDNA multi-exon sequence, note that sequences that overlap a splice site cannot be used as guide sequences. If you only have a cDNA sequence, please BLAST or BLAT your sequence first against the genome, then use the resulting exon from the genome for CRISPOR.<br>This warning also appears if you have selected the wrong or no genome.")
        elif repCount > 0:
            print ("Repeat")
            htmlHelp("At <= 4 mismatches, %d alignments were found in the genome for this sequence, without looking at the PAM sequence around these alignments.<br>This guide is a repeated region, it is too unspecific.<br>Usually, CRISPR cannot be used to target repeats. Also, note that sequences that include long repeats will make the CRISPOR website slow. You can mask repeats with Ns to speed up the search." % repCount)
        else:
            print(otDesc)
            print("<br>")

            # mismatch description, last 12 bp
            print('<small style="color:grey">'+last12Desc+"</small><br>")
            otCount = len(otData)
            print("<br><small>%d off-targets</small>" % otCount)
        print("</td>")

        # links to offtargets
        print("<td><small>")
        if otData!=None:
            if len(otData)>500 and len(guideData)>1:
                otData, cutoff = findOtCutoff(otData)
                if cutoff==None:
                    print("More than 1000 off-targets, showing only top "+str(len(otData)))
                else:
                    print("More than 500 off-targets, showing %d with score &gt;%0.1f " % (len(otData), cutoff))

                htmlHelp("This guide sequence has a high number of off-targets, its use is discouraged.<br>To show all off-targets, paste only the guide sequence into the input sequence box.")

            otLinks = makeOtBrowserLinks(otData, chrom, dbInfo, pamId)

            print("\n".join(otLinks[:3]))
            if len(otLinks)>3:
                cssPamId = pamId.replace("-","minus").replace("+","plus") # +/-: not valid in css
                cssPamId = cssPamId+"More"
                print('<div id="%s" class="otMore" style="display:none; width:100%%">' % cssPamId)
                print("\n".join(otLinks[3:]))

                print('''<a style="float:right;text-decoration:underline" href="%s?batchId=%s&pamId=%s&otPrimers=1" id="%s">''' % (myName, batchId, urllib.parse.quote(pamId), cssPamId))
                print('<strong>Off-target primers</strong></a>')

                print('</div>')

                print('''<a id="%sMoreLink" class="otMoreLink" onclick="showAllOts('%s')">''' % (cssPamId, cssPamId))
                print('show all...</a>')


                print('''<a id="%sLessLink" class="otLessLink" style="display:none" onclick="showLessOts('%s')">''' % (cssPamId, cssPamId))
                print('show less...</a>')



        print("</small></td>")

        print("</tr>")
        count = count+1

    print("</table>")
    printDownloadTableLinks(batchId, addTsv=True)

    printNoEffScoreFoundWarn(effScoresCount, pam)

def linkLocalFiles(listFname):
    """ write a <link> statement for each filename in listFname. Version them via mtime
    (-> browser cache)
    """
    for fname in open(listFname).read().splitlines():
        fname = fname.strip()
        if not isfile(fname):
            fname = join(HTMLDIR, fname)
            if not isfile(fname):
                print("missing: %s<br>" % fname)
                continue
        mTime = str(os.path.getmtime(fname)).split(".")[0] # seconds is enough
        if fname.endswith(".css"):
            #url = fname.replace("/var/www/", "http://tefor.net/")
            print("<link rel='stylesheet' media='screen' type='text/css' href='%s%s?%s'/>" % (HTMLPREFIX, fname, mTime))

def printHeader(batchId, title):
    " print the html header "

    print('''<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">''')
    print("<html><head>")

    if title==None:
        if batchName!="":
            print("""<title>CRISPOR - %s</title>""" % batchName)
        else:
            print("""<title>CRISPOR</title>""")
    else:
        print("""<title>%s</title>""" % title)
    print("""
<meta name='description' content='Design CRISPR guides with off-target and efficiency predictions, for more than 100 genomes.'/>
<meta http-equiv='Content-Type' content='text/html; charset=utf-8' />
<meta property='fb:admins' content='692090743' />
<meta name="google-site-verification" content="OV5GRHyp-xVaCc76rbCuFj-CIizy2Es0K3nN9FbIBig" />
<meta property='og:type' content='website' />
<meta property='og:url' content='http://crispor.gi.ucsc.edu/' />
<meta property='og:image' content='http://crispor.gi.ucsc.edu/image/CRISPOR.png' />
<script src="https://cdn.jsdelivr.net/npm/clipboard@2/dist/clipboard.min.js"></script>

""")

    # load jquery from local copy, not from CDN, for offline use
    print(("""<script src='%sjs/jquery.min.js'></script>
<script src='%sjs/jquery-ui.min.js'></script>
""" % (HTMLPREFIX, HTMLPREFIX)))

    #print('<link rel="stylesheet" href="//fonts.googleapis.com/css?family=Roboto:300,300italic,700,700italic" />')
    #print('<link rel="stylesheet" type="text/css" href="https://cdnjs.cloudflare.com/ajax/libs/normalize/5.0.0/normalize.min.css" />')
    #print('<link rel="stylesheet" href="//cdn.rawgit.com/milligram/milligram/master/dist/milligram.min.css">')

    linkLocalFiles("includes.txt")

    print('<link rel="stylesheet" type="text/css" href="%sstyle/tooltipster.css" />' % HTMLPREFIX)
    print('<link rel="stylesheet" type="text/css" href="%sstyle/tooltipster-shadow.css" />' % HTMLPREFIX)
    print('<link rel="stylesheet"  href="https://cdnjs.cloudflare.com/ajax/libs/chosen/1.6.2/chosen.css" />')
    print('<link rel="stylesheet"  href="https://cdn.jsdelivr.net/npm/source-code-pro@2.38.0/source-code-pro.css" />')

    # the UFD combobox, https://code.google.com/p/ufd/wiki/Usage
    # patched to allow mouse wheel
    # https://code.google.com/p/ufd/issues/detail?id=86&q=mouse%20wheel
    print('<script type="text/javascript" src="%sjs/jquery.ui.ufd.js"></script>' % HTMLPREFIX)
    #print '<link rel="stylesheet" type="text/css" href="%sstyle/ufd-base.css" />' % HTMLPREFIX
    print('<link rel="stylesheet" type="text/css" href="%sstyle/plain.css" />' % HTMLPREFIX)
    print('<link rel="stylesheet" type="text/css"  href="%sstyle/jquery-ui.css" />' % HTMLPREFIX)
    print('<script type="text/javascript" src="js/jquery.tooltipster.min.js"></script>')
    print('<script src="https://cdnjs.cloudflare.com/ajax/libs/chosen/1.6.2/chosen.jquery.min.js"></script>')

    # override the main TEFOR css
    print("""
<style>

select { font-size: 80%; }

body {
   text-align: left;
   /* float: left; */
}
p {
}
ul {
    -webkit-margin-before: 0;
    -webkit-margin-after: 0;
}

.editTable {
    border: 1px solid black;
    background-color: white;
}

mut {
    color: blue;
    background-color: yellow;
}

.editTable th {
    background-color: #F0F0F0;
}

tt { font-size: 90% }
div.contentcentral { text-align: left; float: left}

/* for chosen.js */
.chosen-container { width: 600px }
.chosen-container .chosen-results li.active-result { float: left}

</style>
""")

    # activate tooltipster
   #theme: 'tooltipster-shadow',
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

    $(function () {
       $(".tooltipAuto").tooltip({
       contentAsHtml : true
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
      """)

    # if we're showing all scores, we have very little space, so turn by
    # 90degrees. otherwise, we can afford 45 degrees, which is easier to read.
    #if cgiParams.get("showAllScores", "0")=="1":
    print("""
    -webkit-transform: rotate(-90);
    -moz-transform: rotate(270deg);
    -ms-transform: rotate(270deg);
    -o-transform: rotate(270deg);
    transform: rotate(270deg);
    width: 25px;""")

    #else:
          #print("""
         #-webkit-transform: rotate(-45deg);
         #-moz-transform: rotate(315deg);
         #-ms-transform: rotate(315deg);
         #-o-transform: rotate(315deg);
         #transform: rotate(315deg);
         #width: 25px;""")

    print("""
       }

       th.rotate > div > span {
         /* border-bottom: 1px solid #ccc; */
         padding: 0px 3px;
         white-space: nowrap;
       }
    </style>""")


    print("</head>")

    print('<body id="wrapper">')

def firstFreeLine(lineMasks, y, start, end):
    " recursively search for first free line to place a feature (start, end) "
    #print "first free line called with y", y, "<br>"
    if y>=len(lineMasks):
        return None
    lineMask = lineMasks[y]
    for x in range(start, end):
        #print "checking pos", x, "<br>"
        if lineMask[x]!=0:
            return firstFreeLine(lineMasks, y+1, start, end)
    return y
    #return None

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

        # Cannot use Unicode here: these symbols are not part of the
        # monospace font on some platforms and therefore their width
        # is not the same as the other characters
        #arrNE = u'\u2197'
        #arrSE = u'\u2198'
        arrNE = '/'
        arrSE = '\\'
        #arrNE = u'\u2a3c' # hebrew
        #arrSE = u'\ufb27' # math
        ftSeq = seq[start:end]
        if strand=="+":
            if pamIsFirst:
                label = '%s'%(ftSeq)+'.................%s....%s' % (arrNE, arrSE)
                startFt = start
                endFt = start+len(label)
            else:
                #label = '%s..%s'%(seq[start-3].lower(), ftSeq)
                #label = '---%s'%(ftSeq)
                #label = '&#45;&#45;&#45;%s'%(ftSeq)
                label = '&#8722;&#8722;&#8722;%s'%(ftSeq)
                startFt = start - 3
                endFt = end
        else:
            if pamIsFirst:
                spc1 = "...."
                spc2 = "................."
                labelPrefix = '%s%s%s%s' % (arrSE, spc1, arrNE, spc2)
                label = labelPrefix + ftSeq
                startFt = start - len(labelPrefix)
                endFt = startFt+len(label)
            else:
                #label = '%s..%s'%(ftSeq, seq[end+2].lower())
                label = '%s&#45;&#45;&#45;'%(ftSeq)
                startFt = start
                endFt = end + 3

        #print "feature", strand, start, startFt, endFt, SLOP,"<br>"
        #print "mask", lineMasks[0][startFt:endFt], "<br>"
        y = firstFreeLine(lineMasks, 0, startFt, endFt)
        #print "free line: %s<br>" % y
        if y==None:
            errAbort("not enough space to plot features")

        # fill the current mask
        mask = lineMasks[y]
        maskStart = max(startFt-SLOP, 0)
        maskEnd = min(endFt+SLOP, len(seq))
        #print "mask:", maskStart, maskEnd
        for i in range(maskStart, maskEnd):
            mask[i]=1

        maxY = max(y, maxY)

        pamId = "s%d%s" % (start, strand)
        ft = (startFt, endFt, label, strand, pamId)
        #print "labelLen: %d<br>" % len(label)
        #print "ft: %s<br>" % repr(ft)
        ftsByLine[y].append(ft )
    return ftsByLine, maxY

def writePamFlank(seq, startDict, pam, faFname):
    " write pam flanking sequences to fasta file, optionally with versions where each nucl is removed "
    #print "writing pams to %s<br>" % faFname
    faFh = open(faFname, "w")
    for pamId, pamStart, guideStart, strand, flankSeq, pamSeq, pamPlusSeq in flankSeqIter(seq, startDict, len(pam), True):
        faFh.write(">%s\n%s\n" % (pamId, flankSeq))
    faFh.close()

def runCmd(cmd, ignoreExitCode=False, useShell=True):
    " run shell command, check ret code, replaces BIN and SCRIPTS special variables "
    if useShell:
        cmd = cmd.replace("$BIN", binDir)
        cmd = cmd.replace("$PYTHON", sys.executable)
        cmd = cmd.replace("$SCRIPT", scriptDir)
        cmd = "set -o pipefail; " + cmd
        executable = "/bin/bash"
    else:
        cmd = [x.replace("$BIN", binDir).replace("$PYTHON", sys.executable).replace("$SCRIPT", scriptDir) for x in cmd]
        executable=None

    debug("Running %s" % cmd)
    ret = subprocess.call(cmd, shell=useShell, executable=executable)
    if ret!=0 and not ignoreExitCode:
        if not useShell:
            cmd = " ".join(cmd)
        if commandLineMode:
            logging.error("Error: could not run command %s." % cmd)
            sys.exit(1)
        else:
            print("Server error: could not run command %s, error %d.<p>" % (cmd, ret))
            print("please send us an email, we will fix this error as quickly as possible. %s " % contactEmail)
            raise
            sys.exit(0)

def isAltChrom(chrom):
    """ return true is chrom name looks like it's not on the primary assembly. This is mostly relevant for hg38.
    examples: chr6_*_alt (hg38)
    """
    return chrom.endswith("_alt")

def parseOfftargets(db, batchId, onTargetChrom=""):
    """ parse a bed file with annotataed off target matches from overlapSelect,
    has two name fields, one with the pam position/strand and one with the
    overlapped segment
    return as dict pamId -> editDist -> (chrom, start, end, seq, strand, segType, segName, totalAlnCount, isRep)
    segType is "ex" "int" or "ig" (=intergenic)
    if intergenic, geneNameStr is two genes, split by |

    The isRep flag is true if BWA reported more than one alignment with X0+X1 but didn't report these with the
    XA tag. It means that we can't get the alignments for this sequence from BWA (=repeats).
    """
    # edge case: target is on chr6_alt -> we remove a single off-target with 0 mismatches on chr6
    targetIsAlt = isAltChrom(onTargetChrom)
    # keep track of pamIds already handled for this edge case
    skippedPams = set()
    # ideally we would check if the offtarget falls into the chrom area that gave rise to the alt
    # but that would mean parsing yet another non-small file and this case should be sufficiently rare
    # to not bother 99% of users.

    batchBase = join(batchDir, batchId)
    bedFname = batchBase+".bed.gz"
    # example input:
    # chrIV 9864393 9864410 s41-|-|5|ACTTGACTG|0    chrIV   9864303 9864408 ex:K07F5.16
    # chrIV   9864393 9864410 s41-|-|5|ACTGTAGCTAGCT|9999    chrIV   9864408 9864470 in:K07F5.16
    debug("reading offtargets from %s" % bedFname)

    # first sort into dict (pamId,chrom,start,end,editDist,strand)
    # -> (segType, segName)
    pamData = {}

    #ifh = open(bedFname) # switched to gzip compression in Dec 2018, converted old files with bash script
    ifh = gzip.open(bedFname, "rt")

    for line in ifh:
        fields = line.rstrip("\n").split("\t")
        chrom, start, end, name, segment = fields
        logging.debug("off-target: %s" % name)
        # hg38: ignore alternate chromosomes otherwise the
        # regions on the main chroms look as if they could not be
        # targeted at all with Cas9

        if isAltChrom(chrom):
            logging.debug("skipping off-target: on alt-chromosome")
            continue
        nameFields = name.split("|")
        pamId, strand, editDist, seq = nameFields[:4]
        #print pamId, strand, editDist, seq, chrom, start, end, name, segment, "<br>"

        if targetIsAlt:
            if editDist=='0' and onTargetChrom.split("_")[0]==chrom.split("_")[0] and not pamId in skippedPams:
                logging.debug("altChrom edge case: target is on alt-chrom, skipping a single 0-mismatch off-target on primary chrom")
                skippedPams.add(pamId)
                continue

        isRep = 0
        totalAlnCount = 0
        # for compatibility with old bed files, only parse these fields if they are present
        # note: for some reason, the MIT hitScore was always written to these files
        # However, it's not parsed here and never was. In order to not break the old files
        # I kept it in the files, but am not reading it here
        # these are the different formats until now:
        # seqId+"|"+strand+"|"+editDist+"|"+seq+"|"+str(hitScore) # five fields
        # seqId+"|"+strand+"|"+editDist+"|"+seq+"|"+x1Score+"|"+str(hitScore) # six fields
        # (x1Score was roughly the alnCount, similar enough for practical purposes, fixed in 2019)
        # seqId+"|"+strand+"|"+editDist+"|"+seq+"|"+alnCount+"|"+str(hitScore)+"|"+isRep # seven fields
        if len(nameFields)>5:
            totalAlnCount = int(nameFields[4])
            if len(nameFields)>6:
                isRep = bool(int(nameFields[6]))

        editDist = int(editDist)
        # some gene models include colons
        if ":" in segment:
            segType, segName = segment.split(":", maxsplit=1)
        else:
            segType, segName = "", segment
        start, end = int(start), int(end)
        otKey = (pamId, chrom, start, end, editDist, seq, strand, totalAlnCount, isRep)

        # if an offtarget is in the PAR region, we keep only the chrY off-target
        parNum = isInPar(db, chrom, start, end)
        # keep only matches on chrX
        if parNum is not None and chrom=="chrX":
            logging.debug("off-target on PAR region, skipping")
            continue

        # if a offtarget overlaps an intron/exon or ig/exon boundary it will
        # appear twice; in this case, we only keep the exon offtarget
        if otKey in pamData and segType!="ex":
            logging.debug("skipping off-target: ex/ig boundary edge case")
            continue
        pamData[otKey] = (segType, segName)


    # index by pamId and edit distance
    indexedOts = defaultdict(dict)
    for otKey, otVal in pamData.items():
        pamId, chrom, start, end, editDist, seq, strand, totalAlnCount, isRep = otKey
        segType, segName = otVal
        otTuple = (chrom, start, end, seq, strand, segType, segName, totalAlnCount, isRep)
        indexedOts[pamId].setdefault(editDist, []).append( otTuple )

    return indexedOts

class ConsQueue:
    """ a pseudo job queue that does nothing but report progress to the console """
    def startStep(self, batchId, desc, label):
        logging.info("Progress %s - %s - %s" % (batchId, desc, label))

def annotateBedWithPos(inBed, outBed, genome):
    """
    given an input bed4 and an output bed filename, add an additional column 5 to the bed file
    that is a descriptive text of the chromosome pos (e.g. chr1:1.23 Mbp).
    """
    ofh = gzip.open(outBed, "wt")
    for line in open(inBed):
        chrom, start = line.split("\t")[:2]
        chrom = applyChromAlias(genome, chrom)
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

def findAllGuides(seq, pam):
    startDict, endSet = findAllPams(seq, pam)
    pamInfo = list(flankSeqIter(seq, startDict, len(pam), False))
    return pamInfo

def extractMutScores(scoreDict, pamIds):
    " make a list of the guide-related outcome scores in the order of pamIds "
    res = []
    for pamId in pamIds:
        res.append(scoreDict[pamId][0])
    return res

def calcSaveEffScores(batchId, seq, extSeq, pam, queue):
    """ given a sequence and an extended sequence, get all potential guides
    with pam, extend them to 100mers and score them with various eff. scores.
    Return a
    list of rows [headers, (guideSeq, 100mer, score1, score2, score3,...), ... ]

    Also write the results to a database so they can be retrieved later.

    extSeq can be None, if we were unable to extend the sequence
    """
    seq = seq.upper()
    if extSeq:
        extSeq = extSeq.upper()

    pamInfo = findAllGuides(seq, pam)

    pamIds = []
    guides = []
    longSeqs = []

    for pamId, startPos, guideStart, strand, guideSeq, pamSeq, pamPlusSeq in pamInfo:
        logging.debug("PAM ID: %s - guideSeq %s" % (pamId, guideSeq))
        gStart, gEnd = pamStartToGuideRange(startPos, strand, len(pam))
        longSeq = getExtSeq(seq, gStart, gEnd, strand, 50-GUIDELEN, 50, extSeq) # +-50 bp from the end of the guide
        if longSeq!=None:
            longSeqs.append(longSeq)
            pamIds.append(pamId)
            guides.append(guideSeq+pamSeq)

    if len(longSeqs)>0 and doEffScoring:
        enz = None
        if pamIsCpf1(pam) and not pam=="NGTN":
            enz = "cpf1"
        elif pamIsSaCas9(pam):
            enz = "sacas9"

        # for spcas9, we use the extended list for the calculation
        global scoreNames
        if enz is None:
            scoreNames = allScoreNames

        effScores = crisporEffScores.calcAllScores(longSeqs, enzyme=enz, scoreNames=scoreNames)

        # these are slow algorithms, so store the results for later
        queue.startStep(batchId, "outcome", "Calculating editing outcomes")
        mutScores = crisporEffScores.calcMutSeqs(pamIds, longSeqs, enz, scoreNames=mutScoreNames)
        saveOutcomeData(batchId, mutScores)

        # for output and sorting, it's easier to treat the outcome-derived scores like an efficiency score
        for mutScoreName in mutScoreNames:
            if mutScoreName in mutScores:
                effScores[mutScoreName] = extractMutScores(mutScores[mutScoreName], pamIds)

        # make sure the "N bug" reported by Alberto does never happen again:
        # we must get back as many scores as we have sequences
        for scoreName, scores in effScores.items():
            if len(scores)!=len(longSeqs):
                print("Internal error when calculating score %s" % scoreName)
                assert(False)
    else:
        effScores = {}

    activeScoreNames = list(effScores.keys())

    # reformat to rows, write all scores to file
    assert(len(pamIds)==len(guides)==len(longSeqs))
    rows = []
    for i, (guideId, guide, longSeq) in enumerate(zip(pamIds, guides, longSeqs)):
        row = [guideId, guide, longSeq]
        for scoreName in activeScoreNames:
            scoreList = effScores[scoreName]
            if len(scoreList) > 0:
                row.append(scoreList[i])
            else:
                row.append("noScore?")
        rows.append(row)

    headerRow = ["guideId", "guide", "longSeq"]
    headerRow.extend(activeScoreNames)
    rows.insert(0, headerRow)
    return rows

def writeRow(ofh, row):
    " write list to file as tab-sep row "
    row = [str(x) for x in row]
    ofh.write("\t".join(row))
    ofh.write("\n")

def createBatchEffScoreTable(batchId, queue):
    """ annotate all potential guides with efficiency scores and write to file.
    tab-sep file for easier debugging, no pickling
    """
    outFname = join(batchDir, batchId+".effScores.tab")

    # Todo: why don't we get these from the caller as arguments instead of reading the batch?
    batchInfo = readBatchAsDict(batchId)
    seq = batchInfo["seq"]
    extSeq = batchInfo.get("extSeq") # cannot always extend a sequence, e.g. when no perfect match
    pam = batchInfo["pam"]
    pam = setupPamInfo(pam)
    seq = seq.upper()
    if extSeq:
        extSeq = extSeq.upper()

    guideRows = calcSaveEffScores(batchId, seq, extSeq, pam, queue)
    guideFh = open(outFname, "w")
    for row in guideRows:
        writeRow(guideFh, row)
    guideFh.close()
    logging.info("Wrote eff scores to %s" % guideFh.name)

def readEffScores(batchId):
    " parse eff scores from tab sep file and return as dict pamId -> dict of scoreName -> value "
    effScoreFname = join(batchDir, batchId)+".effScores.tab"

    seqToScores = {}

    if isfile(effScoreFname):
        for row in lineFileNext(open(effScoreFname)):
            scoreDict = {}
            rowDict = row._asdict()
            # the first three fields are the pamId, shortSeq, longSeq, they are not scores
            allScoreNames = row._fields[3:]
            for scoreName in allScoreNames:
                score = rowDict[scoreName]
                if score=="None":
                    score = "NA"
                elif "." in score or "e" in score:
                    score = float(score)
                else:
                    score = int(score)
                scoreDict[scoreName] = score
            seqToScores[row.guideId] = scoreDict

    return seqToScores

def findOfftargetsBwa(queue, batchId, batchBase, faFname, genome, pamDesc, bedFname):
    " align faFname to genome and create matchedBedFname "
    matchesBedFname = batchBase+".matches.bed"
    saFname = batchBase+".sa"
    pam = setupPamInfo(pamDesc)
    pamLen = len(pam)
    genomeDir = genomesDir # make var local, see below

    open(matchesBedFname, "w") # truncate to 0 size

    # increase MAXOCC if there is only a single query, but only in CGI mode
    #if len(parseFasta(open(faFname)))==1 and not commandLineMode:
        #global MAXOCC
        #global maxMMs
        #MAXOCC=max(HIGH_MAXOCC, MAXOCC)
        #maxMMs=HIGH_maxMMs

    maxDiff = maxMMs
    queue.startStep(batchId, "bwa", "Alignment of potential guides, mismatches <= %d" % maxDiff)
    convertMsg = "Converting alignments"
    seqLen = GUIDELEN

    bwaM = MFAC*MAXOCC # -m is queue size in bwa
    cmd = "$BIN/bwa aln -o 0 -m %(bwaM)s -n %(maxDiff)d -k %(maxDiff)d -N -l %(seqLen)d %(genomeDir)s/%(genome)s/%(genome)s.fa %(faFname)s > %(saFname)s" % locals()
    runCmd(cmd)

    queue.startStep(batchId, "saiToBed", convertMsg)
    maxOcc = MAXOCC # create local var from global
    # EXTRACTION OF POSITIONS + CONVERSION + SORT/CLIP
    # the sorting should improve the twoBitToFa runtime
    python = sys.executable
    cmd = "$BIN/bwa samse -n %(maxOcc)d %(genomeDir)s/%(genome)s/%(genome)s.fa %(saFname)s %(faFname)s | $SCRIPT/xa2multi.pl | %(python)s $SCRIPT/samToBed %(pam)s %(seqLen)d | sort -k1,1 -k2,2n | $BIN/bedClip stdin %(genomeDir)s/%(genome)s/%(genome)s.sizes stdout >> %(matchesBedFname)s " % locals()
    runCmd(cmd)

    filtMatchesBedFname = batchBase+".filtMatches.bed"
    queue.startStep(batchId, "filter", "Removing matches without a PAM motif")
    altPats = ",".join(offtargetPams.get(pam, ["na"]))
    bedFnameTmp = bedFname+".tmp"
    altPamMinScore = str(ALTPAMMINSCORE)
    shmFaFname = join("/dev/shm", genome+".fa")

    # EXTRACTION OF SEQUENCES + ANNOTATION - big headache!!
    # twoBitToFa was 15x slower than python's twobitreader, after markd's fix it is better
    # but bedtools uses an fa.idx file and also mmap, so is a LOT faster
    # arguments: guideSeq, mainPat, altPats, altScore, passTotalAlnCount
    if isfile(shmFaFname):
        logging.info("Using bedtools and genome fasta on ramdisk, %s" % shmFaFname)
        cmd = "time bedtools getfasta -s -name -fi %(shmFaFname)s -bed %(matchesBedFname)s -fo /dev/stdout | $SCRIPT/filterFaToBed %(faFname)s %(pam)s %(altPats)s %(altPamMinScore)s > %(filtMatchesBedFname)s" % locals()
    else:
        cmd = "time $BIN/twoBitToFa %(genomeDir)s/%(genome)s/%(genome)s.2bit stdout -bed=%(matchesBedFname)s | %(python)s $SCRIPT/filterFaToBed %(faFname)s %(pam)s %(altPats)s %(altPamMinScore)s > %(filtMatchesBedFname)s" % locals()
    #cmd = "$SCRIPT/twoBitToFaPython %(genomeDir)s/%(genome)s/%(genome)s.2bit %(matchesBedFname)s | $SCRIPT/filterFaToBed %(faFname)s %(pam)s %(altPats)s %(altPamMinScore)s %(maxOcc)d > %(filtMatchesBedFname)s" % locals()
    runCmd(cmd)

    segFname = "%(genomeDir)s/%(genome)s/%(genome)s.segments.bed" % locals()

    # if we have gene model segments, annotate them, otherwise just use the chrom position
    if isfile(segFname):
        queue.startStep(batchId, "genes", "Annotating matches with genes")
        cmd = "cat %(filtMatchesBedFname)s | $BIN/overlapSelect %(segFname)s stdin stdout -mergeOutput -selectFmt=bed -inFmt=bed | cut -f1,2,3,4,8 | gzip > %(bedFnameTmp)s " % locals()
        runCmd(cmd)
    else:
        queue.startStep(batchId, "chromPos", "Annotating matches with chromosome position")
        annotateBedWithPos(filtMatchesBedFname, bedFnameTmp, genome)

    # make sure the final bed file is never in a half-written state,
    # as it is our signal that the job is complete
    shutil.move(bedFnameTmp, bedFname)
    queue.startStep(batchId, "done", "Job completed")

    # remove the temporary files
    tempFnames = [saFname, matchesBedFname, filtMatchesBedFname]
    if not DEBUG:
        for tfn in tempFnames:
            if isfile(tfn):
                os.remove(tfn)
    return bedFname

def makeVariants(seq):
    " generate all possible variants of sequence at 1bp-distance"
    seqs = []
    for i in range(0, len(seq)):
        for l in "ACTG":
            if l==seq[i]:
                continue
            newSeq = seq[:i]+l+seq[i+1:]
            seqs.append((i, seq[i], l, newSeq))
    return seqs

def expandIupac(seq):
    """ expand all IUPAC characters to nucleotides, returns list.
    >>> expandIupac("NY")
    ['GC', 'GT', 'AC', 'AT', 'TC', 'TT', 'CC', 'CT']
    """
    # http://stackoverflow.com/questions/27551921/how-to-extend-ambiguous-dna-sequence
    d = {'A': 'A', 'C': 'C', 'B': 'CGT', 'D': 'AGT', 'G': 'G', \
        'H': 'ACT', 'K': 'GT', 'M': 'AC', 'N': 'GATC', 'S': 'CG', \
        'R': 'AG', 'T': 'T', 'W': 'AT', 'V': 'ACG', 'Y': 'CT', 'X': 'GATC'}
    seqs = []
    for i in product(*[d[j] for j in seq]):
       seqs.append("".join(i))
    return seqs

def writeBowtieSequences(inFaFname, outFname, pamPat):
    """ write the sequence and one-bp-distant-sequences + all possible PAM sequences to outFname
    Return dict querySeqId -> querySeq and a list of all
    possible PAMs, as nucleotide sequences (not IUPAC-patterns)
    """
    ofh = open(outFname, "w")
    outCount = 0
    inCount = 0
    guideSeqs = {} # 20mer guide sequences
    qSeqs = {} # 23mer query sequences for bowtie, produced by expanding guide sequences
    allPamSeqs = expandIupac(pamPat)
    for seqId, seq in parseFastaAsList(open(inFaFname)):
        inCount += 1
        guideSeqs[seqId] = seq
        for pamSeq in allPamSeqs:
            # the input sequence + the PAM
            newSeqId = "%s.%s" % (seqId, pamSeq)
            newFullSeq = seq+pamSeq
            ofh.write(">%s\n%s\n" % (newSeqId, newFullSeq))
            qSeqs[newSeqId] = newFullSeq

            # all one-bp mutations of the input sequence + the PAM
            for nPos, fromNucl, toNucl, newSeq in makeVariants(seq):
                newSeqId = "%s.%s.%d:%s>%s" % (seqId, pamSeq, nPos, fromNucl, toNucl)
                newFullSeq = newSeq+pamSeq
                ofh.write(">%s\n%s\n" % (newSeqId, newFullSeq))
                qSeqs[newSeqId] = newFullSeq
                outCount += 1
    ofh.close()
    logging.debug("Wrote %d variants+expandedPam of %d sequences to %s" % (outCount, inCount, outFname))
    return guideSeqs, qSeqs, allPamSeqs

def applyModifStr(seq, modifStrs, strand):
    """ bowtie: given a list of pos:toNucl>fromNucl and a seq, return the original seq.
    position is 0-based
    >>> applyModifStr("ACAATAAGACATAAACATATCGG", "14:T>A,21:A>G,22:C>G".split(","), "+")
    'ACAATAAGACATAATCATATCAC'
    """
    seq = list(seq)
    for modifStr in modifStrs:
        #logging.debug( modifStr)
        pos, toFromNucl = modifStr.split(":")
        fromNucl, toNucl = toFromNucl.split(">")
        pos = int(pos)
        if strand=="-":
            fromNucl = revComp(fromNucl)
        seq[pos] = fromNucl
    return "".join(seq)

def parseRefout(tmpDir, guideSeqs, pamLen):
    """ parse all .map file in tmpDir and return as list of chrom,start,end,strand,guideSeq,tSeq
    """
    fnames = glob.glob(join(tmpDir, "*.map"))

    # while parsing, make sure we keep only the hit with the lowest number of mismatches
    # to the guide. Saves time when parsing.
    posToHit = {}
    hitBestMismCount = {}
    for fname in fnames:
        for line in open(fname):
           # s20+.17:A>G     -       chr8    26869044        CCAGCACGTGCAAGGCCGGCTTC IIIIIIIIIIIIIIIIIIIIIII 7       4:C>G,13:T>G,15:C>G
           guideIdWithMod, strand, chrom, start, tSeq, weird, someScore, alnModifStr = \
               line.rstrip("\n").split("\t")

           guideId = guideIdWithMod.split(".")[0]
           modifParts = alnModifStr.split(",")
           if modifParts==['']:
               modifParts = []
           mismCount = len(modifParts)
           hitId = (guideId, chrom, start, strand)
           oldMismCount = hitBestMismCount.get(hitId, 9999)
           if mismCount < oldMismCount:
               hit = (mismCount, guideIdWithMod, strand, chrom, start, tSeq, modifParts)
               posToHit[hitId] = hit
               hitBestMismCount[hitId] = mismCount # thanks to github user mbsimonovic

    ret = []
    for guideId, hit in posToHit.items():
           mismCount, guideIdWithMod, strand, chrom, start, tSeq, modifParts = hit
           if strand=="-":
               tSeq = revComp(tSeq)
           guideId = guideIdWithMod.split(".")[0]
           guideSeq = guideSeqs[guideId]
           genomeSeq = applyModifStr(tSeq, modifParts, strand)
           start = int(start)
           bedRow = (guideId, chrom, start, start+GUIDELEN+pamLen, strand, guideSeq, genomeSeq)
           ret.append( bedRow )

    return ret

def getEditDist(str1, str2):
    """ return edit distance between two strings of equal length
    >>> getEditDist("HIHI", "HAHA")
    2
    """
    assert(len(str1)==len(str2))
    str1 = str1.upper()
    str2 = str2.upper()

    editDist = 0
    for c1, c2 in zip(str1, str2):
        if c1!=c2:
            editDist +=1
    return editDist

def findOfftargetsBowtie(queue, batchId, batchBase, faFname, genome, pamPat, bedFname):
    " align guides with pam in faFname to genome and write off-targets to bedFname "
    tmpDir = batchBase+".bowtie.tmp"
    os.mkdir(tmpDir)

    # make sure this directory gets removed, no matter what
    global tmpDirsDelExit
    tmpDirsDelExit.append(tmpDir)
    if not DEBUG:
        atexit.register(delTmpDirs)

    # write out the sequences for bowtie
    queue.startStep(batchId, "seqPrep", "preparing sequences")
    bwFaFname = abspath(join(tmpDir, "bowtieIn.fa"))
    guideSeqs, qSeqs, allPamSeqs = writeBowtieSequences(faFname, bwFaFname, pamPat)

    genomePath =  abspath(join(genomesDir, genome, genome))
    oldCwd = os.getcwd()

    # run bowtie
    queue.startStep(batchId, "bowtie", "aligning with bowtie")
    os.chdir(tmpDir) # bowtie writes to hardcoded output filenames with --refout
    # -v 3 = up to three mismatches
    # -y   = try hard
    # -t   = print time it took
    # -k   = output up to X alignments
    # -m   = do not output any hit if a read has more than X hits
    # --max = write all reads that exceed -m to this file
    # --refout = output in bowtie format, not SAM
    # --maxbts=2000 maximum number of backtracks
    # -p 4 = use four threads
    # --mm = use mmap
    maxOcc = MAXOCC # meaning in BWA: includes any PAM, in bowtie we have the PAM in the input sequence
    cmd = "$BIN/bowtie -e 1000 %(genomePath)s -f %(bwFaFname)s  -v 3 -y -t -k %(maxOcc)d -m %(maxOcc)d dummy --max tooManyHits.txt --mm --refout --maxbts=2000 -p 4" % locals()
    runCmd(cmd)
    os.chdir(oldCwd)

    queue.startStep(batchId, "parse", "parsing alignments")
    pamLen = len(pamPat)
    hits = parseRefout(tmpDir, guideSeqs, pamLen)

    queue.startStep(batchId, "scoreOts", "scoring off-targets")
    # make the list of alternative PAM sequences
    altPats = offtargetPams.get(pamPat, [])
    altPamSeqs = []
    for altPat in altPats:
        altPamSeqs.extend(expandIupac(altPat))

    # iterate over bowtie hits and write to a BED file with scores
    # if the hit looks OK (right PAM + score is high enough)
    tempBedPath = join(tmpDir, "bowtieHits.bed")
    tempFh = open(tempBedPath, "w")

    offTargets = {}
    isSaCas9 = pamIsSaCas9(pamPat)

    for guideIdWithMod, chrom, start, end, strand, _, tSeq in hits:
        guideId = guideIdWithMod.split(".")[0]
        guideSeq = guideSeqs[guideId]
        genomePamSeq = tSeq[-pamLen:]
        logging.debug( "PAM seq: %s of %s" % (genomePamSeq, tSeq))
        if genomePamSeq in altPamSeqs:
            minScore = ALTPAMMINSCORE
        elif genomePamSeq in allPamSeqs:
            minScore = MINSCORE
        else:
            logging.debug("Skipping off-target for %s: %s:%d-%d" % (guideId, chrom, start, end))
            continue

        logging.debug("off-target minScore = %f" % minScore )

        # check if this match passes the off-target score limit
        if pamIsCpf1(pamPat):
            otScore = 0.0
        else:
            tSeqNoPam = tSeq[:-pamLen]

            if isSaCas9:
                otScore = calcSaHitScore(guideSeq, tSeqNoPam)
            else:
                otScore = calcHitScore(guideSeq, tSeqNoPam)

            if otScore < minScore:
                logging.debug("off-target not accepted")
                continue

        editDist = getEditDist(guideSeq, tSeqNoPam)
        guideHitCount = 0
        guideId = guideId.split(".")[0] # full guide ID looks like s33+.0:A>T
        name = guideId+"|"+strand+"|"+str(editDist)+"|"+tSeq+"|"+str(guideHitCount)+"|"+str(otScore)
        row = [chrom, str(start), str(end), name]
        # this way of collecting the features will remove the duplicates
        otKey = (chrom, start, end, strand, guideId)
        logging.debug("off-target key is %s" % str(otKey))
        offTargets[ otKey ] = row

    for rowKey, row in offTargets.items():
        tempFh.write("\t".join(row))
        tempFh.write("\n")

    tempFh.flush()

    # create a tempfile which is moved over upon success
    # makes sure we do not leave behind a half-written file if
    # we crash later
    tmpFd, tmpAnnotOffsPath = tempfile.mkstemp(dir=tmpDir, prefix="annotOfftargets")
    tmpFh = open(tmpAnnotOffsPath, "w")

    # get name of file with genome locus names
    genomeDir = genomesDir # make local var
    segFname = "%(genomeDir)s/%(genome)s/%(genome)s.segments.bed" % locals()

    # annotate with genome locus names
    cmd = "$BIN/overlapSelect %(segFname)s %(tempBedPath)s stdout -mergeOutput -selectFmt=bed -inFmt=bed | cut -f1,2,3,4,8 > %(tmpAnnotOffsPath)s" % locals()
    runCmd(cmd)

    shutil.move(tmpAnnotOffsPath, bedFname)
    queue.startStep(batchId, "done", "Job completed")

    if DEBUG:
        logging.info("debug mode: Not deleting %s" % tmpDir)
    else:
        shutil.rmtree(tmpDir)

def processSubmission(faFname, genome, pamDesc, bedFname, batchBase, batchId, queue):
    """ search fasta file against genome, filter for pam matches and write to bedFName
    optionally write status updates to work queue. Remove faFname.
    """
    batchInfo = readBatchAsDict(batchId)

    if genome=="noGenome":
        posStr = "?"
    elif "batchName" in batchInfo and batchInfo["batchName"].count(":")==2: # chrom:start-end:strand
        posStr = batchInfo["batchName"]
    else:
        queue.startStep(batchId, "bwasw", "Searching genome for one 100% identical match to input sequence")
        posStr = findPerfectMatch(batchId)

    batchInfo["posStr"] = posStr

    if posStr!="?":
        # get a 100bp-extended version of the input seq
        chrom, start, end, strand = parsePos(posStr)
        extSeq = extendAndGetSeq(genome, chrom, start, end, strand, batchInfo["seq"])
        if extSeq is None:
            # this can only happen if there is a 100%-M match but small SNPs in it compared to the input sequence
            # so the extension of the input fails.
            # in this case, we also invalidate the position, as there was no perfect match and the user
            # has to do something to fix it
            batchInfo["posStr"] = "?"
        else:
            logging.debug("100pb-extended seq (len: %d) is: %s" % (len(extSeq), extSeq))
            batchInfo["extSeq"] = extSeq

    # must save the batch again, as otherwise display won't work, we need the position saved
    writeBatchAsDict(batchInfo, batchId)

    if doEffScoring:
        queue.startStep(batchId, "effScores", "Calculating guide efficiency scores")
        createBatchEffScoreTable(batchId, queue)

    if genome=="noGenome":
        # skip the off-target search entirely
        open(bedFname, "w") # create a 0-byte file to signal job completion
        queue.startStep(batchId, "done", "Job completed")
        return

    if useBowtie:
        findOfftargetsBowtie(queue, batchId, batchBase, faFname, genome, pamDesc, bedFname)
    else:
        findOfftargetsBwa(queue, batchId, batchBase, faFname, genome, pamDesc, bedFname)

    if not DEBUG:
        os.remove(faFname)

    return bedFname

def lineFileNext(fh):
    """
        parses tab-sep file with headers as field names
        yields collection.namedtuples
        strips "#"-prefix from header line
    """
    line1 = fh.readline()
    while line1.startswith("##"):
        line1 = fh.readline()
    line1 = line1.strip("\n").strip("#")
    headers = line1.split("\t")
    Record = namedtuple('tsvRec', headers)

    for line in fh:
        line = line.rstrip("\n")
        fields = line.split("\t")
        try:
            rec = Record(*fields)
        except Exception as msg:
            logging.error("Exception occured while parsing line, %s" % msg)
            logging.error("Filename %s" % fh.name)
            logging.error("Line was: %s" % repr(line))
            logging.error("Does number of fields match headers?")
            logging.error("Headers are: %s" % headers)
            #raise Exception("wrong field count in line %s" % line)
            continue
        # convert fields to correct data type
        yield rec

allGenomes = None

def readGenomes():
    " return list of all genomes supported "
    global allGenomes
    if allGenomes:
        return allGenomes
    genomes = {}

    myDir = dirname(__file__)
    genomesDir = join(myDir, "genomes")

    inFnames = []
    globalFname = join(genomesDir, "genomeInfo.all.tab")
    if isfile(globalFname):
        inFnames = [globalFname]
    else:
        for subDir in os.listdir(genomesDir):
            infoFname = join(genomesDir, subDir, "genomeInfo.tab")
            if isfile(infoFname):
                inFnames.append(infoFname)

    for infoFname in inFnames:
        for row in lineFileNext(open(infoFname)):
            # add a note to identify UCSC genomes
            if row.server.startswith("ucsc"):
                addStr="UCSC "
            else:
                addStr = ""
            genomes[row.name] = row.scientificName+" - "+row.genome+" - "+addStr+row.description

    genomes = list(genomes.items())
    genomes.sort(key=operator.itemgetter(1))
    allGenomes = genomes
    return allGenomes

def printOrgDropDown(lastorg, genomes):
    " print the organism drop down box. "
    print('<select id="genomeDropDown" class style="max-width:600px" name="org" tabindex="2">')
    print('<option ')
    if lastorg == "noGenome":
        print('selected ')
    print('value="noGenome">-- No Genome: no specificity, only cleavage efficiency scores (max. len 25kbp)</option>')

    for db, desc in genomes:
        print('<option ')
        if db == lastorg :
            print('selected ')
        print('value="%s">%s</option>' % (db, desc))

    print("</select>")
    #print ('''
      #<script type="text/javascript">
      #$("#genomeDropDown").ufd({maxWidth:350, listWidthFixed:false});
      #</script>''')
    print ('''<br>''')

    print ("""<script>
    $("#genomeDropDown").chosen();
    $(".chosen-choices li").css("background","red");
    </script>
    """)


def printPamDropDown(lastpam):

    print('<select style="float:left" name="pam" tabindex="3">')
    for key,value in pamDesc:
        print('<option ')
        if key == lastpam :
            print('selected ')
        print('value="%s">%s</option>' % (key, value))
    print("</select>")

def printForm(params):
    " print html input form "
    scriptName = basename(__file__)

    genomes = readGenomes()

    haveHuman = False
    for g in genomes:
        if g[0]=="hg19":
             haveHuman = True

    # The returned cookie is available in the os.environ dictionary
    cookies=http.cookies.SimpleCookie(os.environ.get('HTTP_COOKIE'))
    if "lastorg" in cookies and "lastseq" in cookies and "lastpam" in cookies:
       lastorg   = cookies['lastorg'].value
       lastseq   = cookies['lastseq'].value
       lastpam   = cookies['lastpam'].value
    else:
       if not haveHuman:
           global DEFAULTSEQ
           global DEFAULTORG
           DEFAULTSEQ = ALTSEQ
           DEFAULTORG = ALTORG
       lastorg = DEFAULTORG
       lastseq = DEFAULTSEQ
       lastpam = DEFAULTPAM

    # SerialCloner is sending us the sequence via a HTTP get parameter
    if "seq" in params:
        lastseq = params["seq"]
    if "org" in params:
        lastorg = params["org"]

    seqName = ""
    if "seqName" in params:
        seqName = params["seqName"]

    printTeforBodyStart()
    #print('''March 6 2023: Sorry, no CRISPOR on the new UCSC-based-server (with RS3 scores) today. Too many performance problems on the new server. We were able to renew the old server. Please use the <a href="http://37.187.154.234">old server</a> temporarily.''')
    #sys.exit(0)

    print("""
<form id="main-form" method="post" action="%s">

 <div style="text-align:left; margin-left: 10px">
 CRISPOR (<a href="https://academic.oup.com/nar/article/46/W1/W242/4995687">citation</a>) is a program that helps design, evaluate and clone guide sequences for the CRISPR/Cas9 system. <a target=_blank href="/manual/">CRISPOR Manual</a>

<br><i>July 18, 2024: The old server has been retired. The new Python3 server is still lacking the Najm 2018 saCas9 score. See <a href="doc/changes.html">Full list of changes</a></i><br>

 </div>

<div class="windowstep subpanel" style="width:40%%;">
    <div class="substep">
        <div class="title">
            Step 1
        </div>

        Planning a lentiviral gene knockout screen? Use <a href="crispor.py?libDesign=1">CRISPOR Batch</a><br>

        Sequence name (optional): <input type="text" name="name" size="20" value="%s"><br>

        Enter a single genomic sequence, &lt; %d bp, typically an exon
        <img src="%simage/info-small.png" title="CRISPOR conserves the lowercase and uppercase format of your sequence, allowing to highlight sequence features of interest such as ATG or STOP codons.<br>Avoid using cDNA sequences as input, CRISPR guides that straddle splice sites are unlikely to work.<br>You can paste a single >23bp sequence and even multiple sequences, separated by N characters." class="tooltipster">
    <br>
    <small><a href="javascript:clearInput()">Clear Box</a> - </small>
    <small><a href="javascript:resetToExample()">Reset to default</a></small>
    </div>

    <textarea tabindex="1" style="width:100%%" name="seq" rows="12"
              placeholder="Paste here the genomic - not a cDNA - sequence of the exon you want to target. The sequence has to include the PAM site for your enzyme of interest, e.g. NGG. Maximum size %d bp. If you only have a cDNA, please BLAST or BLAT the cDNA first to find the right exon sequence for CRISPOR.">%s</textarea>
      <small>Text case is preserved, e.g. you can mark ATGs with lowercase.<br>Instead of a sequence, you can paste a chromosome range, e.g. chr1:11,130,540-11,130,751</small>
</div>
<div class="windowstep subpanel" style="width:50%%">
    <div class="substep" style="margin-bottom: 1px">
        <div class="title" style="cursor:pointer;" onclick="$('#helpstep2').toggle('fast')">
            Step 2
        </div>
        Select a genome
    </div>
    """% (scriptName, seqName, MAXSEQLEN, HTMLPREFIX, MAXSEQLEN, lastseq))

    printOrgDropDown(lastorg, genomes)
    print("""
    <div id="trackHubNote" style="margin-bottom:5px">
    <small>Note: pre-calculated exonic guides for this species are on the <a id='hgTracksLink' target=_blank href="">UCSC Genome Browser</a>.</small>
    </div>
    """)
    print('<small style="float:left">We have %d genomes, but not yours? Search <a href="https://www.ncbi.nlm.nih.gov/assembly">NCBI assembly</a> and send a GCF_/GCA_ ID to <a href="mailto:%s">CRISPOR support</a>.</small>' % (len(genomes), contactEmail))
    print("""
    </div>
    <div class="windowstep subpanel" style="width:50%%; height:158px">
    <div class="substep">
    <div class="title" style="cursor:pointer;" onclick="$('#helpstep3').toggle('fast')">
        Step 3
        <img src="%simage/info-small.png" title="The most common system uses the NGG PAM recognized by Cas9 from S. <i>pyogenes</i>. The VRER and VQR mutants were described by <a href='http://www.nature.com/nature/journal/vaop/ncurrent/abs/nature14592.html' target='_blank'>Kleinstiver et al</a>, Cas9-HF1 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4851738/'>Kleinstiver 2016</a>, eSpCas1.1 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4714946/'>Slaymaker 2016</a>, Cpf1 by <a href='http://www.cell.com/abstract/S0092-8674(15)01200-3'>Zetsche 2015</a>, SaCas9 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/25830891/'>Ran 2015</a> and KKH-SaCas9 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/26524662/'>Kleinstiver 2015</a>, modified As-Cpf1s by <a href='http://biorxiv.org/content/early/2016/12/04/091611'>Gao et al. 2017</a>." class="tooltipsterInteract">
        </div>
        Select a Protospacer Adjacent Motif (PAM)
    </div>
    """ % HTMLPREFIX)

    printPamDropDown(lastpam)

    print("""<br>See <a target=_blank href="manual/manual.html#enzymes">notes on enzymes</a> in the manual.<br>""")

    print("""
    <div style="width:40%; margin-top: 10px; margin-left:50px; text-align:center; display:block">
    <input type="submit" name="submit" value="SUBMIT" tabindex="4"/>
    </div>
    </div>
    """)
    print("""
<script>
/* set the dropbox to hg19 and paste the example sequence into the input box. */
function resetToExample() {
    $("textarea[name='seq']").val("%s");
    $("#genomeDropDown").val("%s");
    $("select[name='pam']").val("NGG");
    }

/* clear the sequence input box */
function clearInput() {
    $("textarea[name='seq']").val("");
    }

</script>
<script>
    /* hide the track hub note if genome is not hg19 */
    ucscTrackDbs=['hg19', 'hg38', 'rn5', 'mm10', 'mm9', 'ci2', 'danRer7', 'sacCer3', 'dm6'];
    function showHideHubNote() {
        var valSel = $("#genomeDropDown").val();
        if (jQuery.inArray(valSel, ucscTrackDbs)!=-1)
            {
            $("#trackHubNote").css('visibility', 'visible');
            $("#hgTracksLink").attr("href", "http://genome.ucsc.edu/cgi-bin/hgTracks?db="+valSel+"&crispr=show");
            }
        else
            $("#trackHubNote").css('visibility', 'hidden');
    }
    $("#genomeDropDown").on('change', showHideHubNote);
    showHideHubNote();
</script>

</form>
    """ % (DEFAULTSEQ, DEFAULTORG))

def readBatchAsDict(batchId):
    " return contents of batch as a dictionary or None "
    batchBase = join(batchDir, batchId)
    jsonFname = batchBase+".json"
    if isfile(jsonFname):
        params = json.load(open(jsonFname))
    else:
        db = sqlite3.connect(batchArchive)
        c = db.cursor()
        c.execute("select data from jobArchive where id=?", (batchId,))
        data = None
        for row in c.fetchall():
            data = row[0]
        db.close()

        if data is None:
            return None

        jsonStr = gzip.decompress(data).decode("utf8")
        params = json.loads(jsonStr)

    if "batchName" in params:
        global batchName
        batchName = params["batchName"]

    return params

def writeBatchAsDict(batchInfo, batchId):
    batchBase = join(batchDir, batchId)
    tmpFname = batchBase+".json.tmp"

    ofh = open(tmpFname, "w")
    json.dump(batchInfo, ofh)
    ofh.close()

    jsonFname = batchBase+".json"
    os.rename(tmpFname, jsonFname)
    logging.debug("Wrote batch info to %s: %s" % (jsonFname, batchInfo))

def readBatchParams(batchId):
    """ given a batchId, return the genome, the pam, the input sequence and the
    chrom pos and extSeq, a 100bp-extended version of the input sequence.
    Returns None for pos if not found. """
    params = readBatchAsDict(batchId)
    if params != None:
        return params["seq"], params["org"], params["pam"], params.get("posStr"), params.get("extSeq")

    # FROM HERE UP TO END OF FUNCTION: legacy cold for old batches pre-end-2016 (no json files back then)
    # remove in 2017
    batchBase = join(batchDir, batchId)
    inputFaFname = batchBase+".input.fa"
    if not isfile(inputFaFname):
        errAbort('Could not find the batch %s. We cannot keep Crispor runs for more than '
                'a few months. Please resubmit your input sequence via'
            ' <a href="crispor.py">the query input form</a>' % batchId)

    ifh = open(inputFaFname, encoding="utf8")
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

    return inSeq, genome, pamSeq, position, extSeq

def gzipStr(s):
    " compress a string with gzip and return "
    out = StringIO()
    with gzip.GzipFile(fileobj=out, mode="w") as f:
         f.write(s)
    return out.getvalue()

def gunzipStr(s):
    " uncompress a string with gzip and return "
    print(len(s), type(s), dir(s))
    f = gzip.GzipFile(StringIO(s))
    result = f.read()
    f.close()
    return result

def openDbm(dbFname, mode):
    " some distributions don't include the dbm module anymore "
    #import dbm.ndbm
    #dbMod = dbm
    #import dbm.gnu
    #dbMod = gdbm
    #import semidbm
    # lmdbm is faster than everything else: https://pypi.org/project/lmdbm/
    # though semidbm is not bad either
    # Also see leveldb Wiki page
    from lmdbm import Lmdb
    db = Lmdb.open(dbFname+".lmdb", mode)
    return db

def saveOutcomeData(batchId, data):
    """ save outcome data of batch. data is a dictionary with key = score name """
    batchBase = join(batchDir, batchId)
    dbFname = batchBase
    db = openDbm(dbFname, "c")

    #conn = sqlite3.connect(dbFname, "w")
    #c = conn.cursor()
    #c.execute('''CREATE TABLE outcomes (id text PRIMARY KEY, data blob))''' % scoreName)
    #c.commit()

    for scoreName, data in data.items():
        #c.execute("INSERT INTO outcomes values (?, ?)", (scoreName, gzipStr(json.dumps(data))))
        db[scoreName] = zlib.compress(json.dumps(data).encode("utf8"))

    db.close()

    #c.commit()

def readOutcomeData(batchId, scoreName):
    """ open outcome data of batch, key is score name """
    batchBase = join(batchDir, batchId)
    #conn = sqlite3.connect(dbFname, "r")
    #c = conn.cursor()
    #binData = c.execute("SELECT data FROM outcomes where id=?", scoreName)
    #try:
    #    import dbm.ndbm
    #    db = dbm.ndbm.open(batchBase, "r") # dbm always adds .db to the file name
    #except:
    #    # old batches on crispor.org are still using gdbm
    #    dbFname = batchBase+".gdbm"
    #    import dbm.gnu
    #    db = dbm.gnu.open(dbFname, "r")
    db = openDbm(dbFname, "r")
    dbObj = db[scoreName]
    jsonStr = zlib.decompress(dbObj)
    data = json.loads(jsonStr)
    db.close()
    return data

def findAllPams(seq, pam):
    """ find all matches for PAM and return as dict startPos -> strand and a set
    of end positions. The start positions for the negative strand are for the
    rev-complemented PAM
    """
    seq = seq.upper()
    startDict, endSet = findPams(seq, pam, "+", {}, set())
    startDict, endSet = findPams(seq, revComp(pam), "-", startDict, endSet)

    if pam in multiPams:
        for pam2 in multiPams[pam]:
            startDict, endSet = findPams(seq, pam2, "+", startDict, endSet)
            startDict, endSet = findPams(seq, revComp(pam2), "-", startDict, endSet)

    return startDict, endSet

def newBatch(batchName, seq, org, pam):
    """ obtain a batch ID and write seq/org/pam to their files.
    Return batchId.
    """
    batchId = makeTempBase(seq, org, pam, batchName)

    batchData = {}
    batchData["org"] = org
    batchData["pam"] = pam
    batchData["batchName"] = batchName
    batchData["seq"] = seq
    batchData["posStr"] = ""

    writeBatchAsDict(batchData, batchId)
    return batchId

def readDbInfo(org):
    " return a dbInfo object with the columsn in the genomeInfo.tab file "
    myDir = dirname(__file__)
    genomesDir = join(myDir, "genomes")
    infoFname = join(genomesDir, org, "genomeInfo.tab")
    if not isfile(infoFname):
        return None
    dbInfo = next(lineFileNext(open(infoFname)))
    return dbInfo

def printQueryNotFoundNote(dbInfo):
    print("<div class='title'>Query sequence, not found in the selected genome, %s (%s)</div>" % (dbInfo.scientificName, dbInfo.name))
    print("<div class='substep' style='border: 1px black solid; padding:5px; background-color: aliceblue'>")
    print("<strong>Warning:</strong> The query sequence was not found in the selected genome.")
    print("This can be a valid query, e.g. a GFP sequence.<br>")
    print("If not, you might want to check if you selected the right genome for your query sequence.<br>")
    print("Use a tool like <a target=_blank href='http://genome.ucsc.edu/cgi-bin/hgBlat'>BLAT</a> to check if the " \
        "sequence really has a 100% identical match in the target genome.<p>")
    print("When reading the list of guide sequences and off-targets below, bear in mind that in case that the input sequence is really in the genome and just has a few differences, the software will use the first found match as the on-target as it cannot distinguish 0-mismatch off-targets from 0-mismatch on-targets. In this case, the specificity scores of guide sequences are too low. In other words, some guides may be fine, the problem may just be that the on-target is shown as an off-target. <br>")
    print("Because there is no flanking sequence available, the guides in your sequence that are within 50bp of the ends will have no efficiency scores. The efficiency scores will instead be shown as '--'. Include more flanking sequence > 50bp to obtain these scores.")
    print("</div>")

def getOfftargets(seq, org, pamDesc, batchId, startDict, queue):
    """ write guides to fasta and run bwa or use cached results.
    Return name of the BED file with the matches or None if not yet available.
    Write progress status updates to queue object.
    """
    pam = setupPamInfo(pamDesc)
    assert('-' not in pam)

    batchBase = join(batchDir, batchId)
    otBedFname = batchBase+".bed.gz"

    batchInfo = readBatchAsDict(batchId)

    flagFile = batchBase+".running"
    if isfile(flagFile):
       errAbort("This sequence is still being processed. Please wait for ~20 seconds "
           "and try again, e.g. by reloading this page. If you see this message for "
           "more than 2-3 minutes, please send an email to %s. Thanks!" % contactEmail)

    if not batchInfo or not isfile(otBedFname) or commandLineMode or not "posStr" in batchInfo or \
            (batchInfo["posStr"]=="" and not batchInfo["org"]=="noGenome"): # pre-4.8 batches don't have a posStr at all
        # write potential PAM sites to file
        faFname = batchBase+".fa"
        writePamFlank(seq, startDict, pam, faFname)
        if commandLineMode:
            processSubmission(faFname, org, pamDesc, otBedFname, batchBase, batchId, queue)
        else:
            q = JobQueue()
            q.openSqlite()
            ip = os.environ.get("REMOTE_ADDR", "noIp")

            if ip=="195.176.112.240":
                errAbort("IP address blocked.")

            wasOk = q.addJob("search", batchId, "ip=%s,org=%s,pam=%s" % (ip, org, pamDesc))
            if not wasOk:
                print("CRISPOR job %s failed-running..." % batchId)
                pass
            q.close()
            return None

    return otBedFname

def showPamWarning(pam):
    if pamIsCpf1(pam):
        print('<div style="text-align:left; border: 1px solid; background-color: aliceblue; padding: 3px">')
        print("<strong>Note:</strong> You are using the Cpf1 enzyme or related enzyme.")
        print("While there is an efficiency score specificially for Cpf1, there is no off-target ranking algorithm available in the literature, to our knowledge. We use Hsu and CFD scores below for off-target ranking, but they were developed for spCas9. There is not enough data yet to support their usefulness for Cpf1. Contact us for more info if you need to rank Cpf1 off-targets for validation or if you have a dataset that could elucidate this question. We are showing out-of-frame scores, but they are based on micro-homology that assumes a spCas9 cut site, so most likely the out-of-frame scores are not accurate for the staggered cut of Cpf1 either.")
        print('</div>')
    #elif pamIsSaCas9(pam):
        #print '<div style="text-align:left; border: 1px solid; background-color: aliceblue; padding: 3px">'
        #print "<strong>Note:</strong> Your query is using a Cas9 from S. aureus.<br>"
        #print "Please note that while the efficiency scoring was built for saCas9, the off-target ranking below and specificity scores are based on CFD/Hsu models, which were developed for spCas9. The ranking of off-targets could be very inaccurate. If you have a saCas9 off-target dataset, you can contact us for further info, we are only aware of the BLESS dataset by <a href='https://www.nature.com/articles/nature14299' target=_blank>Ran et al. 2015</a>.<br>As for out-of-frame and micro-homology, this model is also based on spCas9, but <a target=_blank href='https://www.nature.com/articles/nature14299'>Ran et al 2015</a> showed that the saCas9 cleavage pattern looks identical to spCas9's, so the OOF micro-homology model should work with saCas9."
        #print '</div>'
    elif not pamIsSpCas9(pam) and not pamIsSaCas9(pam):
        print('<div style="text-align:left; border: 1px solid; background-color: aliceblue; padding: 3px">')
        print("<strong>Warning:</strong> Your query involves a Cas9 that is not from S. Pyogenes and is also not Cpf1 nor saCas9.")
        print("Please bear in mind that specificity and efficiency scores were designed using data with S. Pyogenes Cas9 and will very likely not be applicable to this particular Cas9. There is nothing we can do about this, we are unaware of a published dataset for this enzyme. If you know one, please contact us. Also contact us if you think another one of the existing scoring model would be more appropriate for this enzyme.<br>")
        print('</div>')

    if pam=="NNNNACA":
        printNote("You selected the old version of the CjCas9 PAM. You may want to select the more recent "+
        "PAMs from the menu on the first page, based on the study by "+
        "<a target=_blank href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5473640/'>Kim et al 2017</a>.")

    if pam=="NGN":
        printNote("You have selected the NGN pam for xCas9. While this PAM is documented to work, "+
        "if you read the paper in detail, you will notice that the editing efficiency is much lower. "+
        "For optimal efficiency, consider going back and switching to the 'high-efficiency' xCas9 PAM.")

    if pam=="NGK":
        printNote("You have selected the most efficient PAM for xCas9. You can also select the more general/flexible"+
        " NGN PAM from the menu when you submit your job. If you read the xCas9 paper in detail, you will find "+
        "that NGN is not as efficient though.")

    if pam=="TTTN":
        printWarning("You selected TTTN as the PAM for Cpf1. " +
            "This is not the best PAM. The actual PAM is TTTV, as shown in Fig. 2a of " +
            "<a href='https://www.ncbi.nlm.nih.gov/pubmed/27992409'>Kim HK et al. Nat Meth 2017</a>.<br>")

def showNoGenomeWarning(dbInfo):
    if dbInfo==None:
        printNote('As there is no genome that can be used to get flanking sequence for your sequence, efficiency scores 50bp from the start or the end of your sequence cannot be calculated and are shown as "--". If needed, extend the input sequence and retry.')

def getSeq(db, posStr):
    """
    given a database name and a string with the position as chrom:start-end, return the sequence as
    a string.
    """
    chrom, start, end, strand =  parsePos(posStr)

    if end-start > MAXSEQLEN and db!="noGenome":
        errAbort("Input sequence range too long. Please retry with a sequence range shorter than %d bp." % MAXSEQLEN)
    genomeDir = genomesDir # pull in global var
    twoBitFname = getTwoBitFname(db)
    binPath = join(binDir, "twoBitToFa")

    chromSizes = parseChromSizes(db)
    if chrom not in chromSizes:
        errAbort("Sorry, the chromosome '%s' is not valid in the genome %s. Check upper/lowercase, e.g. for most mammalian genomes, " \
            "it is chrX not chrx, and chr1, not Chr1." % (cgi.escape(chrom), db))
    if start<0 or end<0 or start>chromSizes[chrom] or end>chromSizes[chrom]:
        errAbort("Sorry, the coordinates '%d-%d' are not valid in the genome %s. Coordinates must not be outside chromosome boundaries or less than 0." % (start, end, db))

    cmd = [binPath, twoBitFname, "-seq="+chrom, "-start="+str(start), "-end="+str(end), "stdout"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    seqStr = proc.stdout.read()

    retCode = proc.wait()
    if retCode!=0:
        errAbort("Error on sequence retrieval. This looks like a bug. Please contact us and tell us the input and genome, we will fix this error.")

    # remove fasta header line
    lines = seqStr.decode("utf8").splitlines()
    if len(lines)>0:
        lines.pop(0)
    seq = "".join(lines)
    if len(seq) < 23:
        errAbort("Sorry, the sequence range %s on genome %s is not longer than 23bp. To find a valid CRISPR/Cas9 site, one needs at least a 23bp long sequence." % (db, posStr))

    if strand=="-":
        seq = revComp(seq)

    return seq

def printStatus(batchId, msg):
    " print status, not using any Ajax "
    q = JobQueue()
    q.openSqlite()
    status = q.getStatus(batchId)
    q.close()

    errorState = False

    if "Traceback" in status:
        print("<!--")
        print(status)
        print("-->")
        status = "An error occured during the processing.<br> Please send an email to %s and tell us that the failing batchId was %s.<br>We can usually fix this quickly. Thanks! <br>If you submit the same sequence/genome/name again, it will not be re-run, Crispor will pickup the old error. We will have to reset it before you can resubmit this particular sequence, so you will have to contact us or change the sequence to get a new job into the system." % (contactEmail, batchId)
        errorState = True
    else:
        print('<meta http-equiv="refresh" content="10" >')
        if len(msg)!=0:
            print((msg+"<p>"))
        print("CRISPOR job has been submitted.<p>")

    if status==None:
        status = "Batch completed. Refresh page to show results."

    print(("Job Status: <tt>%s</tt><p>" % status))

    if not errorState:
        print("<p><small>This page will refresh every 10 seconds</small><br>")
        print(("<p><small>If you see this message for longer than 5 minutes, please <a href='mailto:%s'>contact us</a>." % contactEmail))

def readVarDbs(db):
    """ find all possible variant VCFs and return as list of (shortLabel, fname, label, hasAF)
    hasAF = file has the AF field (allele frequency). Means that the UI
    will show the "frequency filter" button.
    """
    # parse the descriptions of the VCF files
    # descriptions are optional
    labelFname = join(genomesDir, db, "vcfDescs.txt")
    ret = []
    if isfile(labelFname):
        for line in open(labelFname):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields)==4:
                shortLabel, fname, desc, hasAF = fields
            else:
                errAbort("not four fields in vcfDescs.txt: %s" % fields)

            fpath = join(genomesDir, db, fname)
            if not isfile(fpath):
                print("Error: Cannot find VCF file %s" % fpath)
                continue
            hasAF = (hasAF=="1")
            ret.append( (shortLabel, fname, desc, hasAF) )
    return ret

def parseVcfInfo(info):
    " parse a VCF info string and return as a dict key=val "
    fs = info.split(";")
    ret = {}
    for field in fs:
        parts = field.split("=", maxsplit=1)
        if len(parts)==2:
            key, val = parts
        elif len(parts)==1:
            key = parts[0]
            val = True
        ret[key] = val
    return ret

def findVariantsInRange(vcfFname, chrom, start, end, strand, minFreq):
    """ find variants that overlap the position.
    varDb is a tuple of label, vcfFname
    return as a dict relative position -> (chrom, pos, refAllele, altAllele, list of info-dicts)
    special position is "label" which is the label of the variant db
    """
    minFreq = float(minFreq)
    seqLen = end-start
    if not isfile(vcfFname):
        errAbort("%s not found" % vcfFname)
    tb = tabix.open(vcfFname)
    chrom = chrom.replace("chr","")
    try:
        records = tb.query(chrom, start+1, end) # VCF is 1-based
    except tabix.TabixError:
        sys.stderr.write("Chromosome in query does not exist in VCF file? chrom: %s, VCF file: %s\n" % (chrom, vcfFname))
        records = []


    varDict = defaultdict(list)
    for rec in records:
        chrom, varPos, varId, refAll, altAllStrList, qual, filterFlag, info = rec[:8]
        infoDict = parseVcfInfo(info)
        altAllList = altAllStrList.split(",")
        if "AF" in infoDict:
            afList = infoDict["AF"].split(",")
        else:
            afList = [None] * len(altAllList)

        for altAll, allFreq in zip(altAllList, afList):
            # 1000 genomes had <CN0> at AAAAATTTTTAAAAATTAGCTGG
            # no idea what this is supposed to represent. issue #7
            if "<" in altAll:
                continue
            if minFreq is not None and allFreq is not None:
                allFreq = float(allFreq)
                if not allFreq > minFreq:
                    continue

            attribs = {}
            #afList = infoDict["AF"].split(",")
            #altAllList = altAll.split(",")
            #newAltAllList = []
            #newAfList = []
            #for af, altAll in zip(afList, altAllList):
                #afNum = float(af)
                #if afNum < minFreq:
                    #continue
                #newAltAllList.append(altAll)
                #newAfList.append(af)
            ##if len(newAltAllList)==0:
                #continue
            #altAll = ",".join(newAltAllList)
            #infoDict["AF"] = ",".join(newAfList)
            if allFreq!=None:
                attribs["freq"] = allFreq
            relPos = int(varPos)-1-start
            if strand=="-":
                relPos = seqLen - relPos - len(refAll)
                refAll = revComp(refAll)
                altAlls = []
                for altAll in altAll.split(","):
                    altAlls.append(revComp(altAll))
                altAll = ",".join(altAlls)
            if varId != ".":
                attribs["varId"] = varId
            varInfo = (chrom, varPos, refAll, altAll, attribs)
            varDict[relPos].append(varInfo)

    return varDict

def showSeqDownloadMenu(batchId):
    " show a little dropdown menu so user can get annotated sequence in genbank format "
    print("""<div style="padding-top:4px"><small>Download for: """)

    htmls = []

    baseUrl = "crispor.py?batchId=%s" % batchId

    myUrl = baseUrl+"&download=serialcloner"
    html = "<a href='%s'>SerialCloner</a> (<a target=_blank href='http://serialbasics.free.fr/Serial_Cloner-Download.html'>free</a>)" % myUrl
    htmls.append(html)

    myUrl = baseUrl+"&download=ape"
    html = '<a href="%s">ApE</a> (<a target=_blank href="http://biologylabs.utah.edu/jorgensen/wayned/ape/">free</a>)' % myUrl
    htmls.append(html)

    myUrl = "http://crispor.tefor.net/crispor.py?batchId=%s&download=genomecompiler" % batchId
    #backUrl = "https://designer.genomecompiler.com/plasmid_iframe?file_url=%s#/plasmid" % urllib.quote(myUrl)
    backUrl = "https://designer.genomecompiler.com/plasmid_iframe?file_url=%s#/plasmid" % urllib.parse.quote(myUrl)
    html = "<a target=_blank href='%s'>GenomeCompiler</a>" % backUrl
    htmls.append(html)

    myUrl = baseUrl+"&download=benchling"
    html = "<a href='%s'>Benchling</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl+"&download=snapgene"
    html = "<a href='%s'>SnapGene</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl+"&download=geneious"
    html = "<a href='%s'>Geneious</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl+"&download=vnti"
    html = "<a href='%s'>Vector NTI</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl+"&download=lasergene"
    html = "<a href='%s'>LaserGene</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl+"&download=genbank"
    html = "<a href='%s'>Genbank</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl+"&download=fasta"
    html = "<a href='%s'>FASTA</a>" % myUrl
    htmls.append(html)

    #html = "<div id='copyLink' data-clipboard-target='#seqAsText'>Copy sequence to clipboard</div>"
    #htmls.append(html)

    print(" - ".join(htmls))

    print("</small></div>")

def mapToGenome(seqStart, seqStrand, pamStart, guideStart, guideStrand):
    if pamIsFirst:
        # thick part = PAM comes first
        chromStart = seqStart+pamStart
        thickStart = seqStart+guideStart
        thickEnd = thickStart+GUIDELEN
        chromEnd = thickEnd
    else:
        chromStart = seqStart+guideStart
        thickStart = chromStart
        thickEnd   = chromStart+GUIDELEN
        chromEnd   = thickEnd

    strands = seqStrand+guideStrand
    chromStrand = "+"
    if strands=='+-' or strands=='-+':
        chromStrand = "-"

    return chromStart, chromEnd, thickStart, thickEnd, chromStrand

def makeCustomTrack(org, chrom, seqStart, seqEnd, seqStrand, guideData, batchId, batchName):
    " create a custom track file for a given batch and return the filename "
    ctDir = join(batchDir, "customTracks")
    if not isdir(ctDir):
        os.makedirs(ctDir)

    ctFname = join(ctDir, batchId+".bed") # temporary bed file
    bbFname = join(ctDir, batchId+".bb")  # bigBed file
    ctFname = join(ctDir, batchId+".txt") # custom track settings
    if isfile(ctFname):
        return bbFname

    seqStart = int(seqStart)
    seqEnd = int(seqEnd)

    rows = []
    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, pamStart, guideStart, guideStrand, pamId, guideSeq, \
            pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, repCount = guideRow

        rgb = hexToRgb(scoreToColor(guideScore)[0])
        chromStart, chromEnd, thickStart, thickEnd, chromStrand = mapToGenome(seqStart, seqStrand, pamStart, guideStart, guideStrand)

        mitScore = str(guideScore)
        fusiScore = str(effScores.get("fusi", -1))
        crisprScanScore = str(effScores.get("crisprScan", -1))
        oofScore = str(effScores.get("oof", -1))

        name = pamId
        bed = [chrom, chromStart, chromEnd, name, mitScore, chromStrand, thickStart, thickEnd, rgb, guideSeq, pamSeq, mitScore, fusiScore, crisprScanScore, oofScore, batchId]
        rows.append(bed)

    # sort and write to file
    rows.sort()
    ofh = open(ctFname, "w")
    for row in rows:
        row = [str(x) for x in row]
        ofh.write("\t".join(row))
        ofh.write("\n")
    ofh.close()

    sizeFname = getSizeFname(org)
    asFname = "crispor.as"
    cmd = ["$BIN/bedToBigBed", "-type=bed9+", "-tab", "-as="+asFname, ctFname, sizeFname, bbFname]
    runCmd(cmd, useShell=False)

    bbUrl = baseUrl+"/%s.bb" % batchId

    ofh = open(ctFname, "w")
    if batchName=="":
        batchName="Results"
    ofh.write("browser position %s:%d-%d\n" % (chrom, seqStart, seqEnd))
    ofh.write('track type=bigBed name="CRISPOR %(batchName)s" description="CRISPOR Results %(batchName)s %(batchId)s" bigDataUrl=%(bbUrl)s itemRgb=On visibility=pack\n' % locals())
    ofh.close()

    #hubFname = join(ctDir, batchId+".txt")
    #ofh = open(hubFname, "w")
    #ofh.write("hub CRISPOR\n")
    #ofh.write("shortLabel CRISPOR %s\n")
    #ofh.write("longLabel CRISPOR batch %s %s\n")
    #ofh.write("genomesFile %s\n" % hubFname)
    #ofh.write("email crispor@tefor.net\n")
    #ofh.write("descriptionUrl http://crispor.org\n")
    #ofh.write("genome %s\n" % org)
    #ofh.write("trackDb %s\n" % hubFname)
    #ofh.write("\n" % hubFname)
    #ofh.write("track crispor%s\n" % batchId)
    #ofh.write("shortLabel M-CAP
    #longLabel M-CAP
    #group genes
    #visibility dense
    #type bigBed 9 +
    ##os.remove(ctFname)

    #return bbFname
    ctUrl = baseUrl+"/%s.txt" % batchId
    return ctUrl

def iterBbLines(bbPath, chrom, start, end, strand):
    " yield bigGenePred rows from bigBed that overlap pos "
    binPath = join(binDir, "bigBedToBed")
    cmd = [binPath, bbPath, "stdout", "-chrom="+chrom, "-start="+str(start), "-end="+str(end)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, encoding="latin1")
    for line in proc.stdout:
        yield line.split("\t")

def trimExonAndFlip(exStart, exEnd, exStrand, seqLen, seqStrand):
    """ Put the exon into the current sequence window:
    - trim exon to the window (0, seqLen), return None if completely outside the view.
    - reverse the exon coordinates if seqStrand=="-"
    """
    if exStart < 0:
        if exEnd < 0:
            # the whole exon is outside the view on the left side
            return None, None, None
        else:
            # truncate the exon to start at 0
            exStart = 0
    if exEnd > seqLen:
        if exStart > seqLen:
            # the whole exon is outside the view on the right side
            return None, None, None
        else:
            # truncate the end
            exEnd = seqLen

    if seqStrand=="-":
        oldExEnd = exEnd
        exEnd = seqLen - exStart
        exStart = seqLen - oldExEnd
        # inputSeq forw and transcript forw -> exon is forw
        # inputSeq forw and transcript rev -> exon is rev
        # inputSeq rev and transcript forw -> exon is rev
        # inputSeq rev and transcript rev -> exon is forw
        if exStrand=="+":
            exStrand = "-"
        else:
            exStrand = "+"

    return exStart, exEnd, exStrand

def getExonInfo(org, geneName, position):
    """ retrieve exon info between position, return format transId -> (exNumber, start, end, exFrame, nextExonFrame, strand)
    - start and end are relative to position!
    - exNumber starts at 0
    - nextExonFrame is relative to strand: for a transcript on the - strand, it is for the more 5' exon
    """
    # bigGenePred format:
    #string chrom;       "Reference sequence chromosome or sca
    #uint   chromStart;  "Start position in chromosome"
    #uint   chromEnd;    "End position in chromosome"
    #string name;        "Name or ID of item, ideally both hum
    #uint score;         "Score (0-1000)"
    #char[1] strand;     "+ or - for strand"
    #uint thickStart;    "Start of where display should be thi
    #uint thickEnd;      "End of where display should be thick
    #uint reserved;       "RGB value (use R,G,B string in inpu
    #int blockCount;     "Number of blocks"
    #int[blockCount] blockSizes; "Comma separated list of bloc
    #int[blockCount] chromStarts; "Start positions relative to
    #string name2;       "Alternative/human readable name"
    #string cdsStartStat; "enum('none','unk','incmpl','cmpl')"
    #string cdsEndStat;   "enum('none','unk','incmpl','cmpl')"
    #int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame
    #string type;        "Transcript type"
    #string geneName;    "Primary identifier for gene"
    #string geneName2;   "Alternative/human readable gene name
    #string geneType;    "Gene type"

    ret = defaultdict(list)
    seqChrom, seqStart, seqEnd, seqStrand = parsePos(position)
    seqLen = seqEnd - seqStart

    fname = join(genomesDir, org, geneName+".bb")

    maxIdLen = 0

    for row in iterBbLines(fname, seqChrom, seqStart, seqEnd, seqStrand):
        chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, reserved, \
            blockCount, blockSizes, blockStarts, name2, cdsStartStat, cdsEndStat, exonFrames, \
            tType, geneName, geneName2, geneType = row

        chromStart = int(chromStart)
        chromEnd = int(chromEnd)
        thickStart = int(thickStart)
        thickEnd = int(thickEnd)

        blockSizes = [int(x) for x in blockSizes.split(",") if x!='']
        blockStarts = [int(x) for x in blockStarts.split(",") if x!='']
        exonFrames = [int(x) for x in exonFrames.split(",") if x!='']
        assert(len(blockSizes)==len(blockStarts)==len(exonFrames))

        symbol = ""
        if geneName2!="":
            symbol = geneName2

        blocks = list(zip(blockSizes, blockStarts, exonFrames))
        for exIdx, (blockSize, blockStart, exonFrame) in enumerate(blocks):
            exChromStart=chromStart+blockStart
            exChromEnd = exChromStart+blockSize
            exStrand = strand

            nextFrame = None
            if strand=="+":
                if exIdx+1 < len(blocks):
                    nextFrame = blocks[exIdx+1][-1]
            else:
                if exIdx > 0:
                    nextFrame = blocks[exIdx-1][-1]
            #print("next frame", nextFrame, "<br>")

            # figure out exon start/end: special case for UTRs: trim down to CDS start/end
            if exChromStart < thickStart < exChromEnd:
                exStart = thickStart-seqStart
                # add the UTR as a special exon
                utrStart, utrEnd, utrStrand  = trimExonAndFlip(exChromStart-seqStart, exStart, exStrand, seqLen, seqStrand)
                if utrStart!=None:
                    ret[(name, symbol)].append((-1, utrStart, utrEnd, -1, nextFrame, utrStrand))
            else:
                exStart = chromStart+blockStart-seqStart

            if exChromStart < thickEnd < exChromEnd:
                exEnd = thickEnd-seqStart
                # add the UTR as a special exon
                utrStart, utrEnd, utrStrand = trimExonAndFlip(exEnd, exChromEnd-seqStart, exStrand, seqLen, seqStrand)
                if utrStart!=None:
                    ret[(name, symbol)].append((-1, utrStart, utrEnd, -1, nextFrame, utrStrand))
            else:
                exEnd = chromStart+blockStart+blockSize-seqStart


            #print chromStart, chromStart+blockSize, seqStart, "<br>"
            exStart, exEnd, exStrand = trimExonAndFlip(exStart, exEnd, exStrand, seqLen, seqStrand)
            if exStart==None:
                # whole exon is outside of view
                continue

            symbol = ""
            if geneName2!="":
                symbol = geneName2
            ret[(name, symbol)].append((exIdx, exStart, exEnd, exonFrame, nextFrame, exStrand))
            maxIdLen = max(maxIdLen, len(name))

    return ret, maxIdLen

def checkOtherArgs(params):
    # check if minFreq was specified
    minFreq = params.get("minFreq", "0.0")
    try:
        minFreq = float(minFreq)
    except ValueError:
        errAbort("minFreq has to be a floating point number")

    varDb = params.get("varDb", None)

    pam = params.get("pam", None)
    org = params.get("org", None)

    if pamIsXCas9(pam) and org=="noGenome":
        errAbort("You selected no genome, so only efficiency scoring is active. "
           "You also selected the enzyme xCas9. "
           "This does not work, since no efficiency score has been published yet for xCas9. "+
           "Please contact us if you think there is a xCas9 scoring model available somewhere. Thanks!")

    return minFreq, varDb

def crisprSearch(params):
    " do crispr off target search and eff. scoring "
    if "org" in params:
        db = params["org"]
        twoBitFname = getTwoBitFname(db)
        if not isfile(twoBitFname) and db!="noGenome":
            errAbort("Sorry, a genome assembly called %s is not on Crispor "\
                "yet or not anymore. "\
                "Please send us an email if you want us to add it." % db)

    # retrieve sequence if not provided
    if "pos" in params and not "seq" in params:
        params["seq"] = getSeq(params["org"], params["pos"])

    if "batchId" in params:
        # if we're getting only the batchId, extract the parameters from the batch
        # this allows a stable link to a batch that is done
        batchId = params["batchId"]
        seq, org, pamDesc, _, _ = readBatchParams(batchId)
        # pamDesc can include additional options, like guidelen and base editor
        # added after the pam, e.g. "NGG-BE1". setupPamInfo(pam) will set the globals
        # based on it
        seq, warnMsg = cleanSeq(seq, org)
    else:
        # this is a new sequence: create a new batch (TODO: and add it to the queue?)
        seq, org, pamDesc = params["seq"], params["org"], params["pam"]
        newBatchName = params.get("name", "")

        # the "seq" parameter can contain a chrom:start-end position instead of the sequence.
        if re.match(" *[a-zA-Z0-9_-]+: *[0-9, ]+ *- *[0-9,]+(:[+-])? *", seq):
            seq = getSeq(params["org"], seq.replace(" ","").replace(",",""))

        seq, warnMsg = cleanSeq(seq, org)

        if len(seq) > MAXSEQLEN2 and (isSlowPam(pamDesc)):
            errAbort("Sorry, but xCas9, SCanis and enCas12a have so many PAM sites that we are restricting "
                " the input sequence length to %d bp at the moment to keep the "
                "web site fast enough. We will revisit this in a few months. Let us know if "
                "you think this is too short." % MAXSEQLEN2, isWarn=True)

        if len(seq) > MAXSEQLEN3 and (pamDesc in verySlowPams):
            errAbort("Sorry, but SpRY has so many PAM sites that we are restricting "
                " the input sequence length to %d bp at the moment to keep the "
                "web site usable. We will revisit this in a few months. Please let us know if "
                "you think this is too short or have other ideas how to handle the issue." % MAXSEQLEN3, isWarn=True)

        batchId = newBatch(newBatchName, seq, org, pamDesc)
        print ("<script>")
        print(('''history.replaceState('crispor.py', document.title, '?batchId=%s');''' % (batchId)))
        print ("</script>")

    pam = setupPamInfo(pamDesc)
    assert("-" not in pam)

    minFreq, varDb = checkOtherArgs(params)

    if len(warnMsg)!=0:
        print(warnMsg+"<p>")

    batchBase = join(batchDir, batchId)

    # read genome info tab file into memory
    dbInfo = readDbInfo(org)

    # search for PAMs
    uppSeq = seq.upper()
    startDict, endSet = findAllPams(uppSeq, pam)
    otDone = getOfftargets(uppSeq, org, pamDesc, batchId, startDict, None)

    if otDone is None:
        # Job has been added to the queue or is not done yet. 
        printStatus(batchId, warnMsg)
        return

    # if we reach this, the batch has been processed
    batchInfo = readBatchAsDict(batchId)
    position = batchInfo.get("posStr") # if there was no match, the posStr key is "?"

    if dbInfo==None:
        print("<div class='title'>No Genome selected, specificity scoring is deactivated</div>")
        print('<div style="text-align:left;"><strong>Note:</strong> There are no predicted off-targets below and all specificity scores are shown in red as their score is 0. <br></div>')
        chrom = ""

    elif position=='?':
        printQueryNotFoundNote(dbInfo)
        chrom = ""
    else:
        genomePosStr = ":".join(position.split(":")[:2])
        chrom, start, end, strand = parsePos(position)
        start = str(int(start)+1)
        chrom = applyChromAlias(org, chrom)
        oneBasedPosition = "%s:%s-%s" % (chrom, start, end)

        print("<div class='title'><em>")
        if batchName!="":
            print(batchName+":")

        ctUrl = None
        #if org in ["hg19", "mm10"]:
            #ctUrl = ctBaseUrl+"/%s.txt" % batchId

        print("%s (%s)</em>, " % (dbInfo.scientificName, dbInfo.name))
        print('<span style="text-decoration:underline">')
        #mouseOver = "link to UCSC,Ensembl or Gbrowse Genome Browser"
        mouseOver = None
        if dbInfo.server=="manual":
            mouseOver = "no genome browser link available for this organism"
        print(makeBrowserLink(dbInfo, genomePosStr, oneBasedPosition, mouseOver, ["tooltipster"], ctUrl=ctUrl)+"</span>, ")
        if strand=="+":
            print(" forward genomic strand")
        else:
            print(" reverse genomic strand")
        print("</div>")
        #print " (link to Genome Browser)</div>"

    otMatches = parseOfftargets(org, batchId, chrom)
    effScores = readEffScores(batchId)
    sortBy = (params.get("sortBy", None))
    guideData, guideScores, hasNotFound, pamIdToSeq = mergeGuideInfo(uppSeq, startDict, pam, otMatches, \
        position, effScores, sortBy, org=org)


    if len(guideScores)==0:
        print("Found no possible guide sequence. Make sure that your input sequence is long enough and contains at least one match to the PAM motif %s." % pam)
        print('<br><a class="neutral" href="crispor.py">')
        print('<div class="button" style="margin-left:auto;margin-right:auto;width:150px;">New Query</div></a>')
        return

    if hasNotFound and not position=="?":
        printNote("At least one of the possible guide sequences was not found "
        "in the genome. If you pasted a cDNA sequence, note that sequences with "
        "score 0, e.g. splice junctions, are not in the genome, only in the cDNA "
        "and are not usable as CRISPR guides. To find the reference genomic exon sequence "
        "for your cDNA (which contains possibly PCR mutations), please use BLAST "
        "or <a href='https://www.genome.ucsc.edu/cgi-bin/hgBlat'>BLAT</a> to find the best match, "
        "copy the exon sequence from the reference genome and paste it into CRISPOR. <br>"
        "This also applies to any sequence that is different from the reference genome, e.g. mouse "
        "strain sequences: you will have to first map these to a reference genome, then enter the "
        "the reference genome sequence, as otherwise CRISPOR cannot be sure where the target is. "
        "If you have a strain where a reference genome is available, you can contact us, "
        "ideally send us the NCBI genome "
        "identifier (GCA_xxx or GCF_xxx), to crispor@tefor.net")

    chrom, start, end, strand = parsePos(position)

    parNum = isInPar(org, chrom, start, end)
    if parNum!=None:
        print(("<div style='text-align:left; background-color: aliceblue; padding:5px; border: 1px solid black'><strong>Note</strong>: The target sequence is in the PAR%s region. The off-targets on chrY's PAR copy have been removed from the off-target search. We treat the PAR regions as a single region, as all guides are assumed to modify both copies.</div>" % parNum))

    # get list of variant databases
    varLabel = None
    varDbs = readVarDbs(org)

    if len(varDbs)>0 and not position=="?":
        if varDb is None:
            varDb = varDbs[0][1]

        # pull out label of the variant database
        varLabel = None
        for shortLabel, varKey, lab, hasAF in varDbs:
            if varKey==varDb:
                varLabel = lab
                break
        if varLabel is None:
            errAbort("variant DB %s was not found in vcfDescs.txt" % varDb)

        vcfFname = join(genomesDir, org, varDb)
        varDict = findVariantsInRange(vcfFname, chrom, start, end, strand, minFreq)
        varDict["label"] = varLabel
        varShortLabel = shortLabel
    else:
        varDict = None
        varLabel = None
        varShortLabel = None

    varHtmls = varDictToHtml(varDict, seq, varShortLabel)
    showSeqAndPams(org, seq, startDict, pam, guideScores, varHtmls, varDbs, varDb, minFreq, position, pamIdToSeq)

    showSeqDownloadMenu(batchId)

    showGuideTable(guideData, pam, otMatches, dbInfo, batchId, org, chrom, varHtmls)


    print('<br><a class="neutral" href="crispor.py">')
    print('<div class="button" style="margin-left:auto;margin-right:auto;width:150px;">New Query</div></a>')

    #makeCustomTrack(org, chrom, start, end, strand, guideData, batchId, batchName)

def printFile(fname):
    if "/" in fname:
        path = fname
    else:
        myDir = dirname(__file__)
        path = "%s/%s" % (myDir, fname)
#
    if not isfile(path):
        print("install error: %s not found" % path)
        return
    print(open(path).read())

def printCrisporBodyStart():
    #print("""<a href='crispor.py'><img style='width:70px' src='%simage/2021-Logo-Do-3.jpg' alt='UCSC Logo'></a>""" % (HTMLPREFIX))
    print("""<div style='margin-top:10px'><a href='crispor.py'>&nbsp;&larr; Back to CRISPOR homepage</a></div>""")
    print('<div id="bd">')
    print('<div class="centralpanel" style="margin-left:0px">')
    print('<div class="subpanel" style="background:transparent;box-shadow:none;">')
    print('<div class="contentcentral" style="margin-left:0px; width:100%; background:none">')

def printTeforBodyStart():
    print("""<div style="float:left">""")
    #print """<a href='http://genome.ucsc.edu'><img style='vertical-align: top; height: 40px' src='%s/image/ucscBioinf.jpg' alt=''></a>""" % (HTMLPREFIX)
    print("""<a href='crispor.py'><img style='width:150px' src='%simage/2021-Logo-Do-3.jpg' alt='UCSC Logo'></a>""" % (HTMLPREFIX))
    print("""<a href='crispor.py'><img style='width:70px' src='%simage/logo_tefor.png' alt=''></a>""" % (HTMLPREFIX))
    print("</div>")

    print('<div id="bd">')
    print('<div class="centralpanel" style="margin-left:0px">')
    print('<div class="subpanel" style="background:transparent;box-shadow:none;">')
    print('<div class="contentcentral" style="margin-left:0px; width:100%; background:none">')

def printTeforBodyEnd():
    print('<div style="clear:both; text-align:center">Version %s - ' % versionStr)
    print('<a target=_blank href="/manual/">Documentation</a>&nbsp; - ')
    print("""<a href='mailto:%s'>Contact us</a> - <a href="downloads/">Downloads/local installation</a> - <a href="https://academic.oup.com/nar/article/46/W1/W242/4995687">Citation</a> - <a href="https://github.com/maximilianh/crisporWebsite/blob/master/LICENSE.txt">License</a></div>""" % (contactEmail))

    print('</div>')
    print ("""
<script>
$('.tooltipster').tooltipster({
    minWidth: 0,
    contentAsHTML: true,
    maxWidth:400,
    arrow: false,
    interactive: true,
    speed : 0
});

$('.tooltipsterInteract').tooltipster({
    minWidth: 0,
    contentAsHTML: true,
    maxWidth:400,
    interactive: true,
    onlyOne: true,
    arrow: false,
    speed : 0
});

    </script> """)

    print('''
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');
  ga('create', 'UA-38389239-1', 'auto');
  ga('send', 'pageview');
</script>
''')

pamIdRe = re.compile(r's([0-9]+)([+-])g?([0-9]*)')

def intToExtPamId(pamId):
    " convert the internal pam Id like s20+ to the external one, like 21Forw "
    pamPos, strand, rest = pamIdRe.match(pamId).groups()
    if strand=="+":
        strDesc = 'forw'
    else:
        strDesc = 'rev'
    guideDesc = str(int(pamPos)+1)+strDesc
    return guideDesc

def concatGuideAndPam(guideSeq, pamSeq, pamPlusSeq=""):
    " return guide+pam or pam+guide, depending on pamIsFirst "
    if pamIsFirst:
        return pamPlusSeq+pamSeq+guideSeq
    else:
        return guideSeq+pamSeq+pamPlusSeq

def makeGuideHeaders():
    " return list of the headers of the guide output file "
    headers = list(tuple(guideHeaders)) # make a copy of the list

    logging.debug("active scoreNames: %s" % scoreNames)
    tableScoreNames = list(tuple(scoreNames))
    if not pamIsFirst:
        tableScoreNames.extend(mutScoreNames)

    for scoreName in tableScoreNames:
        headers.append(scoreDescs[scoreName][0]+"-Score")

    if not pamIsFirst:
        headers.append("GrafEtAlStatus")

    return headers, tableScoreNames

def effScorePass(effScores, minFusi):
    if minFusi is None:
        return True
    if effScores.get("fusi", 999) < minFusi:
        logging.debug("Fusi score Does not pass min fusi filter")
        return False
    if effScores.get("najm", 999) < minFusi:
        logging.debug("Najm score Does not pass min fusi filter")
        return False
    return True

def iterGuideRows(guideData, addHeaders=False, seqId=None, satMutOpt=None, minSpec=None, minFusi=None):
    "yield rows from guide data. Need to know if for Cpf1 or not "
    headers, tableScoreNames = makeGuideHeaders()

    if satMutOpt:
        headers.append("Oligonucleotide")
        headers.append("AdapterHandle+PrimerFw")
        headers.append("AdapterHandle+PrimerRev")
        oligoPrefix, oligoSuffix, primerFwPrefix, primerRevPrefix, batchId, genome, position, ampLen, tm = satMutOpt

        otMatches = parseOfftargets(genome, batchId)

        guideData.sort(key=operator.itemgetter(3)) # sort by position, makes more sense here

    if seqId != None:
        headers.insert(0, "#seqId")
    else:
        headers[0] = "#"+headers[0]

    headers.append("grafType")

    if addHeaders:
        yield headers

    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, startPos, guideStart, strand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, repCount = guideRow
        if minSpec and guideScore < minSpec:
            continue

        if not effScorePass(effScores, minFusi):
            continue

        otCount = 0
        if otData!=None:
            otCount = len(otData)

        guideDesc = intToExtPamId(pamId)

        fullSeq = concatGuideAndPam(guideSeq, pamSeq)
        row = [guideDesc, fullSeq, guideScore, guideCfdScore, otCount, ontargetDesc]

        for scoreName in tableScoreNames:
            row.append(effScores.get(scoreName, "NotEnoughFlankSeq"))

        grafType = crisporEffScores.getGrafType(guideSeq)
        if grafType is None:
            grafType = "GrafOK"
        row.append(grafType)

        if satMutOpt:
            oligo = oligoPrefix+guideSeq+oligoSuffix
            row.append(oligo)

            chrom, start, end, strand, gene, isUnique = findOntargetPos(otMatches, pamId, position)
            lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq, ampRange, flankSeq, addTags = \
                designPrimer(genome, chrom, start, end, strand, 0, batchId, ampLen, tm)

            fwPrimer = lSeq
            if fwPrimer!=None:
                fwPrimer = primerFwPrefix + fwPrimer

            revPrimer = rSeq
            if revPrimer!=None:
                revPrimer = primerRevPrefix + revPrimer

            row.append(fwPrimer)
            row.append(revPrimer)


        row = [str(x) for x in row]
        if seqId != None:
            row.insert(0, seqId)
        yield row

def iterOfftargetRows(guideData, addHeaders=False, skipRepetitive=True, seqId=None):
    " yield bulk offtarget rows for the tab-sep download file "
    otRows = []

    headers = list(offtargetHeaders) # clone list
    if seqId:
        headers.insert(0, "seqId")

    if addHeaders:
        otRows.append(headers)

    skipCount = 0

    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, startPos, guideStart, strand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, \
            ontargetDesc, repCount = guideRow

        if otData!=None:
            otCount = len(otData)
            if otCount > 5000 and skipRepetitive:
                skipCount += otCount
                continue

            for otSeq, mitScore, cfdScore, editDist, pos, gene, alnHtml, inLinkage in otData:
                gene = gene.replace(",", "_").replace(";","-")
                chrom, start, end, strand = parsePos(pos)
                guideDesc = intToExtPamId(pamId)
                mismStr = highlightMismatches(guideSeq, otSeq, len(pamSeq))
                fullSeq = concatGuideAndPam(guideSeq, pamSeq)
                row = [guideDesc, fullSeq, otSeq, mismStr, editDist, mitScore, cfdScore, chrom, start, end, strand, gene]
                if seqId:
                    row.insert(0, seqId)
                row = [str(x) for x in row]
                otRows.append(row)

    if skipCount != 0:
        newRow = [""]*len(headers)
        newRow[0] = "# %d off-targets are not shown: guides with more than 5000 off-targets were considered too repetitive" % skipCount
        otRows.insert(0, newRow)

    return otRows

def writeHtmlTable(rows, outFile):
    " write list of rows to outFile in html format "
    outFile.write('<link rel="stylesheet" href="https://unpkg.com/purecss@1.0.0/build/pure-min.css" integrity="sha384-nn4HPE8lTHyVtfCBi5yW9d20FjT8BJwUXyWZT9InLYax14RDjBj46LmSztkmNP9w" crossorigin="anonymous">\n')
    outFile.write("<table class='pure-table'>\n")
    headDone = False
    for row in rows:
        if headDone:
            tag = "td"
        else:
            tag = "th"

        outFile.write("<tr>\n")
        for field in row:
            outFile.write("<%s>%s</%s>" % (tag, field, tag))
        outFile.write("</tr>\n")
        headDone = True
    outFile.write("</table>\n")

def writeTable(fileFormat, rows, ofh):
    " write table to ofh, currently writes xls files as tsv "
    if fileFormat=="tsv" or fileFormat=="xls" or fileFormat=="csv":
        sep = "\t"
        if fileFormat=="csv":
            sep =","
        for row in rows:
            ofh.write(sep.join(row))
            ofh.write("\n")
        ofh.close()
    else:
        writeHtmlTable(rows, ofh)

def xlsWrite(rows, title, outFile, colWidths, fileFormat, seq, org, pam, position, batchId, optFields=None):
    """ given rows, writes a XLS binary stream to outFile, if xlwt is available
    Otherwise writes a tab-sep file.
    colWidths is a list of widths of columns, in Arial characters.
    """
    if xlwtLoaded and fileFormat in ["xls"]:
        seqStyle = xlwt.easyxf('font: name Courier New')
        charSize = 269 # see http://reliablybroken.com/b/2011/10/widths-heights-with-xlwt-python/
        wb = xlwt.Workbook()
        ws = wb.add_sheet(title)

        ws.write(0, 0, "# Name")
        ws.write(0, 1, batchName)
        ws.write(1, 0, "# Sequence")
        ws.write(1, 1, seq)
        ws.write(3, 0, "# PAM")
        ws.write(3, 1, pam)
        ws.write(2, 0, "# Genome")
        ws.write(2, 1, org)
        ws.write(4, 0, "# Position")
        ws.write(4, 1, position)

        ws.write(5, 0, "# Version")
        #http://stackoverflow.com/questions/4530069/python-how-to-get-a-value-of-datetime-today-that-is-timezone-aware
        FORMAT='%Y-%m-%dT%H:%M:%S%Z'
        dateStr=time.strftime(FORMAT, time.localtime())
        ws.write(5, 1, "CRISPOR %s, %s" % (versionStr, dateStr))

        ws.write(6, 0, "# Results")
        url = "http://crispor.gi.ucsc.edu/crispor.py?batchId=%s" % batchId
        #ws.write(6, 1, xlwt.Formula('HYPERLINK("%s";"Link")' % (url)))
        ws.write(6, 1, url)

        startRow = 7
        curRow = startRow
        if optFields is not None:
            for key, val in optFields.items():
                ws.write(curRow, 0, "# %s" % key)
                ws.write(curRow, 1, val)
                curRow+=1

        skipRows = curRow + 1

        seqCols = [1, 7, 8, 9] # columns with sequences -> fixed width font

        for rowCount, row in enumerate(rows):
            if rowCount==65534-startRow:
                ws.write(rowCount+skipRows, 0, "WARNING: cannot write more than 65535 rows to an Excel file. Switch to .tsv format to get all off-targets.")
                break

            isHeader = False
            if "Id" in row[0]:
                isHeader = True

            for colCount, col in enumerate(row):
                if col.isdigit():
                    col = int(col)
                else:
                    # -0.1 is not a digit, so try to convert to float
                    try:
                        col = float(col)
                    except ValueError:
                        pass
                if colCount in seqCols and not isHeader:
                    ws.write(rowCount+skipRows, colCount, col, seqStyle)
                else:
                    ws.write(rowCount+skipRows, colCount, col)

        # set sizes in characters per column
        for colId, colWidth in enumerate(colWidths):
            ws.col(colId).width = charSize*colWidth

        try:
            outFile.flush() # flush out the Content-type header
            wb.save(outFile.buffer)
        except:
            print("error")

    elif fileFormat=="html":
        writeHtmlTable(rows, outFile)
    else:
        # raw ASCII tsv output mode
        sep = "\t"
        if fileFormat=="csv":
            sep = ","
        for row in rows:
            outFile.write(sep.join(row))
            outFile.write("\n")

    outFile.flush()

def seqToGenbankLines(seq):
    """ chunk sequence string into lines each with six parts of 10bp, return as a list
    >>> seqToGenbankLines("aacacacatggtacacactgactagctagctacgatccagtacgatcgacgtagctatcgatcgatcgatcgactagcta")
    ['aacacacatg gtacacactg actagctagc tacgatccag tacgatcgac gtagctatcg', 'atcgatcgat cgactagcta']
    """
    # first chunk into 10bp parts
    parts = [seq[i:i+10] for i in range(0, len(seq), 10)]

    # put into lines of 6*10 bp
    lines = []
    for i in range(0, len(parts), 6):
        lines.append(" ".join(parts[i:i+6]))
    return lines

def writeLn(fh, line, indent=None, doWrap=True):
    " write line to file, using \r\n "
    if indent==None:
        fh.write(line)
        fh.write("\r\n")
    else:
        if doWrap:
            lineSize = 80-indent
        else:
            lineSize = 10000
        parts = [line[i:i+lineSize] for i in range(0, len(line), lineSize)]
        spacer = "".join(([" "]*indent))
        for p in parts:
            fh.write(spacer+p)
            fh.write("\r\n")

def genbankWrite(batchId, fileFormat, desc, seq, org, position, pam, guideData, ofh):
    " write a description of the current job in genbank format to ofh "
    if fileFormat=="serialcloner":
        # a bug in serial cloner means that we cannot use the linear format
        writeLn(ofh, """LOCUS       %s    %d bp      DNA     circular   1/1/17""" % (desc, len(seq)))
    elif fileFormat=="snapgene":
        writeLn(ofh, "LOCUS       Exported                 239 bp ds-DNA     linear   SYN 22-MAR-2017")
    else:
        writeLn(ofh, """LOCUS       %s    %d bp      DNA     linear   1/1/17""" % (desc, len(seq)))

    batchUrl = "http://crispor.gi.ucsc.edu/crispor.py?batchId=%s" % batchId
    seqDesc1 = """Sequence exported from CRISPOR.org V%s""" % versionStr
    seqDesc2 = "Genome %s, position %s. View full CRISPOR results at %s""" % (org, position, batchUrl)

    if fileFormat in ["ape"]:
        seqDesc2 += " Features indicate PAMs. Click the little triangles next to the features to show scores and guide sequences."

    if fileFormat=="genomecompiler":
        # genomecompiler plasmid viewer can't show more than ~30 characters as the seq definition line
        writeLn(ofh, """DEFINITION  %s""" % desc)
    else:
        writeLn(ofh, """DEFINITION  %s""" % seqDesc1)
        writeLn(ofh, seqDesc2, indent=12)
        writeLn(ofh, """             Export for: %s""" % (fileFormat))

    #writeLn(ofh, """ACCESSION""")
    #writeLn(ofh, """VERSION""")
    writeLn(ofh, """SOURCE      %s""" % org)
    writeLn(ofh, """  ORGANISM  %s""" % org)

    if fileFormat in ["serialcloner"]:
        writeLn(ofh, """COMMENT     Serial Cloner Genbank Format""")
        writeLn(ofh, """COMMENT     SerialCloner_Type=DNA""")
        writeLn(ofh, """COMMENT     SerialCloner_Comments=%s""" % seqDesc1+" "+seqDesc2)
        writeLn(ofh, """COMMENT     SerialCloner_Ends=0,0,,0,""")
    elif fileFormat in ["ape", "genomecompiler", "vnti"]:
        # ape does not show the definition line, only the comment block
        writeLn(ofh, """COMMENT     %s""" % seqDesc1)
        writeLn(ofh, seqDesc2, indent=12)

    if fileFormat=="snapgene":
        # this tells snapgene that the file is really in snapgene format
        writeLn(ofh, """KEYWORDS    snapgene3""")
        writeLn(ofh, """REFERENCE   1  (bases 1 to 239)""")
        writeLn(ofh, """  AUTHORS   .""")
        writeLn(ofh, """  TITLE     Direct Submission""")
        writeLn(ofh, """  JOURNAL   Exported Wednesday, Mar 22, 2017 from SnapGene Viewer 3.3.3""")
        writeLn(ofh, """          http://www.snapgene.com""")

    writeLn(ofh, """FEATURES             Location/Qualifiers""")

    i = 1
    guideData.sort(key=operator.itemgetter(3)) # sort by position

    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, startPos, guideStart, strand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, repCount = guideRow

        guideUrl = batchUrl+"&pamId="+urllib.parse.quote(pamId)+"&pam="+pam

        otCount = 0
        if otData!=None:
            otCount = len(otData)

        fullSeq = concatGuideAndPam(guideSeq, pamSeq)

        if fileFormat in ["geneious", "snapgene"]:
            # this code annotates only the guide sequence, suggested by Alyce Chen
            start = guideStart + 1
            end = start + len(guideSeq) - 1
        else:
            # most viewers don't handle overlaps well. We highlight only the PAM in these cases
            if strand=="+":
                start = startPos + 1
                end   = startPos + len(pamSeq)
            else:
                start = startPos + 1
                end = startPos + len(pamSeq)

        colorHex, colorName = scoreToColor(guideScore)
        guideName = intToExtPamId(pamId)
        urlLink = "full details/primers at %s"  % guideUrl

        if pamIsSpCas9(pam):
            fusiScore = str(effScores.get("fusi", -1))
            crisprScanScore = str(effScores.get("crisprScan", -1))
            descStr = "%s: Spec %s, Eff %s/%s" % (guideName, guideScore, fusiScore, crisprScanScore)
            guideSeqDescSeq = "Guide %s MIT-Spec %s, Eff Doench2016 %s, Eff Mor.-Mat. %s" % (guideSeq, guideScore, fusiScore, crisprScanScore)
            longDesc = "MIT-Specificity score: %s, Efficiency Doench2016 = %s, Efficiency Moreno-Mateos = %s, guide sequence: %s, " \
                "full details/primers at %s" % (guideScore, fusiScore, crisprScanScore, guideSeq, guideUrl)
        else:
            mainEffName = scoreNames[0]
            descStr = "%s: Spec %s, Eff %s" % (guideName, guideScore, effScores.get(mainEffName, -1))

            shortStrList = []
            longStrList = []
            for scoreName in scoreNames:
                scoreVal = str(effScores.get(scoreName, -1))
                shortStrList.append("Eff %s %s" % (scoreName, scoreVal))
                longStrList.append("Efficiency %s = %s" % (scoreDescs[scoreName][0], scoreVal))

            guideSeqDescSeq = "Guide %s MIT-Spec %s, %s" % (guideSeq, guideScore, ", ".join(shortStrList))
            longDesc = "MIT-Specificity score: %s, %s, guide sequence: %s, %s" \
                % (guideScore, ", ".join(longStrList), guideSeq, urlLink)


        featType = "misc_feature"

        if fileFormat in ["genomecompiler"]:
            # genomecompiler does not store colors in the genbank file.
            # we use the feature type to simulate colors
            if colorName=="green":
                featType = "promoter"
            elif colorName=="red":
                featType = "terminator"
            elif colorName=="yellow":
                featType = "misc_feature"

        if strand=="+":
            writeLn(ofh, """     %s    %d..%d""" % (featType, start, end))
        else:
            writeLn(ofh, """     %s    complement(%d..%d)""" % (featType, start, end))

        if fileFormat in ["serialcloner"]:
            # serialcloner has no label
            writeLn(ofh, '''                     /note="%s"''' % descStr)
            writeLn(ofh, """                     /SerialCloner_Color=&h%s""" % colorHex.replace("#", ""))
            writeLn(ofh, """                     /SerialCloner_Show=True""")
            writeLn(ofh, """                     /SerialCloner_Protect=True""")
            writeLn(ofh, """                     /SerialCloner_Arrow=True""")
            writeLn(ofh, """                     /SerialCloner_Desc=%s""" % longDesc)
        elif fileFormat in ["ape"]:
            # ape can use multiple note lines
            writeLn(ofh, '''                     /locus_tag="%s"''' % descStr)
            writeLn(ofh, '''                     /note=MIT Specificity: %s''' % guideScore)
            writeLn(ofh, '''                     /note=Efficiency: Doench2016 %s Mor-Mat. %s''' % (str(effScores["fusi"]), str(effScores["crisprScan"])))
            writeLn(ofh, """                     /ApEinfo_fwdcolor=%s""" % colorHex)
            writeLn(ofh, """                     /ApEinfo_revcolor=%s""" % colorHex)
            writeLn(ofh, """                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0""")
            writeLn(ofh, '''                     /note=Guide %s''' % guideSeq)
        elif fileFormat=='benchling':
            writeLn(ofh, '''                     /note="%s"''' % guideSeqDescSeq)
            writeLn(ofh, """                     /ApEinfo_fwdcolor=%s""" % colorHex)
            writeLn(ofh, """                     /ApEinfo_revcolor=%s""" % colorHex)
            writeLn(ofh, """                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0""")
        elif fileFormat in ['genomecompiler']:
            writeLn(ofh, '''                     /label="%s"''' % guideName)
            writeLn(ofh, '''                     /note="%s"''' % guideSeqDescSeq)
        elif fileFormat in ['snapgene']:
            writeLn(ofh, '''                     /label="%s"''' % guideName)
            writeLn(ofh, '''                     /note="%s"''' % longDesc)
            if strand=="+":
                writeLn(ofh, '''                     /note="color: %s; direction: RIGHT"''' % colorHex)
            else:
                writeLn(ofh, '''                     /note="color: %s; direction: LEFT"''' % colorHex)
        # vector NTI treats attributes as a key-val list
        # lasergene shows all attributes, vector NTI is the most complete way to show all data
        elif fileFormat in ["vnti", "lasergene"]:
            writeLn(ofh, '''/label="%s"''' % guideName, indent=21)
            writeLn(ofh, '''/note="%s"''' % descStr, indent=21)
            writeLn(ofh, '''/MITSpecScore="%s"''' % str(guideScore), indent=21)
            # longDesc = "MIT-Specificity score: %s, Efficiency Doench2016 = %s, Efficiency Moreno-Mateos = %s, guide sequence: %s, full details/primers at %s" % (guideScore, str(effScores["fusi"]), str(effScores["crisprScan"]), guideSeq, guideUrl)
            writeLn(ofh, '''/Doench2016Eff="%s"''' % str(effScores["fusi"]), indent=21)
            writeLn(ofh, '''/Mor-MateosEff="%s"''' % str(effScores["crisprScan"]), indent=21)
            writeLn(ofh, '''/guide_sequence="%s"''' % guideSeq, indent=21)
            writeLn(ofh, '''/url="%s"''' % guideUrl, indent=21)
        elif fileFormat in ["geneious"]:
            writeLn(ofh, '''/label="%s"''' % guideName, indent=21, doWrap=False) # geneious translates \n to spaces, breaks link
            writeLn(ofh, '''/MIT-spec_score="%s"''' % guideScore, indent=21, doWrap=False)
            writeLn(ofh, '''/guide_sequence="%s"''' % guideSeq, indent=21, doWrap=False)
            writeLn(ofh, '''/note="%s"''' % urlLink, indent=21, doWrap=False)
            if "crisprScan" in effScores:
                writeLn(ofh, '''/Efficiency="%s"''' % effScores["crisprScan"], indent=21, doWrap=False)
        else:
            writeLn(ofh, '''/label="%s"''' % guideName, indent=21)
            writeLn(ofh, '''/note="%s"''' % longDesc, indent=21)

        i += 1

    writeLn(ofh, """ORIGIN""")
    i = 1
    for line in seqToGenbankLines(seq):
        writeLn(ofh, """%9d %s""" % (i, line))
        i += 60

    writeLn(ofh, "//")
    ofh.close()

def writeHttpAttachmentHeader(fname, doDownload=True):
    " write the http header for attachments "
    if doDownload:
        print('Content-Disposition: attachment; filename="%s"' % (fname))
    else:
        print('Content-type: text/html')
    print("") # = end of http headers
    sys.stdout.flush()

def buildPoolOptions(barcodeId, custPrefix="", custSuffix=""):
    " return a list of pool settings and a dictionary with pool options "
    barcodeDict = dict(satMutBarcodes)
    barcodeLabel = barcodeDict[barcodeId]

    if int(barcodeId)==0:
        barcodePre, barcodePost = "", ""
    else:
        barcodePre, barcodePost = barcodeLabel.split()[-1].split("/")

    optFields = OrderedDict()
    optFields["Subpool Barcode"] = barcodeLabel

    optFields["Subpool Prefix"] = barcodePre
    optFields["Subpool Suffix"] = barcodePost

    if custPrefix!="" and custSuffix!="":
        optFields["Custom Prefix"] = custPrefix
        optFields["Custom Suffix"] = custSuffix
        optFields["Oligonucl. structure"] = "Subpool Prefix + Custom Prefix + sgRNA + Custom Suffix + Subpool Suffix"

        fullPrefix = barcodePre+custPrefix
        fullSuffix = barcodePost+custSuffix
    else:

        prePrimer = "GGAAAGG"
        pre = "ACGAAACACCG"
        post = "GTTTTAGAGCTAGAAATA"
        postPrimer = "GCAAGTTAAAATAAGGC"
        fullPrefix = barcodePre+prePrimer+pre
        fullSuffix = post+postPrimer+barcodePost

        optFields["pLentiGuidePre primer"] = prePrimer
        optFields["pLentiGuidePre"] = pre

        optFields["pLentiGuidePost"] = post
        optFields["pLentiGuidePost primer"] = postPrimer

        optFields["Oligonucl. structure"] = "Subpool Prefix + pLentiGuidePre Primer + pLentiGuidePre + sgRNA + pLentiGuide Post + pLentiGuidePost Primer + Subpool Suffix"

    primerFwPrefix = 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
    primerRevPrefix = 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'

    satMutOpt = [fullPrefix, fullSuffix, primerFwPrefix, primerRevPrefix]

    return satMutOpt, optFields

def writeSatMutFile(barcodeId, ampLen, tm, batchId, minSpec, minFusi, fileFormat, outFh):
    " write saturating mutagenesis table of all guides, scores, primers and amplicons around them to outFh "
    seq, org, pam, position, guideData = readBatchAndGuides(batchId)

    satMutOpt, optFields = buildPoolOptions(barcodeId)
    satMutOpt.extend( [batchId, org, position, ampLen, tm] )

    guideRows = iterGuideRows(guideData, addHeaders=True, satMutOpt=satMutOpt, minSpec=minSpec, minFusi=minFusi)
    xlsWrite(guideRows, "guides", outFh, [20,28,10,10,10,10,10,60,21,21], fileFormat, seq, org, pam, position, batchId, optFields=optFields)

def readBatchAndGuides(batchId):
    " parse the input file, the batchId-json file and the offtargets and link everything together "
    seq, org, pam, position, extSeq = readBatchParams(batchId)
    chrom, _, _, _ = parsePos(position)
    pam = setupPamInfo(pam)
    uppSeq = seq.upper()

    startDict, endSet = findAllPams(uppSeq, pam)

    effScoreFname = join(batchDir, batchId+".effScores.tab")

    otMatches = parseOfftargets(org, batchId, chrom)
    effScores = readEffScores(batchId)
    guideData, guideScores, hasNotFound, pamIdToSeq = mergeGuideInfo(uppSeq, startDict, pam, otMatches, position, effScores, org=org)
    return seq, org, pam, position, guideData

def writeOntargetAmpliconFile(outType, batchId, ampLen, tm, ofh, fileFormat="tsv", minSpec=0, minFusi=0):
    """ design primers with approx ampLen and tm around each guide's target.
    outType can be "primers" or "amplicons"
    """
    inSeq, db, pamPat, position, extSeq = readBatchParams(batchId)
    chrom, _, _, _ = parsePos(position)
    otMatches = parseOfftargets(db, batchId, chrom)

    startDict, endSet = findAllPams(inSeq, pamPat)
    pamSeqs = list(flankSeqIter(inSeq, startDict, len(pamPat), True))

    allEffScores = readEffScores(batchId)
    guideData, guideScores, hasNotFound, pamIdToSeq = mergeGuideInfo(inSeq, startDict, pamPat, otMatches, position, allEffScores, sortBy="pos", org=db)

    if outType=="primers":
        headers = ["#guideId", "forwardPrimer", "leftPrimerTm", "revPrimer", "revPrimerTm", "ampliconSequence", "guideSequence"]
    else:
        headers = ["#guideId", "ampliconSequence", "guideSequence"]

    rows = []
    rows.append(headers)
    #ofh.write("\t".join(headers))
    #ofh.write("\n")

    #for pamId, pamStart, guideStart, strand, guideSeq, pamSeq, pamPlusSeq in pamSeqs:
    for guideScore, guideCfdScore, effScores, startPos, guideStart, strand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, \
            ontargetDesc, repCount in guideData:

        if guideScore < minSpec:
            continue
        if not effScorePass(effScores, minFusi):
            continue

        chrom, start, end, strand, gene, isUnique = findOntargetPos(otMatches, pamId, position)
        #effScores = allEffScores.get(pamId, None)

        note = ""
        if not isUnique:
            note = "warning: guide has no unique match in genome"

        lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq, ampRange, flankSeq, addTags = \
            designPrimer(db, chrom, start, end, strand, 0, batchId, ampLen, tm)

        pamName = intToExtPamId(pamId)
        if outType=="primers":
            row = [pamName, lSeq, lTm, rSeq, rTm, targetSeq, guideSeq]
        else:
            row = [pamName, targetSeq, guideSeq]

        row = [str(x) for x in row]
        rows.append(row)

    writeTable(fileFormat, rows, ofh)

def writeTargetSeqs(guideData, ofh, minSpec=None, minFusi=None):
    " write the guide sequences and their pam to ofh "
    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, startPos, guideStart, strand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, repCount = guideRow
        if minSpec and guideScore < minSpec:
            continue
        if minFusi and effScores["fusi"] < minFusi:
            continue

        fullSeq = concatGuideAndPam(guideSeq, pamSeq)
        row = [fullSeq]
        ofh.write("\t".join(row))
        ofh.write("\n")
    ofh.close()

def writeTargetLocs(position, guideData, ofh, fileFormat, minSpec=None, minFusi=None):
    " write the guide locations and their sequences to ofh, in the format for crisprSurf"
    seqChrom, seqStart, seqEnd, seqStrand = parsePos(position)

    rows = []
    header = ["Chr","Start","Stop","sgRNA_Sequence","Strand","sgRNA_Type"]
    rows.append(header)

    for guideRow in guideData:
        guideScore, guideCfdScore, effScores, startPos, guideStart, guideStrand, pamId, \
            guideSeq, pamSeq, otData, otDesc, last12Desc, mutEnzymes, ontargetDesc, totalAlnCount = guideRow
        if minSpec and guideScore < minSpec:
            continue
        if not effScorePass(effScores, minFusi):
            continue

        chromStart, chromEnd, _, _, chromStrand = mapToGenome(seqStart, seqStrand, startPos, guideStart, guideStrand)

        row = [seqChrom, chromStart, chromEnd, chromStrand, guideSeq, 'observation']
        row = [str(x) for x in row]
        rows.append(row)

    # crisprsurf takes only csv files
    if fileFormat=="tsv":
        fileFormat="csv"
    writeTable(fileFormat, rows, ofh)

def fastaWrite(seqId, seq, fh, width=80):
    """ output fasta seq to file object, break to 80 char width """
    fh.write(">"+seqId+"\n")
    if len(seq)>width:
        last = 0
        for l in range(width,len(seq),width):
            fh.write(seq[last:l])
            fh.write("\n")
            last = l
        fh.write(seq[last:len(seq)])
    else:
        fh.write(seq)
    fh.write("\n")

def downloadFile(params):
    " "
    global scoreNames

    batchId = params["batchId"]
    seq, org, pam, position, guideData = readBatchAndGuides(batchId)

    if batchName!="":
        queryDesc = batchName+"_"
    else:
        queryDesc = ""

    if position=="?":
        queryDesc += org+"-unknownLoc"
    else:
        queryDesc += org+"-"+position.strip(":+-").replace(":","-")
        #print org, position, queryDesc

    fileType = params["download"]
    # the assistant cannot set the argument as it uses multi-submit buttons
    if fileType=="useGet":
        for key in params:
            if key.startswith("get-"):
                fileType = key.split("-")[1]
                break

    fileFormat = params.get("format", "tsv")
    if not fileFormat in ["tsv", "xls", "html"]:
        errAbort("Invalid fileFormat argument")

    doDownload = True
    if fileFormat=="html":
        doDownload = False

    if fileType=="guides":
        if cgiParams.get("showAllScores", "0")=="1":
            if not pamIsCpf1(pam) and not pamIsSaCas9(pam):
                scoreNames = allScoreNames
        writeHttpAttachmentHeader("guides_%s.%s" % (queryDesc, fileFormat), doDownload)
        xlsWrite(iterGuideRows(guideData, addHeaders=True), "guides", sys.stdout, [9,28,10,10], fileFormat, seq, org, pam, position, batchId)

    elif fileType=="offtargets":
        writeHttpAttachmentHeader("offtargets_%s.%s" % (queryDesc, fileFormat), doDownload)
        skipRepetitive = (fileFormat=="xls")
        otRows = list(iterOfftargetRows(guideData, addHeaders=True, skipRepetitive=skipRepetitive))
        doReverse = (not pamIsFirst)
        otRows.sort(key=operator.itemgetter(4), reverse=doReverse)
        xlsWrite(otRows, "offtargets", sys.stdout, [9,28,28,5], fileFormat, seq, org, pam, position, batchId)

    elif fileType=="targetSeqs":
        writeHttpAttachmentHeader("targetSeqs_%s.txt" % (queryDesc), doDownload)
        minSpec = cgiGetNum(params, "minSpec", 0)
        minFusi = cgiGetNum(params, "minFusi", 0)
        writeTargetSeqs(guideData, sys.stdout, minSpec=minSpec, minFusi=minFusi)

    elif fileType=="targetLocs":
        writeHttpAttachmentHeader("targetLocs_%s.csv" % (queryDesc), doDownload)
        minSpec = cgiGetNum(params, "minSpec", 0)
        minFusi = cgiGetNum(params, "minFusi", 0)
        writeTargetLocs(position, guideData, sys.stdout, fileFormat, minSpec=minSpec, minFusi=minFusi)

    elif fileType=="amplicons":
        # write amplicons of all off-targets for a single guide
        fname = makeCrispressoFname(batchName, batchId)
        writeHttpAttachmentHeader(fname, doDownload)
        pamId = cgiGetStr(params, "pamId")
        writeAmpliconFile(params, batchId, pamId, sys.stdout)

    elif fileType=="ontargetAmplicons":
        # design primers around all targets in input sequence
        writeHttpAttachmentHeader("ontargetAmplicons_%s.tsv" % (queryDesc), doDownload)
        ampLen = cgiGetNum(params, "ampLen", 140)
        tm = cgiGetNum(params, "tm", 60)
        minSpec = cgiGetNum(params, "minSpec", 0)
        minFusi = cgiGetNum(params, "minFusi", 0)
        writeOntargetAmpliconFile("amplicons", batchId, ampLen, tm, sys.stdout, fileFormat, minSpec, minFusi)

    elif fileType=="ontargetPrimers":
        writeHttpAttachmentHeader("ontargetPrimers_%s.tsv" % (queryDesc), doDownload)
        ampLen = cgiGetNum(params, "ampLen", 140)
        tm = cgiGetNum(params, "tm", 60)
        minSpec = cgiGetNum(params, "minSpec", 0)
        minFusi = cgiGetNum(params, "minFusi", 0)
        writeOntargetAmpliconFile("primers", batchId, ampLen, tm, sys.stdout, fileFormat, minSpec, minFusi)

    elif fileType=="satMut":
        fileName = "satMutOligos-%s.%s" % (queryDesc, fileFormat)

        barcodeId = cgiGetNum(params, "barcode", 0)
        barcodeDict = dict(satMutBarcodes)
        if not barcodeId in barcodeDict:
            errAbort("'barcodeId' parameter is not a valid index in our barcode table.")

        ampLen = cgiGetNum(params, "ampLen", 140)
        tm = cgiGetNum(params, "tm", 60)

        writeHttpAttachmentHeader(fileName, doDownload)
        minSpec = cgiGetNum(params, "minSpec", 0)
        minFusi = cgiGetNum(params, "minFusi", 0)
        writeSatMutFile(barcodeId, ampLen, tm, batchId, minSpec, minFusi, fileFormat, sys.stdout)

    elif fileType in ["serialcloner", "ape", "genomecompiler", "fasta", "benchling", "snapgene", "genbank", "vnti", "lasergene", "geneious"]:
        fileFormat = params['download']
        ext = "gb"
        if fileFormat=="serialcloner":
            ext = "xdna"
        elif fileFormat=="ape":
            ext = "ape"
        elif fileFormat=="snapgene":
            ext = "dna"
        elif fileFormat=="fasta":
            ext = "fa"
        elif fileFormat=="vnti":
            ext = "gb"
        fileName = "crispor_%s-%s.%s" % (queryDesc, fileFormat, ext)

        if fileFormat!="genomecompiler":
            writeHttpAttachmentHeader(fileName)
        else:
            print("Content-type: text/plain\n")

        if fileFormat!="fasta":
            genbankWrite(batchId, fileFormat, queryDesc, seq, org, position, pam, guideData, sys.stdout)
        else:
            fastaWrite("crispor-"+queryDesc, seq, sys.stdout)

    else:
        errAbort("invalid value for download parameter, fileType=%s" % fileType)

def makeCrispressoFname(batchName, batchId):
    fnameDesc = ["crisporAmplicons"]
    if batchName!="":
        fnameDesc.append(batchName)
    fnameDesc.append(batchId)
    fname = "_".join(fnameDesc)+".txt"
    return fname

def designOfftargetPrimers(inSeq, db, pam, position, extSeq, pamId, ampLen, tm, otMatches):
    " return a list of off-target primers sorted by CFD score "
    targetChrom, targetStart, targetEnd, strand = parsePos(position)
    chromSizes = parseChromSizes(db)
    pam = setupPamInfo(pam)

    guideSeq, pamSeq, pamPlusSeq, guideSeqWPam, guideStrand, guideSeqHtml, guideStart, guideEnd \
        = findGuideSeq(inSeq, pam, pamId)

    # get the coords
    coords = []
    names = []
    nameToOtScoreSeq = {}
    flankSize = 1000
    for mismCount, otMatchRows in otMatches.items():
        for otMatch in otMatchRows:
            chrom, start, end, otSeq, strand, segType, segName, totalAlnCount, fromXaTag = otMatch
            if chrom==targetChrom and start>=targetStart and end<=targetEnd:
                prefix = "ontarget_"
            else:
                prefix = ""
            segDesc = segTypeConv.get(segType, "") # some genomes do not have descriptions
            name = "%(prefix)smm%(mismCount)d_%(segDesc)s_%(segName)s_%(chrom)s_%(start)d" % locals()
            if start-flankSize < 0 or end+flankSize > chromSizes[chrom]:
                print(("Cannot design primer for %s, too close to chromosome boundary" % name))
            else:
                coords.append( (chrom, start-flankSize, end+flankSize, name) )
                otScore = calcCfdScore(guideSeq, otSeq)
                nameToOtScoreSeq[name] = (otScore, otSeq)

    # coords -> sequences
    flankSeqs = getGenomeSeqsBin(db, coords, doRepeatMask=True)
    targetSeqs = [(x[3], x[6]) for x in flankSeqs] # strip coords, keep name+seq
    nameToSeq = dict(targetSeqs)

    # sequences -> primers
    # check the input parameters: ampLen, tm
    ampMin = ampLen-110
    ampRange = "%d-%d" % (ampMin, ampLen)

    primers = runPrimer3(targetSeqs, flankSize, GUIDELEN+len(pamSeq), ampRange, tm, {})

    # sort primers by CFD off-target score
    scoredPrimers = []
    for name, primerInfo in primers.items():
        score, otSeq = nameToOtScoreSeq[name]
        scoredPrimers.append( (score, name, primerInfo) )
    scoredPrimers.sort(reverse=True)

    return scoredPrimers, nameToSeq, nameToOtScoreSeq, guideSeqHtml

def makeCrispressoOfftargetRows(scoredPrimers, nameToSeq, nameToOtScoreSeq):
    " yield the crispresso off-target rows "
    for score, seqName, primerInfo in scoredPrimers:
        flankSeq = nameToSeq[seqName]
        score, otSeq = nameToOtScoreSeq[seqName]
        lSeq, lTm, lPos, rSeq, rTm, rPos = primerInfo
        if lSeq is None:
            continue
        ampSeq = flankSeq[lPos:rPos+1]
        row = [seqName, ampSeq, otSeq, "NA", "NA"]
        yield row

def writeAmpliconFile(params, batchId, pamId, outFh):
    " create the table of off-target amplicons for crispressoPooled and write to outFh "
    ampLen = cgiGetNum(params, "ampLen", 140)
    tm = cgiGetNum(params, "tm", 60)

    inSeq, db, pam, position, extSeq = readBatchParams(batchId)

    pamOtMatches = parseOfftargets(db, batchId)
    otMatches = pamOtMatches[pamId]

    scoredPrimers, nameToSeq, nameToOtScoreSeq, guideSeqHtml = \
        designOfftargetPrimers(inSeq, db, pam, position, extSeq, pamId, ampLen, tm, otMatches)

    for row in makeCrispressoOfftargetRows(scoredPrimers, nameToSeq, nameToOtScoreSeq):
        outFh.write("\t".join(row))
        outFh.write("\n")

def otPrimerPage(params):
    " print a page with all primers for all off-targets "
    # initialize everything
    batchId, pamId = params["batchId"], params["pamId"]

    ampLen = cgiGetNum(params, "ampLen", 140)
    tm = cgiGetNum(params, "tm", 60)

    inSeq, db, pam, position, extSeq = readBatchParams(batchId)
    pamOtMatches = parseOfftargets(db, batchId)
    otMatches = pamOtMatches[pamId]

    scoredPrimers, nameToSeq, nameToOtScoreSeq, guideSeqHtml = \
        designOfftargetPrimers(inSeq, db, pam, position, extSeq, pamId, ampLen, tm, otMatches)

    # primers -> table rows
    primerTable = []
    for otScore, seqName, primerInfo in sorted(scoredPrimers, reverse=True):
        if otScore is None:
            otScore = "noOfftargetScore"
        else:
            otScore = "%.3f" % otScore
            otScore = otScore[:4]
        lSeq, lTm, lPos, rSeq, rTm, rPos = primerInfo
        if batchName!="":
            baseName = batchName+"_"+seqName
        else:
            baseName = seqName

        # Nextera handle sequences, from Matt Canver
        prefixForw="<b>TCGTCGGCAGCGTC</b>"
        prefixRev="<b>GTCTCGTGGGCTCGG</b>"

        if lSeq is None:
            primerTable.append( (baseName+"_F", "Primer3: not found at this Tm", "N.d.", otScore) )
        else:
            primerTable.append( (baseName+"_F", prefixForw+lSeq, lTm[:4], otScore) )

        if rSeq is None:
            primerTable.append( (baseName+"_R", "Primer3: not found at this Tm", "N.d.", otScore) )
        else:
            primerTable.append( (baseName+"_R", prefixRev+rSeq, rTm[:4], otScore) )

    prefix = ""
    if batchName!="":
        prefix = "<strong>%s :</strong>" % batchName

    printBackLink()
    print("Contents:<br>")
    print("<ul>")
    print("<li><a href='#primers'>PCR Primers</a>")
    print("<li><a href='#ampliconSeqs'>Amplicon sequences</a></li>")
    print("<li><a href='#crispresso'>Crispresso support</a></li>")
    print("</ul>")
    print("<hr>")

    print(("<h2 id='primers'>%sPCR primers for off-targets of %s</h2>" % (prefix, guideSeqHtml)))
    print("<p>In the table below, Illumina Nextera Handle sequences have been added and highlighted in bold. Primers for the on-target have been added for convenience. The table below is sorted by the CFD off-target score. Sites with very low CFD scores &lt; 0.02 are unlikely to be cleaved, see our study <a href='http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2'>Haeussler et al. 2016</a>, Figure 2.</p>")
    print("<p>In the protocol by Matthew Canver, Harvard, two PCRs are run: one PCR to amplify the potential off-target, then a second PCR to extend the handles with Illumina barcodes. Please <a href='downloads/prot/canverProtocol.pdf'>click here</a> to download the protocol. Alternatively, you can have a look at <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3988262/#S8'>Fu et al, 2014</a>.</p>")
    print("<p>If a primer was not found, the reason is usually that the region around the off-target is too repetitive. "
        "To avoid unspecific primers, all repeats are masked for the primer design (not for off-target search). "
        "If you think that we should change the parameters here or should use different primer3 settings, please let "
        "us know.</p>")

    if pamId not in pamOtMatches:
        errAbort("pamId %s not valid" % pamId)

    print(('''<form id="paramForm" action="%s" method="GET">''' % basename(__file__)))

    printAmpLenAndTm(ampLen, tm)
    printHiddenFields(params, {"batchId":batchId, "pamId":pamId, "otPrimers":"1"})
    print("""<input type="submit" name="submit" value="Update">""")
    print("</form>")
    print("<p>")


    printPrimerTable(primerTable, withTm=True, withScore=True)

    print("<h2 id='ampliconSeqs'>Off-target amplicon sequences with primers</h2>")
    print("<p>These only list off-targets that have primers in the table above. Primers underlined, off-targets in bold.</p>")

    print("<table>")
    for score, seqName, primerInfo in scoredPrimers:
        print("<tr>")
        flankSeq = nameToSeq[seqName]
        lSeq, lTm, lPos, rSeq, rTm, rPos = primerInfo
        if lSeq is None:
            continue
        ampSeq = flankSeq[lPos:rPos+1]
        ulCoords = [(0, len(lSeq)), (rPos-lPos-len(rSeq), rPos-lPos+1)]
        boldCoords = [(1000-lPos, 1000-lPos+GUIDELEN+len(pam) )]
        print(("<td>%s</td> " % seqName))
        print("<td><tt>")
        print(markupSeq(ampSeq, ulCoords, boldCoords))
        print("</tt></td>")
        print("</tr>")
    print("<table>")

    #print("<h2>Input file for <a href='http://crispresso.rocks/'>Crispresso</a></h2>")
    print("<h2 id='crispresso'>Input file for Crispresso</h2>")
    print("<p><a href='http://crispresso.rocks/'>Crispresso</a>, written by Luca Pinello, is a software package to quantify the Cas9-induced mutations on off- or on-targets.</p>")

    downloadUrl = cgiGetSelfUrl({"download":"amplicons"})
    print(("<p><a href='%s'>Click here</a> to download an amplicon input file for Crispresso. For each off-target, it includes the off-target name, its PCR amplicon and the guide sequence. Keep a copy of this file.</p>" % downloadUrl))
    print("After sequencing, run CRISPRessoPooled. The tool will map the reads to the amplicons and analyse the mutations:<br>")
    fname = makeCrispressoFname(batchName, batchId)
    print(("<tt>CRISPRessoPooled -r1 Reads1.fastq.gz -r2 Reads2.fastq.gz -f %s --name MY_EXPERIMENT</tt></p>" % fname))

def printBackLink():
    " print a link back to the main batch page "
    newParams = {}
    newParams["batchId"] = cgiParams.get('batchId', "")
    if "batchName" in cgiParams:
        newParams["batchName"] = cgiParams["batchName"]
    paramStrs = ["%s=%s" % (key, val) for key, val in newParams.items()]
    paramStr = "&".join(paramStrs)
    url = basename(__file__)+"?"+paramStr

    print(("<p><a href='%s'>&larr; return to the list of all guides</a></p>" % url))

def microHomPage(params):
    " show the Bae et al microhomology sequences "
    printBackLink()
    batchId, pamId = params["batchId"], params["pamId"]
    inSeq, db, pam, position, extSeq = readBatchParams(batchId)
    pam = setupPamInfo(pam)

    guideSeq, pamSeq, pamPlusSeq, guideSeqWPam, guideStrand, guideSeqHtml, guideStart, guideEnd \
        = findGuideSeq(inSeq, pam, pamId)

    targetChrom, targetStart, targetEnd, strand = parsePos(position)
    #gStart = int(pamId.replace("s", "").replace("+","").replace("-",""))
    #gEnd = gStart+GUIDELEN
    strand = pamId[-1]

    scoreCode = params["showMh"]

    print("<h2>DNA break repair outcome predictions</h2>")
    print("<strong>Guide:</strong> %s<p>" % guideSeqHtml)

    #mutScores = crisporEffScores.calcMutSeqs(pamIds, longSeqs, enz, scoreNames=[scoreCode])
    outcome = readOutcomeData(batchId, scoreCode)[pamId]
    # extend guide to get +- 50 bp around it
    longSeq = getExtSeq(inSeq, guideStart, guideEnd, strand, 50-GUIDELEN, 50, extSeq=extSeq)
    if scoreCode=="oof":
        oof, mhSeqs = outcome
        mhSeqs.sort(reverse=True)
        print("""The following table lists possible deletions. """)

    elif scoreCode=="lindel":
        fsProb, mhSeqs = outcome
        print("""The following table lists possible deletions and insertions, scored by the <i>Lindel</i> repair model.<br> """)

    else:
        errAbort("Invalid score code")

    print("""Each sequence below represents the context around the guide's target, with deleted nucleotides shown as "-".<br>""")

    if scoreCode=="oof":
        print("""Sequences  are ranked by micro-homology score. This score is correlated to the likelihood of finding a particular deletion, as predicted by the method of <a target=_blank href='http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html'>Bae et al. 2014</a>.<p>""")
        scoreName = "mh- Score"

    elif scoreCode=="lindel":
        print("""Sequences are ranked by their probability, as predicted by the method of <a target=_blank href='https://academic.oup.com/nar/article/47/15/7989/5511473'>Chen et al. 2019</a>.<p>""")
        scoreName = "Probability"

    else:
        errAbort("Unknown score code")

    print("<table>")
    print("<tr>")
    print("<th>%s</th>" % scoreName)
    print("<th>Sequence</th>")
    print("<th>Effect</th>")
    print("</th>")

    if scoreCode=="oof":
        print("<tr>")
        print("<td></td>")
        print("<td><tt>%s</tt></td>" % longSeq[10:-10].upper())
        print("<td>(Wild-type sequence)</td>")
        print("</tr>")

    for row in mhSeqs:
        score, seq = row[:2]
        print("<tr>")
        if score==0:
            print("<td></td>")
        else:
            if scoreCode=="oof":
                print("<td>%d</td>" % score)
            else:
                print("<td>%.2f%%</td>" % score)

        print("<td><tt>%s</tt></td>" % seq)
        if scoreCode=="oof":
            delCount = seq.count("-")
            if delCount % 3 == 0:
                resStr = "no frameshift"
            else:
                resStr = "frameshift"
            print("<td><tt>%d bp deleted &rarr; %s</tt></td>" % (delCount, resStr))
        else:
            if score==0.0:
                print("<td>wild-type</td>") # for lindel
            else:
                print("<td>%s</td>" % (row[2]))
        print("</tr>")
    print("</table>")

def printSatMutPage(params):
    " special form for sat mutagenesis downloads, for our Nat Prot paper "
    defBarcode = "1"
    barcode = params.get("barcode", defBarcode)
    batchId = params["batchId"]

    batchInfo = readBatchAsDict(batchId)
    pam = batchInfo["pam"]

    printBackLink()

    print("<h3>Oligonucleotides for a Lentiviral Saturating Mutagenesis Screen</h3>")
    print("""
    <p>This page allows you to download all files for the ordering and the analysis of your oligonucleotide pool.<br>""")
    #<br> You can use Excel filters or command line programs like AWK to filter the list of oligos below
    #, e.g. to include only guides with a specificy score > 50.</p>

    print(("""<p><strong>1 - Subpool barcodes:</strong> The minimum order is
usually several thousand guides and it is cheaper to order more oligos. To reduce the cost per oligo, subsets can be
tagged with a "barcode" (unrelated to Illumina sequencing index barcodes) so they
can be selectively amplified from the pool.<br>
    <form id="paramForm" action="%(action)s" method="GET">
    Subpool barcode:
    """ % {"action": basename(__file__)}))

    printDropDown("barcode", satMutBarcodes, barcode)
    print("</p>")

    print("<p><strong>2 - PCR primers:</strong> Output file C includes two primers per target, for cleavage analysis with high-throughput sequencing. Primers are prefixed with Illumina Adapters. The forward primer prefix is TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG and the reverse prefix is GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG. These prefixes are already added to the primers in the Excel/Tab-sep tables.<br>")
    ampLen = "150"
    tm = "60"
    printAmpLenAndTm(ampLen, tm)
    print("</p>")

    print("<strong>3 - Filters</strong><br>Minimum specificity score: ")
    print("""<input id="minSpec" type="text" size="5" name="minSpec" value="30" /><br>""")
    print("""<small>A minimal value of 30 will remove only the guides in repeated regions. For screens, many researchers do not care a lot about off-targets. Increase this threshold if you want to more aggressively remove guides with many predicted off-targets.</small>""")
    print("<br>")

    effScoreName = "fusi"
    if pamIsSaCas9(pam):
        effScoreName = "najm"
    scoreLabel = scoreDescs.get(effScoreName, ["- Error: NO KNOWN SCORE FOR THIS ENZYME - "])[0]

    print(("Minimum %s efficiency score: " % scoreLabel))
    print("""<input id="minFusi" type="text" size="5" name="minFusi" value="10" /><br>""")
    print("""<small>A minimal value of 10 will only remove the least efficient guides. Increase this if want to enrich more for predicted high efficiency guides.</small>""")
    print("<br>")

    outTypes = [
        ("satMut", "A - saturating mutagenesis screen oligonucleotides"),
        ("targetSeqs", "B1 - Selection Analysis: guide sequences (CrispressoCount, .txt format)"),
        ("targetLocs", "B2 - Selection Analysis: guide locations (CrisprSurf, .csv format)"),
        ("ontargetPrimers", "C - On-target sequencing: PCR primers, one pair for each guide target"),
        ("ontargetAmplicons", "D - Cleavage Validation/Analysis: PCR Amplicons and guide sequence (CrispressoPooled)")
    ]

    formats = [
        ("html", "Do not download, display as webpage to get an idea"),
        ("xls", "Default: Excel for A and text for B/C/D"),
        ("tsv", "For programmers: always text"),
    ]

    print("<p><strong>4 - File format:</strong> Excel tables include a header with information how the oligonucleotides were constructed.<br>Text files have no header and are easier to process with other software.<br>Cleavage analysis files for Crispresso need to be in text format and are therefore not available as Excel files.<br>")

    print("File format:")
    printDropDown("format", formats, "html")
    print("<br>")

    printHiddenFields(params, {"satMut":None, "download":"useGet"})

    #print("The output table consists of the guide target sequences with their scores, the oligonucleotides to order and the possible sequencing primers for targeted PCR amplification of the target.")

    print("<p><strong>Output files:</strong><ul>")
    print("<li>A - Saturating Mutagenesis Oligonucleotides: the oligonucleotides to order from your Custom Oligonucleotide Array Supplier")
    print("<li>B - Selection Analysis: for every oligo from A, its sequence and genome location. For CrisprSurf quantification</li>")
    print("<li>C - On-target sequencing primers: one forward and one reverse primer for every target from B in the input sequence</li>")
    print("<li>D - Cleavage Analysis: for each pair of primers from C, a table with the PCR amplicon and the guide. For CrispressoPooled, to analyze DNA cleavage induced by this guide</ul>")

    print("You can click the four buttons below and save all four files.")

    print("</p>")
    print("""<p>""")
    print("""<input id="submitForm" style="" type="submit" name="get-satMut" value="Get A - oligo list">""")
    print("""<input id="submitForm" style="" type="submit" name="get-targetLocs" value="Get B - oligo locations">""")
    print("""<input id="submitForm" style="" type="submit" name="get-ontargetPrimers" value="Get C - target primers">""")
    print("""<input id="submitForm" style="" type="submit" name="get-ontargetAmplicons" value="Get D - target amplicons">""")

    print("</form>")

def printLibForm(params):
    """
    """
    sampleGenes = "PITX2\nMTOR\nTP53\nABO\n3661\nNM_134261"
    print(("""<form id="libForm" action="%s" method="POST">""" % basename(__file__)))

    #print("Organism: ")
    #printDropDown("org", [("human", "Human"),("mouse", "Mouse")], "human")
    #print("<br>")

    url = "crispor.py"
    print(("<p><a href='%s'>&larr; return to the CRISPOR main page</a></p>" % url))

    print("<h2>CRISPOR Batch Gene Targeting Assistant: Paste a list of genes to download a list of guides</h2>")

    print("Note: if you are planning a saturating mutagenesis screen, e.g. of a non-coding sequence, this is not the right tool. Submit your sequence on the normal <a href='crispor.py'>CRISPOR</a> page, then use the link 'Saturating mutagenesis' at the top of the guide table to get oligonucleotides of all guides in the input sequence.")

    print("Both tools are best used in conjunction with our wet-lab <a href='https://www.nature.com/articles/nprot.2018.005'>protocol</a> (<a href='http://biorxiv.org/content/early/2017/04/07/125245'>preprint</a>)<p>")

    print("<strong>Lentiviral screen library:</strong>")
    printDropDown("libName", libLabels, "geckov2")
    print("<br>")

    print("<strong>Number of guides per gene:</strong> ")
    printDropDown("guideCount", [(1,1), (2,2), (3,3), (4,4), (5,5), (6,6)], "3")
    print("<br>")

    print("<strong>Number of non-targeting control guides (max: 1000): </strong>")
    print("""<input id="ctrlCount" type="text" size="5" name="ctrlCount" value="10" />""")
    print("<br>")

    print("<strong>Subpool barcode: </strong>")
    printDropDown("barcode", satMutBarcodes, "1")
    print("<br>")

    print("<strong>Custom oligo prefix and suffix : </strong>")
    print("""<input id="custPrefix" type="text" size="20" name="custPrefix" value="" />""")

    #print("<strong>Custom oligo suffix (empty for defaults): </strong>")
    print("""<input id="custSuffix" type="text" size="20" name="custSuffix" value="" />""")
    print("both empty = defaults for cloning into pLentiGuide-Puro<p>")

    print(("""
    Enter a list of gene symbols, Entrez Gene IDs or Refseq IDs, one per line (case-insensitive):<br>
    <small>Type 'all' below to get all guides in the library and no gene filtering.</small>
    <small>You may need to <a href='https://discover.nci.nih.gov/matchminer/MatchMinerLookup.jsp'>convert old symbols</a>.</small>
    <small>To convert from UniProt IDs, try <a href="http://www.uniprot.org/uploadlists/">the UniProt converter.</a></small>
    <textarea tabindex="1" style="width:100%%" name="geneIds" rows="25" \
    placeholder="Paste a list of gene symbols here">%s</textarea><p>""" % sampleGenes))

    print("""<input id="submitGenes" type="submit" name="submit" value="Submit">""")
    #print('<input type="hidden" name="libDesign" value="1">')
    print("""</form>""")

def isInPar(db, chrom, start, end):
    """ return None if not in PAR or "1" or "2" if genome is hg19 or hg38 and chrom:start-end is in a PAR1/2 region """
    if db not in ("hg19", "hg38"):
        return None
    if not chrom in ("chrX", "chrY"):
        return None

    # all coordinates are from https://en.wikipedia.org/wiki/Pseudoautosomal_region
    # and look like they're 1-based
    if db=="hg38":
        if chrom=="chrX":
            if start >= 10001 and end < 2781479:
                return "1"
            if start >= 155701383 and end < 156030895:
                return "2"
        elif chrom=="chrY":
            if start >= 10001 and end < 2781479:
                return "1"
            if start >= 56887903 and end < 57217415:
                return "2"
    elif db=="hg19":
        if chrom=="chrX":
            if start >= 60001 and end < 2699520:
                return "1"
            if start >= 154931044 and end < 155260560:
                return "2"
        elif chrom=="chrY":
            if start >= 10001 and end < 2649520:
                return "1"
            if start >= 59034050 and end < 59363566:
                return "2"

    return None

def getControls(org):
    " return controls as a list "
    dbFname = "screenData/%s_controls.sqlite" % (org)
    conn = sqlite3.Connection(dbFname, timeout=10)
    cur = conn.cursor()
    sql = "SELECT guideSeq from guides"
    try:
        cur.execute(sql)
    except sqlite3.OperationalError:
        errAbort("Cannot open the file %s" % dbFname)

    rows = cur.fetchall()
    guideList = []
    for guideSeq in rows:
        guideList.append(guideSeq[0])
    return list(set(guideList))

def getLibGuides(org, libName, geneIdStr):
    " return a dict with geneId -> list of guide sequences "
    dbFname = "screenData/%s.sqlite" % (libName)
    conn = sqlite3.Connection(dbFname, timeout=10)
    cur = conn.cursor()

    guideData = defaultdict(list) # geneId -> guideRows
    orderedGenes = [] # ordered list of guideIds in guideData
    inSearchIds = geneIdStr.split("\n")
    notFoundGenes = []
    for inId in inSearchIds:
        inId = inId.strip().split()[0]
        if inId.lower()=="all":
            sql = "SELECT geneSym, entrezId, refseqId, guideSeq, pam from guides"
            try:
                cur.execute(sql)
            except sqlite3.OperationalError:
                errAbort("Cannot open the file %s" % dbFname)
        else:
            if inId.isdigit():
                searchField = "entrezId"
            elif inId.startswith("NM_"):
                searchField = "refseqId"
            else:
                searchField = "geneSym"

            sql = "SELECT geneSym, entrezId, refseqId, guideSeq, pam from guides WHERE %s = ? COLLATE NOCASE" % searchField
            try:
                cur.execute(sql, (inId,))
            except sqlite3.OperationalError:
                errAbort("Cannot open the file %s" % dbFname)

        rows = cur.fetchall()

        if len(rows)==0:
            notFoundGenes.append(inId)
            continue

        # reformat to dict gene -> seq
        doneGenes = set()
        for geneSym, entrezId, refseqId, guideSeq, pam in rows:
            guideData[geneSym].append( (entrezId, refseqId, guideSeq, pam) )
            if geneSym not in doneGenes:
                orderedGenes.append(geneSym)
                doneGenes.add(geneSym)

    return guideData, orderedGenes, notFoundGenes

def createGuideTable(lentiTemp, geneIdStr, guideCount, org, libName, barcodeId, controlCount, custPrefix, custSuffix):
    " write a table with lenti lib guides and oligos to temp/lenti/ and return its filename"

    keyStr = geneIdStr+str(guideCount)+org+libName+str(barcodeId)+str(controlCount)+custPrefix+custSuffix
    hasher = hashlib.sha1(keyStr.encode("latin1"))
    digest = hasher.digest()[0:20]
    lentiJobId = base64.urlsafe_b64encode(digest)
    lentiJobId = lentiJobId.decode("latin1").translate(transTab)

    outFname = join(lentiTemp, lentiJobId+".tsv")

    if isfile(outFname) and not DEBUG:
        return outFname

    ofh = open(outFname, "w")

    ofh.write("## CRISPOR %s library extractor\n" % versionStr)
    ofh.write("## CRISPOR JobID\t%s\n" % lentiJobId)
    ofh.write("## Organism\t%s\n" % org)
    ofh.write("## Library\t%s\n" % libName)
    ofh.write("## Number of guides per gene\t%d\n" % guideCount)

    barcodeDict = dict(satMutBarcodes)
    satMutOpt, optFields = buildPoolOptions(barcodeId, custPrefix, custSuffix)

    oligoPrefix, oligoSuffix = satMutOpt[:2]

    for label, val in optFields.items():
        ofh.write("## %s\t%s\n" % (label, val))

    guideData, orderedGenes, notFoundGenes = getLibGuides(org, libName, geneIdStr)
    ofh.write("## Not found genes\t%s\n" % ", ".join(notFoundGenes))

    ofh.write("guideId\tentrezId\trefseqId\tguideSeq\toligoSeq\tpam\n")

    for geneSym in orderedGenes:
        guideRows = guideData[geneSym]
        for guideId, (entrezId, refseqId, guideSeq, pam) in enumerate(guideRows):
            row = ["%s_%d" % (geneSym, guideId+1), str(entrezId), refseqId,
                guideSeq, "%s%s%s" % (oligoPrefix, guideSeq, oligoSuffix), pam]
            ofh.write("\t".join(row))
            ofh.write("\n")

            if guideId+1 >= guideCount:
                break

    ctrls = getControls(org)
    for ctrlId, guideSeq in enumerate(ctrls[:controlCount]):
        row = ["control_%d" % (ctrlId), "", "", guideSeq, "%s%s%s" % (oligoPrefix, guideSeq, oligoSuffix), ""]
        ofh.write("\t".join(row))
        ofh.write("\n")

    ofh.close()

    return ofh.name

def printLibGuides(params):
    """
    print a table with the guides we have on file for a list of gene IDs
    """
    lentiTemp = join(batchDir, "lenti")
    if not isdir(lentiTemp):
        os.makedirs(lentiTemp)

    geneIdStr = cgiGetStr(params, "geneIds", "").strip()
    guideCount = cgiGetNum(params, "guideCount", 3)
    #org = cgiGetStr(params, "org", "human")
    libName = cgiGetStr(params, "libName", "geckov2")
    barcodeId = cgiGetNum(params, "barcode", 1)
    org = libName.split("_")[0]
    controlCount = cgiGetNum(params, "ctrlCount", 10)
    custPrefix = cgiGetStr(params, "custPrefix", "")
    custSuffix = cgiGetStr(params, "custSuffix", "")

    # need to check libName, as it's used to open a file
    validNames = dict(libLabels)
    if libName not in validNames:
        errAbort("Invalid library name")

    tabFname = createGuideTable(lentiTemp, geneIdStr, guideCount, org, libName, barcodeId, controlCount, custPrefix, custSuffix)

    for line in open(tabFname):
        if "Not found genes" in line:
            notFoundGenes = line.strip("\n").split("\t")[1].split(",")

    url = "crispor.py?libDesign=1"
    print(("<p><a href='%s'>&larr; return to the CRISPOR Batch input page</a></p>" % url))

    #print("Organism: %s<br>" % org)
    libLabel = dict(libLabels)[libName].split("(")[0] # strip the 'recommended' note
    print(("<strong>Library:</strong> %s<br>" % libLabel))
    print(("<strong>Number of guides per gene:</strong> %d<br>" % guideCount))
    print(("<strong>Number of non-targeting controls:</strong> %d (all controls are from the GeCKOV2 library)<br> " % controlCount))

    barcodeDict = dict(satMutBarcodes)
    satMutOpt, optFields = buildPoolOptions(barcodeId, custPrefix, custSuffix)
    oligoPrefix, oligoSuffix = satMutOpt[:2]

    for label, val in optFields.items():
        print(("<strong>%s:</strong> %s<br>\n" % (label, val)))
    print("<p>")

    print("For details on these sequences, see our <a href='https://www.nature.com/articles/nprot.2018.005'>protocol</a> (<a href='https://www.biorxiv.org/content/early/2017/09/11/125245'>preprint</a>).<p>")

    if len(notFoundGenes)!=0 and notFoundGenes!=[""]:
        print(("Input gene identifiers that were not found: %s<p>" % ",".join(notFoundGenes)))

    if guideCount>4 and "gecko" not in libName:
        print(("Note: you asked for %d guides per gene but this library includes only four guides per gene, so the maximum number of guides per genes below is four.<p>" % guideCount))

    print(("<a href='%s'>Download table</a><p>" % relpath(tabFname, dirname(abspath(__file__)))))

    print('<table class="libTable">')
    print("<tr><th style='width:10em'>ID of guide</th>")
    print("<th style='width:10em'>Target Entrez ID</th>")
    print("<th style='width:10em'>Target Refseq ID</th>")
    print("<th style='width:14em'>Guide RNA<br>(click to show in CRISPOR)</th>")
    print(("<th style='width:40em'>Full oligonucleotide including guide RNA<br><small>%s</small></th>" % optFields["Oligonucl. structure"]))
    print("</tr>")

    genomeDbs = {
        "human" : "hg19",
        "mouse" : "mm10"
    }
    genomeDb = genomeDbs.get(org)

    for row in lineFileNext(open(tabFname)):
        print('<tr>')
        print(('<td>%s</td>' % (row.guideId)))
        print(('<td>%s</td>' % (row.entrezId)))
        print(('<td>%s</td>' % (row.refseqId)))
        if row.pam!="":
            print(('<td><tt><a target=_blank href="crispor.py?org=%s&seq=%s&pam=NGG">%s</a></tt></td>' % (genomeDb, row.guideSeq+row.pam, row.guideSeq)))
        else:
            print(('<td><tt>%s</tt></td>' % row.guideSeq))
        print(('<td><tt>%s</tt></td>' % (row.oligoSeq)))
        print('</tr>')
    print('</table>')

def printBody(params):
    " main dispatcher function "

    # TODO: first, if batchId is the only parameters,
    # we will check for a flag file to see if the job is running and output the status file if it is.

    global doCfdFix
    if "fixCfd" in params:
        doCfdFix = True

    if "batchId" in params and not "satMut" in params:
        printCrisporBodyStart()
        if "pamId" in params:
            if "pam" in params:
                primerDetailsPage(params)
            elif "otPrimers" in params:
                otPrimerPage(params)
            elif "showMh" in params:
                microHomPage(params)
            else:
                errAbort("Unrecognized CGI parameters.")
        else:
            crisprSearch(params)

    elif "satMut" in params:
        printCrisporBodyStart()
        printSatMutPage(params)

    elif ("seq" in params or "pos" in params) and "org" in params and "pam" in params:
        printCrisporBodyStart()
        crisprSearch(params)
    elif "libDesign" in params:
        printCrisporBodyStart()
        printLibForm(params)
    elif "geneIds" in params:
        printCrisporBodyStart()
        printLibGuides(params)
    else:
        printForm(params)

def iterParseBoulder(tmpOutFname):
    " parse a boulder IO style file, as output by Primer3 "
    data = {}
    for line in open(tmpOutFname):
        if line == "=\n":
            if data!={}:
                yield data
                data = {}
                continue
        key, val = line.rstrip("\n").split("=")
        data[key] = val
    if data!={}:
        yield data

def runPrimer3(seqs, targetStart, targetLen, prodSizeRange, tm, addTags):
    """
    input is a list of (seqId, seq)
    return dict seqId -> (primerseq1, tm1, pos1, primerseq2, tm2, pos2)
    addTags is a dict tag -> value
    """
    p3OutFh = makeTempFile("primer3Out", ".txt")

    minTm = tm-3
    maxTm = tm+3
    optTm = tm


    p3InFh = makeTempFile("primer3In", ".txt")
    # values from https://www.ncbi.nlm.nih.gov/tools/primer-blast/
    # MAX_END_STABILITY is strange but it seems to be set that way by NCBI
    primer3ConfigDir = abspath(join(baseDir, "bin", "src", "primer3-2.3.6", "src", "primer3_config"))
    # PRIMER_MAX_POLY suggested by Yueh-Chiang.Hu@cchmc.org

    # PRIMER_MAX_POLY_X=4 suggested by Chiang.Hu@cchmc.org
    conf = """PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_MIN_TM=%(minTm)d
PRIMER_OPT_TM=%(optTm)d
PRIMER_MAX_TM=%(maxTm)d
PRIMER_MAX_POLY_X=4
PRIMER_MIN_SIZE=18
PRIMER_OPT_SIZE=20
PRIMER_MAX_SIZE=25
PRIMER_MAX_END_STABILITY=9
PRIMER_MAX_POLY=4
PRIMER_PRODUCT_SIZE_RANGE=%(prodSizeRange)s
SEQUENCE_TARGET=%(targetStart)s,%(targetLen)s
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%(primer3ConfigDir)s/
""" % locals()

    for seqId, seq in seqs:
        seqId = seqId.replace("=", " ")
        p3InFh.write("SEQUENCE_ID=%s\n" % seqId)
        p3InFh.write("SEQUENCE_TEMPLATE=%s\n" % seq)
        p3InFh.write(conf)
        for key, val in addTags.items():
            p3InFh.write("%s=%s\n" % (key, val))
        p3InFh.write("=\n")

    p3InFh.flush()

    #shutil.copy(p3InFh.name, "/tmp/primer3.tmp")

    binName = join(binDir, "primer3_core")
    cmdLine = binName+" %s > %s" % (p3InFh.name, p3OutFh.name)
    runCmd(cmdLine, ignoreExitCode=True)

    ret = {}
    #print "<br>".join(open(p3OutFh.name).read().splitlines())
    for p3 in iterParseBoulder(p3OutFh.name):
        if "PRIMER_LEFT_NUM_RETURNED" not in p3 or "PRIMER_RIGHT_NUM_RETURNED" not in p3 or \
            p3["PRIMER_LEFT_NUM_RETURNED"]=="0" or p3["PRIMER_RIGHT_NUM_RETURNED"]=="0":
        #if "PRIMER_LEFT_0_SEQUENCE" not in p3 or "PRIMER_RIGHT_0_SEQUENCE" not in p3:
            #errAbort("No primer found for this Tm. Please contact us if you think there should be a primer found.")
            ret[seqId] = None, None, None, None, None, None
            continue

        seqId = p3["SEQUENCE_ID"]
        seq1 = p3["PRIMER_LEFT_0_SEQUENCE"]
        seq2 = p3["PRIMER_RIGHT_0_SEQUENCE"]
        tm1 = p3["PRIMER_LEFT_0_TM"]
        tm2 = p3["PRIMER_RIGHT_0_TM"]
        pos1 = int(p3["PRIMER_LEFT_0"].split(",")[0])
        pos2 = int(p3["PRIMER_RIGHT_0"].split(",")[0])
        ret[seqId] = seq1, tm1, pos1, seq2, tm2, pos2

    # make sure that unsuccessful sequences are not missing from the response
    for seqId, seq in seqs:
        if seqId not in ret:
            ret[seqId] = None, None, None, None, None, None

    return ret

def parseFastaAsList(fileObj):
    " parse a fasta file, return list (id, seq) "
    seqs = []
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n").rstrip("\r")
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
        line = line.rstrip("\n").rstrip("\r")
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

def findPerfectMatch(batchId):
    """ find best match for input sequence from batchId in genome and return as
    a string chrom:start-end:strand or "?" if not found "
    """
    if skipAlign:
        return "?"

    batchInfo = readBatchAsDict(batchId)
    seq = batchInfo["seq"]
    genome = batchInfo["org"]

    # write seq to tmp file
    #tmpFaFh = tempfile.NamedTemporaryFile(dir=TEMPDIR, prefix="crisporBestMatch", suffix=".fa")
    tmpFaFh = makeTempFile(prefix="bwaswInput", suffix=".fa")
    tmpFaFh.write(">tmp\n%s" % seq)
    tmpFaFh.flush()
    logging.debug("seq: %s" % open(tmpFaFh.name).read())
    faFname = tmpFaFh.name

    # create temp SAM file
    #tmpSamFh = tempfile.NamedTemporaryFile(dir=TEMPDIR, prefix="crisporBestMatch", suffix=".sam")
    tmpSamFh = makeTempFile(prefix="bwaswOutput", suffix=".sam")
    samFname = tmpSamFh.name

    bwaIndexPath = abspath(join(genomesDir, genome, genome+".fa"))
    remoteAddr = pipes.quote(os.environ.get("REMOTE_ADDR", "noIp")) # so I can figure out is someone is hammering the server
    cmd = "true %(batchId)s %(remoteAddr)s && $BIN/bwa bwasw -b 100 -q 100 -T 20 %(bwaIndexPath)s %(faFname)s > %(samFname)s" % locals()
    runCmd(cmd)

    chrom, start, end = None, None, None
    logging.debug("Parsing SAM file %s" % samFname)
    matchByChrom = defaultdict(list)
    for l in open(samFname):
        if l.startswith("@"):
            continue
        l = l.rstrip("\n")
        fs = l.split("\t")
        logging.debug("SAM input-line: %s" % repr(fs))
        qName, flag, rName, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = fs[:11]
        logging.debug("qName=%s, flag=%s, rName=%s, pos=%s, mapq=%s, cigar=%s" % \
            (qName, flag, rName, pos, mapq, cigar))
        if (int(flag) and 2) == 2:
            strand = "-"
        else:
            strand = "+"
        if not re.compile("[0-9]*").match(cigar):
            logging.debug("CIGAR is not number")
            continue
        if cigar=="*":
            logging.debug("CIGAR is *")
            continue
            #errAbort("Sequence not found in genome. Are you sure you have pasted the correct sequence and also selected the right genome?")
        # Todo: why do we get soft clipped sequences from BWA? repeats?
        if "S" in cigar:
            logging.debug("match found, but soft-clipped, cigar: %s" % cigar)
            continue
        cleanCigar = cigar.replace("M","")
        if not cleanCigar.isdigit():
            logging.debug("match found, but cigar string was %s" % cigar)
            continue
        matchLen = int(cleanCigar)
        chrom, start, end =  rName, int(pos)-1, int(pos)-1+matchLen # SAM is 1-based
        assert("|" not in chrom) # We do not allow '|' in chrom name. I use this char to sep. info fields in BED.
        matchByChrom[chrom].append( (chrom, start, end, strand) )

    # second pass, to handle the PAR matches properly
    matches = []
    for chrom, matchList in matchByChrom.items():
        if isInPar(genome, chrom, start, end) is not None:
            # only keep matches on chrY
            if chrom=="chrX" and "chrY" in matchByChrom:
                logging.debug("In PAR region, so skipping chrX")
                continue
        for m in matchList:
            matches.append(m)

    # delete the temp files
    tmpSamFh.close()
    tmpFaFh.close()

    if len(matches)==0:
        return "?"

    nonAltMatches = [x for x in matches if not isAltChrom(x[0])]
    if len(nonAltMatches)!=0:
        bestMatch = nonAltMatches[0]
    else:
        bestMatch = matches[0]

    logging.debug("Found %d best matches, %d on non-alts. matches: %s, best match %s" % (len(matches), len(nonAltMatches), matches, bestMatch))
    return "%s:%d-%d:%s" % (bestMatch)

def maskLowercase(seq):
    " replace all lowercase letters with 'N' "
    newSeq = []
    for c in seq:
        if c.islower():
            newSeq.append("N")
        else:
            newSeq.append(c)
    return "".join(newSeq)

def getGenomeSeqsBin(genome, coordList, doRepeatMask=False):
    """ use the binary twoBitToFa to get seqs:
    coordList has format (chrom, start, end, name, score, strand)
    returns list (chrom, start, end, name, seq)
    """
    twoBitPath = join(genomesDir, "%(genome)s/%(genome)s.2bit" % locals())
    twoBitPath = abspath(twoBitPath)

    #bedFh = makeTempFile("getGenomeSeqs", ".bed")
    bedFh = open("/var/tmp/primer3In8w4woj1a.txt", "w")
    faFh = makeTempFile("getGenomeSeqs", ".fa")

    for row in coordList:
        line = "\t".join([str(x) for x in row])
        bedFh.write(line)
        bedFh.write("\n")

    bedFh.flush()
    bedFh.close()

    #print("output")
    #for line in open(bedFh.name):
        #print line
    #print("output2")

    cmd = ["$BIN/twoBitToFa","-bed="+bedFh.name, twoBitPath, faFh.name]
    runCmd(cmd, useShell=False)

    faSeqs = parseFastaAsList(open(faFh.name))

    assert(len(faSeqs)==len(coordList))

    seqs = []
    for coordTuple, seqTuple in zip(coordList, faSeqs):
        seqId, seq = seqTuple
        if len(coordTuple)==4:
            chrom, start, end, name = coordTuple
            strand = "+"
            score = "0"
        else:
            chrom, start, end, name, score, strand = coordTuple

        #seq = tbf[chrom][start:end]
        if strand=="-":
            seq = revComp(seq)
        if doRepeatMask:
            seq = maskLowercase(seq)
        seqs.append( (chrom, start, end, name, score, strand, seq) )
    return seqs

def getGenomeSeqs(genome, coordList, doRepeatMask=False):
    """ DOES NOT WORK ON PYTHON3 ANYMORE - NOT USED ANYMORE, as long twobitreader has not been updated

    return dict of genome sequences,
    coordList has format (chrom, start, end, name, score, strand)
    returns list (chrom, start, end, name, seq)
    """
    twoBitPath = join(genomesDir, "%(genome)s/%(genome)s.2bit" % locals())
    twoBitPath = abspath(twoBitPath)

    tbf = twobitreader.TwoBitFile(twoBitPath)
    seqs = []
    for coordTuple in coordList:
        if len(coordTuple)==4:
            chrom, start, end, name = coordTuple
            strand = "+"
            score = "0"
        else:
            chrom, start, end, name, score, strand = coordTuple

        seq = tbf[chrom][start:end]
        if strand=="-":
            seq = revComp(seq)
        if doRepeatMask:
            seq = maskLowercase(seq)
        seqs.append( (chrom, start, end, name, score, strand, seq) )
    return seqs

def designPrimer(
    genome,
    chrom,
    start,
    end,
    strand,
    guideStart,
    batchId,
    ampLen,
    tm,
    hdrDist = None,
):
    " create primer for region around chrom:start-end, write output to batch "
    " returns (leftPrimerSeq, lTm, lPos, rightPrimerSeq, rTm, rPos, amplified sequence)"

    # get 1kbp of flanking sequence
    flankSeq = getFlankSeq(genome, chrom, start, end)

    targetStart, targetLen, ampRange, addTags = getTargetForPrimerDesign(guideStart, ampLen, hdrDist, strand)

    primers = runPrimer3([("seq1", flankSeq)], targetStart, targetLen, ampRange, tm, addTags)
    lSeq, lTm, lPos, rSeq, rTm, rPos = list(primers.values())[0]

    if lSeq==None or rSeq==None:
        return None, None, None, None, None, None, flankSeq, ampRange, flankSeq, addTags

    targetSeq = flankSeq[lPos:rPos+1]
    if 'N' in targetSeq:
        # Get the flankSeq again but without the Ns... see maskLowercase
        flankSeq = getFlankSeq(genome, chrom, start, end, doRepeatMask=False)
        targetSeq = flankSeq[lPos:rPos+1]

    return lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq, ampRange, flankSeq, addTags


def getFlankSeq(genome, chrom, start, end, doRepeatMask=True):
    flankStart = start - 1000
    flankEnd = end + 1000

    chromSizes = parseChromSizes(genome)

    if flankStart < 0 or flankEnd > chromSizes[chrom]:
        errAbort("Not enough space on genome sequence to design primer. Need at least 1kbp on each side of the input sequence to design primers. Please design primers manually, choose a more recent genome assembly with longer contig sequences or paste a shorter input sequence (e.g. just the guide sequence alone with the PAM). Still questions? Email %s" % contactEmail)

    # get 1kbp of flanking sequence
    flankSeq = getGenomeSeqsBin(
        genome,
        [(chrom, flankStart, flankEnd, "seq")],
        doRepeatMask=doRepeatMask
    )[0][-1]
    return flankSeq


def getTargetForPrimerDesign(guideStart, ampLen, hdrDist, strand):
    if hdrDist is not None:
        # Cut for Cas9 is always between 3 and 4 nt from the PAM
        if strand == '+': # '+' here is guide strand relative to the genome
            cutAt = guideStart + 17
        else:
            cutAt = guideStart + 6
        # hdrDist is defined as the distance to the right of the cut to the insert
        insertAt = cutAt - hdrDist

        homologyArmLen = 55
        repairBuffer = 50
        targetStart = 1000 + insertAt - homologyArmLen - repairBuffer
        targetLen = homologyArmLen * 2 + repairBuffer * 2  # ==210
        # Assuming 250bp paired reads and 90bp insert
        ampRange = '250-310'
    else:
        # try to get some good heuristics for the primer placement
        # primers must not overlap the target but also not be
        # too far away, especially when using next-gen sequencing
        targetStart = 1000+guideStart
        targetLen = 23
        if ampLen<=250:
            ampDist = 50
        else:
            ampDist = 150
        ampRange = str(ampLen-ampDist)+"-"+str((ampLen))

    # for long reads: make sure that the primers are not too close to the target
    addTags = {}
    if ampLen >= 600:
        addTags = {
            "SEQUENCE_EXCLUDED_REGION" : "%d,%d" % (targetStart-150, targetLen+300),
            #"SEQUENCE_EXCLUDED_REGION" : "0,1500"
        }

    if hdrDist is not None:
        # Primer products should be as short as possible while staying at or above 250bp.
        # See also prodSizeRange and PRIMER_PRODUCT_SIZE_RANGE.
        addTags.update({
            # The optimum size for the PCR product.
            'PRIMER_PRODUCT_OPT_SIZE': ampRange.split('-')[0], # 250
            # Penalty weight for products longer than PRIMER_PRODUCT_OPT_SIZE.
            'PRIMER_PAIR_WT_PRODUCT_SIZE_GT': 1,
        })

    return targetStart, targetLen, ampRange, addTags


def markupSeq(seq, ulPosList, boldPosList, annots = {}):
    """ return seq as html with some parts underlined or in bold.
    annots is a dict with (start,end) -> dict with keys like "color"
    """
    annotStarts = {}
    annotEnds = defaultdict(set)
    for (start, end), aDict in annots.items():
        annotStarts[start] = aDict
        aDict["end"] = end

    ulStarts = set([x[0] for x in ulPosList])
    ulEnds = set([x[1] for x in ulPosList])
    boldStarts = set([x[0] for x in boldPosList])
    boldEnds = set([x[1] for x in boldPosList])
    ret = []
    openAnnots = defaultdict(int) # current number of open spans, per cssString
    openTags = set()
    for i, nucl in enumerate(seq):
        if i in annotEnds:
            for tagStr in annotEnds[i]:
                if tagStr in openAnnots:
                    openAnnots[tagStr]-=1
                    if openAnnots[tagStr]==0:
                        ret.append("</span>")
                        del openAnnots[tagStr]

        if i in annotStarts:
            aDict = annotStarts[i]
            cssParts = []
            for key, val in aDict["css"].items():
                cssParts.append("%s:%s" % (key, val))
            cssStr = ";".join(cssParts)
            tagStr = "<span style='%s'>" % cssStr
            if not tagStr in openAnnots:
                ret.append(tagStr)
            openAnnots[tagStr]+=1
            annotEnds[aDict["end"]].add(tagStr)

        if i in ulStarts:
            ret.append("<u>")
            openTags.add("u")
        if i in ulEnds:
            ret.append("</u>")
            if "u" in openTags:
                openTags.remove("u")
        if i in boldStarts:
            ret.append("<b>")
            openTags.add("b")
        if i in boldEnds:
            ret.append("</b>")
            if "strong" in openTags:
                openTags.remove("b")
        ret.append(nucl)
        if (i+1) % 80==0:
            ret.append("<br>")
    for tag in openTags:
        ret.append("</%s>" % tag)
    return "".join(ret)
    #return seq[:start]+"<u>"+seq[start:end]+"</u>"+seq[end:]

def makeHelperPrimers(guideName, guideSeq, plasmid, pam):
    " return dict with various names -> primer for primer page "
    primers = defaultdict(list)
    guideRnaFw = "guideRna%sT7sense" % guideName
    guideRnaRv = "guideRna%sT7antisense" % guideName

    # T7 plasmids
    if guideSeq.lower().startswith("g"):
        primers["T7"].append((guideRnaFw, "TAG<b>%s</b>" % guideSeq))
        primers["T7"].append((guideRnaRv, "AAAC<b>%s</b>" % revComp(guideSeq[1:])))
    else:
        primers["T7"].append((guideRnaFw, "TAGG<b>%s</b>" % guideSeq))
        primers["T7"].append((guideRnaRv, "AAAC<b>%s</b>" % revComp(guideSeq)))

    # T7 in-vitro
    prefix = ""
    if not guideSeq.lower().startswith("g"):
        prefix = "G"
    #specPrimer = "TAATACGACTCACTATA%s<b>%s</b>GTTTTAGAGCTAGAAATAGCAAG" % (prefix, guideSeq)

    if pamIsSaCas9(pam):
        # See http://dx.doi.org/10.1016/j.cell.2015.08.007
        specPrimer = "GAAATTAATACGACTCACTATA%s<b>%s</b>GTTTTAGTACTCTGGAAACAGAA" % (prefix, guideSeq)
        tracrPrimer = "GTTTTAGTACTCTGGAAACAGAATCTACTAAAACAAGGCAAAATGCCGTGTTTATCTCGTCAACTTGTTGGCGAGATTTTTTT"
    else:
        specPrimer = "GAAATTAATACGACTCACTATA%s<b>%s</b>GTTTTAGAGCTAGAAATAGCAAG" % (prefix, guideSeq)
        tracrPrimer = "AAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC"

    primers["T7iv"].append(("guideRNA%sT7crTarget" % guideName, specPrimer))
    primers["T7iv"].append(("guideRNAallT7common (constant primer used for all guide RNAs)", tracrPrimer))

    # geneart primers
    primers["geneArt"].append(("guideRNA%sGeneArtFw" % guideName, "TACGACTCACTATAG<b>"+guideSeq+"</b>"))
    primers["geneArt"].append(("guideRNA%sGeneArtRev" % guideName, "TTCTAGCTCTAAAAC<b>"+revComp(guideSeq)+"</b>"))

    # U6 - mammalian cells
    u6Prefix = ""
    if not guideSeq.lower().startswith("g"):
        u6Prefix = "gN20-"
    fwName = "%sguideRNA%sU6sense" % (u6Prefix, guideName)
    revName = "%sguideRNA%sU6antisense" % (u6Prefix, guideName)

    if pamIsCpf1(pam):
        primers["mammCells"].append((fwName, "AGAT<b>%s</b>" % guideSeq))
        primers["mammCells"].append((revName, "AAAA<b>%s</b>" % revComp(guideSeq)))
    else:
        if guideSeq.lower().startswith("g"):
            addGPrefix = ""
            addCSuffix = ""
        else:
            addGPrefix = "<u>G</u>"
            addCSuffix = "<u>C</u>"
            primers["mammCellsNote"] = True
        guideSeqTrunc = guideSeq[1:]


        u6FwPrefix, u6RwPrefix, u6Suffix = addGenePlasmidInfo[plasmid][:3]
        for pi in addGenePlasmids:
            if pi[0]==plasmid:
                plasmidLabel = pi[1][1]
        primers["mammCells"].append((fwName+plasmidLabel, "%s%s<b>%s</b>%s" % (u6FwPrefix, addGPrefix, guideSeq, u6Suffix)))
        primers["mammCells"].append((revName+plasmidLabel, "%s<b>%s</b>%s%s" % (u6RwPrefix, revComp(guideSeq), addCSuffix, u6Suffix)))

        if not guideSeq.lower().startswith("g"):
            primers["mammCells19"].append(("gN19-"+fwName+plasmidLabel, "%s%s<b>%s</b>%s" % (u6FwPrefix, addGPrefix, guideSeqTrunc, u6Suffix)))
            primers["mammCells19"].append(("gN19-"+revName+plasmidLabel, "%s<b>%s</b>%s%s" % (u6RwPrefix, revComp(guideSeqTrunc), addCSuffix, u6Suffix)))


        #if guideSeq.lower().startswith("g"):
            #primers["mammCells"].append((fwName, "ACACC<b>%s</b>G" % guideSeq))
            #primers["mammCells"].append((revName, "AAAAC<b>%s</b>G" % revComp(guideSeq)))
        #else:
            #primers["mammCells"].append((fwName, "ACACC<u>G</u><b>%s</b>G" % guideSeq))
            #primers["mammCells"].append((revName, "AAAAC<b>%s</b><u>C</u>G" % revComp(guideSeq)))
            #primers["mammCellsNote"] = True

    # add the prefix to all names
    newPrimers = {}
    for key, primList in primers.items():
        if key.endswith("Note"):
            newPrimers[key] = primList
            continue

        newList = []
        for name, seq in primList:
            if batchName!="":
                newName = batchName+"_"+name
            else:
                newName = name
            newList.append( (newName, seq) )
        newPrimers[key] = newList
    return newPrimers

def printPrimerTableAll(primers):
    print('<table class="primerTable">')
    for key, primerList in primers.items():
        if key.endswith("Note"):
            continue
        for name, seq in primerList:
            name = name.split()[0]
            print('<tr>')
            print("<td>%s</td>" % name)
            print("<td><tt>%s</tt></td>" % seq)
            print("</tr>")
    print("</table>")

def printPrimerTable(primerList, onRows=False, withTm=False, withScore=False, seqType="Primer"):
    " given a list of (name, seq) tuples, print a table "
    print('<table class="primerTable">')
    print(("<tr><th>Name</th><th>%s Sequence</th>" % seqType))
    if withTm:
        print("<th>Tm</th>")
    if withScore:
        print("<th>CFD Score</th>")
    print("</tr>")
    for row in primerList:
        name, seq = row[:2]
        if onRows:
            print("<tr><td>%s</td></tr>" % name)
            print("<tr><td><tt>%s</tt></td></tr>" % seq)
        else:
            print('<tr>')
            print("<td>%s</td>" % name)
            print("<td><tt>%s</tt></td>" % seq)
            if withTm:
                tm = row[2]
                print("<td>%s</td>" % tm)
            if withScore:
                score = row[3]
                print("<td>%s</td>" % score)
            print("</tr>")
    print("</table>")

#def printDropboxLink():
    #print """
    #<a id="dropbox-link"></a>
    #<script type="text/javascript" src="https://www.dropbox.com/static/api/2/dropins.js" id="dropboxjs" data-app-key="lp7njij87kjrcxv"></script>");
    #<script type="text/javascript">
    #options = {
        #success: function(files) {
           #$('textarea[name=hgct_customText]').val(files[0].link);
           #alert('Here is the file link: ' + files[0].link);
        #},
        #cancel: function() {},
        #linkType: 'direct',
        #multiselect: true,
    #};
    #var button = Dropbox.createChooseButton(options);
    #document.getElementById('dropbox-link').appendChild(button);
    #</script>
    #"""

def mergeParamDicts(params, changeParams):
    """ changeParams is a dict that can override elements in params.
    if value==None in changeParams, the whole element will get removed.
    if onlyParams is set, only copy over the keys in onlyParams (a list)
    """
    newParams = {}
    newParams.update(params)
    newParams.update(changeParams)
    for key, val in changeParams.items():
        if val==None:
            del newParams[key]
    return newParams

def printHiddenFields(params, changeParams):
    " "
    newParams = mergeParamDicts(params, changeParams)
    for key, val in newParams.items():
        print(('<input type="hidden" name="%s" value="%s">' % (key, val)))

def cgiGetSelfUrl(changeParams, anchor=None, onlyParams=None):
    """ create a URL to myself with different parameters than what we have.
    changeParams is a dict with the new arguments.
    onlyParams is an optional list of CGI variables that should be exported.
    """
    if onlyParams:
        cgiSubs = {}
        for x in onlyParams:
            if x not in cgiParams:
                continue
            cgiSubs[x] = cgiParams.get(x)
        newParams = mergeParamDicts(cgiSubs, changeParams)
    else:
        newParams = mergeParamDicts(cgiParams, changeParams)
    paramStrs = ["%s=%s" % (key, urllib.parse.quote(val)) for key, val in newParams.items()]
    paramStr = "&".join(paramStrs)
    url = basename(__file__)+"?"+paramStr
    if anchor is not None:
        url += "#"+anchor
    return url

def printDropDown(name, nameValList, default, onChange=None, style=None):
    """ print a dropdown box and set a default """
    addStr = ""
    if onChange is not None:
        addStr = """ onchange="%s" """ % onChange
    addStr2 = ""
    if style is not None:
        addStr2 = """ style='%s' """ % style

    print(('<select id="dropdown" name="%s"%s%s>' % (name, addStr,addStr2)))
    for name, desc in nameValList:
        name = str(name)
        addString = ""
        if default is not None and str(name)==str(default):
            addString = ' selected="selected"'
        print(('   <option value="%s"%s>%s</option>' % (name, addString, desc)))
    print('</select>')

def findGuideSeq(inSeq, pam, pamId):
    """ given the input sequence and the pamId, return the guide sequence,
    the sequence with the pam and its strand.
    """
    startDict, endSet = findAllPams(inSeq, pam)
    pamInfo = list(flankSeqIter(inSeq, startDict, len(pam), False))
    for guidePamId, pamStart, guideStart, guideStrand, guideSeq, pamSeq, pamPlusSeq in pamInfo:
        if guidePamId!=pamId:
            continue

        guideSeqWPam = concatGuideAndPam(guideSeq,pamSeq)
        # prettify guideSeqWPam to highlight the PAM
        if pamIsFirst:
            guideSeqHtml = "<i>%s</i> %s" % \
                (guideSeqWPam[:len(pam)].upper(), guideSeqWPam[len(pam):].upper())
        else:
            guideSeqHtml = "%s <i>%s</i>" % \
                (guideSeqWPam[:-len(pam)].upper(), guideSeqWPam[-len(pam):].upper())

        guideEnd = guideStart + GUIDELEN
        return guideSeq, pamSeq, pamPlusSeq, guideSeqWPam , guideStrand, guideSeqHtml, \
                guideStart, guideEnd
    errAbort("pamId %s not found? This is a bug." % pamId)

def findOntargetPos(otMatches, pamId, position, absentOk=False):
    " find position of guide sequence in genome at MM0 "
    if pamId not in otMatches or 0 not in otMatches[pamId]:
        if absentOk:
            return None, None, None, None, None, False
        else:
            errAbort("No perfect match found for guide sequence in the genome. Cannot design primers for a non-matching guide sequence.<p>Are you sure you have selected the right genome? <p> If you have selected the right genome and entered a cDNA as the query sequence, please note that sequences that overlap a splice site are not part of the genome and cannot be used as guide sequences.")

    matchList = otMatches[pamId][0] # = get all matches with 0 mismatches
    isUnique = True
    if len(matchList)>1:
        # this only gets executed if there are multiple exact matches for the target
        # we iterate over all offs and try find the correct one.
        targetChrom, targetStart, targetEnd, strand = parsePos(position)

        filtMatch = None
        # search for off-target that is the on-target
        for match in matchList:
            # example match: ('scaffold_1', 578, 601, 'TATTGGATTGGTCCAATCGTTGG', '-', 'ex', 'GSADVT00000001001', 293)
            chrom, start, end = match[:3]
            if chrom==targetChrom and start>=targetStart and end<=targetEnd:
                filtMatch = match
                break

        if filtMatch is None:
            errAbort("Multiple matches for this guide, but no single match is within the target sequence? This can happen when your input is not part of the genome, in which case bulk-primer design makes no sense, as far as we can tell. If you think this would make sense for you, contact us at %s." % contactEmail)

        chrom, start, end = filtMatch[:3]
        isUnique = False

        matchList = [filtMatch]

    global batchName
    batchName = batchName.replace(" ", "_")

    chrom, start, end, seq, strand, segType, segName, alnCount, hasXa = matchList[0]
    start = int(start)
    end = int(end)
    return chrom, start, end, strand, segName, isUnique

def printAmpLenAndTm(ampLen, tm):
    " print form fields for amplicon length and TM "
    print ("Maximum amplicon length:")
    dropDownSizes = [
        ("100", "100 bp - for >= 75bp paired reads"),
        ("150", "150 bp - for >= 100bp paired reads "),
        ("200", "200 bp - for >= 125bp paired reads"),
        ("300", "300 bp - for >= 200bp paired reads"),
        ("400", "400 bp - for >= 250bp paired reads"),
        ("500", "500 bp - for >= 300bp paired reads"),
        ("600", "600 bp - for Sanger reads"),
        ("800", "800 bp - for Sanger reads")
    ]

    printDropDown("ampLen", dropDownSizes, ampLen, onChange="""$('#submitPcrForm').click()""")

    print ("&nbsp;&nbsp;&nbsp; Primer Tm:")
    tmList = [
        ("56", "56 deg."),
        ("57", "57 deg."),
        ("58", "58 deg."),
        ("59", "59 deg."),
        ("60", "60 deg."),
        ("61", "61 deg."),
        ("62", "62 deg."),
        ("63", "63 deg."),
        ("64", "64 deg."),
        ("65", "65 deg."),
        ("66", "66 deg."),
        ("67", "67 deg."),
        ("68", "68 deg."),
        ("70", "70 deg."),
    ]
    printDropDown("tm", tmList, tm, onChange="""$('#submitPcrForm').click()""")

def printValidationPcrSection(batchId, genome, pamId, position, params,
        guideStart, guideEnd, primerGuideName, guideSeq):
    " print the PCR section of the primer page "

    # check the input parameters: ampLen, tm
    ampLen = params.get("ampLen", "400")
    if not ampLen.isdigit():
        errAbort("ampLen parameter must be a number")
    ampLen = int(ampLen)

    tm = params.get("tm", "60")
    if not tm.isdigit():
        errAbort("tm parameter must be a number")
    tm = int(tm)

    # Optional, for HDR knock-in experiments
    hdrDist = params.get("hdrDist", None)
    if hdrDist == '':
        hdrDist = None
    if hdrDist is not None:
        try:
            hdrDist = int(hdrDist)
        except (ValueError, TypeError):
            errAbort("hdrDist parameter must be a number")

    print("<h2 id='ontargetPcr'>PCR to amplify the on-target site</h2>")

    otMatches = parseOfftargets(genome, batchId)
    chrom, start, end, strand, gene, isUnique = findOntargetPos(otMatches, pamId, position, absentOk=True)
    if chrom is None:
        print("Not found in genome, cannot design primer")
        return None, None, None

    if not isUnique:
        print("<div class='substep' style='border: 1px black solid; padding:5px; background-color: aliceblue'>")
        print(("<strong>Warning</strong>: Found multiple perfect matches for this guide sequence in the genome. For the PCR, we are using the on-target match in the input sequence %s:%d-%d (gene: %s), but this guide will not be specific. Is this a polyploid organism? Try selecting another guide sequence or email %s to discuss your strategy or modifications to this software.<p>" % (chrom, start+1, end, gene, contactEmail)))
        print("</div>")

    lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq, ampRange, flankSeq, addTags = \
        designPrimer(genome, chrom, start, end, strand, 0, batchId, ampLen, tm, hdrDist)

    primerPosList = []
    if lSeq is not None:
        primerPosList.append( (0, len(lSeq)) )
        primerPosList.append( ( (len(targetSeq)-len(rSeq)), len(targetSeq) ) )

    guideStartOnTarget = None
    guideEndOnTarget = None
    targetHtml = ""
    if lPos!=None:
        guideStartOnTarget = 1000-lPos
        guideEndOnTarget = guideStartOnTarget+GUIDELEN+PAMLEN
        annots = defaultdict(dict)
        annots[(guideStartOnTarget, guideEndOnTarget)]["css"] = {"font-weight":"bold", "background-color" : "yellow"}
        targetHtml = markupSeq(targetSeq, primerPosList, [], annots)

    allPrimersFound = True

    if batchName=="":
        primerPrefix = ""
    else:
        primerPrefix = batchName+"_"

    print("Use these primers to amplify a genomic fragment around the on-target site:<br>")

    print('<table class="primerTable">')
    print('<tr>')
    print("<td>%sOntargetGuideRna%sLeft</td>" % (primerPrefix, primerGuideName))

    if lSeq is not None:
        print("<td>%s</td>" % (lSeq))
    else:
        allPrimersFound = False
        print("<td>Not found</td>")

    print("<td>Tm %s</td>" % (lTm))
    print("</tr><tr>")
    print("<td>%sOntargetGuideRna%sRight</td>" % (primerPrefix, primerGuideName))

    if rSeq is not None:
        print("<td>%s</td>" % (rSeq))
    else:
        allPrimersFound = False
        print("<td>Not found</td>")

    print("<td>Tm %s</td>" % (rTm))
    print('</tr></table><p>')

    print("<h3>Genome fragment with validation primers (underlined) and guide sequence (yellow)</h3>")

    print(("""<form id="ampLenForm" action="%s" method="GET">""" %
        basename(__file__)))

    printAmpLenAndTm(ampLen, tm)

    printHiddenFields(params, {"ampLen":None, "tm":None})

    print("""<input id="submitPcrForm" style="display:none" type="submit" name="submit" value="submit">""")
    print('</form><p>')

    if strand=="-":
        print("Your guide sequence is on the reverse strand relative to the genome sequence, so it is reverse complemented in the sequence below.<p>")

    if not chrom.startswith("ch"):
        chromLong = "chr"+chrom
    else:
        chromLong = chrom

    print('''<div style='word-wrap: break-word; word-break: break-all;'>''')
    if allPrimersFound:
        print("<strong>Genomic sequence %s:%d-%d including primers, genomic forward strand:</strong>" % (chromLong, start, end))
        print("<br><tt>%s</tt><br>" % (targetHtml))
    else:
        print("<strong>Warning: No primers were found at this Tm, please design them manually e.g. with NCBI PrimerBlast.</strong><br>")
    print("<br><tt>%s</tt><br>" % (targetHtml))

    print('''</div>''')
    if rPos is not None:
        print("<strong>Sequence length:</strong> %d<p>" % (rPos-lPos))

    # primer3 may have used some special tags, add them to the description as a string
    p3ArgStr = ""
    if len(addTags)!=0:
        primer3Tags = []
        for key, val in addTags.items():
            primer3Tags.append("%s=%s" % (key, val))
        p3ArgStr = ", ".join(primer3Tags)

    print('<small>Method: Primer3.2 with default settings, target length %s bp, %s</small>' % (ampRange, p3ArgStr))

    return targetSeq, guideStartOnTarget, guideEndOnTarget

def printEnzymeSection(mutEnzymes, targetSeq, guideSeqWPam, guideStartOnTarget, guideEndOnTarget):
    " print the section about restriction enzymes in the target seq "
    print("<h2 id='restrSites'>Restriction Sites for PCR product validation</h2>")

    print("Cas9 induces mutations, usually 3bp 5' of the PAM site.")
    print("If a mutation is induced, then it is very likely that one of the following enzymes no longer cuts your PCR product amplified from the mutant sequence.")
    print("For each restriction enzyme, the guide sequence with the restriction site underlined is shown below.<p>")

    print("<table>")
    print("<tr>")
    print("<th>Enzyme</th><th>Pattern</th><th>Guide with Restriction Site</th><th>Suppliers</th>")
    print("</tr>")
    allSitePos = set()
    patList = []
    for (enzName, pattern, suppliers), posList in mutEnzymes.items():
        print("<tr>")
        patList.append((enzName, pattern))
        print("<td>%s</td><td>%s</td>" % (enzName, pattern))

        print("<td><tt>")
        guideHtmls = []
        for start, end in posList:
            annots = defaultdict(dict)
            annots[(start, end)]["css"] = {"font-weight":"bold"}
            guideHtmls.append( markupSeq(guideSeqWPam.upper(), [], [], annots))
            #print guideSeqWPam
            #if strand=="-":
                #allSitePos.add( (guideEnd-end, guideEnd-start) )
            #else:
            allSitePos.add( (guideStartOnTarget, guideEndOnTarget) )


        print(", ".join(guideHtmls))
        print("</tt></td>")

        supplNames = [rebaseSuppliers.get(x, x) for x in suppliers]
        print("<td>%s</td>" % ", ".join(sorted(supplNames)))
        print("</tr>")
        #print "<br>"
    print("</table>")

    print("<h3>All restriction enzyme sites on the amplicon sequence</h3>")
    print("Restriction sites are shown in yellow, the guide sequence is highlighted in bold. Use this schema to check if the sites are unique enough to give separate bands on a gel:<p>")

    for enzName, pat in patList:
        annots = defaultdict(dict)
        fragLens = []
        lastEnd = 0
        for pos in findPat(targetSeq, pat):
            start = pos
            end = pos+len(pat)
            annots[(start, end)]["css"] = {"background-color":"yellow"}

            fragLens.append( str(start - lastEnd)+"bp" )
            lastEnd = end
        #print markupSeq(targetSeq, [], [(guideStart, guideEnd)], annots)
        #annots[(guideStart, guideEnd)]["css"] = {"font-weight":"bold"}
        fragLens.append( str(len(targetSeq) - lastEnd)+"bp" )

        print(("<div style='margin-top: 6px'>Enzyme: <strong>%s</strong>, Site: %s, Restriction fragment lengths: %s<br>" % (enzName, pat, ", ".join(fragLens))))
        print("<tt>")
        print(markupSeq(targetSeq, [], [(guideStartOnTarget, guideEndOnTarget)], annots))
        print("</tt><br></div>")

def printCloningSection(batchId, primerGuideName, guideSeq, params, pam):
    " print the cloning/expression section of the primer page "
    print("<h2 id='cloning'>Cloning and expression of guide RNA</h2>")

    plasmid = params.get("plasmid", addGenePlasmids[0][0])
    plasmidToName = dict(addGenePlasmids)
    if plasmid not in plasmidToName:
        errAbort("Invalid value for the parameter 'plasmid'")

    primers = makeHelperPrimers(primerGuideName, guideSeq, plasmid, pam)

    # T7 from plasmids
    if not pamIsCpf1(pam):
        print("<h3 id='t7plasmid'>T7 <i>in vitro</i> expression from a plasmid</h3>")
        print('To produce guide RNA by in vitro transcription with T7 RNA polymerase, the guide RNA sequence can be cloned into a variety of plasmids (see <a href="http://addgene.org/crispr/empty-grna-vectors/">AddGene website</a>).<br>')

        print("For the guide sequence %s, the following primers should be ordered for cloning into the BsaI-digested plasmid <a href='https://www.addgene.org/42250/'>DR274</a> generated by the Joung lab.<p>" % guideSeq)

    printPrimerTable(primers["T7"])

    # T7 from primers, in vitro
    print("<h3 id='t7oligo'>T7 <i>in vitro</i> expression from overlapping oligonucleotides</h3>")

    primerType = "spCas9"
    if pamIsSaCas9(pam):
        primerType = "saCas9"

    if not pamIsSaCas9(pam) and not pamIsSpCas9(pam):
        printWarning("The T7 primer information below only applies to %s enzymes. Make sure that you understand the importance of the scaffold sequence in the primers below (both primers depend on the scaffold sequence) and how it is specific to the tracrRNA sequence of your enzyme. Contact us if in doubt, we can quickly adapt the primer for this particular enzyme. " % primerType)

    print("For %s, template for <i>in vitro</i> synthesis of guide RNA with T7 RNA polymerase can be prepared by annealing and primer extension of the following primers:<p>" % primerType)

    printPrimerTable(primers["T7iv"])

    print("""T7 RNA polymerase starts transcription most efficiently if the first two nucleotides to be transcribed are GG. A common recommendation is to add the prefix GG- if our guide does not start with G (5'-N20-(NGG)-3'), to add G- if your guide starts with a single G (5'-GN19-(NGG)-3') and to not add anything if your guide starts with GG already (5'-GGN18-(NGG)-3').<p>""")

    print('One protocol for template preparation from oligonucleotides and in-vitro transcription can be found in <a href="http://www.cell.com/cell-reports/abstract/S2211-1247(13)00312-4">Bassett et al. Cell Rep 2013</a>. We also provide our own <a href="downloads/prot/sgRnaSynthProtocol.pdf">optimized protocol</a> for T7 guide expression.<p>')

    print('''<a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4038517/?report=classic">Gagnon et al. PLoS ONE 2014</a> prefixed guides with GG to ensure high efficiency in vitro transcription by T7 RNA polymerase. It has been shown by other authors that the 5' nucleotides of the guide have little or no role in target specificity and it is therefore generally accepted that prefixing guides with GG should not affect activity.<br>''')

    print('However, in our lab, we found that in vitro transcription with T7 RNA polymerase is efficient enough when the sequence starts with a single G rather than with GG. This took some optimization of the reaction conditions including using large amounts of template DNA and running reactions overnight. <a href="downloads/prot/sgRnaSynthProtocol.pdf">Click here</a> to download our optimized protocol for T7 guide expression.<br>')

    print('Do not use G-prefixing with high-fidelity Cas9 Variants like HF1 and eSpCas9 1.1 when this adds a mismatch in the genome as the efficiency will most likely be very low.')

    # T7 from primers, in vitro
    print("<h3 id='geneArt'>T7 <i>in vitro</i> expression with the GeneArt kit</h3>")
    print("Use these two primers for the Invitrogen GeneArt kit:<p>")

    printPrimerTable(primers["geneArt"])

    # MAMMALIAN CELLS
    print("<h3 id='u6plasmid'>U6 expression from an Addgene plasmid</h3>")
    if "tttt" in guideSeq.lower():
        print("The guide sequence %s contains the motif TTTT, which terminates RNA polymerase. This guide sequence cannot be transcribed in mammalian cells." % guideSeq)
    else:
        print("The guide sequence %s does not contain the motif TTTT, which terminates RNA polymerase, so it can be transcribed in mammalian cells." % guideSeq)

        print("<br>")
        if not pamIsCpf1(pam):
            print(("""<p><form style="margin-bottom: 0px" id="plasmidForm" action="%s#u6plasmid" method="GET">""" %
                basename(__file__)))
            print("Select your Addgene plasmid: ")

            # we need a separate form here (not PCR form), as the target anchor
            # to jump to after a submit is different
            plasmidNames = [(x,y) for x,(y,z) in addGenePlasmids]
            printDropDown("plasmid", plasmidNames, plasmid, onChange="""$('#submitPlasmidForm').click()""")
            printHiddenFields(cgiParams, {"plasmid":None, "submit":None})
            print("""<input id="submitPlasmidForm" style="display:none" type="submit" name="submit" value="submit">""")
            print("""</form></p>""")

            print(("<p>To clone the guide into <i><a href='https://www.addgene.org/%s/'>%s</a></i>, use these primers:" % (plasmid, plasmidToName[plasmid][0])))
        else:
            print("To express guide RNA for Cpf1 in mammalian cells, two plasmids are available. To clone the guide RNA sequence into the plasmids <a href='https://www.addgene.org/78956/'>pU6-As-crRNA</a> or <a href='https://www.addgene.org/78957/'>pU6-Lb-crRNA</a>, where guide RNA expression is driven by a human U6 promoter, the following primers should be used :")

        print("<br>")
        if "mammCellsNote" in primers:
            print("<p><strong>Note:</strong> Efficient transcription from the U6 promoter requires a 5' G. This is not the case for this guide. Several options are possible, you can either add an additional G- prefix to the N20 guide sequence, called  gN20 guides here, or replace the first with a G and create a gN19 guide. For users of HF1 and eSpCas9: G- prefixing with the high-fidelity variants may reduce efficiency, as it introduces a mismatch.</p>")

            print("<strong>Primers for gN20 guides:</strong>")
            printPrimerTable(primers["mammCells"])

            print("<p><strong>Primers for gN19 guides:</strong><br>")
            print("<a href='https://www.nature.com/articles/s41551-019-0505-1'>Kim et al 2020</a>. showed that changing "
            "the first nucleotide to 'G' is slightly more efficient.</p>")

            printPrimerTable(primers["mammCells19"])
        else:
            printPrimerTable(primers["mammCells"])


        _, _, _, enzyme, protoUrl = addGenePlasmidInfo[plasmid]
        print(("<p>The plasmid has to be digested with: <i>%s</i><br>" % enzyme))
        print(("<a href='%s'>Click here</a> to download the cloning protocol for <i>%s</i>" % (protoUrl, plasmidToName[plasmid][0])))

    if pamIsCpf1(pam):
        print("<h3 id='ciona'>Direct PCR for <i>C. intestinalis</i></h3>")
        print ("""A method only used at the moment by <i>Ciona intestinalis</i> (alias <i>Ciona robusta</i>) labs. The DNA construct is assembled during the PCR reaction; expression cassettes are generated with One-Step Overlap PCR (OSO-PCR) <a href="http://www.sciencedirect.com/science/article/pii/S0012160616306455">Gandhi et al., Dev Bio 2016</a> (<a href="http://biorxiv.org/content/early/2017/01/01/041632">preprint</a>) following <a href="downloads/prot/cionaProtocol.pdf">this protocol</a>. The resulting unpurified PCR product can be directly electroporated or injected into Ciona eggs.<br>""")
        if batchName!="":
            primerStart = batchName
        else:
            primerStart = "sg"
        ciPrimers = [
            (batchName+".%s.sgF" % primerGuideName, "g<b>"+guideSeq[1:]+"</b>gtttaagagctatgctggaaacag"),
            (batchName+".%s.U6R" % primerGuideName, "<b>"+revComp(guideSeq[1:])+"</b>catctataccatcggatgccttc")
        ]
        printPrimerTable(ciPrimers)

    print("<h3 id='gibson'>Lentiviral vectors: cloning with Gibson assembly</h3>")
    print ("""Order the following oligonucleotide to clone with Gibson assembly into the vector <a href='https://www.addgene.org/52963/'>pLentiGuide-puro</a>. See the <a href="https://www.nature.com/articles/nprot.2018.005">protocol by Matt Canver</a>.<br>""")
    print ("""To clone with restriction enzymes into this vector, see the section <a href="#u6plasmid">U6 expression from an AddGene plasmid</a> and choose pLentiGuide-puro from the list of AddGene plasmids.<br>""")

    satMutUrl = cgiGetSelfUrl({"satMut":"1"}, onlyParams=["batchId"])
    print(("""If you use lentiviral vectors, you may be interested in our tools for <a href="%s">saturating mutagenesis</a>"""
    """ and for <a href="crispor.py?libDesign=1">gene knockout libraries</a>."""  % satMutUrl))

    lentiPrimers = [
        ("batchOligo%s" % primerGuideName, "<i>GGAAAGGACGAAACACCG</i>"+guideSeq+"<i>GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC</i>"),
    ]

    printPrimerTable(lentiPrimers, seqType="Oligonucleotide")

    print("<h3 id='primerSummary'>Summary of main cloning/expression primers</h3>")
    printPrimerTableAll(primers)

def primerDetailsPage(params):
    """ create primers with primer3 around site identified by pamId in batch
    with batchId. Output primers as html
    """
    # retrieve batch information
    batchId, pamId, pam = params["batchId"], params["pamId"], params["pam"]
    pam = setupPamInfo(pam)

    inSeq, genome, pamSeq, position, extSeq = readBatchParams(batchId)
    seqLen = len(inSeq)
    batchBase = join(batchDir, batchId)

    guideSeq, pamSeq, pamPlusSeq, guideSeqWPam, guideStrand, guideSeqHtml, guideStart, guideEnd \
        = findGuideSeq(inSeq, pam, pamId)

    # search for restriction enzymes that overlap the mutation site
    allEnzymes = readRestrEnzymes()
    mutEnzymes = matchRestrEnz(allEnzymes, guideSeq.upper(), pamSeq.upper(), pamPlusSeq, pam)

    # create a more human readable name of this guide
    guidePos = int(pamId.strip("s+-"))+1
    guideStrand = pamId[-1]
    if guideStrand=="+":
        primerGuideName = str(guidePos)+"fw"
    else:
        primerGuideName = str(guidePos)+"rv"

    # primer helper
    print("""
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

        table.libTable td {
            border-width: 1px;
            table-layout: fixed;
            border-collapse: collapse;
        }
        table.libTable td {
            border-color: #DDDDDD;
        }
    </style>
    """)

    # output the page header
    print('''<div style='width: 80%; margin-left:10%; margin-right:10%; text-align:left;'>''')
    printBackLink()
    print("<h2>")
    if batchName!="":
        print(batchName+":")
    print("Guide sequence: %s</h2>" % (guideSeqHtml))

    print("Contents:<br>")
    print("<ul>")
    print("<li><a href='#cloning'>Cloning or expression of guide RNA</a>")
    print("<ul><li><a href='#t7plasmid'>T7 <i>in vitro</i> expression from a plasmid</a></li></ul>")
    print("<ul><li><a href='#t7oligo'>T7 <i>in vitro</i> expression from overlapping oligonucleotides</a></li></ul>")
    print("<ul><li><a href='#geneArt'>T7 expression with the GeneArt kit</a></li></ul>")
    print("<ul><li><a href='#u6plasmid'>U6 expression from an Addgene plasmid</a></li></ul>")
    print("<ul><li><a href='#ciona'>Direct PCR for <i>C. intestinalis</i></a></li></ul>")
    print("<ul><li><a href='#gibson'>Lentiviral vectors: Cloning with Gibson assembly</a></li></ul>")
    print("<ul><li><a href='#primerSummary'>Summary of main cloning/expression primers</a></li></ul>")
    print("<li><a href='#ontargetPcr'>PCR to amplify the on-target site</a></li>")
    if len(mutEnzymes)!=0:
        print("<li><a href='#restrSites'>Restriction sites for PCR validation</a></li>")
    print("<li><a href='#offtargetPcr'>PCR to amplify off-target sites</a></li>")
    print("<li><a href='#donorGuide'>Guide mutations that minimize on-target activity</a></li>")
    print("<li><a href='#satMut'>Saturating mutagenesis using all guides</a></li>")
    print("</ul>")
    print("<hr>")

    printCloningSection(batchId, primerGuideName, guideSeq, params, pam)
    print("<hr>")

    targetSeq, guideStartOnTarget, guideEndOnTarget = printValidationPcrSection(batchId, genome, pamId, position, params, \
        guideStart, guideEnd, primerGuideName, guideSeq)
    print("<hr>")

    if len(mutEnzymes)!=0 and targetSeq is not None:
        printEnzymeSection(mutEnzymes, targetSeq, guideSeqWPam, guideStartOnTarget, guideEndOnTarget)
        print("<hr>")

    print("<h2 id='offtargetPcr'>PCR to amplify off-target sites</h2>")
    offtUrl = cgiGetSelfUrl({"otPrimers":"1"}, onlyParams=["batchId", "pamId"])
    print(("<p>Primers for all off-targets can be downloaded from the <a href='%s'>Off-target PCR</a> page.</p>" % offtUrl))

    print("<h2 id='donorGuide'>BETA: Guide mutations to minimize on-target activity</h2>")
    doDonorGuide = cgiGetStr(params, "donorGuide", "off")
    if doDonorGuide=="on":
        #print("Guide sequence without PAM is: <tt>%s</tt><p>" % guideSeq)
        #inSeq, genome, pamSeq, position, extSeq = readBatchParams(batchId)
        seq, org, pam, position, guideData = readBatchAndGuides(batchId)
        for guideRow in guideData:
            guideScore, guideCfdScore, effScores, startPos, guideStart, \
            strand, rowPamId, guideSeq, pamSeq, otData, otDesc, \
            last12Desc, mutEnzymes, ontargetDesc, repCount = guideRow
            if rowPamId != pamId:
                continue

            cfdSums = defaultdict(float)
            for otSeq, score, cfdScore, editDist, pos, gene, \
                alnHtml, inLinkage in otData:
                #print("offtarget=%s, cfd=%f, " % (otSeq, cfdScore))
                lastThree = otSeq[:-len(pamSeq)][-3:]
                #print("last 3 bp=%s<br>" % lastThree)
                cfdSums[lastThree]+=cfdScore

            sumList = [(y,x) for (x,y) in list(cfdSums.items())]
            sumList.sort()
            print("Sums of CFD scores for all possible mutations of the last three nucleotides of the guide:<br>")
            for cfdSum, seq in sumList:
                print(("<li>%s: %f" % (seq, cfdSum)))

    else:
        url = cgiGetSelfUrl({"donorGuide":"on"})
        print(("<a href='%s#donorGuide'>Click here to list mutated guides " \
            "sorted by off-targetactivity</a>" % url))

    print("<h2 id='satMut'>Saturating mutagenesis using all guides</h2>")
    satMutUrl = cgiGetSelfUrl({"satMut":"1"}, onlyParams=["batchId"])
    print(("<p>Oligonucleotides of all guides for pooled cloning into a lentiviral vector can be downloaded from the <a href='%s'>Saturating mutagenesis page</a>.</p>" % satMutUrl))

    print("<hr>")

    print('</div>')

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

    guideSeq = "GAGTCCGAGCAGAAGAAGAA"
    for seq, expScore in testRes2.items():
        score = calcHitScore(guideSeq, seq)

def parseArgs():
    " parse command line options into args and options "
    parser = optparse.OptionParser("""usage: %prog [options] org inFile guideOutFile

Command line interface for the Crispor tool.

    org          = genome identifier, like hg19 or ensHumSap
    inFile       = Fasta or BED input file
    guideOutFile = tab-sep file, one row per guide

    If many guides have to scored: Add GGG to them to make them valid
    guides, separate these sequences by N characters and supply as a single
    fasta sequence, a few dozen to ~100 per file. This is faster than providing a multi-FASTA file
    or providing a BED file.

    Examples:
    crispor.py hg38 regions.bed scoreGuides.tsv
    crispor.py mm10 exons.fa scoreGuides.tsv -o offtargets.tsv
    """)

    parser.add_option("-d", "--debug", dest="debug", \
        action="store_true", help="show debug messages, do not delete temp directory")
    parser.add_option("-t", "--test", dest="test", \
        action="store_true", help="run internal tests")
    pamNames = (",".join([x for x,y in pamDesc]))
    parser.add_option("-p", "--pam", dest="pam", \
        action="store", help="PAM-motif to use, default %default. TTTN triggers special Cpf1 behavior: the PAM is assumed to be 5' of the guide. Common PAMs are: " + pamNames, default="NGG")
    parser.add_option("-o", "--offtargets", dest="offtargetFname", \
        action="store", help="write offtarget info to this filename")
    parser.add_option("-m", "--maxOcc", dest="maxOcc", \
        action="store", type="int", help="MAXOCC parameter, guides with more matches are not even processed, default %default", default=MAXOCC)
    parser.add_option("", "--mm", dest="mismatches", \
        action="store", type="int", help="maximum number of mismatches, default %default", default=4)
    parser.add_option("", "--guideLen", dest="guideLen", type="int", \
        action="store", help="Lenght of the guide. Default is: 21 for PAM=NNGRRT/NNNRRT, 23 for Cpf1, 20 otherwise. Note: 19bp guides are less efficient")
    parser.add_option("", "--bowtie", dest="bowtie", \
        action="store_true", help="new: use bowtie as the aligner. Careful: misses off-targets. Do not use.")
    parser.add_option("", "--skipAlign", dest="skipAlign", \
        action="store_true", help="Assume that the input is not in the genome: do not align the input sequence. The on-target will be a random match with 0 mismatches. Switches off efficiency scoring as there is no sequence context.")
    parser.add_option("", "--noEffScores", dest="noEffScores", \
        action="store_true", help="do not calculate any efficiency score")
    parser.add_option("", "--effScores", dest="effScores", \
        action="store", help="calculate only these efficiency scores. Comma-sep list. Possible values, per enzyme: %s" % crisporEffScores.possibleScores)
    parser.add_option("", "--minAltPamScore", dest="minAltPamScore", \
        action="store", type="float", help="minimum MIT off-target score for alternative PAMs, default %default", \
        default=ALTPAMMINSCORE)
    parser.add_option("", "--worker", dest="worker", \
        action="store_true", help="Run as worker process for web server: watches job queue and runs jobs")
    #parser.add_option("", "--user", dest="user", \
        #action="store", help="for the --worker option: switch to this user at program start")
    parser.add_option("", "--clear", dest="clear", \
        action="store_true", help="clear the worker job table and exit")
    parser.add_option("-g", "--genomeDir", dest="genomeDir", \
        action="store", help="directory with genomes, default %default", default=genomesDir)
    parser.add_option("", "--tempDir", dest="tempDir", \
        action="store", help="temp directory for command line. If not specified, remove all temp files on exit")
    parser.add_option("", "--ampliconDir", dest="ampDir", \
        action="store", help="For each guide, write a file with the off-target amplicons to this directory. Filename is <seqId>_<pamId>.txt. See repeat masking note under --satMutDir.")
    parser.add_option("", "--satMutDir", dest="satMutDir", \
        action="store", help="write saturating mutagenesis files to this directory, one per input sequence: ontargetAmplicons.tsv, satMutOligos.tsv, ontargetPrimers.tsv and targetSeqs.txt. Repeats are masked when designing primers, so not all guides may have primers (flagged with 'None').")
    #parser.add_option("", "--ontargetPrimers", dest="ontargetPrimers", \
        #action="store", help="write on-target primers and amplicons, one per guide, to this tab-sep file")
    parser.add_option("", "--ampLen", dest="ampLen", type="int", default=140, \
        action="store", help="for --ampliconDir/--satMutDir: amplicon length, default %default")
    parser.add_option("", "--ampTm", dest="tm", type="int", default=60, \
        action="store", help="for --ampliconDir/--satMutDir: Tm for PCR, default %default")
    #parser.add_option("", "--notInGenome", dest="notInGenome", \
        #action="store_true", help="Input is an articial sequence: do not try to find the input sequence in the genome, assume that the first perfect match of every guide is the on-target")

    (options, args) = parser.parse_args()

    if len(args)==0 and not options.test and not options.worker and not options.clear:
        parser.print_help()
        sys.exit(0)

    if options.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    return args, options

def delBatchDir():
    " called at program exit, for command line mode "
    delTmpDirs() # first remove any subdirs
    if not isdir(batchDir):
        return
    logging.debug("Deleting dir %s" % batchDir)
    fnames = glob.glob(join(batchDir, "*"))
    if len(fnames)>50:
        raise Exception("cowardly refusing to remove many temp files")
    for fname in fnames:
        os.remove(fname)
    os.removedirs(batchDir)

tmpDirsDelExit = []

def delTmpDirs():
    " signal handler at program exit, to remove registered tmp dirs "
    global tmpDirsDelExit
    logging.debug("Removing tmpDirs: %s" % ",".join(tmpDirsDelExit))
    for tmpDir in tmpDirsDelExit:
        if isdir(tmpDir):
            shutil.rmtree(tmpDir)
    tmpDirsDelExit = []

def runQueueWorker():
    " in an infinite loop, take jobs from the job queue in jobs.db and finish them "
    #if userName!=None:
        #uid =  pwd.getpwnam(userName)[2]
        #os.setuid(uid)

    try:
       # Store the Fork PID
       pid = os.fork()

       if pid > 0:
         print('PID: %d' % pid)
         os._exit(0)

    except OSError as error:
      print('Unable to fork. Error: %d (%s)' % (error.errno, error.strerror))
      os._exit(1)

    print(("%s: Worker daemon started. Waiting for jobs." % datetime.ctime(datetime.now())))
    sys.stdout.flush()

    q = JobQueue()
    q.openSqlite()
    print(("Job queue: %s" % JOBQUEUEDB))

    global doAbort
    doAbort = False

    while True:
        if q.waitCount()==0:
            #q.dump()
            time.sleep(1+random.random()/10)
            continue

        if isfile("/tmp/stopCrispor"):
            logging.info("Stop signal received")
            sys.exit(0)

        jobType, batchId, paramStr = q.popJob()
        if jobType=="search":
            print("found job")
            jobError = False
            try:
                seq, org, pamDesc, position, extSeq = readBatchParams(batchId)
                pam = setupPamInfo(pamDesc) # pamDesc includes info on guidelen, etc
                assert('-' not in pam)
                uppSeq = seq.upper()
                startDict, endSet = findAllPams(uppSeq, pam)
                print("searching for offtargets:  ", seq, org, pam, position)
                getOfftargets(uppSeq, org, pamDesc, batchId, startDict, q)
            except:
                exStr = traceback.format_exc()
                print(" - WORKER CRASHED WITH EXCEPTION -")
                print(exStr)
                q.startStep(batchId, "crash", exStr.replace("\n", "///"))
                jobError = True

            if not jobError:
                q.jobDone(batchId)
        elif jobType is None:
            logging.debug("No job")
            time.sleep(1)
        else:
            #raise Exception()
            logging.error("Illegal jobtype: %s - %s. Marking as done." % (jobType, batchId))
            q.jobDone(batchId)
    q.close()

def clearQueue():
    " empty the job queue "
    q = JobQueue()
    q.openSqlite()
    q.clearJobs()
    q.close()
    print(("Worker queue %s, table queue, now empty" % JOBQUEUEDB))

def handleOptions(options):
    " set global vars based on options "
    if options.test:
        runTests()
        import doctest
        doctest.testmod()
        sys.exit(0)

    #if options.ajax:
        #sendStatus(options.ajax)

    if options.clear:
        clearQueue()
        sys.exit(0)

    if options.debug:
        global DEBUG
        DEBUG = True

    if options.noEffScores:
        global doEffScoring
        doEffScoring = False

    # handle the alignment/filtering options
    if options.maxOcc != None:
        global MAXOCC
        MAXOCC = options.maxOcc
        #HIGH_MAXOCC = options.maxOcc

    if options.minAltPamScore!=None:
        global ALTPAMMINSCORE
        ALTPAMMINSCORE = options.minAltPamScore

    if options.mismatches:
        global maxMMs
        maxMMs = options.mismatches

    if options.bowtie:
        global useBowtie
        useBowtie = True

    if options.pam:
        pam = setupPamInfo(options.pam)

    # show all scores in command line mode output files
    global scoreNames
    if options.effScores:
        scoreNames = options.effScores.split(",")
    else:
        scoreNames = allScoreNames
    logging.debug("Active efficiency scores are: %s" % scoreNames)

    # this comes after setupPamInfo, so it overwrites the defaults
    if options.guideLen:
        global GUIDELEN
        GUIDELEN=options.guideLen
        logging.info("Overriding guide length with %d bp as set on command line" % GUIDELEN)

def parseBed(inFname):
    " parse bed  file and return rows "
    rows = []
    for line in open(inFname):
        row = line.rstrip("\n").split()
        name = "noName"
        score = "0"
        strand = "+"

        if len(row)==3:
            chrom, start, end = row
        elif len(row)==4:
            chrom, start, end, name = row
        else:
            chrom, start, end, name, score, strand = row[:6]

        start, end = int(start), int(end)
        assert(strand in ".+-")
        rows.append( (chrom, start, end, name, score, strand) )

    return rows


def mainCommandLine():
    " main entry if called from command line "
    global commandLineMode
    commandLineMode = True

    args, options = parseArgs()

    if options.worker:
        runQueueWorker()
        sys.exit(0)

    handleOptions(options)
    org, inSeqFname, outGuideFname = args

    global skipAlign
    global doEffScoring
    skipAlign = False
    if options.skipAlign:
        skipAlign = True
        doEffScoring = False

    # different genomes directory?
    if options.genomeDir != None:
        global genomesDir
        genomesDir = options.genomeDir

    # get sequence
    if inSeqFname.endswith(".bed"):
        regions = parseBed(inSeqFname)
        seqList = getGenomeSeqsBin(org, regions)
        seqs = {}
        for chrom, start, end, name, score, strand, seq in seqList:
            seqId = "%s:%d-%d:%s" % (chrom, start, end, strand)
            seqs[seqId] = seq
    else:
        seqs = parseFasta(open(inSeqFname))

    # make a directory for the temp files
    # and put it into a global variable, so all functions will use it
    global batchDir

    if options.tempDir:
        batchDir = options.tempDir
        if not isdir(batchDir):
            errAbort("%s is not a directory or does not exist." % batchDir)
    else:
        batchDir = tempfile.mkdtemp(dir=TEMPDIR, prefix="crispor")
        logging.debug("Created directory %s" % batchDir)
        if options.debug:
            logging.info("debug-mode, temporary directory %s will not be deleted" % batchDir)
        else:
            atexit.register(delBatchDir)

    if options.ampDir and not isdir(options.ampDir):
        errAbort("%s does not exist" % options.ampDir)

    # putting multiple sequences into the input file is possible
    # but very inefficient. Rather separate them with a stretch of 10 Ns
    # as explained in the docs
    guideFh = None
    offtargetFh = None
    for seqId, seq in seqs.items():
        seq = seq.upper()
        logging.info(" * running on sequence '%s', guideLen=%d, seqLen=%d" % (seqId, GUIDELEN, len(seq)))
        # get the other parameters and write to a new batch
        seq = seq.upper()
        pamPat = options.pam
        pamPat = setupPamInfo(pamPat)
        batchId = newBatch(seqId, seq, org, pamPat)
        logging.debug("Temporary output directory: %s/%s" % (batchDir, batchId))

        #if position=="?":
            #logging.error("no match found for sequence %s in genome %s" % (inSeqFname, org))

        startDict, endSet = findAllPams(seq, pamPat)

        getOfftargets(seq, org, pamPat, batchId, startDict, ConsQueue())

        batchInfo = readBatchAsDict(batchId)
        position = batchInfo["posStr"]

        otMatches = parseOfftargets(org, batchId)

        # Special batch primer / Crispresso mode
        if options.ampDir:
            pamSeqs = list(flankSeqIter(seq, startDict, len(pamPat), True))
            for pamId, pamStart, guideStart, strand, guideSeq, pamSeq, pamPlusSeq in pamSeqs:
                cPath = join(options.ampDir, "crispresso_%s_%s.txt" % (seqId, pamId))
                logging.info("Writing Crispresso table for seq %s, PAM %s to %s" % (seqId, pamId, cPath))
                ampLen = options.ampLen
                tm = options.tm

                seq, org, pam, position, extSeq = readBatchParams(batchId)
                scoredPrimers, nameToSeq, nameToOtScoreSeq, guideSeqHtml = \
                    designOfftargetPrimers(seq, org, pamPat, position, extSeq, pamId, ampLen, tm, otMatches[pamId])

                cFh = open(cPath, "w")
                for row in makeCrispressoOfftargetRows(scoredPrimers, nameToSeq, nameToOtScoreSeq):
                    cFh.write("\t".join(row))
                    cFh.write("\n")

        # special saturation mutagenesis mode
        if options.satMutDir:
            seq, org, pam, position, guideData = readBatchAndGuides(batchId)
            satMutFname = join(options.satMutDir, seqId+"_satMutOligos.tsv")
            smFh = open(satMutFname, "w")
            logging.info("Writing saturating mutagenesis oligos to %s" % satMutFname)
            writeSatMutFile(0, options.ampLen, options.tm, batchId, None, None, "tsv", smFh)

            primerFname = join(options.satMutDir, seqId+"_ontargetPrimers.tsv")
            pFh = open(primerFname, "w")
            logging.info("Writing primers to %s" % primerFname)
            writeOntargetAmpliconFile("primers", batchId, options.ampLen, options.tm, pFh)

            ampFname = join(options.satMutDir, seqId+"_ontargetAmplicons.tsv")
            ampFh = open(ampFname, "w")
            logging.info("Writing amplicons to %s" % ampFname)
            writeOntargetAmpliconFile("amplicons", batchId, options.ampLen, options.tm, ampFh)

            guideFname = join(options.satMutDir, seqId+"_targetSeqs.tsv")
            gFh = open(guideFname, "w")
            logging.info("Writing guide sequences to %s" % guideFname)
            writeTargetSeqs(guideData, gFh)

        if not doEffScoring:
            effScores = {}
        else:
            effScores = readEffScores(batchId)
        logging.debug("Got efficiency scores: %s" % effScores)

        guideData, guideScores, hasNotFound, pamIdToSeq = \
            mergeGuideInfo(seq, startDict, pamPat, otMatches, position, effScores, org=org)

        # write guide headers
        if guideFh is None:
            guideFh = open(join(batchDir, "guideInfo.tab"), "w")
            guideHeaders, _ = makeGuideHeaders()
            guideHeaders.insert(0, "#seqId")
            guideFh.write("\t".join(guideHeaders)+"\n")

        # write offtarget headers
        if options.offtargetFname and offtargetFh is None:
            offtargetFh = open(join(batchDir, "offtargetInfo.tab"), "w")
            offtargetHeaders.insert(0, "seqId")
            offtargetFh.write("\t".join(offtargetHeaders)+"\n")

        for row in iterGuideRows(guideData, seqId=seqId):
            guideFh.write("\t".join(row))
            guideFh.write("\n")

        if options.offtargetFname:
            for row in iterOfftargetRows(guideData, seqId=seqId, skipRepetitive=False):
                offtargetFh.write("\t".join(row))
                offtargetFh.write("\n")

    guideFh.close()
    shutil.move(guideFh.name, outGuideFname)
    logging.info("guide info written to %s" % outGuideFname)

    if options.offtargetFname:
        offtargetFh.close()
        shutil.move(offtargetFh.name, options.offtargetFname)
        logging.info("off-target info written to %s" % options.offtargetFname)

    if not options.debug and not options.tempDir:
       shutil.rmtree(batchDir)

def sendStatus(batchId):
    " send batch status as json "
    q = JobQueue()
    q.openSqlite()
    status = q.getStatus(batchId)
    q.close()

    if status==None:
        d = {"status":status}
    elif "Traceback" in status:
        d = {"status" : "An error occured. Please send an email to %s and tell us that the failing batchId was %s. We can usually fix this quickly. Thanks!" % (contactEmail, batchId)}
    else:
        d = {"status":status}
    print(json.dumps(d))

def cleanJobs():
    """ look for flag file cleanJobs in current dir. If present, remove jobs.db.
    this is the only way to remove the file, as the jobs.db file is owned by apache
    """
    if isfile("cleanJobs"):
        os.remove(JOBQUEUEDB)
        os.remove("cleanJobs")


def printAssistant(params):
    " "
    print("<h4>What is the aim of your experiment?</h4>")
    print('<form action="crispor.py" name="main" method="get">')
    print('<input type=radio checked name="expType" value="ko">Knock-out of a gene<p>')
    print('<input type=radio name="expType" value="ki">Knock-in of a sequence at a genome position<p>')
    print('<input type=hidden name="assist" value="1">')
    print('<input type=submit name="submit" value="submit">')
    print('</form>')

def mainCgi():
    " main entry if called from apache "
    # XX I need a throttling system
    ip = os.environ.get("REMOTE_ADDR", "noIp")
    if ip=="18.141.51.207" or ip=="80.11.166.200":
        print("Content-type: text/html\n")
        print("Your IP address is hammering crispor and has brought down the server for dozens of other users.")
        print("Please contact me at maxh@ucsc.edu.")
        sys.exit(0)

    # XX IS THE SCRIPT SYMLINKED ? XX
    if os.getcwd()!="/var/www/crispor":
        # only activate stackdumps if running on a development machine
        import cgitb
        cgitb.enable()

    # make all output files world-writable. Useful so we can clean the tempfiles
    os.umask(000)
    cleanJobs()

    # parse incoming parameters and clean them
    params = cgiGetParams()
    batchId = None

    #print "Content-type: text/html\n"
    if "batchId" in params and "download" in params:
        downloadFile(params)
        return

    if "ajaxStatus" in params and "batchId" in params:
        sendStatus(params["batchId"])
        return

    # save seq/org/pam into a cookie, if they were provided
    if "seq" in params and "org" in params and "pam" in params:
        seq, org, pam = params["seq"], params["org"], params["pam"]
        seq, warnMsg = cleanSeq(seq, org)
        saveSeqOrgPamToCookies(seq, org, pam)

    # print headers
    if "downloadCrispresso" not in params:
        print("Content-type: text/html\n")
        # errAbort must know if it has to print this line again
        global contentLineDone
        contentLineDone = True

        title = "CRISPOR"
        if "org" in params:
            title = "CRISPOR: "+params["org"]

        printHeader(batchId, title)

    if "assist" in params:
        printTeforBodyStart()
        printAssistant(params)
        printTeforBodyEnd()
        return

    printBody(params)     # main dispatcher, branches based on the params dictionary

    printTeforBodyEnd()

    # some keywords for google searches
    print("""<div style='display:none'>CRISPR/Cas9 Guide Designer for chordate
    vertebrate ecdysozoans lophotrochozoans protostomes spongi corals plants
    butterflies metazoans genomes fruitflies insects nematodes mammals.
    </div>""")
    print("</body></html>")

def main():
    # detect if running under apache or not
    #print ("Content-type: text/html\n")
    if 'REQUEST_METHOD' in os.environ and sys.argv[-1] != "--worker":
        mainCgi()
    else:
        mainCommandLine()

if __name__=="__main__":
    main()
