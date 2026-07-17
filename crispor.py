#!/data/www/crispor/venv/bin/python3
# if you do not want the hardcoded PATH above, delete this line and the one above to use the default Python3 interpreter
#!/usr/bin/env python3
# I know that this line looks unprofessional to you, but modifying the PATH on a shared Apache webserver is not obvious.

# the tefor crispr tool
# can be run as a CGI or from the command line

# python std library
import subprocess, tempfile, optparse, logging, atexit, glob, shutil, signal, pdb
import http.cookies, time, sys, cgi, re, random, platform, os, pipes, html
import hashlib, base64, string, logging, operator, urllib.request, urllib.parse, urllib.error, time
import traceback, json, pwd, gzip, zlib
import math, difflib

from io import StringIO
from collections import defaultdict, namedtuple
from datetime import datetime
from itertools import product
from os.path import abspath, basename, dirname, isdir, isfile, join, relpath

try:
    # prefer the pip package, it's more up-to-date than the native package
    import pysqlite3 as sqlite3

    SQLITEERROR = pysqlite3.dbapi2.OperationalError
except:
    import sqlite3

    SQLITEERROR = sqlite3.OperationalError

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import (
        OrderedDict,
    )  # python2.6 users: run 'sudo pip install ordereddict'

# for matplotlib, improves "import" performance
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-cache"

# try to load external dependencies
# we're going into great lengths to create a readable error message
needModules = set(["pytabix", "twobitreader", "pandas", "matplotlib", "scipy"])
try:
    import tabix  # if not found, install with 'pip install pytabix'

    needModules.remove("pytabix")
except:
    pass

try:
    import twobitreader  # if not found, install with 'pip install twobitreader'

    needModules.remove("twobitreader")
except:
    pass

try:
    import pandas  # required by doench2016 score. install with 'pip install pandas'

    needModules.remove("pandas")
    import scipy  # required by doench2016 score. install with 'pip install scipy'

    needModules.remove("scipy")
    import matplotlib  # required by doench2016 score. install with 'pip install matplotlib'

    needModules.remove("matplotlib")
    import numpy  # required by doench2016 score. install with 'pip install numpy'

    needModules.remove("numpy")
except:
    pass

# Fix for Azimuth model unpickling with newer numpy versions (>=1.16.1)
# Moved out of the try/except block to ensure it is applied if numpy is loaded
try:
    from numpy.random import mtrand

    if hasattr(mtrand, "__randomstate_ctor"):

        def __randomstate_ctor_patch(bit_generator_name="MT19937", seed=None):
            try:
                from numpy.random import MT19937, RandomState

                if seed is not None:
                    return RandomState(MT19937(seed))
                return RandomState(MT19937())
            except ImportError:
                from numpy.random import RandomState

                if seed is not None:
                    return RandomState(seed)
                return RandomState()

        mtrand.__randomstate_ctor = __randomstate_ctor_patch
except Exception:
    pass

if len(needModules) != 0:
    print("Content-type: text/html\n")
    print(("Python interpreter path: %s<p>" % sys.executable))
    print(("These python modules were not found: %s<p>" % ",".join(needModules)))
    print(
        "To install all requirements in one line, run: sudo pip install biopython numpy scikit-learn==0.16.1 pandas twobitreader<p>"
    )
    sys.exit(0)

# our own eff scoring library
import crisporEffScores
from subserversConf import SUBSERVERS

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
try:
    import MySQLdb

    mysqldbLoaded = True
except:
    mysqldbLoaded = False

# version of crispor
versionStr = "5.4 test"

# Current release note

releaseNote = "April 2026: test version for the new knock-out and knock-in modes."
# contact email
contactEmail = "crispor@tefor.net"

# url to this server
ctBaseUrl = "http://crispor-max.tefor.net/temp/customTracks"

# write debug output to stdout
DEBUG = False
# DEBUG = True

# use bowtie for off-target search?
useBowtie = False

# calculate the efficienc scores?
doEffScoring = True

# system-wide temporary directory
# TEMPDIR = os.environ.get("TMPDIR", "/var/tmp")
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
HTMLPREFIX = ""
# alternative directory on local disk where image/, style/ and js/ are located
HTMLDIR = "/usr/local/apache/htdocs/crispor/"

# directory of crispor.py
baseDir = dirname(__file__)

# filename of this script, usually crispor.py
myName = basename(__file__)

# the segments.bed files use abbreviated genomic region names
segTypeConv = {"ex": "exon", "in": "intron", "ig": "intergenic"}

# directory for processed batches of offtargets ("cache" of bwa results)
batchDir = join(baseDir, "temp")

# sqlite3 db with gzipped old json batch files, to avoid hitting the ext4 inode limits
batchArchive = "/data/crisporJobArchive.db"

# the file where the sqlite job queue is stored
# JOBQUEUEDB = join(TEMPDIR, "crisporJobs.db") # TEMPDIR is mapped away for security reasons under Redhat/Centos for CGIs
JOBQUEUEDB = "/data/www/temp/crisporJobs.db"

# alternatively: connection info for mysql
jobQueueMysqlConn = {"socket": None, "host": None, "user": None, "password": None}

# directory for platform-independent scripts (e.g. Heng Li's perl SAM parser)
scriptDir = join(baseDir, "scripts")

# directory for helper binaries (e.g. BWA)
# system() is one of 'Linux', 'Darwin', 'Windows', machine() is one of 'x86_64', 'arm64', 'aarch64'
os_name = platform.system()  # 'Linux', 'Darwin', 'Windows'
arch = platform.machine()  # 'x86_64', 'arm64', 'aarch64'
binDir = abspath(join(baseDir, "bin", platform.system() + "-" + platform.machine()))

# directory for genomes
genomesDir = join(baseDir, "genomes")

DEFAULTORG = "hg19"
DEFAULTSEQ = "cttcctttgtccccaatctgggcgcgcgccggcgccccctggcggcctaaggactcggcgcgccggaagtggccagggcgggggcgacctcggctcacagcgcgcccggctattctcgcagctcaccatgGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCC"

DEFAULTKISEQ = "cctgcgaactagtcggtggctcgggcgccggcggggagctgctcggcggcggacagtgtaATGTTGGGTGGGAGTGCGGGACGCCTCAAAATGTCTTCCAGTGGCACCCTCAGCAACTA"
DEFAULTINSERTSEQ = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTA"
DEFAULTINSERT = (
    DEFAULTKISEQ[0:63].lower() + DEFAULTINSERTSEQ + DEFAULTKISEQ[63:].lower()
)
DEFAULTDEL = DEFAULTKISEQ[0:45].lower() + DEFAULTKISEQ[60:].lower()
DEFAULTSUBST = DEFAULTKISEQ[0:60].lower() + "G" + DEFAULTKISEQ[61:].lower()
DEFAULTREPL = DEFAULTKISEQ[0:60].lower() + "TAG" + DEFAULTKISEQ[63:].lower()


# used if hg19 is not available
ALTORG = "sacCer3"
ALTSEQ = "ATTCTACTTTTCAACAATAATACATAAACatattggcttgtggtagCAACACTATCATGGTATCACTAACGTAAAAGTTCCTCAATATTGCAATTTGCTTGAACGGATGCTATTTCAGAATATTTCGTACTTACACAGGCCATACATTAGAATAATATGTCACATCACTGTCGTAACACTCT"

pamDesc = [
    ("NGG", "20bp-NGG - Sp Cas9, SpCas9-HF1, eSpCas9 1.1"),
    ("NNG", "20bp-NNG - Cas9 S. canis"),
    ("NGN", "20bp-NGN - SpG"),
    ("NNGT", "20bp-NNGT - Cas9 S. canis - high efficiency PAM, recommended"),
    ("NAA", "20bp-NAA - iSpyMacCas9"),
    (
        "TTN",
        "TTN-23bp - hfCas12Max - as recommended by Synthego",
    ),  # Casey Jowdy by email
    ("TNN", "TNN-23bp - hfCas12Max, broader PAM, as recommended by Synthego"),  #
    ("NGG-22", "NGG-22bp - eSpOT-ON (ePsCas9), as recommended by Synthego"),
    ("NNGRRT", "21bp-NNG(A/G)(A/G)T - Cas9 S. Aureus"),
    ("NNGRRT-20", "20bp-NNG(A/G)(A/G)T - Cas9 S. Aureus with 20bp-guides"),
    ("NGK", "20bp-NG(G/T) - xCas9, recommended PAM, see notes"),
    ("NRTH", "20bp-N(A/G)T(A/C/T), Sp Cas9 engineered PAM"),
    ("NRRH", "20bp-N(A/G)(A/G)(A/C/T), Sp Cas9 engineered PAM"),
    ("NRCH", "20bp-N(A/G)C(A/C/T), Sp Cas9 engineered PAM"),
    # ('NGN','20bp-NGN or GA(A/T) - xCas9 (low efficiency, not recommended)'),
    # ('NGG-BE1', '20bp-NGG - BaseEditor1, modifies C->T'),
    ("NNNRRT", "21bp-NNN(A/G)(A/G)T - KKH SaCas9"),
    ("NNNRRT-20", "20bp-NNN(A/G)(A/G)T - KKH SaCas9 with 20bp-guides"),
    ("NGA", "20bp-NGA - Cas9 S. Pyogenes mutant VQR"),
    ("NNNNCC", "24bp-NNNNCC - Nme2Cas9"),
    ("NGCG", "20bp-NGCG - Cas9 S. Pyogenes mutant VRER"),
    ("NNAGAA", "20bp-NNAGAA - Cas9 S. Thermophilus"),
    ("NGGNG", "20bp-NGGNG - Cas9 S. Thermophilus"),
    ("NNNNGMTT", "20bp-NNNNG(A/C)TT - Cas9 N. Meningitidis"),
    ("NNNNACA", "20bp-NNNNACA - Cas9 Campylobacter jejuni, original PAM"),
    ("NNNNRYAC", "22bp-NNNNRYAC - Cas9 Campylobacter jejuni, revised PAM"),
    ("NNNVRYAC", "22bp-NNNVRYAC - Cas9 Campylobacter jejuni, opt. efficiency"),
    ("TTCN", "TTCN-20bp - CasX"),
    ("TTTV", "TTT(A/C/G)-23bp - Cas12a (Cpf1)  - recommended, 23bp guides"),
    ("TTTV-21", "TTT(A/C/G)-21bp - Cas12a (Cpf1) - 21bp guides recommended by IDT"),
    ("TTTN", "TTTN-23bp - Cas12a (Cpf1) - low efficiency"),
    ("ATTN", "ATTN-23bp - BhCas12b v4"),
    ("NGTN", "NGTN-23bp - ShCAST/AcCAST, Strecker et al, Science 2019"),
    ("TYCV", "T(C/T)C(A/C/G)-23bp - TYCV As-Cpf1 K607R"),
    ("TATV", "TAT(A/C/G)-23bp - TATV As-Cpf1 K548V"),
    ("TTTA", "TTTA-23bp - TTTA LbCpf1"),
    ("TCTA", "TCTA-23bp - TCTA LbCpf1"),
    ("TCCA", "TCCA-23bp - TCCA LbCpf1"),
    ("CCCA", "CCCA-23bp - CCCA LbCpf1"),
    ("GGTT", "GGTT-23bp - CCCA LbCpf1"),
    ("YTTV", "YTTV-20bp - MAD7 Nuclease, Lui, Schiel, Maksimova et al, CRISPR J 2020"),
    (
        "TTYN",
        "TTYN- or VTTV- or TRTV-23bp - enCas12a E174R/S542R/K548R - Kleinstiver et al Nat Biot 2019",
    ),
    ("NNNNCNAA", "20bp-NNNNCNAA - Thermo Cas9 - Walker et al, Metab Eng Comm 2020"),
    (
        "NNN",
        "20bp-NNN - SpRY, Walton et al Science 2020",
    ),  # https://science.sciencemag.org/content/368/6488/290.abstract
    ("NRN", "20bp-NRN - SpRY (high efficiency PAM)"),
    ("NYN", "20bp-NYN - SpRY (low efficiency PAM)"),
    # ('VTTV','(A/C)TT(A/C)-23bp - enCas12a S542R - Kleinstiver et al Nat Biot 2019'),
    # ('TRTV','T(A/G)T(A/C)-23bp - enCas12a K548R - Kleinstiver et al Nat Biot 2019'),
]

# list of base editors
beDesc = [("NGG-BE1", "20bp-NGG - BaseEditor1, modifies C->T")]

# Ideally, pass pams in the batch parameters (to allow the user to select pams) ?
multiPamDict = {
    "20bp-NGG": (["NGG"], "20bp-NGG - Sp Cas9, SpCas9-HF1, eSpCas9 1.1"),
    "commercial": (
        ["NGG", "TTN", "TTTV-21", "NNGRRT"],
        "PAMs from commercially available nucleases including only non-engineered PAMs",
    ),
    "pamless": (
        [
            "TNN",
        ],
        "PAMs from commercially available nucleases including lower specificity engineered PAMs",
    ),
    "plasmid": (
        [],
        "PAMs from nucleases with expression plasmid available from Addgene (to be added)",
    ),
}


# dict of mutations possible with base editos
# tuples are (fromNucl, toNucl)
# [(on + strand), (on - strand)]: enzyme
possibleEdits = {"CBE": [("C", "T"), ("G", "A")],
                 "ABE": [("A", "G"), ("T", "C")],
                 "CGBE": [("C", "G"), ("G", "C")]}

# List of Base editor, with their respective models and editing windows
# work in progress
allBeModels = {
    "ABE": [
      {"tool": "DeepBe",     "model": "DeepNG-BE_8e",  "win": (2, 11)},
      {"tool": "DeepBe",     "model": "DeepNG-BE_17m", "win": (2, 11)},
      {"tool": "ForecastBe", "model": "ABE",           "win": (4, 9)}
    ],
    "CBE":  [
      # results from Kim et al show a window of 2-13 for SsAPOBEC3B, but the model seems to predict for positions 2-9 only ?
      # window is hardcoded to be seq[6:13] for Ss and seq[7:11] for YE1
      {"tool": "DeepBe",     "model": "DeepNG-BE_Ss",  "win": (2, 9)},
      {"tool": "DeepBe",     "model": "DeepNG-BE_YE1", "win": (3, 8)},
      {"tool": "ForecastBe", "model": "CBE",           "win": (3, 10)}  # 3-10
      ],
    "CGBE": [
      {"tool": "DeepBe",     "model": "DeepNG-BE_mini",  "win": (3, 8)},
      {"tool": "DeepBe",     "model": "DeepNG-BE_CGBE1", "win": (3, 8)},
      {"tool": "DeepBe",     "model": "DeepNG-BE_Bi", "win": (3, 8)}
      ]
    }

pamVariantModels = {
        "NGG": 1,  # PAM_variant_SpCas9_model.h5
        # 1: "VRQR",  # PAM_variant_VRQR_model.h5
        # "NGN": 2,  # PAM_variant_NG_model.h5
        "NRRH": 3,  # PAM_variant_NRRH_model.h5
        "NRTH": 4,  # PAM_variant_NRTH_model.h5
        "NRCH": 5,  # PAM_variant_NRCH_model.h5
        "NGN": 6,  # PAM_variant_SpG_model.h5
        "NRN": 7,  # PAM_variant_SpRY_model.h5
        "NNGT": 8,  # PAM_variant_Sc++_model.h5
        }

# mapping of Base editing model to their respective enzyme
modelToEnzyme = {
        "ABE":
        ("SpCas9-ABE", "Generalist Adenine Base Editor model trained on data from. See <a href='https://doi.org/10.1093/nar/gkac161' target='blank'>Pallaseni et al. 2022</a>"),
        "DeepNG-BE_8e":
        ("SpCas9-ABE8e(V106W)", "Model trained on data from SpCas9 fused with ABE8e(V106W). See <a href='https://doi.org/10.1038/s41587-023-01792-x' target='blank'>Kim et et al. 2023</a>"),
        "DeepNG-BE_17m":
        ("SpCas9-ABE8.17-m+V106W", "Model trained on data from SpCas9 fused with ABE8e(V106W). See <a href='https://doi.org/10.1038/s41587-023-01792-x' target='blank'>Kim et et al. 2023</a>"),
        "CBE":
        ("SpCas9-CBE", "Generalist Cytosine Base Editor model trained on data from. See <a href='https://doi.org/10.1093/nar/gkac161' target='blank'>Pallaseni et al. 2022</a>"),
        "DeepNG-BE_Ss":
        ("SpCas9-SsAPOBEC3B", "Model trained on data from SpCas9 fused with ABE8e(V106W). See <a href='https://doi.org/10.1038/s41587-023-01792-x' target='blank'>Kim et et al. 2023</a>"),
        "DeepNG-BE_YE1":
        ("SpCas9-YE1-BE4max", "Model trained on data from SpCas9 fused with ABE8e(V106W). See <a href='https://doi.org/10.1038/s41587-023-01792-x' target='blank'>Kim et et al. 2023</a>"),
        "DeepNG-BE_mini":
        ("SpCas9-miniCGBE1", "Model trained on data from SpCas9 fused with ABE8e(V106W). See <a href='https://doi.org/10.1038/s41587-023-01792-x' target='blank'>Kim et et al. 2023</a>"),
        "DeepNG-BE_CGBE1":
        ("SpCas9-CGBE1", "Model trained on data from SpCas9 fused with ABE8e(V106W). See <a href='https://doi.org/10.1038/s41587-023-01792-x' target='blank'>Kim et et al. 2023</a>"),
        "DeepNG-BE_Bi":
        ("SpCas9-APOBEC-nCas9-Ung", "Model trained on data from SpCas9 fused with ABE8e(V106W). See <a href='https://doi.org/10.1038/s41587-023-01792-x' target='blank'>Kim et et al. 2023</a>")
        }

DEFAULTPAM = "NGG"

# the default base editor modification window
DEFAULTBEWIN = "4-9"

# for some PAMs, there are alternative main PAMs. These are also shown on the main sequence panel
multiPams = {
    # "NGN" : ["GAW"],
    "TTYN": ["VTTV", "TRTV"]
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
    "NGG": ["NAG", "NGA"],
    # "NGN" : ["GAW"],
    "NGK": ["GAW"],
    "NGA": ["NGG"],
    "NNGRRT": ["NNGRRN"],
    "TTTV": ["TTTN"],
    "ATTN": ["TTTN", "GTTN"],
    "TTYN": ["VTTV", "TRTV"],
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
maxMMs = 4

# maximum number of occurences in the genome to get flagged as repeats.
# This is used in bwa samse, when converting the same file
# and for warnings in the table output.
MAXOCC = 60000

# the BWA queue size is 2M by default. We derive the queue size from MAXOCC
MFAC = 2000000 / MAXOCC

# the length of the guide sequence, set by setupPamInfo
GUIDELEN = None
# length of the PAM sequence
PAMLEN = None

# the name of the base editor, if any. This is the flag to activate
# baseEditor mode in the UI
baseEditor = None

# input sequences are extended by X basepairs so we can calculate the efficiency scores
# and can better design primers
FLANKLEN = 100

# the name of the currently processed batch, assigned only once
# in readBatchParams and only for json-type batches
batchName = ""

# are we doing a Cpf1 run?
# this variable changes almost all processing and
# has to be set on program start, as soon as we know
# the PAM we're running on
pamIsFirst = None
saCas9Mode = False

# Highly-sensitive mode (not for CLI mode):
# MAXOCC is increased in processSubmission() and in the html UI if only one
# guide seq is run
# Also, the number of allowed mismatches is increased to 5 instead of 4
# HIGH_MAXOCC=600000
# HIGH_maxMMs=5

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
koCas9ScoreNames = ["rs3", "EVA", "crisprScan"]
cas9ScoreNames = ["fusi", "crisprScan", "rs3", "EVA"]
allScoreNames = [
    "fusi",
    "chariRank",
    "ssc",
    "wuCrispr",
    "doench",
    "wang",
    "crisprScan",
    "ccTop",
    "rs3",
    "EVA",
]

mutScoreNames = []
spCas9MutScoreNames = ["oof", "lindel"]  # lindel is only added for spCas9
otherMutScoreNames = ["oof"]  # lindel is only added for spCas9

cpf1ScoreNames = ["seqDeepCpf1"]

saCas9ScoreNames = ["najm"]

# to make the CFD more comparable to the MIT score, Nicholas Parkinson suggests to multiply it with 100.
# can be switched on with the URL argument fixCfd=1
doCfdFix = False

# how many digits shall we show for each score? default is 0
scoreDigits = {
    "ssc": 1,
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
    (
        "61591",
        (
            "pX601-AAV-CMV::NLS-SaCas9-NLS-3xHA-bGHpA;U6::BsaI-sgRNA (Zhang lab)",
            "pX601",
        ),
    ),
    ("61592", ("pX600-AAV-CMV::NLS-SaCas9-NLS-3xHA-bGHpA (Zhang lab)", "pX600")),
    (
        "61593",
        (
            "pX602-AAV-TBG::NLS-SaCas9-NLS-HA-OLLAS-bGHpA;U6::BsaI-sgRNA (Zhang lab)",
            "pX602",
        ),
    ),
    ("65779", ("VVT1 (Joung lab)", "VVT1")),
]

# list of AddGene primer 5' and 3' extensions, one for each AddGene plasmid
# format: prefixFw, prefixRw, u6-G-suffix, restriction enzyme, link to protocol
addGenePlasmidInfo = {
    "43860": (
        "ACACC",
        "AAAAC",
        "G",
        "BsmBI",
        "https://www.addgene.org/static/data/plasmids/43/43860/43860-attachment_T35tt6ebKxov.pdf",
    ),
    "49330": (
        "TTC",
        "AAC",
        "",
        "Bsp QI",
        "http://bio.biologists.org/content/3/1/42#sec-9",
    ),
    "42230": (
        "CACC",
        "AAAC",
        "",
        "Bbs1",
        "https://www.addgene.org/static/data/plasmids/52/52961/52961-attachment_B3xTwla0bkYD.pdf",
    ),
    "52961": (
        "CACC",
        "AAAC",
        "",
        "BsmBI",
        "https://www.addgene.org/static/data/plasmids/52/52961/52961-attachment_B3xTwla0bkYD.pdf",
    ),
    "61591": (
        "CACC",
        "AAAC",
        "",
        "BsaI",
        "https://www.addgene.org/static/data/plasmids/61/61591/61591-attachment_it03kn5x5O6E.pdf",
    ),
    "61592": (
        "CACC",
        "AAAC",
        "",
        "BsaI",
        "https://www.addgene.org/static/data/plasmids/61/61592/61592-attachment_iAbvIKnbqNRO.pdf",
    ),
    "61593": (
        "CACC",
        "AAAC",
        "",
        "BsaI",
        "https://www.addgene.org/static/data/plasmids/61/61592/61592-attachment_iAbvIKnbqNRO.pdf",
    ),
    "65779": (
        "CACC",
        "AAAC",
        "",
        "BsmBI (aka Esp3l)",
        "https://www.addgene.org/static/data/plasmids/65/65779/65779-attachment_G8oNyvV6pA78.pdf",
    ),
    "52963": (
        "CACC",
        "AAAC",
        "",
        "BsmBI (aka Esp3l)",
        "https://www.addgene.org/static/data/plasmids/52/52963/52963-attachment_IPB7ZL_hJcbm.pdf",
    ),
}

# the barcodes for subpool tagging for oligo pool tables
satMutBarcodes = [
    (0, "No Subpool barcode"),
    (1, "Subpool 1: CGGGTTCCGT/GCTTAGAATAGAA"),
    (2, "Subpool 2: GTTTATCGGGC/ACTTACTGTACC"),
    (3, "Subpool 3: ACCGATGTTGAC/CTCGTAATAGC"),
    (4, "Subpool 4: GAGGTCTTTCATGC/CACAACATA"),
    (5, "Subpool 5: TATCCCGTGAAGCT/TTCGGTTAA"),
    (6, "Subpool 6: TAGTAGTTCAGACGC/ATGTACCC"),
    (7, "Subpool 7: GGATGCATGATCTAG/CATCAAGC"),
    (8, "Subpool 8: ATGAGGACGAATCT/CACCTAAAG"),
    (9, "Subpool 9: GGTAGGCACG/TAAACTTAGAACC"),
    (10, "Subpool 10: AGTCATGATTCAG/GTTGCAAGTCTAG"),
]

# Restriction enzyme supplier codes
rebaseSuppliers = {
    "B": "Life Technologies",
    "C": "Minotech",
    "E": "Agilent",
    "I": "SibEnzyme",
    "J": "Nippon Gene",
    "K": "Takara",
    "M": "Roche",
    "N": "NEB",
    "O": "Toyobo",
    "Q": "Molecular Biology Resources",
    "R": "Promega",
    "S": "Sigma",
    "V": "Vivantis",
    "X": "EURx",
    "Y": "SinaClon BioScience",
}

# tag, linker and marker sequences

taggingSeqs = {
    "eGFP": "GTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG",
    "Streptavidin": "TGGAGCCACCCGCAGTTCGAAAAA",
    "(GGGGS)x2": "GGTGGTGGTGGTTCTGGTGGTGGTGGTTCT",
    "GGGGS": "GGAGGCGGCGGCAGC",
    "GSGGG": "GGATCCGGCGGAGGA",
    "GSG-2A": "GGATCCGGA",
    "XTEN": "AGCGGCAGCGAGACTCCCGGGACCTCAGAGTCCGCCACACCCGAAAGT",
    "Blast": "ATGGCCAAGCCTTTGTCTCAAGAAGAATCCACCCTCATTGAAAGAGCAACGGCTACAATCAACAGCATCCCCATCTCTGAAGACTACAGCGTCGCCAGCGCAGCTCTCTCTAGCGACGGCCGCATCTTCACTGGTGTCAATGTATATCATTTTACTGGGGGACCTTGTGCAGAACTCGTGGTGCTGGGCACTGCTGCTGCTGCGGCAGCTGGCAACCTGACTTGTATCGTCGCGATCGGAAATGAGAACAGGGGCATCTTGAGCCCCTGCGGACGGTGCCGACAGGTGCTTCTCGATCTGCATCCTGGGATCAAAGCCATAGTGAAGGACAGTGATGGACAGCCGACGGCAGTTGGGATTCGTGAATTGCTGCCCTCTGGTTATGTGTGGGAGGGC",
    "Puro": "ATGACCGAGTACAAGCCCACGGTGCGCCTCGCCACCCGCGACGACGTCCCCAGGGCCGTACGCACCCTCGCCGCCGCGTTCGCCGACTACCCCGCCACGCGCCACACCGTCGATCCGGACCGCCACATCGAGCGGGTCACCGAGCTGCAAGAACTCTTCCTCACGCGCGTCGGGCTCGACATCGGCAAGGTGTGGGTCGCGGACGACGGCGCCGCGGTGGCGGTCTGGACCACGCCGGAGAGCGTCGAAGCGGGGGCGGTGTTCGCCGAGATCGGCCCGCGCATGGCCGAGTTGAGCGGTTCCCGGCTGGCCGCGCAGCAACAGATGGAAGGCCTCCTGGCGCCGCACCGGCCCAAGGAGCCCGCGTGGTTCCTGGCCACCGTCGGCGTCTCGCCCGACCACCAGGGCAAGGGTCTGGGCAGCGCCGTCGTGCTCCCCGGAGTGGAGGCGGCCGAGCGCGCCGGGGTGCCCGCCTTCCTGGAGACCTCCGCGCCCCGCAACCTCCCCTTCTACGAGCGGCTCGGCTTCACCGTCACCGCCGACGTCGAGGTGCCCGAAGGACCGCGCACCTGGTGCATGACCCGCAAGCCCGGTGCC",
    "Zeo": "ATGGCCAAGTTGACCAGTGCCGTTCCGGTGCTCACCGCGCGCGACGTCGCCGGAGCGGTCGAGTTCTGGACCGACCGGCTCGGGTTCAGCCGGGACTTCGTGGAGGACGACTTCGCCGGTGTGGTCCGGGACGACGTGACCCTGTTCATCAGCGCGGTCCAGGACCAGGTGGTGCCGGACAACACCCTGGCCTGGGTGTGGGTGCGCGGCCTGGACGAGCTGTACGCCGAGTGGTCGGAGGTCGTGTCCACGAACTTCCGGGACGCCTCCGGGCCGGCCATGACCGAGATCGGCGAGCAGCCGTGGGGGCGGGAGTTCGCCCTGCGCGACCCGGCCGGCAACTGCGTGCACTTCGTGGCCGAGGAGCAGGAC",
    "moxGFP": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTGAACGGCCACAAGTTCTCCGTGCGGGGCGAGGGCGAGGGCGATGCCACCAACGGCAAGCTGACCCTGAAGTTCATCAGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGAGCTTCTCCCGCTACCCCGACCACATGAAGCGCCACGACTTCTTCAAGAGCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTCCTTCAAGGACGACGGCACCTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTTCAACTCCCACAACGTCTATATCACCGCCGACAAGCAGAAGAACGGCATCAAGGCCAACTTCAAGATCCGCCATAACGTGGAGGACGGCTCCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGTCCACCCAGTCCAAGCTGTCCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCACGGCATGGACGAGCTGTACAAG",
    "2A ribosomal skipping peptide": "GCTACTAACTTCAGCCTGCTGAAGCAGGCCGGAGACGTGGAGGAGAACCCTGGACCT",
    "EF1α promoter": "GGGCAGAGCGCACATCGCCCACAGTCCCCGAGAAGTTGGGGGGAGGGGTCGGCAATTGAACCGGTGCCTAGAGAAGGTGGCGCGGGGTAAACTGGGAAAGTGATGTCGTGTACTGGCTCCGCCTTTTTCCCGAGGGTGGGGGAGAACCGTATATAAGTGCAGTAGTCGCCGTGAACGTTCTTTTTCGCAACGGGTTTGCCGCCAGAACACAG",
    "EF1α intron": "GTAAGTGCCGTGTGTGGTTCCCGCGGGCCTGGCCTCTTTACGGGTTATGGCCCTTGCGTGCCTTGAATTACTTCCACCTGGCTGCAGTACGTGATTCTTGATCCCGAGCTTCGGGTTGGAAGTGGGTGGGAGAGTTCGAGGCCTTGCGCTTAAGGAGCCCCTTCGCCTCGTGCTTGAGTTGAGGCCTGGCCTGGGCGCTGGGGCCGCCGCGTGCGAATCTGGTGGCACCTTCGCGCCTGTCTCGCTGCTTTCGATAAGTCTCTAGCCATTTAAAATTTTTGATGACCTGCTGCGACGCTTTTTTTCTGGCAAGATAGTCTTGTAAATGCGGGCCAAGATCTGCACACTGGTATTTCGGTTTTTGGGGCCGCGGGCGGCGACGGGGCCCGTGCGTCCCAGCGCACATGTTCGGCGAGGCGGGGCCTGCGAGCGCGGCCACCGAGAATCGGACGGGGGTAGTCTCAAGCTGGCCGGCCTGCTCTGGTGCCTGGCCTCGCGCCGCCGTGTATCGCCCCGCCCTGGGCGGCAAGGCTGGCCCGGTCGGCACCAGTTGCGTGAGCGGAAAGATGGCCGCTTCCCGGCCCTGCTGCAGGGAGCTCAAAATGGAGGACGCGGCGCTCGGGAGAGCGGGCGGGTGAGTCACCCACACAAAGGAAAAGGGCCTTTCCGTCCTCAGCCGTCGCTTCATGTGACTCCACGGAGTACCGGGCGCCGTCCAGGCACCTCGATTAGTTCTCGAGCTTTTGGAGTACGTCGTCTTTAGGTTGGGGGGAGGGGTTTTATGCGATGGAGTTTCCCCACACTGAGTGGGTGGAGACTGAAGTTAGGCCAGCTTGGCACTTGATGTAATTCTCCTTGGAATTTGCCCTTTTTGAGTTTGGATCTTGGTTCATTCTCAAGCCTCAGACAGTGGTTCAAAGTTTTTTTCTTCCATTTCAG",
    "Cre": "ATGAGCAATTTACTGACCGTACACCAAAATTTGCCTGCATTACCGGTCGATGCAACGAGTGATGAGGTTCGCAAGAACCTGATGGACATGTTCAGGGATCGCCAGGCGTTTTCTGAGCATACCTGGAAAATGCTTCTGTCCGTTTGCCGGTCGTGGGCGGCATGGTGCAAGTTGAATAACCGGAAATGGTTTCCCGCAGAACCTGAAGATGTTCGCGATTATCTTCTATATCTTCAGGCGCGCGGTCTGGCAGTAAAAACTATCCAGCAACATTTGGGCCAGCTAAACATGCTTCATCGTCGGTCCGGGCTGCCACGACCAAGTGACAGCAATGCTGTTTCACTGGTTATGCGGCGGATCCGAAAAGAAAACGTTGATGCCGGTGAACGTGCAAAACAGGCTCTAGCGTTCGAACGCACTGATTTCGACCAGGTTCGTTCACTCATGGAAAATAGCGATCGCTGCCAGGATATACGTAATCTGGCATTTCTGGGGATTGCTTATAACACCCTGTTACGTATAGCCGAAATTGCCAGGATCAGGGTTAAAGATATCTCACGTACTGACGGTGGGAGAATGTTAATCCATATTGGCAGAACGAAAACGCTGGTTAGCACCGCAGGTGTAGAGAAGGCACTTAGCCTGGGGGTAACTAAACTGGTCGAGCGATGGATTTCcGTCTCTGGTGTAGCTGATGATCCGAATAACTACCTGTTTTGCCGGGTCAGAAAAAATGGTGTTGCCGCGCCATCTGCCACCAGCCAGCTATCAACTCGCGCCCTGGAAGGGATTTTTGAAGCAACTCATCGATTGATTTACGGCGCTAAGGATGACTCTGGTCAGAGATACCTGGCCTGGTCTGGACACAGTGCCCGTGTCGGAGCCGCGCGAGATATGGCCCGCGCTGGAGTTTCAATACCGGAGATCATGCAAGCTGGTGGCTGGACCAATGTAAATATTGTCATGAACTATATCCGTAACCTGGATAGTGAAACAGGGGCAATGGTGCGCCTGCTGGAAGATGGCGAC",
    "mStrayGold": "ATGGTGTCTACAGGCGAGGAACTGTTTACCGGCGTGGTGCCCTTCAAGTTCCAGCTGAAGGGCACCATCAACGGCAAGAGCTTCACCGTGGAAGGCGAGGGCGAGGGCAATAGCCACGAGGGCAGCCACAAAGGCAAGTATGTGTGCACCAGCGGCAAACTGCCAATGTCTTGGGCCGCCCTGGGAACTAGCTTCGGCTATGGCATGAAATACTATACCAAGTACCCCAGCGGCCTGAAAAACTGGTTCCACGAGGTGATGCCTGAGGGCTTCACCTACGACAGACACATCCAGTACAAGGGCGACGGCAGCATCCACGCCAAGCACCAGCACTTCATGAAGAACGGCACCTACCACAACATCGTGGAGTTCACCGGCCAGGACTTCAAGGAGAACAGCCCCGTGCTGACCGGCGACATGGACGTGAGCCTGCCCAACGAGGTGCAGCACATCCCCATTGATGACGGCGTGGAGTGCACAGTGACCCTGCAGTACCCTCTGCTGAGCGACGAAAGCAAGTGCGTGGAAGCCTATCAGAACACCATCATCAAGCCCCTGCACAATCAGCCAGCCCCCGATGTGCCATTTCACTGGATCAGAAAGCAGTACACCCAGAGCAAGGACGACACCGAGGAGAGAGACCACATCATCCAGAGCGAGACCCTGGAGGCCCACCTG",
    "mNeon": "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCTCTCTCCCAGCGACACATGAGTTACACATCTTTGGCTCCATCAACGGTGTGGACTTTGACATGGTGGGTCAGGGCACCGGCAATCCAAATGATGGTTATGAGGAGTTAAACCTGAAGTCCACCAAGGGTGACCTCCAGTTCTCCCCCTGGATTCTGGTCCCTCATATCGGGTATGGCTTCCATCAGTACCTGCCCTACCCTGACGGGATGTCGCCTTTCCAGGCCGCCATGGTAGATGGCTCCGGATACCAAGTCCATCGCACAATGCAGTTTGAAGATGGTGCCTCCCTTACTGTTAACTACCGCTACACCTACGAGGGAAGCCACATCAAAGGAGAGGCCCAGGTGAAGGGGACTGGTTTCCCTGCTGACGGTCCTGTGATGACCAACTCGCTGACCGCTGCGGACTGGTGCAGGTCGAAGAAGACTTACCCCAACGACAAAACCATCATCAGTACCTTTAAGTGGAGTTACACCACTGGAAATGGCAAGCGCTACCGGAGCACTGCGCGGACCACCTACACCTTTGCCAAGCCAATGGCGGCTAACTATCTGAAGAACCAGCCGATGTACGTGTTCCGTAAGACGGAGCTCAAGCACTCCAAGACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATGTGATGGGCATGGACGAGCTGTACAAG",
    "mScarlet": "ATGGTGAGCAAGGGCGAGGCAGTGATCAAGGAGTTCATGCGGTTCAAGGTGCACATGGAGGGCTCCATGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCTCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAGGGCCTTCACCAAGCACCCCGCCGACATCCCCGACTACTATAAGCAGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGCCGTGACCGTGACCCAGGACACCTCCCTGGAGGACGGCACCCTGATCTACAAGGTGAAGCTCCGCGGCACCAACTTCCCTCCTGACGGCCCCGTAATGCAGAAGAAGACAATGGGCTGGGAAGCGTCCACCGAGCGGTTGTACCCCGAGGACGGCGTGCTGAAGGGCGACATTAAGATGGCCCTGCGCCTGAAGGACGGCGGCCGCTACCTGGCGGACTTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGATGCCCGGCGCCTACAACGTCGACCGCAAGTTGGACATCACCTCCCACAACGAGGACTACACCGTGGTGGAACAGTACGAACGCTCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG",
    "mCherry": "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAG",
    "sTagRFP": "ATGGTGTCTAAGGGCGAGGAACTGATTAAGGAGAATATGCACATGAAGCTGTACATGGAGGGCACCGTGAACAACCACCACTTCAAATGCACCTCCGAGGGCGAAGGCAAGCCCTACGAGGGCACCCAGACCATGAGAATCAAGGTGGTCGAGGGCGGCCCTCTCCCCTTCGCCTTCGACATCCTGGCTACCAGCTTCATGTACGGCAGCAGAACCTTCATCAACCACACCCAGGGCATCCCCGATTTCTTTAAGCAGTCCTTCCCAGAGGGCTTCACATGGGAGAGAGTCACCACATACGAGGACGGGGGCGTGCTGACCGCCACCCAGGACACCAGCCTCCAGGACGGCTGCCTCATCTACAACGTCAAGATCAGAGGGGTGAACTTCCCATCCAACGGCCCTGTGATGCAGAAGAAAACACTCGGCTGGGAGGCCAACACCGAGATGCTGTACCCCGCTGACGGCGGCCTGGAAGGCAGAACCGtCATGGCCCTGAAGCTCGTGGGCGGGGGCCACCTGATCTGCAACTTCAAGACCACATACAGGTCCAAGAAACCCGCTAAGAACCTGAAGATGCCCGGAGTGTACTATGTGGACCACAGACTGGAGAGAATCAAGGAGGCCGACAAAGAGACATACGTCGAGCAGCACGAGGTGGCTGTGGCCAGATACTGCGACCTCCCTAGCAAACTGGGCCACAAGCTGAACGGCATGGACGAGCTGTACAAG",
    "miRFP670nano3": "ATGGCAAACCTGGACAAGATGCTGAACACCACCGTGACCGAGGTGCGCAAGTTCCTGCAAGCAGACAGAGTGTGCGTGTTCAAGTTCGAGGAAGATTACTCCGGCACCGTCAGCCACGAAGCCGTGGACGACAGATGGATTAGCATCCTGAAAACCCAGGTGCAGGACAGATACTTCATGGAAACCAGAGGCGAGGAATACGTCCACGGCAGATACCAGGCCATCGCCGACATCTACACAGCCAATCTGGTCGAGTGCTACAGAGACCTGCTGATCGAGTTTCAGGTGCGGGCCATTCTGGCTGTCCCCATCCTGCAAGGCAAGAAGCTGTGGGGCCTGCTGGTGGCCCACCAACTGGCCGGCCCTCGGGAGTGGCAGACCTGGGAAATCGACTTCCTGAAACAGCAAGCCGTGGTGATGGGCATCGCCATCCAGCAGAGC",
    "miniTurbo": "ATGGCGATCCCGCTGCTGAACGCTAAACAGATTCTGGGACAGCTGGACGGCGGGAGCGTGGCAGTCCTGCCTGTGGTCGACTCCACCAATCAGTACCTGCTGGATCGAATCGGCGAGCTGAAGAGTGGGGATGCTTGCATTGCAGAATATCAGCAGGCAGGGAGAGGAAGCAGAGGGAGGAAATGGTTCTCTCCTTTTGGAGCTAACCTGTACCTGAGTATGTTTTGGCGCCTGAAGCGGGGACCAGCAGCAATCGGCCTGGGCCCGGTCATCGGAATTGTCATGGCAGAAGCGCTGCGAAAGCTGGGAGCAGACAAGGTGCGAGTCAAATGGCCCAATGACCTGTATCTGCAGGATAGAAAGCTGGCAGGCATCCTGGTGGAGCTGGCCGGAATAACAGGCGATGCTGCACAGATCGTCATTGGCGCCGGgatTAACGTGGCTATGAGGCGCGTGGAGGAAAGCGTGGTCAATCAGGGCTGGATCACACTGCAGGAAGCAGGGATTAACCTGGACAGGAATACTCTGGCCGCTATGCTGATCCGAGAGCTGCGGGCAGCCCTGGAACTGTTCGAGCAGGAAGGCCTGGCTCCATATCTGTCACGGTGGGAGAAGCTGGATAACTTCATCAATAGACCCGTGAAGCTGATCATTGGGGACAAAGAGATTTTCGGGATTAGCCGGGGGATTGATAAACAGGGAGCCCTGCTGCTGGAACAGGACGGAGTTATCAAACCCTGGATGGGCGGAGAAATCAGTCTGCGGTCTGCCGAAAAG",
    "ultraID": "ATGTTCAAGAACCTGATCTGGCTGAAGGAGGTGGACAGCACCCAGGAGAGACTGAAGGAGTGGAACGTGTCCTACGGCACCGCCCTGGTGGCCGACAGACAGACCAAGGGCAGAGGCGGCCCCGGCAGAAAGTGGCTGAGCCAGGAGGGCGGCCTGTACTTCAGCTTCCTGCTGAACCCCAAGGAGTTCGAGAACCTGCTGCAGCTGCCCCTGGTGCTGGGCCTGAGCGTGAGCGAGGCCCTGGAGGAGATCACCGAGATCCCCTTCAGCCTGAAGTGGCCCAACGACGTGTACTTCCAGGAGAAGAAGGTGAGCGGCGTGCTGTGCGAGCTGAGCAAGGACAAGCTGATCGTGGGCATCGGCATCAACGTGAACCAGAGAGAGATCCCCGAGGAGATCAAGGACAGAGCCACCACCCTGTACGAGATCACCGGCAAGGACTGGGACAGAAAGGAGGTGCTGCTGAAGGTGCTGAAGAGAATCAGCGAGAACCTGAAGAAGTTCAAGGAGAAG",
    "dTAG": "GGAGTGCAGGTGGAAACCATCTCCCCAGGAGACGGGCGCACCTTCCCCAAGCGCGGCCAGACCTGCGTGGTGCACTACACCGGGATGCTTGAAGATGGAAAGAAAGTTGATTCCAGCCGGGACAGAAACAAGCCCTTTAAGTTTATGCTAGGCAAGCAGGAGGTGATCCGAGGCTGGGAAGAAGGGGTTGCCCAGATGAGTGTGGGTCAGAGAGCCAAACTGACTATATCTCCAGATTATGCCTATGGTGCCACTGGGCACCCAGGCATCATCCCACCACATGCCACTCTCGTCTTCGATGTGGAGCTTCTAAAACTGGAA",
    "3XFLAG": "GATTACAAGGATGACGACGATAAGGACTATAAGGACGATGATGACAAGGACTACAAAGATGATGACGATAAA",
    "3XHA": "TACCCATACGATGTTCCAGATTACGCTTACCCCTACGACGTGCCTGATTATGCCTACCCATACGATGTGCCAGACTATGCC",
    "V5": "GGCAAGCCCATCCCCAACCCCCTGCTGGGCCTGGACAGCACC",
    "lox66": "ATAACTTCGTATAGCATACATTATACGAACGGTA",
    "lox71": "TACCGTTCGTATAGCATACATTATACGAAGTTAT",
    "loxP": "ATAACTTCGTATAGCATACATTATACGAAGTTAT",
    "SBP": "GACGAGAAGACCACTGGTTGGCGAGGTGGACACGTTGTTGAAGGACTGGCTGGGGAACTTGAACAACTTCGTGCACGACTGGAGCATCACCCACAAGGTCAACGTGAACCA",
    "3Flag": "GACTACAAAGACCATGACGGTGATTATAAAGATCATGACATCGATTACAAGGATGACGATGACAAG",
    "SBP3Flag": "GACGAGAAGACCACTGGTTGGCGAGGTGGACACGTTGTTGAAGGACTGGCTGGGGAACTTGAACAACTTCGTGCACGACTGGAGCATCACCCACAAGGTCAACGTGAACCAGGAGGCAGCGACTACAAAGACCATGACGGTGATTATAAAGATCATGACATCGATTACAAGGATGACGATGACAAG",
    "3FlagSBP": "GGAGGCAGCGACGAGAAGACCACTGGTTGGCGAGGTGGACACGTTGTTGAAGGACTGGCTGGGGAACTTGAACAACTTCGTGCACGACTGGAGCATCACCCACAAGGTCAACGTGAACCA",
    "none": "",
}

tagToColor = {
    "eGFP": "#66ff33",
    "Streptavidin": "#bfbfbf",
    "(GGGGS)x2": "#ffff66",
    "GGGGS": "#ffff66",
    "GSGGG": "#ffff66",
    "XTEN": "#ffff66",
    "Blast": "#bfbfbf",
    "Puro": "#bfbfbf",
    "Zeo": "#bfbfbf",
    "moxGFP": "#66ff33",
    "mScarlet": "#ff3300",
    "2A ribosomal skipping peptide": "#cc99ff",
    "EF1α promoter": "#cc99ff",
    "mStrayGold": "#66ff33",
    "mNeon": "#99ff66",
    "moxGFP": "#66ff33",
    "mScarlet": "#ff3300",
    "mCherry": "#ff3300",
    "sTagRFP": "#9932CC",
    "miRFP670nano3": "#ff3300",
    "miniTurbo": "#bfbfbf",
    "ultraID": "#bfbfbf",
    "dTAG": "#bfbfbf",
    "FLAG": "#ccccff",
    "3XFLAG": "#ccccff",
    "3Flag": "#bfbfbf",
    "SBP3Flag": "#bfbfbf",
    "3FlagSBP": "#bfbfbf",
    "SBP": "#ccccff",
    "HA": "#bfbfbf",
    "V5": "#bfbfbf",
    "loxP": "#ffff66",
    "lox66": "#ffff66",
    "lox71": "#ffff66",
    "none": "#bfbfbf",
}
# labels and descriptions of eff. scores
scoreDescs = {
    "doench": (
        "Doench '14",
        "Range: 0-100. Linear regression model trained on 880 guides transfected into human MOLM13/NB4/TF1 cells (three genes) and mouse cells (six genes). Delivery: lentivirus. The Fusi score can be considered an updated version this score, as their training data overlaps a lot. See <a target='_blank' href='http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html'>Doench et al.</a>",
    ),
    "wuCrispr": (
        "Wu-Crispr",
        "Range 0-100. Aka 'Wong score'. SVM model trained on previously published data. The aim is to identify only a subset of efficient guides, many guides will have a score of 0. Takes into account RNA structure. See <a target='_blank' href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0784-0'>Wong et al., Gen Biol 2015</a>",
    ),
    "ssc": (
        "Xu",
        "Range ~ -2 - +2. Aka 'SSC score'. Linear regression model trained on data from &gt;1000 genes in human KBM7/HL60 cells (Wang et al) and mouse (Koike-Yusa et al.). Delivery: lentivirus. Ranges mostly -2 to +2. See <a target='_blank' href='http://genome.cshlp.org/content/early/2015/06/10/gr.191452.115'>Xu et al.</a>",
    ),
    "crisprScan": [
        "Moreno-Mateos",
        "Also called 'CrisprScan'. Range: mostly 0-100. Linear regression model, trained on data from 1000 guides on &gt;100 genes, from zebrafish 1-cell stage embryos injected with mRNA. See <a target=_blank href='http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3543.html'>Moreno-Mateos et al.</a>. Recommended for guides transcribed <i>in-vitro</i> (T7 promoter). Click to sort by this score. Note that under 'Show all scores', you can find a Doench2016 model trained on Zebrafish scores, Azimuth in-vitro, which should be slightly better than this model for zebrafish.",
    ],
    "wang": (
        "Wang",
        "Range: 0-100. SVM model trained on human cell culture data on guides from &gt;1000 genes. The Xu score can be considered an updated version of this score, as the training data overlaps a lot. Delivery: lentivirus. See <a target='_blank' href='http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3972032/'>Wang et al.</a>",
    ),
    "chariRank": (
        "Chari",
        "Range: 0-100. Support Vector Machine, converted to rank-percent, trained on data from 1235 guides targeting sequences that were also transfected with a lentivirus into human 293T cells. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v12/n9/abs/nmeth.3473.html'>Chari et al.</a>",
    ),
    "fusi": (
        "Doench '16",
        "Aka the 'Fusi-Score', since V4.4 using the version 'Azimuth', scores are slightly different than before April 2018 but very similar (click 'show all' to see the old scores). Range: 0-100. Boosted Regression Tree model, trained on data produced by Doench et al (881 guides, MOLM13/NB4/TF1 cells + unpublished additional data). Delivery: lentivirus. See <a target='_blank' href='http://biorxiv.org/content/early/2015/06/26/021568'>Fusi et al. 2015</a> and <a target='_blank' href='http://www.nature.com/nbt/journal/v34/n2/full/nbt.3437.html'>Doench et al. 2016</a> and <a target=_blank href='https://crispr.ml/'>crispr.ml</a>. Recommended for guides expressed in cells (U6 promoter). Click to sort the table by this score.",
    ),
    "fusiOld": (
        "OldDoench '16",
        "The original implementation of the Doench 2016 score, as received from John Doench. The scores are similar, but not exactly identical to the 'Azimuth' version of the Doench 2016 model that is currently the default on this site, since Apr 2018.",
    ),
    "rs3": (
        "Doench-RS3",
        "The Doench Rule Set 3 (RS3) score (-200-+200). Similar to the Doench 2014 and Doench 2016/Fusi/Azimuth score, but updated and more accurate. See <a href='https://www.nature.com/articles/s41467-022-33024-2' target=_blank>. Scores shown are multiplied with 100 for easier display. RS3 is configured here to use the Hsu-TRACR sequence.",
    ),
    "najm": (
        "Najm 2018",
        "A modified version of the Doench 2016 score ('Azimuth'), by Mudra Hegde for S. aureus Cas9. Range 0-100. See <a target=_blank href='https://www.nature.com/articles/nbt.4048'>Najm et al 2018</a>.",
    ),
    "ccTop": ("CCTop", "The efficiency score used by CCTop, called 'crisprRank'."),
    "aziInVitro": (
        "Azimuth in-vitro",
        "The Doench 2016 model trained on the Moreno-Mateos zebrafish data. Unpublished model, gratefully provided by J. Listgarden. This should be better than Moreno-Mateos, but we have not found the time to evaluate it yet.",
    ),
    "housden": (
        "Housden",
        "Range: ~ 1-10. Weight matrix model trained on data from Drosophila mRNA injections. See <a target='_blank' href='http://stke.sciencemag.org/content/8/393/rs9.long'>Housden et al.</a>",
    ),
    "proxGc": ("ProxGCCount", "Number of GCs in the last 4pb before the PAM"),
    "seqDeepCpf1": (
        "DeepCpf1",
        "Range: ~ 0-100. Convolutional Neural Network trained on ~20k Cpf1 lentiviral guide results. This is the score without DNAse information, 'Seq-DeepCpf1' in the paper. See <a target='_blank' href='https://www.nature.com/articles/nbt.4061'>Kim et al. 2018</a>",
    ),
    "oof": (
        "Out-of-Frame",
        "Range: 0-100. Out-of-Frame score, only for deletions. Predicts the percentage of clones that will carry out-of-frame deletions, based on the micro-homology in the sequence flanking the target site. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html'>Bae et al. 2014</a>. Click the score to show the predicted deletions.",
    ),
    "lindel": (
        "Lindel",
        "Wei Chen Frameshift ratio (0-100). Predicts probability of a frameshift caused by any type of insertion or deletion. See <a href='https://academic.oup.com/nar/article/47/15/7989/5511473'>Wei Chen et al, Bioinf 2018</a>. Click the score to see the most likely deletions and insertions.",
    ),
    "EVA": (
        "EVA score",
        "EVA score, predicts the activity of synthetic guides (range 0-100). This score is derived from a linear model, based on data from human cell lines edited with synthetic gRNAs. See <a href='https://doi.org/10.1038/s41467-025-59947-0'>Riesenberg et al. 2025</a>",
    ),
}

# the headers for the guide and offtarget output files
guideHeaders = [
    "guideId",
    "targetSeq",
    "mitSpecScore",
    "cfdSpecScore",
    "offtargetCount",
    "targetGenomeGeneLocus",
]
offtargetHeaders = [
    "guideId",
    "guideSeq",
    "offtargetSeq",
    "mismatchPos",
    "mismatchCount",
    "mitOfftargetScore",
    "cfdOfftargetScore",
    "chrom",
    "start",
    "end",
    "strand",
    "locusDesc",
]

# library descriptions
libLabels = [
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4486245/
    ("human_brunello", "Human, Brunello, Doench Nat Bio 2016 (recommended)"),
    ("human_avana", "Human, Avana, Doench Nat Bio 2016"),
    ("human_geckov2", "Human, GeCKO V2, Sanjana Nat Meth 2014"),
    ("mouse_brie", "Mouse, Brie, Doench Nat Bio 2016 (recommended)"),
    ("mouse_geckov2", "Mouse, GeCKO V2, Sanjana Nat Meth 2014"),
    ("mouse_asiago", "Mouse, Asiago, Doench Nat Bio 2016"),
]

# a file crispor.conf in the directory of the script allows to override any global variable
myDir = dirname(__file__)
confPath = join(myDir, "crispor.conf")
if isfile(confPath):
    exec(open(confPath).read())
    # execfile(confPath)

cgiParams = None

# ====== END GLOBALS ============


def setupPamInfo(pam):
    "modify a few globals based on the current pam"
    global GUIDELEN

    global pamIsFirst
    global addGenePlasmids
    global PAMLEN
    global scoreNames
    global baseEditor
    global saCas9Mode
    global mutScoreNames
    global isSpg

    baseEditor = None
    saCas9Mode = False
    isSpg = False

    PAMLEN = len(pam)

    GUIDELEN = 20
    pamIsFirst = False
    scoreNames = cas9ScoreNames

    pamOpt = None
    if "-" in pam:
        pam, pamOpt = pam.split("-")
        if pamOpt == "BE1":
            baseEditor = "BE1"
        elif pamOpt == "spg":
            isSpg = True

    if "custom" in pam:
        logging.debug("switching on custom PAM mode")
        pamstr = [elements for elements in pam.split(".")]
        pam = pamstr[
            0
        ]  # may break subsequent calls of setupPamInfo(pam) (but not setupPamInfo(pamDesc)) ?
        PAMLEN = len(pam)
        GUIDELEN = int(pamstr[2])
        ezType = pamstr[1]
        if ezType == "Cas12a":
            pamIsFirst = True
            scoreNames = cpf1ScoreNames
        elif ezType in ["CBE", "ABE"]:
            baseEditor = True

    elif pamIsCasX(pam):
        logging.debug("switching on CasX mode, guide length is 20bp")
        GUIDELEN = 20
        pamIsFirst = True
        scoreNames = cpf1ScoreNames
    elif pamIsCas12max(pam):
        logging.debug("switching on hfCas12max mode, guide length is 20bp")
        GUIDELEN = 20
        pamIsFirst = True
        scoreNames = cpf1ScoreNames
    elif pamIsCpf1(pam):
        logging.debug("switching on Cpf1 mode, guide length is 23bp")
        GUIDELEN = 23
        pamIsFirst = True
        scoreNames = cpf1ScoreNames
        # if pamOpt:
        # GUIDELEN=int(pamOpt)
    elif pam == "NGTN":
        logging.debug("switching on Cpf1 mode for ShCAST, guide length is 23bp")
        GUIDELEN = 23
        pamIsFirst = True
    elif pam == "NNNNRYAC" or pam == "NNNVRYAC":
        GUIDELEN = 22
    elif pam == "NNGRRT" or pam == "NNNRRT":
        logging.debug("switching on S. aureus mode, guide length is 21bp")
        addGenePlasmids = addGenePlasmidsAureus
        GUIDELEN = 21
        # if pamOpt=="20":
        # GUIDELEN=20
        saCas9Mode = True
        scoreNames = saCas9ScoreNames
    elif pam == "NNNNCC":
        GUIDELEN = 24
    elif pam == "NGG-22":
        GUIDELEN = 22
    else:
        GUIDELEN = 20

    if pamOpt and pamOpt.isnumeric():
        GUIDELEN = int(pamOpt)

    if (GUIDELEN == 20 or GUIDELEN == 22) and pam == "NGG":
        mutScoreNames = spCas9MutScoreNames
    else:
        mutScoreNames = otherMutScoreNames

    logging.debug(
        "Enzyme info: pam=%s, guideLen=%d, pamIsFirst=%s, saCas9Mode=%s"
        % (pam, GUIDELEN, pamIsFirst, saCas9Mode)
    )

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
        "CREATE TABLE IF NOT EXISTS %s "
        "("
        "  jobType text,"  # either "index" or "search"
        "  jobId text %s,"  # unique identifier
        "  paramStr text,"  # parameters for jobs, like db, options, etc.
        "  isRunning int DEFAULT 0,"  # indicates steps have started, done jobs are moved to doneJobs table
        "  stepName text,"  # currently step, internal step name for timings
        "  stepLabel text,"  # current step, human-readable status of job, for UI
        "  lastUpdate float,"  # time of last update
        "  stepTimes text,"  # comma-sep list of whole msecs, one per step
        "  startTime text "  # date+time when job was put into queue
        ")"
    )

    # def __init__(self):
    # " no inheritance needed here "
    # self.openSqlite(JOBQUEUEDB)

    def openSqlite(self, dbName=JOBQUEUEDB):
        self.dbName = dbName
        # isolation_level=None = autocommit mode: we manage transactions explicitly
        # timeout=30: wait up to 30s for locks (multiple daemons + CGI share this DB)
        self.conn = sqlite3.connect(dbName, timeout=30, isolation_level=None)
        # WAL mode allows concurrent readers + writer, essential for CGI + daemon(s)
        result = self.conn.execute("PRAGMA journal_mode=WAL;").fetchone()
        if result[0] != "wal":
            logging.warn(
                "Could not set WAL mode, journal_mode is: %s (is another connection open with the old mode?)"
                % result[0]
            )
        # self.conn.set_trace_callback(print) # for debugging: print all sql statements
        self._chmodJobDb()

        try:
            self.conn.execute(self._queueDef % ("queue", "PRIMARY KEY"))
        except SQLITEERROR as ex:
            errAbort("cannot open the sqlite jobs file %s: %s" % (JOBQUEUEDB, ex))

    def _chmodJobDb(self):
        # umask is not respected by sqlite, bug http://www.mail-archive.com/sqlite-users@sqlite.org/msg59080.html
        try:
            os.chmod(JOBQUEUEDB, 0o666)
        except OSError:
            # if the file was created by other job, we can't chmod, as we're the CGI. Just silently ignore this
            pass

    def addJob(self, jobType, jobId, paramStr):
        "create a new job, returns False if not successful"
        self._chmodJobDb()

        sql = (
            "INSERT INTO queue (jobType, jobId, isRunning, lastUpdate, "
            "stepTimes, paramStr, stepName, stepLabel, startTime) VALUES (:jobType, :jobId, :isRunning, :lastUpdate, "
            ":stepTimes, :paramStr, :stepName, :stepLabel, :startTime)"
        )
        now = "%.3f" % time.time()
        values = {
            "jobType": jobType,
            "jobId": jobId,
            "isRunning": 0,
            "lastUpdate": now,
            "stepTimes": "",
            "paramStr": paramStr,
            "stepName": "wait",
            "stepLabel": "Waiting",
            "startTime": now,
        }
        try:
            # in autocommit mode, this INSERT commits immediately
            self.conn.execute(sql, values)
            return True
        except sqlite3.IntegrityError:
            # job already in queue (e.g. resubmit, or daemon restart) - that's fine
            return True
        except SQLITEERROR:
            errAbort(
                "Cannot open DB file %s. Please contact %s"
                % (self.dbName, contactEmail)
            )

    def getStatus(self, jobId):
        "return current job status label or None if job is not in queue"
        sql = "SELECT stepLabel FROM queue WHERE jobId=?"
        try:
            rows = self.conn.execute(sql, (jobId,)).fetchmany(1)
            if len(rows) == 0:
                status = None
            else:
                status = rows[0][0]
        except (StopIteration, IndexError):
            logging.debug("getStatus: job %s not found" % jobId)
            status = None
        return status

    def dump(self):
        "for debugging, write the whole queue table to stdout"
        sql = "SELECT * FROM queue"
        for row in self.conn.execute(sql):
            print("\t".join([str(x) for x in row]))

    def jobInfo(self, jobId, isDone=False):
        "for debugging, return all job info as a tuple"
        print("job info<br>")
        if isDone:
            sql = "SELECT * FROM doneJobs WHERE jobId=?"
        else:
            sql = "SELECT * FROM queue WHERE jobId=?"
        try:
            row = next(self.conn.execute(sql, (jobId,)))
        except StopIteration:
            return []
        return row

    def startStep(self, jobId, newName, newLabel):
        "start a new step. Update lastUpdate, status and stepTime"
        try:
            self.conn.execute("BEGIN IMMEDIATE")

            sql = "SELECT lastUpdate, stepTimes, stepName FROM queue WHERE jobId=?"
            logging.debug(sql)
            rows = self.conn.execute(sql, (jobId,)).fetchmany(1)
            if len(rows) == 0:
                logging.error("startStep: no row for jobId %s" % jobId)
                self.conn.commit()
                return

            lastTime, timeStr, lastStep = rows[0]
            lastTime = float(lastTime)

            # append a string in format "stepName:milliSecs" to the timeStr
            now = time.time()
            timeDiff = "%d" % int((1000.0 * (now - lastTime)))
            newTimeStr = timeStr + "%s=%s" % (lastStep, timeDiff) + ","
            sql = "UPDATE queue SET lastUpdate=?, stepName=?, stepLabel=?, stepTimes=?, isRunning=? WHERE jobId=?"
            self.conn.execute(sql, (now, newName, newLabel, newTimeStr, 1, jobId))

            self.conn.commit()
        except:
            self.conn.rollback()
            raise

    def jobDone(self, jobId):
        "remove the job from the queue and add it to the queue log"
        print("job done<br>")
        try:
            self.conn.execute("BEGIN IMMEDIATE")

            sql = "SELECT * FROM queue WHERE jobId=?"
            row = self.conn.execute(sql, (jobId,)).fetchone()
            if row is None:
                logging.warn("jobDone - job %s has been removed already" % jobId)
                self.conn.commit()
                return

            sql = "DELETE FROM queue WHERE jobId=?"
            self.conn.execute(sql, (jobId,))

            self.conn.commit()
        except:
            self.conn.rollback()
            raise

        # good to have a log file of the old jobs
        with open(
            "doneJobs.tsv", "a"
        ) as ofh:  # if this triggers an error: run 'touch doneJobs.tsv && chmod a+rw doneJobs.tsv' in the crispor dir.
            row = [str(x) for x in row]
            line = "\t".join(row)
            ofh.write(line)
            ofh.write("\n")

    def waitCount(self):
        "return number of waiting jobs"
        sql = "SELECT count(*) FROM queue WHERE isRunning=0"
        return self.conn.execute(sql).fetchone()[0]

    def popJob(self):
        "return (jobType, jobId, params) of first waiting job and set it to running state"
        print("pop job<br>")
        try:
            self.conn.execute("BEGIN IMMEDIATE")

            sql = "SELECT jobType, jobId, paramStr FROM queue WHERE isRunning=0 ORDER BY lastUpdate LIMIT 1"
            row = self.conn.execute(sql).fetchone()
            if row is None:
                logging.debug("popJob: no waiting jobs")
                self.conn.commit()
                return None, None, None

            jobType, jobId, paramStr = row

            sql = "UPDATE queue SET isRunning=1 where jobId=?"
            self.conn.execute(sql, (jobId,))

            self.conn.commit()
        except:
            self.conn.rollback()
            raise

        return jobType, jobId, paramStr

    def clearJobs(self):
        "clear the job table, removing running jobs, too"
        self.conn.execute("DELETE from queue")

    def close(self):
        " "
        self.conn.close()


# ====== FUNCTIONS =====
contentLineDone = False

# the queue workers should be able to never abort
doAbort = True


def getTwoBitFname(db):
    "return the name of the twoBit file for a genome"
    # at UCSC, try to use local disk, if possible
    locPath = join("/scratch", "data", db, db + ".2bit")
    if isfile(locPath):
        return locPath
    path = join(genomesDir, db, db + ".2bit")
    return path


def errAbort(msg, isWarn=False):
    "print err msg and exit"
    if commandLineMode:
        raise Exception(msg)

    if not contentLineDone:
        print("Content-type: text/html\n")

    print(
        '<div style="position: absolute; padding: 10px; left: 100; top: 100; border: 10px solid black; background-color: white; text-align:left; width: 800px; font-size: 18px">'
    )

    if isWarn:
        print("<strong>Warning:</strong><p> ")
    else:
        print("<strong>Error:</strong><p> ")

    print((msg + "<p>"))
    print(
        (
            "If you think this is a bug or you have any other suggestions, please do not hesitate to contact us %s<p>"
            % contactEmail
        )
    )
    if isWarn:
        print("In the email, please also send us the full URL of the page.")
    else:
        print(
            "Please also send us the full URL of the page where you see the error. Thanks!"
        )
    print("</div>")

    if doAbort:
        sys.exit(0)  # cgi must not exit with 1


# allow only dashes, digits, characters, underscores and colons in the CGI parameters
# and +
notOkChars = re.compile(r"[^+a-zA-Z0-9/:\n\r_. -]")


def checkVal(key, inStr):
    """remove special characters from input string, to protect against injection attacks"""
    if key != "geneIds":
        if len(inStr) > 10000:
            errAbort("input parameter %s is too long" % key)
    else:
        if len(inStr) > 100000:
            errAbort(
                "Pasting more than tens of thousands of gene IDs makes little sense. Copy/paste error?"
            )

    matchObj = notOkChars.search(inStr)
    if matchObj != None:
        errAbort(
            "input parameter %s contains an invalid character %s (ASCII %d)"
            % (key, repr(matchObj.group()), ord(matchObj.group()))
        )
    return inStr


def cgiGetParams():
    "get CGI parameters and return as dict"
    form = cgi.FieldStorage()
    global cgiParams
    cgiParams = {}

    # parameters are:
    # "pamId", "batchId", "pam", "seq", "org", "download", "sortBy", "format", "ajax
    for key in list(form.keys()):
        val = form.getfirst(key)
        if val != None:
            # "seq" is cleaned by cleanSeq later
            val = urllib.parse.unquote(val)
            if key not in [
                "seq",
                "name",
                "koGeneId",
                "customseq",
                "insertSeq",
                "globEffScore",
                "linkerseq",
                "tagseq",
                "markerseq",
                "qTag",
                "expressionSeq",
                "exonInfo",
                "guideInfo",
                "geneModel",
                "tagNames",
                "revGuideInfo",
                "fwGuideInfo"
            ]:
                checkVal(key, val)
            cgiParams[key] = val

    if "pam" in cgiParams:
        legalChars = set("ACTGNHMKRYVBE120345/-")
        illegalChars = set(cgiParams["pam"]) - legalChars
        if len(illegalChars) != 0:
            errAbort(
                "Illegal character in PAM-sequence. Only %s are allowed."
                + "".join(legalChars)
            )

    if "batchId" in cgiParams:
        batchId = cgiParams["batchId"]
        if not batchId.isalnum() or len(batchId) > 30:
            errAbort("Invalid batchId")

    return cgiParams


def cgiGetStr(params, argName, default=None):
    val = params.get(argName, None)
    if val == None and default == None:
        errAbort("'%s' parameter must be specified" % argName)
    if val == None:
        return default
    return val


def cgiGetNum(params, argName, default):
    "get CGI parameter which must be a number"
    val = params.get(argName, None)
    if val == None:
        return default
    if not val.isdigit():
        errAbort("'%s' parameter must be a number" % argName)
    val = int(val)
    return val


transTab = str.maketrans("-=/+_", "abcde")


def makeTempBase(seq, org, pam, batchName):
    "create the base name of temp files using a hash function and some prettyfication"
    hasher = hashlib.sha1(
        seq.encode("latin1")
        + org.encode("latin1")
        + pam.encode("latin1")
        + batchName.encode("latin1")
    )
    shortHash = hasher.digest()[0:20]
    batchId = (
        base64.urlsafe_b64encode(shortHash).decode("latin1").translate(transTab)[:20]
    )
    return batchId


def makeTempFile(prefix, suffix):
    "return a temporary file that is deleted upon exit, unless DEBUG is set"
    if DEBUG:
        fname = join("/tmp", prefix + suffix)
        fh = open(fname, "wt")
    else:
        fh = tempfile.NamedTemporaryFile(
            mode="wt", dir=TEMPDIR, prefix="primer3In", suffix=".txt"
        )
    return fh


def pamIsCpf1(pam):
    "if you change this, also change bin/filterFaToBed and bin/samToBed!!!"
    return pam in [
        "TNN",
        "TTN",
        "TTTN",
        "TYCV",
        "TATV",
        "TTTV",
        "TTTR",
        "ATTN",
        "TTTA",
        "TCTA",
        "TCCA",
        "CCCA",
        "YTTV",
        "TTYN",
    ]


def pamIsCas12max(pam):
    return pam in ["TNN", "TTN"]


def pamIsCasX(pam):
    "if you change this, also change bin/filterFaToBed and bin/samToBed!!!"
    return pam in ["TTCN"]


def pamIsSaCas9(pam):
    "only used for notes and efficiency scores, unlike its Cpf1 cousin function"
    return pam.split("-")[0] in ["NNGRRT", "NNNRRT"]


def isSlowPam(pam):
    "do not allow input sequences > 500 bp"
    if pamIsXCas9(pam) or pam == "TTYN" or pam == "NNG" or pam == "TNN":
        return True
    else:
        return False


def pamIsXCas9(pam):
    " "
    return pam in ["NGK", "NGN"]


def pamIsSpCas9(pam):
    "only used for notes and efficiency scores, unlike its Cpf1 cousin function"
    return pam in ["NGG", "NGA", "NGCG"]


def saveSeqOrgPamToCookies(seq, org, pam, koMethod, multipam, expType):
    "create a cookie with seq, org and pam and print it"
    cookies = http.cookies.SimpleCookie()
    expires = 365 * 24 * 60 * 60
    if seq is not None:
        if len(seq) < 3000:
            cookies["lastseq"] = seq
        else:
            cookies["lastseq"] = (
                "(last sequence was too long, could not be saved in Internet Browser cookie)"
            )

    if expType == "ko":
        cookies["lastKOorg"] = org
        cookies["lastKOorg"]["expires"] = expires
        cookies["lastKOpam"] = pam
        cookies["lastKOpam"]["expires"] = expires
        cookies["lastKOmethod"] = koMethod
        cookies["lastKOmethod"]["expires"] = expires

    elif expType == "ki":

        cookies["lastKIorg"] = org
        cookies["lastKIorg"]["expires"] = expires
        cookies["lastKIpam"] = multipam
        cookies["lastKIpam"]["expires"] = expires

    else:
        if "lastseq" in cookies:
            cookies["lastseq"]["expires"] = expires
        cookies["lastorg"] = org
        cookies["lastorg"]["expires"] = expires
        cookies["lastpam"] = pam
        cookies["lastpam"]["expires"] = expires

    print(cookies)


def debug(msg):
    if commandLineMode:
        logging.debug(msg)
    elif DEBUG:
        sys.stderr.write(str(msg) + "\n")
        sys.stderr.write("<br>\n")


def gcContent(seq):
    "return GC content as a float"
    c = 0
    for x in seq:
        if x in ["G", "C"]:
            c += 1
    return float(c) / len(seq)


def getFreeEnergy(seq, temperature=37):
    "calculate the minimum free energy of a given RNA sequence"
    # temporary solution, should it use the viennarna python package instead ?

    progDir = binDir
    cmd = "echo %s | %s/RNAfold -T %s --noPS" % (seq.lower(), progDir, temperature)
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, encoding="utf8")
    vienna = proc.stdout.read()
    proc.wait()  # useless in this case ?
    viennaouts = [out for out in vienna.split(" ")]
    deltaG = viennaouts.pop()
    # structure (unused for now, needs to be displayed with a monospaced font)
    # seqStructure = ''.join(viennaouts)
    if proc.returncode != 0:
        errAbort("Could not run '%s'. Return code %s" % (cmd, str(proc.returncode)))
        print("ERR")
    else:
        return float(deltaG.strip().replace(")", "").replace("(", ""))


def showSecondaryStructure(params, donorSeq=None):
    "Using ViennaRNA, displays MFE structure as a plot"

    guideSeq = params["guideSeq"]
    batchId = params["batchId"]
    pamId = params["pamId"]
    addSeq = params.get("addSeq")

    if addSeq and len(addSeq) < 100:
        legalBases = ["A", "T", "G", "C", "U"]
        for base in addSeq:
            if base not in legalBases:
                doAdd = False
                strSeq = guideSeq
                break
        else:
            strSeq = guideSeq + addSeq
    else:
        strSeq = guideSeq

    if donorSeq:
        strSeq = donorSeq

    temperature = params.get("temperature", "37")
    freeEnergy = getFreeEnergy(strSeq, temperature)

    progDir = binDir
    tmpdir = tempfile.mkdtemp(dir=batchDir)

    if not donorSeq:

        printBackLink()

        print(
            """<div style="margin-left:35%; margin-right:25%; margin-bottom:10%;"> """
        )
        print(
            """<div class="title">Spacer sequence of the guide : %s</div>""" % guideSeq
        )
        if addSeq and strSeq != guideSeq:
            print("<p>3' sequence : %s</p>" % addSeq)
        print(
            """<p>Here is shown the predicted structure of the spacer sequence of the guide. This information is used to calculate the EVA activity score. Structures with a minimum free energy lower than -3.6 kcal/mol are considered detrimental to the activity of the guide.<br> You can also add a tracrRNA sequence in 3' to check the predicted structure of the guide and potential issues with guide RNA folding.</p>
        """
        )
        print(
            """<p>Free energy of this structure : <b>%s kcal/mol</b></p>""" % freeEnergy
        )

        print("<form>")

        print(
            """
        <input type="hidden" name="batchId" value="%(batchId)s"/>
        <input type="hidden" name="guideSeq" value="%(guideSeq)s"/>
        <input type="hidden" name="pamId" value="%(pamId)s"/>
        """
            % locals()
        )

        print(
            """
            <div style="display: flex; align-items: center; gap: 4px; margin-bottom: 45px;">
                <input type="range" id="temperature" name="temperature" value="%(temperature)s" min="5" max="75" style="vertical-align:middle; width:15%%;" oninput="this.nextElementSibling.value = this.value"/>
                at<output>%(temperature)s</output> &#8451
                <textarea cols=50 rows=2 style="margin-left: 12px;" name="addSeq" placeholder="Add a tracrRNA sequence here to check the predicted structure of the guide."></textarea>
                <button style="width: 60px; height: 24px; justify-self: center; margin-left: 12px;" type="submit" name="submit" value="SUBMIT"><small>update</small></button>
            </div>
            """
            % locals()
        )
    else:
        print(
            """<p style="font-weight: 600;">Free enegry : %d kcal/mol (at %s &#8451)</p>"""
            % (freeEnergy, temperature)
        )
    try:

        RNAfoldCmd = os.path.join(progDir, "RNAfold")
        RNAplotCmd = os.path.join(progDir, "RNAplot")

        cmd = "echo %s | %s -p -T %s | %s -o svg" % (
            strSeq.upper(),
            RNAfoldCmd,
            temperature,
            RNAplotCmd,
        )

        subprocess.run(cmd, shell=True, cwd=tmpdir, capture_output=True)

        files = os.listdir(tmpdir)
        for file in files:
            if file.endswith(".svg"):
                rnaPath = join(tmpdir, file)
                rnaPlot = rnaPath

        with open(rnaPlot, "r") as f:
            svgPlot = f.read()

        print("""<div> %s </div>""" % svgPlot)

    finally:
        # Cleanup
        shutil.rmtree(tmpdir)

    print("</div>")
    print("</form>")


def findHomopolymers(seq, basesCount):
    """
    find homopolymers of a given length for the specified bases, in the input sequence
    returns their positions as a list of tuples (base and position as slice)
    basesCount is a dict {'base1': length1, 'base2': length2}
    """

    homoPolymers = []
    n = 0
    start = -1
    lastBase = None

    for pos, char in enumerate(seq.upper()):
        if char == lastBase:
            n += 1
        else:
            if lastBase and lastBase in basesCount and n > basesCount[lastBase]:
                homoPolymers.append((start, pos))
            lastBase = char
            start = pos
            n = 1

    if lastBase and lastBase in basesCount and n > basesCount[lastBase]:
        homoPolymers.append(start, len(seq))

    if homoPolymers:
        return homoPolymers


def findPat(seq, pat):
    """yield positions where pat matches seq, stupid brute force search"""
    seq = seq.upper()
    pat = pat.upper()
    patLen = len(pat)
    for i in range(0, len(seq) - patLen + 1):
        subseq = seq[i : i + patLen]
        if patMatch(subseq, pat):
            yield i


def rndSeq(seqLen):
    "return random seq"
    seq = []
    alf = "ACTG"
    for i in range(0, seqLen):
        seq.append(alf[random.randint(0, 3)])
    return "".join(seq)


def cleanSeq(seq, db, split="nosplit"):
    """remove fasta header, check seq for illegal chars and return (filtered
    seq, user message) special value "random" returns a random sequence.
    split (values : split, nosplit) authorizes "/" in seq
    """
    # print repr(seq)
    if seq.startswith("random"):
        seq = rndSeq(800)
    lines = seq.strip().splitlines()
    # print "<br>"
    # print "before fasta cleaning", "|".join(lines)
    if len(lines) > 0 and lines[0].startswith(">"):
        line1 = lines.pop(0)
    # print "<br>"
    # print "after fasta cleaning", "|".join(lines)
    # print "<br>"
    if split == "split":
        legalChars = "actgACTGNn/"
    else:
        legalChars = "actgACTGNn"
    newSeq = []
    nCount = 0
    for l in lines:
        if len(l) == 0:
            continue
        for c in l:
            if c not in legalChars:
                nCount += 1
            else:
                newSeq.append(c)
    seq = "".join(newSeq)

    msgs = []
    tooLongHint = """
    Please split your input sequence into shorter sequences or use
    the <a href='downloads/'>stand-alone version</a> on your own Linux or Mac server to process longer sequences in batch.<br>
    """

    if len(seq) > MAXSEQLEN and db != "noGenome":
        errMsg = (
            "<strong>Sorry, this tool cannot handle sequences longer than %d bp</strong><br>"
            % (MAXSEQLEN)
        )
        errAbort(errMsg + tooLongHint)
    if len(seq) > MAXSEQLEN_NOGENOME and db == "noGenome":
        errMsg = (
            "<strong>Sorry, this tool cannot handle sequences longer than %d bp when using the 'No Genome' option.</strong><br>"
            % (MAXSEQLEN_NOGENOME)
        )
        errAbort(errMsg + tooLongHint)

    if nCount != 0:
        msgs.append(
            "Sequence contained %d non-ACTGN letters. They were removed." % nCount
        )

    return seq, "<br>".join(msgs)


revTbl = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "M": "K",
    "H": "D",
    "D": "H",
    "K": "M",
    "R": "Y",
    "Y": "R",
    "g": "c",
    "a": "t",
    "c": "g",
    "t": "a",
    "n": "n",
    "h": "d",
    "d": "h",
    "V": "B",
    "v": "b",
    "B": "V",
    "b": "v",
    "W": "W",
    "w": "w",
}


def revComp(seq):
    "rev-comp a dna sequence with UIPAC characters"
    newSeq = []
    for c in reversed(seq):
        newSeq.append(revTbl[c])
    return "".join(newSeq)


def docTestInit(isCpf1, guideLen):
    global pamIsFirst
    global GUIDELEN
    pamIsFirst = isCpf1
    GUIDELEN = guideLen


def findPams(seq, pam, strand, startDict, endSet, exonId=None):
    """return two values: dict with pos -> strand of PAM and set of end positions of PAMs
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
    assert pamIsFirst is not None

    if pamIsFirst:
        maxPosPlus = len(seq) - (GUIDELEN + len(pam))
        minPosMinus = GUIDELEN
    else:
        # -------------------
        #          OKOKOKOKOK
        minPosPlus = GUIDELEN
        # -------------------
        # OKOKOKOKOK
        maxPosMinus = len(seq) - (GUIDELEN + len(pam))

        # search for guides (but not PAMs) in the flanking sequences (if it exists)
        # exonId contains no information, only checks if the function is called in multiseq mode
        if exonId is not None:
            minPosMinus = (
                GUIDELEN + 6
            )  # distance to cut site + splice site (= flank sequence)
            maxPosPlus = len(seq) - minPosMinus

    # print "new search", seq, pam, "minPosPlus=",minPosPlus, "guideLen=", GUIDELEN, "<br>"
    for start in findPat(seq, pam):
        if pamIsFirst:
            # need enough flanking seq on one side
            # return("Cpf1 mode found", start,"<br>")
            if strand == "+" and start > maxPosPlus:
                continue
            if strand == "-" and start < minPosMinus:
                continue
        else:
            # return("non-Cpf1 mode found", start,"<br>")
            if strand == "+":
                if start < minPosPlus:
                    continue
                # prevent searching for PAMs in the flanking sequences
                if exonId is not None and start > maxPosPlus:
                    continue

            if strand == "-":
                if start > maxPosMinus:
                    continue
                if exonId is not None and start < minPosMinus:
                    continue

        # print "match", strand, start, end, "<br>"
        startDict[start] = strand
        end = start + len(pam)
        endSet.add(end)
    return startDict, endSet


def rulerString(maxLen):
    "return line with positions every 10 chars"
    texts = []
    for i in range(0, maxLen, 10):
        numStr = str(i)
        texts.append(numStr)
        spacer = "".join([" "] * (10 - len(numStr)))
        texts.append(spacer)
    return "".join(texts)


def varDictToHtml(varDict, seq, varShortLabel):
    "make a list of one html string per position in the sequence"
    if varDict is None:
        return None

    varHtmls = []
    for i in range(0, len(seq)):
        if not i in varDict:
            varHtmls.append(".")
        else:
            varHooverLines = []
            showStar = False  # show a star if change is non-simple SNP
            varInfos = varDict[i]
            for chrom, pos, refAll, altAll, infoDict in varInfos:
                varHooverLines.append(
                    "%s: %s &rarr; %s<br>" % (varShortLabel, refAll, altAll)
                )
                if "freq" in infoDict:
                    varHooverLines.append(
                        "&nbsp;<b>Freq:</b> %s<br>" % infoDict["freq"]
                    )
                # if "dbg" in infoDict:
                # varHooverLines.append("%s<br>" % infoDict["dbg"])
                if "varId" in infoDict:
                    varHooverLines.append("&nbsp;<b>ID:</b> %s<br>" % infoDict["varId"])
                if "ExAC" in varDict["label"]:
                    endPos = int(pos) + len(refAll)
                    varHooverLines.append(
                        '&nbsp;<a target=_blank href="http://exac.broadinstitute.org/region/%s-%s-%d">ExAC Browser</a><br>'
                        % (chrom, pos, endPos)
                    )

                if len(refAll) != 1 or len(altAll) != 1:
                    showStar = True

            if len(varInfos) != 1:
                showStar = True

            varDesc = "".join(varHooverLines)

            if showStar:
                dispChar = "*"
            else:
                dispChar = altAll
            varHtmls.append(
                "<u class='tooltipsterInteract' title='%s'>%s</u>" % (varDesc, dispChar)
            )
    return varHtmls


def cssClassesFromSeq(guideSeq, suffix=""):
    "The CSS class of guide row and links in seq viewer depend on the first nucl of guide"
    classNames = ["guideRow"]
    if guideSeq[0].upper() != "G":
        classNames.append("guideRowNoPrefixG" + suffix)
    if not guideSeq.startswith("GG"):
        classNames.append("guideRowNoPrefixGG" + suffix)
    if guideSeq[0].upper() != "A":
        classNames.append("guideRowNoPrefixA" + suffix)
    classStr = " ".join(classNames)
    return classStr


def buildCodonTable(key="codon"):
    "from http://www.petercollingridge.co.uk/tutorials/bioinformatics/codon-table/"
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    codon_table = dict(list(zip(codons, amino_acids)))

    # returns a dict {amino_acid : [codons]} instead
    codon_lists = {}
    for codon in codon_table:
        aa = codon_table[codon]
        codon_lists.setdefault(aa, []).append(codon)

    if key == "aa":
        return codon_lists
    else:
        return codon_table


def buildOneToThree():
    "return one-letter -> three-letter conversion table for amino acids"
    oneToThree = {
        "C": "Cys",
        "D": "Asp",
        "S": "Ser",
        "Q": "Gln",
        "K": "Lys",
        "I": "Ile",
        "P": "Pro",
        "T": "Thr",
        "F": "Phe",
        "N": "Asn",
        "G": "Gly",
        "H": "His",
        "L": "Leu",
        "R": "Arg",
        "W": "Trp",
        "A": "Ala",
        "V": "Val",
        "E": "Glu",
        "Y": "Tyr",
        "M": "Met",
        "U": "Sec",
        "*": "Stop",
        "X": "Stop",  # is this really used like that?
        "Z": "Glx",  # special case: asparagine or aspartic acid
        "B": "Asx",  # special case: glutamine or glutamic acid
    }
    return oneToThree


def getSpliceSites(seq, exStart, exEnd):
    """
    Return a list of tuples [(base, start, end)] of splice site bases in a given exon.
    Exon coordinates are replative to the sequence
    """

    spliceBases = []

    if exStart > 2:
        spliceAcc = seq[exStart - 2: exStart]
    else:
        spliceAcc = ""

    if exEnd < len(seq) - 2:
        spliceDon = seq[exEnd: exEnd + 2]
    else:
        spliceDon = ""

    if spliceAcc == "AG":
        for pos in range(exStart - 2, exStart):
            spliceBases.append((seq[pos], pos))

    if spliceDon == "GT":
        for pos in range(exEnd, exEnd + 2):
            spliceBases.append((seq[pos], pos))

    return spliceBases


def makeExonLines(exonInfo, seq, selTransId, koMethod=None, editInfo=None, loadStopGuides=None):
    """create text that draws exons, input is transId -> (exonNumber, exStart, exEnd, exFrame).
    returns a list of (transId (=label), symbol (=mouseover), ASCII-line)"""
    lines = []
    # maxLabelLen = 0
    codonTable = buildCodonTable()
    # oneToThree = buildOneToThree()
    seqLen = len(seq)
    seq = seq.upper()

    if editInfo:
        # need to add ABE / CBE instead of pam (which is just the pam motif)
        pam, editData = editInfo

        # {codon, (position in codon, editor)}
        possibleStops = {"TGG": (1, "CBE"),
                         "CAG": (0, "CBE"),
                         "CGA": (0, "CBE"),
                         "CAA": (0, "CBE"),
                         "TCG": (1, "CGBE"),
                         "TAC": (2, "CGBE"),
                         "TCA": (1, "CGBE"),
                         "TGC": (2, "CGBE")
                         }

        # CGBE can be used for C / T, but CBE is more efficient
        possibleSpliceEdits = {"C": "CBE",
                               "G": "CBE",
                               "A": "ABE",
                               "T": "ABE"
                               }
    else:
        editData = None

    # list of guideIds that result in STOP codons
    if koMethod is not None and editInfo is not None and loadStopGuides is None:
        stopGuides = {}
    else:
        stopGuides = loadStopGuides

    for (transId, symbol), exRows in exonInfo.items():
        if selTransId != "allTrans" and transId != selTransId:
            continue
        line = [" "] * seqLen
        mouseOvers = {}  # position -> mouseOver-text or null, for end of mouse over
        # codingExNum = 0
        for exIdx, (
            exNum,
            exStart,
            exEnd,
            exFrame,
            oldExFrame,
            nextFrame,
            exStrand,
        ) in enumerate(exRows):
            if exFrame == -1:
                for i in range(exStart, exEnd):
                    line[i] = "="
                exonLabel = "noncoding"
                if (exEnd - exStart) > len(exonLabel) + 4:
                    # center the exon label on the exon
                    mid = exStart + int((exEnd - exStart) * 0.5)
                    halfLen = int(len(exonLabel) * 0.5)
                    labStart = mid - halfLen
                    for i in range(0, len(exonLabel)):
                        line[labStart + i] = exonLabel[i]
                    line[labStart - 1] = " "
                    line[labStart + len(exonLabel)] = " "
            else:
                # codingExNum += 1
                if koMethod == "stop" and editData is not None:

                    # add splice sites to the positions to edit
                    spliceBases = getSpliceSites(seq, exStart, exEnd)
                    # print(spliceBases)
                    for spliceBase, pos in spliceBases:
                        isStop, newStopGuides = checkStopCodons(spliceBase, pos, editData, possibleSpliceEdits, stopGuides, splice=True)
                        # print(stopGuides)
                        if len(newStopGuides) > 0:
                            stopGuides.update(newStopGuides)

                if transId == "manual annotation":
                    exonDesc = (
                        "manually annotated coding sequence<br>start phase %s"
                        % oldExFrame
                    )
                else:
                    exonDesc = (
                        "gene %s<br>transcript %s<br>exon %d<br>start phase %s"
                        % (
                            symbol,
                            transId,
                            exNum + 1,
                            oldExFrame,
                        )
                    )
                    if nextFrame is not None:
                        exonDesc += "<br>end phase %s" % nextFrame
                        # oldExFrame is the phase of the whole exon (before trimming to the sequence window)
                        # if (oldExFrame + nextFrame) % 3 == 0:
                        if oldExFrame == nextFrame:
                            exonDesc += (
                                "<br>Removing the exon retains the reading frame"
                            )
                        else:
                            exonDesc += (
                                "<br>Removing the exon will destroy the reading frame"
                            )

                mouseOvers[exStart] = exonDesc
                mouseOvers[exEnd] = None

                # exFrame is the number of bases already accounted for at the end of the previous exon
                exOffset = (3 - exFrame) % 3

                for i in range(exStart, exStart + exOffset):
                    line[i] = "-"

                # iterate in reverse for exons on the opposite strand
                if exStrand == "-":
                    for i in reversed(range(exStart, exEnd - exOffset, 3)):
                        codon = revComp(seq[i: i + 3])
                        if len(codon) == 3:
                            shortAa = codonTable[codon]
                            longAa = "[[" + shortAa
                        else:
                            longAa = "-"
                        for j in range(0, len(longAa)):
                            line[i + j] = longAa[j]
                else:
                    for i in range(exStart + exOffset, exEnd, 3):
                        codon = seq[i: i + 3]

                        if len(codon) == 3:
                            shortAa = codonTable[codon]
                            longAa = shortAa + "]]"
                        else:
                            # codon is split by splice site
                            longAa = "-"
                        for j in range(0, len(longAa)):
                            # in knock-out mode, highlight methionines in green
                            if (
                                koMethod in ["frameshift", "stop"]
                                and codon.upper() == "ATG"
                            ):
                                line[i + j] = (
                                    """<span style="background-color: rgba(0, 255, 0, 0.5); height: 15px; display: inline-block;">%s</span>"""
                                    % longAa[j]
                                )

                            # in base editing mode, these codons can be transformed into STOP codons
                            elif (
                                koMethod == "stop"
                                and codon.upper() in possibleStops
                                and (editData is not None or loadStopGuides is not None)
                            ):
                                isStop = False
                                if loadStopGuides is not None:

                                    # use pre-computed stop guides
                                    for pos, ez in stopGuides.values():
                                        if pos in range(i, i + 3):
                                            isStop = True
                                else:
                                    isStop, newStopGuides = checkStopCodons(
                                        codon.upper(),
                                        i,
                                        editData,
                                        possibleStops,
                                        stopGuides,
                                    )

                                    if len(newStopGuides) > 0:
                                        stopGuides.update(newStopGuides)

                                if isStop:
                                    line[i + j] = (
                                        """<span style="background-color: rgba(255, 51, 0, 0.5); height: 15px; display: inline-block;">%s</span>"""
                                        % longAa[j]
                                    )
                                else:
                                    line[i + j] = longAa[j]
                            else:
                                line[i + j] = longAa[j]

        # now merge the mouse overs as span tags into the ASCII line
        newLine = []
        for pos, char in enumerate(line):
            if pos in mouseOvers:
                overString = mouseOvers[pos]
                if overString is not None:
                    newLine.append(
                        "<span class='tooltipsterInteract' title='%s'>" % overString
                    )
                if char == "-" and selTransId != "manual annotation":
                    newLine.append(
                        "<span class='tooltipsterInteract' title='This codon goes over a splice site. The nucleotides of the split codon are not translated to amino acids but shown as dashes.'>-</span>"
                    )
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
        newLineStr = newLineStr.replace("]", "<lc>&gt;</lc>")

        lines.append((symbol, transId, newLineStr))
        # maxLabelLen = max(maxLabelLen, len(transId))

    if koMethod is not None:
        return lines, stopGuides
    else:
        return lines


def enzymeCoversMutPos(enzyme, mutPos):
    """returns True if at least one model of the given base editor enzyme has an
    editing window that covers mutPos (the edit position relative to the guide).
    Used to reject stop guides whose edit falls outside every model window of the
    chosen enzyme: the selection editData is built with the union window of all
    enzymes (e.g. CBE's DeepNG-BE_Ss 2-12), which is wider than e.g. the CGBE
    windows (all 4-9)."""
    return any(m["win"][0] <= mutPos < m["win"][1] for m in allBeModels[enzyme])


def checkStopCodons(feature, featurePos, editData, possibleEdits, stopGuides, splice=False):
    """
    returns True if a guide in editData can be used to change the codon into a STOP codon
    also returns the list of guides for which isStop is True
    if splice is False, feature is a codon and possibleEdits is a dict of codons that can be changed to a STOP
    if splice is True, feature is a splice site and possibleEdits is a dict of possible substitutions
    """

    # need to do this before calling showExonAndPams and showGuideTable to filter possible guides early
    # but this function needs all the PAMs to be computed..

    isStop = False
    newStopGuides = {}

    for seqIdx, editDict in editData.items():
        # convert to int when editData is loaded from json
        seqIdx = int(seqIdx)
        for editPos, editData in editDict.items():
            editPos = int(editPos)
            # editBase = str(list(editData.keys())[0])
            # print(feature, featurePos, editPos, possibleEdits.keys(), "<br>")
            if splice and editPos == featurePos and feature in possibleEdits.keys():
                # print("SPLICE!")
                enzyme = possibleEdits[feature]
                newStopGuides = {
                    tpl[0]: (editPos, enzyme)
                    for tpl in list(editData.values())[0]
                    if enzymeCoversMutPos(enzyme, tpl[3])
                }
                isStop = len(newStopGuides) > 0

            # may add a double check : mutated codon = revCodonTable[*]
            elif not splice and editPos == featurePos + possibleEdits[feature][0]:
                enzyme = possibleEdits[feature][1]
                newStopGuides = {
                    tpl[0]: (editPos, enzyme)
                    for tpl in list(editData.values())[0]
                    if enzymeCoversMutPos(enzyme, tpl[3])
                }
                isStop = len(newStopGuides) > 0

            # print(editPos, editData, "<br>")

    return isStop, newStopGuides


def getGeneModels(org):
    "read possible gene models for org and return as list (name, desc) or None if no gene models"
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
        name = baseName.split(".")[0]
        desc = geneDescs.get(baseName, name)
        ret.append((name, desc))
    return ret


def getSelGeneModel(org, noGenes=True, manual=False, noAllTrans=False):
    "return (list of (name, desc) of models, selected gene model name)"
    geneModels = getGeneModels(org)
    selGeneModel = None
    selTransId = None

    if manual:
        if geneModels:
            geneModels.append(("manual", "manual annotation"))
        else:
            geneModels = [("manual", "manual annotation")]

    if geneModels:
        if noGenes:
            selGeneModel = cgiParams.get("geneModelSelection", geneModels[0][0])
            geneModels.insert(0, ("noGenes", "Do not show"))
        else:
            selGeneModel = cgiParams.get("geneModelSelection", geneModels[0][0])
            if selGeneModel == "None":
                selGeneModel = geneModels[0][0]
        possNames = [x for x, y in geneModels]
        if selGeneModel not in possNames:
            errAbort(
                "The gene model name specified with the argument geneModel is invalid"
            )

        selTransId = cgiParams.get("selTransId", None)

        # in the donor design page, allTrans needs to be returned as None
        if noAllTrans and selTransId in ["allTrans", "None"]:
            selTransId = None

    return geneModels, selGeneModel, selTransId


def resolveAnnotationParams(org, seq, position):
    """Return a dict of annotation params that correspond to
    with the currently selected geneModel. Stale keys are dropped."""
    selModel = cgiParams.get("geneModelSelection")
    if selModel in (None, "None", "noGenes"):
        return {}

    if selModel == "manual":
        return {
            "manualExStart": cgiParams.get("manualExStart", "1"),
            "manualExEnd": cgiParams.get("manualExEnd", str(len(seq))),
            "manualExFrame": cgiParams.get("manualExFrame", "0"),
        }

    # real gene model: validate selTransId against this model's transcripts
    exonInfo, _ = getExonInfo(org, selModel, position)
    validIds = {"allTrans"} | {tId for tId, _ in exonInfo.keys()}
    transId = cgiParams.get("selTransId")
    if transId not in validIds:
        transId = "allTrans"
    return {"geneModelSelection": selModel, "selTransId": transId}


def printSeqForCopy(seq):
    "print a hidden text area so we can copy the sequence to the clipboard"
    print('<input id="seqAsText" type="text" style="display:none">')
    print(seq)
    print("</input>")


def calcBeScoresServer(seq, guideSeq, pamSeq, pamId, extGuideStart, extGuideEnd, pamStrand, enzyme, extSeq=None, pam=None):
    """
    sends the data to sub-servers running each Base edting models
    enzyme can be CBE, ABE or CGBE
    """

    effs = []
    outcomes = []

    # DeepBE models : for now get scores for all available SpCas9 models
    deepBeSubmodels = {
            "ABE": ["DeepNG-BE_17m", "DeepNG-BE_8e"],
            "CBE": ["DeepNG-BE_Ss", "DeepNG-BE_YE1"],
            "CGBE": ["DeepNG-BE_mini", "DeepNG-BE_CGBE1", "DeepNG-BE_Bi"]
            }

    forecastSubmodels = {
            "ABE": ["ABE"],
            "CBE": ["CBE"]
            }

    # get the PAM variant model for the current PAM (1 = NGG)
    # pam is used in KO mode
    if pam:
        pamPat = pam
    elif "." in pamId:
        pamPat = pamId.split('.')[0]
    else:
        pamPat = "NG"

    pamVariant = pamVariantModels.get(pamPat, 1)

    selDeepBeModels = deepBeSubmodels[enzyme]
    selForecastModels = forecastSubmodels.get(enzyme)
    models = {"ForecastBe": selForecastModels, "DeepBe": selDeepBeModels}
    
    logging.info("MODELS : %s" % models)

    # no CGBE model in FORECasT
    if enzyme == "CGBE":
        del models["ForecastBe"]

    # if the sequence can't be extended use extSeq to extend it
    # if extSeq is not available, add Ts to get a 30bp sequence as a last resort (can't add Ns)
    # or don't calculate the scores ?
    ext3, ext5 = ("", "")
    if extGuideStart < 0:
        if extSeq:
            ext5 = extSeq[100 + extGuideStart: 100]
        else:
            ext5 = ''.join(["T" for i in range(abs(extGuideStart))])
        extGuideStart = 0

    if extGuideEnd > len(seq):
        if extSeq:
            ext3 = extSeq[100 + len(seq): 100 + len(seq) + (extGuideEnd - len(seq))]
        else:
            ext3 = ''.join(["C" for i in range(extGuideEnd - len(seq))])
        extGuideEnd = len(seq)

    extGuideSeq = ext5 + seq[extGuideStart:extGuideEnd] + ext3

    # print(N5, N3, extGuideSeq, len(extGuideSeq), "<br>")

    if pamStrand == "-":
        extGuideSeq = revComp(extGuideSeq)

    # [enzyme, submodel (DeepBE only), input sequence]
    inData = [enzyme, None, pamVariant, extGuideSeq.upper()]

    for model, submodels in models.items():

        if (
                len(pamSeq) != 3 or
                (len(pamSeq) == 3 and pamSeq.upper()[1:3] != "GG")
                ) and model == "ForecastBe":
            # no models for alternative PAMs with FORECasT-BE
            continue

        if submodels is None:

            modelStr = model
            modelOut = callSubServer("run%s" % model, inData)

            if modelOut.get("status") == "processed":

                outcome = modelOut["outcome"]

                if pamStrand == "-":
                    newOutcome = []
                    for outSeq, outEff in modelOut["outcome"]:
                        revOutSeq = revComp(outSeq)
                        newOutcome.append([revOutSeq, outEff])
                    outcome = newOutcome

                outcomes.append((model, outcome))
                effs.append((model, modelOut["eff"]))

        else:
            for submodel in submodels:
                modelStr = "%s - %s" % (model, submodel)
                inData[1] = submodel
                modelOut = callSubServer("run%s" % model, inData)

                if modelOut.get("status") == "processed":
                    # invert lower / uppercase
                    # should do this in runDeepBE.py, but it crashes ?

                    outcome = []
                    for outSeq, outEff in modelOut["outcome"]:
                        if model == "DeepBe":
                            newOutSeq = "".join(
                                [
                                    base.upper() if base.islower() else base.lower()
                                    for base in outSeq
                                ]
                            )
                        else:
                            newOutSeq = outSeq
                        if pamStrand == "-":
                            newOutSeq = revComp(newOutSeq)

                        outcome.append([newOutSeq, outEff])
                    outcomes.append((modelStr, outcome))
                    effs.append((modelStr, modelOut["eff"]))

        # debug message
        if modelOut.get("error") is not None:
            print("""<p>%s crashed with the following error : %s. <br>
            Here is the traceback : %s<br>""" % (model, modelOut.get("error"), modelOut.get("trace")))

    return effs, outcomes


def filterEditInfoToSubst(eiDict, substPamIds):
    """Given an editInfo dict (nucl -> list of guide tuples), keep only the
    guides whose pamId is in substPamIds and drop any nucl whose list ends up
    empty. Bystander edits sharing a substituting pamId are preserved.
    Returns the filtered dict, or None if nothing remains."""
    filtered = {}
    for nucl, guideData in eiDict.items():
        keep = [tpl for tpl in guideData if tpl[0] in substPamIds]
        if keep:
            filtered[nucl] = keep
    return filtered or None


def makeEditLines(
    seq,
    pamSeqs,
    winStart,
    winEnd,
    guideScores,
    exonId=0,
    stopGuides=None,
    substInfo=None,
    batchId=None,
    enzyme="CBE",
    extSeq=None,
    loadJson=False,
    pam=None
):
    """
    Create the lines that show the possible baseEditor edits.
    Can be called to generate JSON data of potential edits,
    or draw the edit lines from existing JSON data (if loadJson is True).
    """

    editInfos = []
    for i in range(0, len(seq)):
        editInfos.append(defaultdict(list))

    # list of pam IDs that result in a substitution in KI mode
    if substInfo is not None:
        insertIdx, insertSeq = substInfo
        substPamIds = []

    upSeq = seq.upper()

    # placeholder
    for pamId, pamStart, guideStart, strand, guideSeq, pamSeq, pamPlusSeq in pamSeqs:

        if loadJson:
            break

        # in knock-out mode, discard guides that don't introduce a STOP codon
        # in knock-in mode, discard substitutions that can't be introduced with BE
        if (
                (stopGuides is not None and pamId not in stopGuides.keys())
                or (substInfo is not None and enzyme is None)
                ):

            continue

        if substInfo is not None and enzyme is None:
            continue

        # should remove this
        if guideScores is not None:
            specScore = guideScores[pamId]
        else:
            specScore = 0

        if strand == "+":
            fromPos = guideStart + winStart
            toPos = guideStart + winEnd
            # dict of possible edits {enzyme: (fromNucl, toNucl)}
            editsTo = {ez: ezList[0] for ez, ezList in possibleEdits.items()}

        else:
            guideEnd = guideStart + GUIDELEN
            fromPos = guideEnd - winEnd
            toPos = guideEnd - winStart
            editsTo = {ez: ezList[1] for ez, ezList in possibleEdits.items()}

        fromNucl, toNucl = editsTo[enzyme]

        for pos in range(fromPos, toPos):
            if pos >= len(seq) or pos < 0:
                continue

            doScore = False

            # only keep pamIds that introduce the desired substitution in KI mode
            if (
                substInfo is not None
                and pos == insertIdx
                and toNucl == insertSeq.upper()
            ):
                substPamIds.append(pamId)
                doScore = True

            # position of mutated nucl on forw strand guide
            if strand == "+":
                mutPos = pos - guideStart
                # extended guide position, to be used by DeepBaseEditor and CRISPRonBE
                # 30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp)
                extGuideStart = guideStart - 4
                extGuideEnd = pamStart + len(pamSeq) + 3
            else:
                mutPos = GUIDELEN - (pos - guideStart) - 1
                extGuideStart = pamStart - 3
                extGuideEnd = guideStart + GUIDELEN + 4

            if upSeq[pos] == fromNucl:

                # in KO mode, calculate the scores only in the second call
                # Should move this outside of the loop and pass the list of guides as input
                # for now, the speed is OK
                if doScore or stopGuides is not None:

                    if stopGuides is not None:
                        # this only works because CBE and CGBE uses the same fromNucl
                        # and STOPs can only be introduced by these enzymes
                        enzyme = stopGuides[pamId][1]
                        toNucl = editsTo[enzyme][1]

                    effs, outcomes = calcBeScoresServer(
                        seq, guideSeq, pamSeq, pamId, extGuideStart, extGuideEnd, strand, enzyme, extSeq=extSeq, pam=pam
                    )

                    # list of models that can be used to mutate this position
                    if stopGuides is not None:
                        stopEnzymes = ["CBE", "CGBE"]
                        possibleModels = {
                            "%s - %s" % (m["tool"], m["model"])
                            for ez in stopEnzymes for m in allBeModels[ez]
                            if m["win"][0] <= mutPos < m["win"][1]
                        }

                    else:
                        possibleModels = {
                            "%s - %s" % (m["tool"], m["model"])
                            for m in allBeModels[enzyme]
                            if m["win"][0] <= mutPos < m["win"][1]
                        }

                    # discard the models that can't result in an edit at this position
                    effs = [(modelStr, eff) for modelStr, eff in effs if modelStr in possibleModels]
                    outcomes = [(modelStr, outcome) for modelStr, outcome in outcomes if modelStr in possibleModels]

                else:
                    effs, outcomes = 0, []

                editInfos[pos][toNucl].append(
                    (pamId, guideSeq, pamSeq, mutPos, effs, outcomes, specScore)
                )

    # if stopGuides is not None or substInfo is not None:
    #     closeBeModels(models)

    # print(editInfos)
    altNucls = ["A", "T", "G", "C"]

    editLabels = []
    for an in altNucls:
        editLabels.append("Edits to " + an)

    editLines = []
    for i in range(0, len(altNucls)):
        editLines.append([" "] * len(seq))

    # rearrange into lines of text + JSON
    if loadJson:
        # the function was colled from the results page : load jsonData from its file
        batchBase = join(batchDir, batchId)
        editFname = batchBase + ".editData.json"
        jsonData = json.load(open(editFname))
        editItems = jsonData.get(str(exonId), {}).items()
    else:
        jsonData = defaultdict(dict)
        editItems = enumerate(editInfos)

    for pos, eiDict in editItems:
        pos = int(pos)
        if not eiDict:
            continue

        # for substitutions, keep only PAMs that result in a substitution
        # (bystander edits sharing the same pamId are preserved).
        # in loadJson mode the stored data is already filtered.
        if substInfo is not None and not loadJson:
            eiDict = filterEditInfoToSubst(eiDict, substPamIds)
            if eiDict is None:
                continue

        # in knock-out mode, separate the data from different sequences
        if not loadJson:
            jsonData[exonId][pos] = dict(eiDict)
        for nucl, guideData in eiDict.items():

            stopPamIds = [editList[0] for editList in guideData]

            # print(pos, stopGuides, "<br>")
            if stopGuides:
                # positions that can introduce a knock-out for in the current exon
                stopGuidePos = []
                for pamId, editTpl in stopGuides.items():
                    stopPos, ez = editTpl
                    if exonId != int(pamId.split('.')[0]) or pamId not in stopPamIds:
                        continue
                    stopGuidePos.append(editTpl[0])
            else:
                stopGuidePos = []

            # print(stopGuidePos, exonId, "<br>")

            if stopGuides is not None and pos in stopGuidePos:
                style = "color: rgb(255, 51, 0); font-weight: bold;"
            elif substInfo is not None and pos == insertIdx:
                style = "color: rgb(255, 127, 4); font-weight: bold;"
            else:
                style = "color: rgb(255, 127, 4);"

            mOverLink = "<d pos=%d exonId=%d>%s</d>" % (pos, exonId, nucl)
            # bystander edits in grey
            if (
                stopGuides is not None
                and pos not in stopGuidePos
                or substInfo is not None
                and pos != insertIdx
            ):
                style = "color: rgb(102, 102, 102)"
                mOverLink = "<d>%s</d>" % nucl

            pamIdLink = guideData[0][0]
            yPos = altNucls.index(nucl)
            editLines[yPos][pos] = (
                """<a href="#%s" name="editBase" style="%s" onclick="showBeTable('beTable')">%s</a>"""
                % (pamIdLink, style, mOverLink)
            )

    ret = []
    for label, lineChars in zip(editLabels, editLines):
        # don't display a line with no edits
        if len([char for char in lineChars if char != " "]) > 0:
            ret.append((label, None, "".join(lineChars)))

    return ret, jsonData


def makeEditGuideLines(editData):
    "draws base editing outcomes on the sequence viewer"
    pamEditData = buildEditData(editData)

    editGuideLines = []

    editSpan = "<span style='background-color: rgba(255, 255, 0, 0.35)'>"
    spanEnd = "</span>"

    for pamId, editInfo in pamEditData.items():
        pamInfo = pamId.split('.')[1]
        pamPos = int(pamInfo[1:-1])
        pamStrand = pamInfo[-1]

        for editTpl in editInfo:
            _, _, effs, allOutcomes = editTpl

            for model, outcomes in allOutcomes:
                if pamStrand == "+":
                    # substract the length of the sequence until PAM + 5' extension
                    spaceRange = range(pamPos - (24 - 4))
                else:
                    spaceRange = range(pamPos - 3)
                space = ''.join([" " for i in spaceRange])

                editGuideLines.append((model, None, ""))
                for i, (outcome, freq) in enumerate(outcomes):
                    freq = round(100 * freq, 2)
                    if freq == 0 or i > 5:
                        continue
                    if pamStrand == "+":
                        outcome = outcome[4:27]
                    else:
                        outcome = outcome[3:26]

                    guide = ''.join([editSpan + base + spanEnd if base.isupper() else base for base in outcome])
                    editGuideLines.append(("freq : %s %%" % freq, None, space + guide))

    return editGuideLines


def makePamLines(
    lines,
    maxY,
    pamIdToSeq,
    guideScores,
    linkToTable=True,
    pamIdToSuppInfo=None,
    pamWindow=None,
    otherPam=None,
    insertIdx=None,
    bePamIds=None
):

    if bePamIds:
        # the function is printed in printKiSteps()
        tableLinkFunc = """onclick="showBeTable('beTable')" """
    elif pamWindow is not None and otherPam is None:
        # KI / HDR PAMs
        tableLinkFunc = """onclick="showBeTable('hdrTable')" """
    else:
        tableLinkFunc = ""

    for y in range(0, maxY + 1):
        texts = []
        lastEnd = 0
        for start, end, name, strand, pamId in lines[y]:

            if bePamIds and pamId not in bePamIds:
                continue

            if otherPam is None:
                guideSeq = pamIdToSeq.get(pamId)
            else:
                guideSeq = None
                # get a simpler version of insertDistance for supplementary PAMs to be filtered in KI mode
                insertDistance = abs(start - insertIdx)
                if pamWindow and insertDistance > pamWindow:
                    continue
            if pamIdToSuppInfo:
                pamIdInfo = pamIdToSuppInfo.get(pamId)
                if pamIdInfo:
                    doRecoding = pamIdInfo[0]
                    insertDistance = pamIdInfo[1]
                    if pamWindow and insertDistance > pamWindow:
                        continue
                    if doRecoding is True:
                        recodeStr = ""
                    else:
                        recodeStr = """, 0 0 10px #33ccff,
                        0 0 20px #33ccff,
                        0 0 30px #33ccff"""
                else:
                    recodeStr = ""
            else:
                recodeStr = ""
            if guideSeq is None and otherPam is None:
                # when there is an N in the guide, the PAM is valid, but the guide is not
                continue

            if otherPam is None:
                classStr = cssClassesFromSeq(guideSeq, suffix="Seq")
                score = guideScores[pamId]
                mouseOver = ""
            else:
                classStr = "tooltipsterInteract"
                mouseOver = """
                title='The scores for the guide sequence corresponding to this PAM are not calculated yet. Click on "New search with the selected PAMs" to show their scores.'
                """
                score = ""

            spacer = "".join([" "] * ((start - lastEnd)))
            lastEnd = end
            texts.append(spacer)
            # XX How can this happen for non-Cpf1 enzymes? Can this ever happen?
            if score is None and not pamIsFirst and otherPam is None:
                continue
            if otherPam is None:
                color = scoreToColor(score)[0]
            else:
                color = "#666666"

            if linkToTable and otherPam is None:
                texts.append(
                    """<a class='%s' %s style="text-shadow: 1px 1px 1px #bbb%s; color: %s " id="list%s" href="#%s" %s>"""
                    % (classStr, mouseOver, recodeStr, color, pamId, pamId, tableLinkFunc)
                )
            else:
                texts.append(
                    """<span class='%s' %s style="text-shadow: 1px 1px 1px #bbb%s; color: %s" id="list%s">"""
                    % (classStr, mouseOver, recodeStr, color, pamId)
                )

            texts.append(name)
            if linkToTable:
                texts.append("</a>")
            else:
                texts.append("</span>")
        yield ("", None, "".join(texts))


def getBeWin(winVal):
    "return (start, end) of base editor window given CGI variable"
    fs = winVal.split("-")
    if len(fs) != 2:
        errAbort("parameter beWin must contain only one dash")
    start = fs[0].strip()
    end = fs[1].strip()
    if not start.isdigit() or not end.isdigit():
        errAbort("parameter beWin must be two dash-separated numbers")
    start = int(start)
    end = int(end)

    return start, end


def printLines(lines, labelLen):
    "print list of (label, string) such that label is at least labelLen characters long"
    for label, mouseOver, line in lines:
        if mouseOver is not None:
            labelStr = (
                '<span class="tooltipsterInteract" title="{:s}">{:'
                + str(labelLen)
                + "s} </span>"
            ).format(label, mouseOver)
        else:
            labelStr = ("{:" + str(labelLen) + "s} ").format(label)

        print((labelStr), end=" ")
        print(line)


def getMaxLen(lines):
    "given a list of tuples where first element is the label, return the longest label len"
    maxLen = 0
    for l in lines:
        label = l[0]
        maxLen = max(maxLen, len(label))
    return maxLen


def printJson(name, obj):

    print("<script %s>")
    print((name), end=" ")
    print(("="), end=" ")
    print((json.dumps(obj)))
    print("</script>")


def showExonAndPams(
    batchId,
    org,
    seq,
    startDict,
    pam,
    guideScores,
    varHtmls,
    varDbs,
    varDb,
    minFreq,
    position,
    pamIdToSeq,
    exonId,
    koMethod,
    browserlink,
    selGeneModel=None,
    selTransId=None,
    exonSelect=None,
    stopGuides=None,
    allEditData=None
):
    pamSeqs = list(flankSeqIter(seq, startDict, len(pam), True, exonId=exonId))
    if koMethod == "splicing":
        if exonSelect.isnumeric():
            originalExon = int(exonSelect)
        if exonId % 2 == 0:
            spliceType = "acceptor"
            if not exonSelect.isnumeric():
                originalExon = (exonId + 1) // 2
        else:
            spliceType = "donor"
            if not exonSelect.isnumeric():
                originalExon = exonId // 2

    if koMethod == "stop":
        extendPos = GUIDELEN + 6
    else:
        extendPos = GUIDELEN - 6

    # control the display of splicing donor / acceptor sites relative to their exons (0-based id)
    if koMethod == "splicing":
        htmlExonId = originalExon
    else:
        htmlExonId = exonId

    if koMethod in ["frameshift"]:
        if exonId == 0:
            exonDisplay = "block"
        else:
            exonDisplay = "none"
    else:
        exonDisplay = "block"

    _, _, _, strand = parsePos(position)
    # don't display the exons where no PAMs were found

    if len(pamSeqs) == 0:
        if koMethod == "splicing":
            if exonSelect.isnumeric():
                exonNumberText = int(exonSelect) + 1
            else:
                exonNumberText = originalExon + 1

            print(
                """
                <p style="display: %s" class="exonGroup%s" name="exonDisplay" > No guide sequences were found in the splicing %s site of exon %s (%s)</p>
                """
                % (exonDisplay, htmlExonId, spliceType, exonNumberText, browserlink)
            )
        else:
            print(
                """
                <p style="display: %s" class="exonGroup%s" name="exonDisplay" > No guide sequences were found at position %s</p>
                 """
                % (exonDisplay, htmlExonId, browserlink)
            )

    # for the region upstream of the TSS, scroll to the left by default
    if koMethod == "excision" and exonId == 0:
        print(
            """
        <script>
        window.addEventListener('DOMContentLoaded', function() {
          const div = document.getElementById('exonPamSeq%s');
            div.scrollLeft = div.scrollWidth - div.clientWidth;
            });
        </script> """
            % exonId
        )
    lines, maxY = distrOnLines(seq.upper(), startDict, len(pam), pam, exonId)
    posLabel = "Position"
    seqLabel = "Sequence"
    varLabel = "Variants"

    # pamLines is empty
    # print(lines, maxY, pamIdToSeq, guideScores)
    pamLines = list(makePamLines(lines, maxY, pamIdToSeq, guideScores))
    labelLen = max(len(seqLabel), len(posLabel), len(varLabel), getMaxLen(pamLines))

    if koMethod == "splicing":
        if exonId % 2 == 0:
            spliceLabel = "Splicing acceptor site"
        else:
            spliceLabel = "Splicing donor site"
        labelLen = max(labelLen, len(spliceLabel))

    exonLines = []

    # only set when base editing is active (stop mode); guard for other methods
    editInfo = None

    # stopGuides is the batch-wide dict (all exons); narrow it to the guides
    # whose PAM lies in the current exon for per-exon counts and display
    exonStopGuides = {}

    if baseEditor and koMethod == "stop" and stopGuides is not None:
        exonPamIds = {tpl[0] for tpl in pamSeqs}
        exonStopGuides = {
            pamId: val for pamId, val in stopGuides.items() if pamId in exonPamIds
        }
        beWinStart, beWinEnd = getBeWin(cgiParams.get("beWin", DEFAULTBEWIN))
        editLines, jsonData = makeEditLines(
            seq,
            pamSeqs,
            beWinStart,
            beWinEnd,
            guideScores,
            exonId,
            stopGuides=stopGuides,
            batchId=batchId,
            loadJson=True
        )
        lines, maxY = distrOnLines(
            seq.upper(), startDict, len(pam), pam, exonId, stopGuides=stopGuides
        )

        editInfo = (pam, jsonData)
        pamLines = list(makePamLines(lines, maxY, pamIdToSeq, guideScores))
        labelLen = max(labelLen, getMaxLen(editLines))

    if selGeneModel is not None:
        exonInfo, maxTransIdLen = getExonInfo(
            org, selGeneModel, position, extendPos=extendPos
        )

        labelLen = max(labelLen, maxTransIdLen)
        exonLines, _ = makeExonLines(
            exonInfo, seq, selTransId, koMethod, editInfo=editInfo, loadStopGuides=exonStopGuides
        )

    # if the last exon on the first third of the coding sequence was extended 14bp inside a coding sequence,
    # avoid flagging this region as an intron
    if selGeneModel and selTransId:
        for exonDesc, exonCoords in exonInfo.items():
            transId, sym = exonDesc
            if transId != selTransId:
                continue
            for exFrame, exStart, exEnd, exFrame, _, _, _ in exonCoords:
                if exFrame >= 0 and (len(seq) - exEnd) < extendPos:
                    seq = seq[:exStart] + seq[exStart:exEnd].upper() + seq[exEnd:]
                    break

    if baseEditor or varDb:
        print(
            (
                """<form style="display:inline" id="paramForm" action="%s" method="GET">"""
                % basename(__file__)
            )
        )

    print(
        """ <div name="exonDisplay" class="exonGroup%s" style="display:%s;"> """
        % (htmlExonId, exonDisplay)
    )
    print(""" <div class="substep" """)
    print('<a id="seqStart"></a>')
    if koMethod in ["frameshift", "stop"]:
        exonLen = len("".join(base for base in seq if base.isupper()))

        # only count the guides that can introduce a stop codon with base editing
        if koMethod == "stop":
            guidesCount = len(exonStopGuides)
        else:
            guidesCount = len(guideScores)

        exonNumberText = exonId + 1

        print(
            "The target coding exon %d (%s) which is %d bp long, and %d bp flanking intron or UTR sequences are shown. They contain %d possible guide sequences.<br>"
            % (exonNumberText, browserlink, exonLen, extendPos, guidesCount)
        )
    elif koMethod == "excision":
        if exonId == 0:
            print(
                " the region %s bp upstream of the TSS (%s) contains %d possible guide sequences.<br>"
                % (len(seq), browserlink, len(guideScores))
            )
        else:
            print(
                " the region %s bp downstream of the TES (%s) contains %d possible guide sequences.<br>"
                % (len(seq), browserlink, len(guideScores))
            )

    elif koMethod == "promoter":
        if exonId == 0:
            print(
                " the region upstream (%s) of the selected promoter contains %d possible guide sequences.<br>"
                % (browserlink, len(guideScores))
            )
        else:
            print(
                " the region downstream (%s) of the selected promoter contains %d possible guide sequences.<br>"
                % (browserlink, len(guideScores))
            )

    elif koMethod == "splicing":
        exonNumberText = originalExon + 1
        if exonId % 2 == 0:
            print(
                " the region around the splicing acceptor site of exon %s (%s) contains %d possible guide sequences.<br>"
                % (exonNumberText, browserlink, len(guideScores))
            )
        else:
            print(
                " the region around the splicing donor site of exon %s (%s) contains %d possible guide sequences.<br>"
                % (exonNumberText, browserlink, len(guideScores))
            )

    print("</div>")

    print(
        """<div class="blueHighlight" name="exonPamSeq" id="exonPamSeq%s" style="text-align: left; overflow-x:scroll; width:98%%; height:7vw; background:#DDDDDD; border-style: solid; border-width: 1px">"""
        % exonId
    )

    print(
        """<pre style="font-family: Source Code Pro; font-size: 80%; display:inline; line-height: 0.95em; text-align:left">"""
    )
    print(("{:" + str(labelLen) + "s} ").format(posLabel), end=" ")
    print(rulerString(len(seq)))

    if varHtmls is not None:
        print(("{:" + str(labelLen) + "s} ").format(varLabel), end=" ")
        print("".join(varHtmls))

    if koMethod == "splicing":
        spliceGap = 17
        print(("{:" + str(labelLen) + "s} ").format(spliceLabel), end=" ")
        # splicing donor site
        if exonId % 2 == 0 and strand == "+" or exonId % 2 != 0 and strand == "-":
            print(
                "".join([" " for i in range(spliceGap)]),
                "".join(["-" for i in range(6)]),
            )
        else:
            print(
                "".join([" " for i in range(spliceGap + 6)]),
                "".join(["-" for i in range(6)]),
            )

    print(("{:" + str(labelLen) + "s} ").format(seqLabel), end=" ")
    # don't display intronic sequences
    # maskedSeq = ''.join([base if base.isupper() else "." for base in seq])
    print(seq)

    printLines(exonLines, labelLen)

    if baseEditor:
        printLines(editLines, labelLen)

    printLines(pamLines, labelLen)

    print("</pre><br>")

    print("""</div>""")
    print("""</div>""")


def showSeqAndPams(
    org,
    seq,
    startDict,
    pam,
    guideScores,
    varHtmls,
    varDbs,
    varDb,
    minFreq,
    position,
    pamIdToSeq,
    browserLinkHtml=None,
    multiPamInfo=None,
    linkToTable=True,
    pamIdToSuppInfo=None,
    pamWindow=None,
    otherPam=None,
    strand=None,
    clippedSeq=None,
    geneId=None,
    useBaseEditor=False,
    extSeq=None,
    editData=None,
    batchId=None,
    noPerfectMatch=None,
):
    "show the sequence and the PAM sites underneath in a sequence viewer"

    global baseEditor

    if multiPamInfo is None:
        pamSeqs = list(flankSeqIter(seq, startDict, len(pam), True))
        lines, maxY = distrOnLines(seq.upper(), startDict, len(pam), pam)
        pamLines = list(makePamLines(lines, maxY, pamIdToSeq, guideScores, linkToTable))
        substInfo = None
        # kiType = None
        # need to get the enzyme from the PAM
        # if baseEditor:
        #    enzyme = "CBE"

    else:
        multipam = multiPamInfo[0]
        # pamList is the list of all PAMs to be used (HDR + BE)
        pamList = multiPamInfo[1]
        insertIdx = multiPamInfo[2]
        kiType = multiPamInfo[3]
        insertSeq = multiPamInfo[4]
        pamSeqs = []
        allPamLines = []

        # used in makeEditLines
        if kiType == "substitution" and useBaseEditor:
            editTpl = (seq[insertIdx].upper(), insertSeq.upper())

            substInfo = (insertIdx, insertSeq)
            enzyme = None

            for ez, editList in possibleEdits.items():
                fw, rev = editList
                if editTpl == fw or editTpl == rev:
                    enzyme = ez
        else:

            substInfo = None
            enzyme = None

        # HDR PAMs : use the list of PAMs selected from the menu
        for pamFullName in multiPamDict[multipam][0]:

            pam = setupPamInfo(pamFullName)
            startDict, endSet = findAllPams(seq.upper(), pam)
            singlePamSeqs = list(
                flankSeqIter(
                    seq, startDict, len(pam), True, exonId=None, pamFullName=pamFullName
                )
            )

            pamSeqs.extend(singlePamSeqs)

            # get lines for the current pam
            currentPamLines = getPamLines(
                seq.upper(), startDict, len(pam), pam, pamFullName=pamFullName
            )
            allPamLines.extend(currentPamLines)

        # layout all lines
        lines, maxY = layoutPamLines(allPamLines, len(seq))
        pamLines = list(
            makePamLines(
                lines,
                maxY,
                pamIdToSeq,
                guideScores,
                linkToTable,
                pamIdToSuppInfo=pamIdToSuppInfo,
                pamWindow=pamWindow,
                otherPam=None,
            )
        )

        # in KI mode, show BE guides separately
        if editData:

            bePamIds = set(buildEditData(editData).keys())
            bePamSeqs = []
            allBePamLines = []

            for pamFullName in pamList:

                pam = setupPamInfo(pamFullName)
                startDict, endSet = findAllPams(seq.upper(), pam)
                beSinglePamSeqs = list(
                    flankSeqIter(
                        seq, startDict, len(pam), True, exonId=None, pamFullName=pamFullName
                    )
                )

                bePamSeqs.extend(beSinglePamSeqs)

                # get lines for the current pam
                beCurrentPamLines = getPamLines(
                    seq.upper(), startDict, len(pam), pam, pamFullName=pamFullName
                )
                allBePamLines.extend(beCurrentPamLines)

            if useBaseEditor:
                # set the global here
                baseEditor = True

            # layout all lines
            beLines, maxY = layoutPamLines(allBePamLines, len(seq))
            bePamLines = list(
                makePamLines(
                    beLines,
                    maxY,
                    pamIdToSeq,
                    guideScores,
                    linkToTable,
                    pamIdToSuppInfo=pamIdToSuppInfo,
                    pamWindow=pamWindow,
                    otherPam=None,
                    bePamIds=bePamIds
                )
            )
            # discard empty pam lines
            bePamLines = [pamTpl for pamTpl in bePamLines if len(pamTpl[2]) > 0]

        # in KI mode, display other PAMs on the sequence if the option was selected
        if otherPam:
            otherPamList = multiPamDict[otherPam][0]
            pamFullName = setupPamInfo(otherPam)
            otherPamLines = []
            otherPamSeqs = []
            for pamFullName in otherPamList:
                pam = setupPamInfo(pamFullName)
                startDict, endSet = findAllPams(seq.upper(), pam)

                singlePamSeqs = list(
                    flankSeqIter(
                        seq,
                        startDict,
                        len(pam),
                        True,
                        exonId=None,
                        pamFullName=pamFullName,
                    )
                )
                otherPamSeqs.extend(singlePamSeqs)
                # get lines for the current pam
                currentPamLines = getPamLines(
                    seq.upper(), startDict, len(pam), pam, pamFullName=pamFullName
                )
                otherPamLines.extend(currentPamLines)

            otherLines, maxY = layoutPamLines(otherPamLines, len(seq))
            otherPamLinesHtml = list(
                makePamLines(
                    otherLines,
                    maxY,
                    None,
                    None,
                    pamWindow=pamWindow,
                    otherPam=otherPam,
                    insertIdx=insertIdx,
                )
            )

    posLabel = "Position"
    varLabel = "Variants"
    if multiPamInfo:
        if kiType == "insertion":
            insertLabel = "Insertion site"
        else:
            insertLabel = kiType
    else:
        insertLabel = ""
    if otherPam:
        otherPamLabel = "Added PAM motifs"
    else:
        otherPamLabel = ""

    seqLabel = "Sequence"
    exonLabelLen = 0
    editLines = []
    exonLines = []

    # selGeneModel = None
    # geneModels = None
    # pre-select the geneId
    if geneId and cgiParams.get("geneModelSelection") is None:
        geneModels, selGeneModel, selTransId = getSelGeneModel(org, manual=True)
        found = False
        for model, modelStr in geneModels:
            exonInfo, _ = getExonInfo(org, model, position)
            for tId, sym in list(exonInfo.keys()):
                if tId == geneId:
                    selGeneModel = model
                    selTransId = tId
                    found = True
                    break
            if found:
                break

    geneModels, selGeneModel, selTransId = getSelGeneModel(org, manual=True)

    # in only show manual annotation to avoid showing wrong sequences
    if noPerfectMatch:
        geneModels = [("manual", "manual annotation")]

    if baseEditor:
        # beWinStart, beWinEnd = getBeWin(cgiParams.get("beWin", DEFAULTBEWIN))
        enzList = allBeModels[enzyme]
        beWinStart = min([enzDict["win"][0] for enzDict in enzList])
        beWinEnd = max([enzDict["win"][1] for enzDict in enzList])

        editLines, jsonData = makeEditLines(
            seq, pamSeqs, beWinStart, beWinEnd, guideScores, substInfo=substInfo, enzyme=enzyme, extSeq=extSeq, batchId=batchId, loadJson=True
        )

        # to show base editing outcomes on the sequence viewer (unused for now)
        # editGuideLines = makeEditGuideLines(editData)

        printJson("editData", jsonData)

    labelLen = max(
        len(varLabel),
        len(insertLabel),
        len(seqLabel),
        len(posLabel),
        len(otherPamLabel),
        getMaxLen(pamLines),
    )

    if selGeneModel is not None:
        exonInfo, maxTransIdLen = getExonInfo(org, selGeneModel, position)
        labelLen = max(labelLen, maxTransIdLen)
    if baseEditor:
        editDetailsLabel = "Base editing outcomes"
        labelLen = max(labelLen, getMaxLen(editLines), len(editDetailsLabel))
    if selGeneModel:
        labelLen = max(labelLen, exonLabelLen)

    if multiPamInfo is not None:
        print("""<details id="results2" open style="margin-bottom: 12px;">""")
        print(
            """<summary style="font-weight: bold; font-size: 20px; margin-top: 24px; margin-bottom: 12px;">Legend and general information about the results</summary>"""
        )
    else:
        print("<div class='substep' style='margin-top: 24px;'>")
        print('<a id="seqStart"></a>')

    print(
        "Your input sequence is %d bp long. It contains %d possible guide sequences.<br>"
        % (len(seq), len(guideScores))
    )

    if multiPamInfo:

        if clippedSeq is True:
            print(
                """<p>The input sequence was trimmed to 100bp on each side of the edit site to save computation time.<br>
                  Guides outside this range are too distant to result in a successful knock-in experiment.</p>"""
            )

        if len(pamList) == 1:
            print(
                """Shown below are ther PAM sites and the expected cleavage position located -3bp 5' of the PAM site, for PAM %s.<br>
                """
                % "".join(pamList)
            )
        else:
            print(
                """Shown below are ther PAM sites and the expected cleavage position, for the following PAMs : %s.<br>
                """
                % ", ".join(pamList)
            )
        print(
            "PAMs highlighted in blue corresponds to guides that overlap the position of the edit.<br>"
        )
        print(
            """Colors <span style="color:#32cd32; text-shadow: 1px 1px 1px #bbb">green</span>, <span style="color:#ffff00; text-shadow: 1px 1px 1px #888">yellow</span> and <span style="text-shadow: 1px 1px 1px #f01; color:#aa0014">red</span> indicate high, medium and low specificity of the PAM's guide sequence in the genome.<p>"""
        )

    elif not pamIsFirst:
        print(
            "Shown below are their PAM sites and the expected cleavage position located -3bp 5' of the PAM site.<br>"
        )
        print(
            """Colors <span style="color:#32cd32; text-shadow: 1px 1px 1px #bbb">green</span>, <span style="color:#ffff00; text-shadow: 1px 1px 1px #888">yellow</span> and <span style="text-shadow: 1px 1px 1px #f01; color:#aa0014">red</span> indicate high, medium and low specificity of the PAM's guide sequence in the genome.<p>"""
        )
        print(
            "Click on a match for the PAM %s below to show its %d bp-long guide sequence.<br>"
            % (pam, GUIDELEN)
        )

    else:
        print(
            "Click on a match for the PAM %s below to show its %d bp-long guide sequence.<br>"
            % (pam, GUIDELEN)
        )

    if multiPamInfo:
        print("</details>")

    print(
        "(Need help? Look at the <a target=_blank href='manual/#annotseq'>CRISPOR manual</a>)<br>"
    )

    if baseEditor or varDb or selGeneModel or multiPamInfo:
        print(
            (
                """<form style="display:inline" id="paramForm" action="%s" method="GET">"""
                % basename(__file__)
            )
        )

    if geneModels:
        # selGeneModel/selTransId may come from defaults computed above
        # (geneId match, first model for multiPamInfo) rather than from a
        # user form submit, in which case cgiParams doesn't yet hold them.
        # Seed cgiParams so resolveAnnotationParams sees the same selection
        # the form is about to render.
        if selGeneModel and selGeneModel != "noGenes":
            cgiParams["geneModelSelection"] = selGeneModel
        if selTransId:
            cgiParams["selTransId"] = selTransId

        # drop stale annotation params from a previous submit and replace them
        # with a set that is consistent with the currently selected gene model
        annotParams = resolveAnnotationParams(org, seq, position)
        for key in ("selTransId", "manualExStart", "manualExEnd", "manualExFrame"):
            cgiParams.pop(key, None)
        cgiParams.update(annotParams)
        if selGeneModel == "noGenes":
            cgiParams.pop("geneModelSelection", None)

        print("""<details id="results3" open autocomplete="off">""")
        print(
            """<summary style="font-weight: bold; font-size: 20px; margin-top: 24px; margin-bottom: 12px;">Annotation of the coding sequence</summary>"""
        )
        print(
            """<div>Select a gene model and a transcript below to display the translated coding sequences.<br>
              If there is no available annotation or if it is incomplete, you can manually annotate the coding sequence by selecting "manual annotation".<br>"""
        )
        if multiPamInfo:
            print(
                "Once you selected a transcript, it will be used as a model to recode the donor DNA, if needed.<br>"
            )
        print("Gene Models:")
        printDropDown(
            "geneModelSelection", geneModels, selGeneModel, style="width:20em"
        )

        if selGeneModel == "manual":
            manualExStart = annotParams["manualExStart"]
            if not manualExStart.isdigit():
                manualExStart = 0
            else:
                manualExStart = int(manualExStart) if int(manualExStart) > 0 else 0

            manualExEnd = annotParams["manualExEnd"]
            if not manualExEnd.isdigit():
                manualExEnd = len(seq)
            else:
                manualExEnd = (
                    int(manualExEnd) if int(manualExEnd) < len(seq) else len(seq)
                )

            manualExFrame = int(annotParams["manualExFrame"])

            # write clamped values back so the donor-DNA URL uses the same values shown in the form
            cgiParams["manualExStart"] = str(manualExStart)
            cgiParams["manualExEnd"] = str(manualExEnd)
            cgiParams["manualExFrame"] = str(manualExFrame)

            print(
                """
                coding start:&nbsp<input name="manualExStart" value=%s style="width: 30px;"/>&nbsp
                coding end:&nbsp<input name="manualExEnd" value=%s style="width: 30px;"/>&nbsp
                <span class="tooltipsterInteract" title="the reading frame is relative to the coding start position">reading frame</span>:&nbsp<select name="manualExFrame">&nbsp
            """
                % (manualExStart, manualExEnd)
            )

            # exFrame is the number of bases already present in the codon (not the shift of the phase)
            for exFrame in [0, 2, 1]:
                if exFrame == 2:
                    frameStr = "1"
                elif exFrame == 1:
                    frameStr = "2"
                else:
                    frameStr = str(exFrame)

                if exFrame == manualExFrame:
                    selected = "selected"
                else:
                    selected = ""
                print(
                    """<option %s value=%d>%s</option>"""
                    % (selected, exFrame, frameStr)
                )
            print("</select>&nbsp")

            exonInfo = {
                ("manual annotation", "manual annotation"): [
                    (
                        1,
                        manualExStart,
                        manualExEnd,
                        manualExFrame,
                        manualExFrame,
                        0,
                        "+",
                    )
                ]
            }

            selTransId = "manual annotation"
            labelLen = max(labelLen, len(selTransId))

        if selGeneModel not in ["noGenes", "manual"]:
            selTransId = annotParams["selTransId"]
            print("Transcript:")
            # only one transcript per symbol in refSeq select
            if selGeneModel != "refSeqSelect":
                transIdInfo = [("allTrans", "All Transcripts")]
            else:
                transIdInfo = []
            for transId, sym in list(exonInfo.keys()):
                transIdInfo.append((transId, sym + " / " + transId))
            printDropDown("selTransId", transIdInfo, selTransId, style="width:20em")

        if selGeneModel != "noGenes":
            exonLines = makeExonLines(exonInfo, seq, selTransId)

        print(
            """<input style="height:18px;margin:0px;font-size:10px;line-height:normal" type="submit" name="submit" value="Update"><br>"""
        )
        if selTransId:
            if "ENST" in selTransId:
                print(
                    """Link to ENSEMBL : <a target="blank" href="https://www.ensembl.org/Multi/Search/Results?q=%s;site=ensembl;page=1">%s</a><br>"""
                    % (selTransId.split("_")[0], selTransId)
                )
            elif "NM" in selTransId or "NR" in selTransId:
                print(
                    """Link to NCBI : <a target="blank"  href="https://www.ncbi.nlm.nih.gov/nuccore/%s/">%s</a><br>"""
                    % (selTransId, selTransId)
                )

        print("""</div><br>""")

    elif multiPamInfo:
        for key in (
            "manualExStart",
            "manualExEnd",
            "manualExFrame",
            "geneModelSelection",
            "selTransId",
        ):
            cgiParams.pop(key, None)
    print("</details>")

    if baseEditor and enzyme is not None:
        print("""<details id="results4" open autocomplete="off">""")
        print(
            """<summary style="font-weight: bold; font-size: 20px; margin-top: 24px; margin-bottom: 12px;">Base editing information</summary>"""
        )
        if multiPamInfo is not None:
            print(
                """<p>Base editing can be used to introduce this substitution.<br>
                  Using base editing bypasses the need for a double strand break and a donor DNA, which may be useful for therapeutic applications.</p>"""
            )
        print(
            """<p>Show below the sequence are the possible edits, using this base editor with the selected modification window.<br>"""
        )
        if multiPamInfo is not None:
            print(
                """
                  <ul>
                        <li>Edits in orange correspond to the substitution.</li>
                        <li>Edits in grey correspond to "bystander" edits (i.e, additional edits that can occur when using the same guide).</li>
                  </ul>
                  """
            )

        print(
            """Hover on an edit to show the corresponding guides and their predicted efficiencies and outcome sequences, or click on it to navigate to the table.<br>"""
        )

        # input te change the base editor modification window
        # unsued now, editing data is pre-calculated by the workers.
        # may use this to show the potential edits in case of new enzymes with extended windows ?

        '''
        print("Base Editor modification window:")
        selBeWin = "%s-%s" % (beWinStart, beWinEnd)
        print(("""<input type="text" name="beWin" size="10" value="%s">""" % selBeWin))
        print(
            """<input style="height:18px;margin:0px;font-size:10px;line-height:normal" type="submit" name="submit" value="Update">"""
        )
        '''

        print("</p>")
        print("</details>")

    if varDb is not None:
        print("Variant database:")
        varDbList = [(b, c) for a, b, c, d in varDbs]  # only keep fname+label
        printDropDown("varDb", varDbList, varDb)

        if minFreq == 0.0:
            minFreq = "0.0"
        else:
            minFreq = str(minFreq)

        # pull out the hasAF field for this varDb
        varDbHasAF = False
        for shortLabel, fname, desc, hasAF in varDbs:
            if fname == varDb:
                varDbHasAF = hasAF
                break

        if varDbHasAF:
            print("""&nbsp; Min. frequency: """)
            print(
                ("""<input type="text" name="minFreq" size="8" value="%s">""" % minFreq)
            )
        print(
            """<input style="height:18px;margin:0px;font-size:10px;line-height:normal" type="submit" name="submit" value="Update">"""
        )
        print(
            (
                "<small style='margin-left:30px'><a href='mailto:%s'>Missing a variant database? We can add it.</a></small>"
                % contactEmail
            )
        )

    if position == "?":
        print(
            "<small style=''>Input sequence not in genome, cannot show genome variants.</small>"
        )
    elif varDb is None:
        print(
            (
                "<small style=''><a href='mailto:%s'>Suggest a genome variants database to show on this page</a></small>"
                % contactEmail
            )
        )

    if multiPamInfo is None:
        print(
            "</div>"
        )  # close <div class='substep'> opened on the multiPamInfo-is-None branch above; without this the rest of the page nests inside .substep and `div.substep table td { border:none }` strips borders from the otTable body
    if browserLinkHtml:
        print(browserLinkHtml)
    print(
        """<div class="blueHighlight" style="text-align: left; overflow-x:scroll; width:98%; background:#e6e6e6; border-style: solid; border-width: 1px">"""
    )

    print(
        """<pre style="font-family: Source Code Pro; font-size: 80%; display:inline; line-height: 0.95em; text-align:left">"""
    )
    print(("{:" + str(labelLen) + "s} ").format(posLabel), end=" ")
    print(rulerString(len(seq)))

    if multiPamInfo:
        # rangeChar was used to highlight a n bp region up/dowstream of the edit site
        if kiType == "deletion":
            delRange = "".join(["x" for i in range(len(insertSeq))])
            insertChar = "\%s/" % delRange
            leftSpace = "".join([" " for i in range(insertIdx - 1)])
        elif kiType in ["substitution", "replacement"]:
            insertChar = insertSeq
            leftSpace = "".join([" " for i in range(insertIdx)])
        else:
            insertChar = "\/"
            leftSpace = "".join([" " for i in range(insertIdx - 1)])

        print(("{:" + str(labelLen) + "s} ").format(insertLabel), end=" ")
        print(leftSpace + insertChar)

    if varHtmls is not None:
        print(("{:" + str(labelLen) + "s} ").format(varLabel), end=" ")
        print("".join(varHtmls))

    print(("{:" + str(labelLen) + "s} ").format(seqLabel), end=" ")
    print(seq)

    printLines(exonLines, labelLen)

    if baseEditor:
        printLines(editLines, labelLen)
        # print("<details><summary>%s</summary>" % editDetailsLabel)
        # printLines(editGuideLines, labelLen)
        # print("</details>")

        print("<details open><summary>Show/hide HDR PAMs</summary>")

    printLines(pamLines, labelLen)

    if baseEditor:
        print("</details>")
        print("<details open><summary>Show/hide BE PAMs</summary>")
        printLines(bePamLines, labelLen)
        print("</details>")

    if otherPam:
        print(("{:" + str(labelLen) + "s} ").format(otherPamLabel), end=" ")
        print("".join(["- - " for i in range(0, int(len(seq) / 4))]))
        printLines(otherPamLinesHtml, labelLen)

    print("</pre><br>")

    print("""</div>""")

    # printSeqForCopy(seq)
    if pamIsCas12max(pam):
        print(
            '<div style="line-height: 1.0; padding-top: 5px; font-size: 15px">Cpf1 has a staggered site: cleavage occurs between the 14th and 16th base on the non-targeted strand (indicated by "\\" in the schema above). Cleavage mostly occurs after the 24rd base on the targeted strand (indicated by "/" in the schema above). See on <a target=_blank href="https://www.synthego.com/products/nuclease/hfcas12max-hifi">Synthego</a></div>'
        )
    elif pamIsCpf1(pam):
        print(
            '<div style="line-height: 1.0; padding-top: 5px; font-size: 15px">Cpf1 has a staggered site: cleavage occurs usually - but not always - after the 18th base on the non-targeted strand which has the TTTV PAM motif (indicated by "\\" in the schema above). Cleavage mostly occurs after the 23rd base on the targeted strand which has the AAAN motif (indicated by "/" in the schema above). See <a target=_blank href="http://www.sciencedirect.com/science/article/pii/S0092867415012003">Zetsche et al 2015</a>, in particular <a target=_blank href="http://www.sciencedirect.com/science?_ob=MiamiCaptionURL&_method=retrieve&_eid=1-s2.0-S0092867415012003&_image=1-s2.0-S0092867415012003-gr3.jpg&_cid=272196&_explode=defaultEXP_LIST&_idxType=defaultREF_WORK_INDEX_TYPE&_alpha=defaultALPHA&_ba=&_rdoc=1&_fmt=FULL&_issn=00928674&_pii=S0092867415012003&md5=11771263f3e390e444320cacbcfae323">Fig 3</a>.</div>'
        )
    elif pamIsCasX(pam):
        print(
            '<div style="line-height: 1.0; padding-top: 5px; font-size: 15px">We have no description yet on how exactly the CasX cleavage looks like. Please contact crispor@tefor.net if you have an idea how to describe the cleavage site.</div>'
        )


def iterOneDelSeqs(seq):
    """given a seq, create versions with each bp removed. Avoid duplicates
    yields (delPos, seq)
    >>> list(iterOneDelSeqs("AATGG"))
    [(0, 'ATGG'), (2, 'AAGG'), (3, 'AATG')]
    """
    doneSeqs = set()
    for i in range(0, len(seq)):
        delSeq = seq[:i] + seq[i + 1 :]
        if delSeq not in doneSeqs:
            yield i, delSeq
        doneSeqs.add(delSeq)


def flankSeqIter(seq, startDict, pamLen, doFilterNs, exonId=None, pamFullName=None):
    """given a seq and dictionary of pamPos -> strand and the length of the pamSite
    yield tuples of (name, pamStart, guideStart, strand, flankSeq, pamSeq)

    flankSeq is the guide sequence (=flanking the PAM).

    if doFilterNs is set, will not return any sequences that contain an N character
    pamPlusSeq are the 5bp after the PAM. If not enough space, pamPlusSeq is None
    """

    startList = sorted(startDict.keys())
    for pamStart in startList:
        strand = startDict[pamStart]

        pamPlusSeq = None
        if pamIsFirst:  # Cpf1: get the sequence to the right of the PAM
            if strand == "+":
                guideStart = pamStart + pamLen
                flankSeq = seq[guideStart : guideStart + GUIDELEN]
                pamSeq = seq[pamStart : pamStart + pamLen]
                if pamStart - pamPlusLen >= 0:
                    pamPlusSeq = seq[pamStart - pamPlusLen : pamStart]
            else:  # strand is minus
                guideStart = pamStart - GUIDELEN
                flankSeq = revComp(seq[guideStart:pamStart])
                pamSeq = revComp(seq[pamStart : pamStart + pamLen])
                if pamStart + pamLen + pamPlusLen < len(seq):
                    pamPlusSeq = revComp(
                        seq[pamStart + pamLen : pamStart + pamLen + pamPlusLen]
                    )
        else:  # common case: get the sequence on the left side of the PAM
            if strand == "+":
                guideStart = pamStart - GUIDELEN
                flankSeq = seq[guideStart:pamStart]
                pamSeq = seq[pamStart : pamStart + pamLen]
                if pamStart + pamLen + pamPlusLen < len(seq):
                    pamPlusSeq = seq[pamStart + pamLen : pamStart + pamLen + pamPlusLen]
            else:  # strand is minus
                guideStart = pamStart + pamLen
                flankSeq = revComp(seq[guideStart : guideStart + GUIDELEN])
                pamSeq = revComp(seq[pamStart : pamStart + pamLen])
                if pamStart - pamPlusLen >= 0:
                    pamPlusSeq = revComp(seq[pamStart - pamPlusLen : pamStart])

        if "N" in flankSeq and doFilterNs:
            continue
        if exonId is not None:
            pamId = "%d.s%d%s" % (exonId, pamStart, strand)
        elif pamFullName is not None:
            pamId = "%s.s%d%s" % (
                pamFullName,
                pamStart,
                strand,
            )  # use parameter "pamPrefix" for both instead
        else:
            pamId = "s%d%s" % (pamStart, strand)

        yield pamId, pamStart, guideStart, strand, flankSeq, pamSeq, pamPlusSeq


def makeBrowserLink(dbInfo, pos, text, title, cssClasses, ctUrl=None, returnUrl=False):
    "return link to genome browser (ucsc or ensembl) at pos, with given text"
    if dbInfo is None:
        errAbort(
            "Your batchID relates to a genome that is not present anymore. You will have to change the version of the site. Or contact us and send us the full URL of this page."
        )

    if dbInfo.server.startswith("Ensembl"):
        baseUrl = "www.ensembl.org"
        urlLabel = "Ensembl"

        # link back to archive, if possible
        if dbInfo.description.startswith("Ensembl "):
            ensVersion = dbInfo.description.split()[1]
            if ensVersion.isdigit():
                baseUrl = "e%s.ensembl.org" % ensVersion

        elif dbInfo.server == "EnsemblPlants":
            baseUrl = "plants.ensembl.org"
        elif dbInfo.server == "EnsemblMetazoa":
            baseUrl = "metazoa.ensembl.org"
        elif dbInfo.server == "EnsemblProtists":
            baseUrl = "protists.ensembl.org"
        org = dbInfo.scientificName.replace(" ", "_")
        pos = pos.replace(":+", "").replace(":-", "")  # remove the strand
        url = "http://%s/%s/Location/View?r=%s" % (baseUrl, org, pos)
    elif (
        dbInfo.server == "ucsc"
        or dbInfo.name.startswith("GCA_")
        or dbInfo.name.startswith("GCF_")
    ):
        urlLabel = "UCSC"
        if len(pos) > 0 and pos[0].isdigit():
            pos = "chr" + pos
        # remove the strand
        pos = pos.replace(":+", "").replace(":-", "")
        url = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s" % (
            dbInfo.name,
            pos,
        )
        if ctUrl is not None:
            url += "&hgt.customText=%s" % ctUrl
    # some limited support for gbrowse
    elif dbInfo.server.startswith("http://"):
        urlLabel = "GBrowse"
        chrom, start, end, strand = parsePos(pos)
        start = start + 1
        url = "%s/?name=%s:%d..%d" % (dbInfo.server, chrom, start, end)
    else:
        chrom, start, end, strand = parsePos(pos)
        if chrom is not None and chrom.startswith("NC_"):
            start = start + 1
            url = (
                "https://www.ncbi.nlm.nih.gov/nuccore/%s?report=graph&log$=seqview&v=%d-%d"
                % (chrom, start, end)
            )
            urlLabel = "NCBI "

        else:
            # return "unknown genome browser server %s, please email services@tefor.net" % dbInfo.server
            urlLabel = None
            url = "javascript:void(0)"

    classStr = ""
    if len(cssClasses) != 0:
        classStr = ' class="%s"' % (" ".join(cssClasses))

    if title is None:
        if urlLabel != None:
            title = "Link to %s Genome Browser" % urlLabel
        else:
            title = "No Genome Browser link available yet for this organism"

    if returnUrl:
        return url
    else:
        return """<a title="%s"%s target="_blank" href="%s">%s</a>""" % (
            title,
            classStr,
            url,
            text,
        )


def highlightMismatches(guide, offTarget, pamLen):
    "return a string that marks mismatches between guide and offtarget with *"
    if pamLen != 0:
        if pamIsFirst:
            offTarget = offTarget[pamLen:]
        else:
            offTarget = offTarget[:-pamLen]
    assert len(guide) == len(offTarget)

    s = []
    for x, y in zip(guide, offTarget):
        if x == y:
            s.append(".")
        else:
            s.append("*")
    return "".join(s)


def parseNewAlias(ifh):
    "part of parseAlias(): IGV-compatible format: first is UCSC, all other columns are aliases"
    toUcsc = {}
    for line in ifh:
        if line.startswith("#"):
            continue
        row = line.rstrip("\n").split("\t")
        for i in range(1, len(row)):
            toUcsc[row[i]] = row[0]
    return toUcsc


def parseAlias(fname):
    "parse tsv file with at least two columns, orig chrom name and new chrom name. copied from chromToUcsc script from the UCSC tools."
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
    "if chrom is in chromAlias, return the human-readable name"
    global chromAlias
    if chromAlias == -1:  # == chromAlias file not present
        return chrom
    elif chromAlias is None:
        chromAliasFname = join("genomes", db, db + ".chromAlias.txt")
        if not isfile(chromAliasFname):
            chromAlias = -1
            return chrom
        else:
            chromAlias = parseAlias(chromAliasFname)

    return chromAlias.get(chrom, chrom)


def makeAlnStr(org, seq1, seq2, pam, mitScore, cfdScore, posStr, chromDist):
    """given two strings of equal length, return a html-formatted string of several lines
    that show the two sequences and a line that highlights where they differ
    """
    lines = [[], [], []]
    last12MmCount = 0
    inLinkage = False
    hlSeed = False

    if pamIsSpCas9(pam):
        hlSeed = True

    if pamIsFirst:
        lines[0].append("<i>" + seq1[: len(pam)] + "</i> ")
        lines[1].append("<i>" + seq2[: len(pam)] + "</i> ")
        lines[2].append("".join([" "] * (len(pam) + 1)))

    if pamIsFirst:
        guideStart = len(pam)
        guideEnd = len(seq1)
    else:
        guideStart = 0
        guideEnd = len(seq1) - len(pam)

    for i in range(guideStart, guideEnd):
        if hlSeed and i == 10:
            lines[1].append("<u>")

        if seq1[i] == seq2[i]:
            lines[0].append(seq1[i])
            lines[1].append(seq2[i])
            lines[2].append(" ")
        else:
            lines[0].append("<b>%s</b>" % seq1[i])

            lines[1].append("<b>%s</b>" % seq2[i])

            lines[2].append("*")
            if i > 7:
                last12MmCount += 1

        if hlSeed and i == guideEnd - 1:
            lines[1].append("</u>")

    if not pamIsFirst:
        lines[0].append(" <i>" + seq1[-len(pam) :] + "</i>")
        lines[1].append(" <i>" + seq2[-len(pam) :] + "</i>")
    lines = ["".join(l) for l in lines]

    chrom, chromPos, strand = posStr.split(":")
    chrom = applyChromAlias(org, chrom)
    posStr = ":".join((chrom, chromPos, strand))

    if len(posStr) > 1 and posStr[0].isdigit():
        posStr = "chr" + posStr

    htmlText1 = (
        "<small><pre>guide:      %s<br>off-target: %s<br>            %s</pre>"
        % (lines[0], lines[1], lines[2])
    )

    if pamIsCpf1(pam) or pamIsCasX(pam):
        htmlText2 = "Cpf1/CasX: No off-target scores available</small>"
    elif saCas9Mode:
        htmlText2 = "SaCas9 Tycko Score: %s" % mitScore
    else:
        if cfdScore == None:
            cfdStr = "Cannot calculate CFD score on non-ACTG characters"
        else:
            cfdStr = "%f" % cfdScore

        htmlText2 = (
            "CFD Off-target score: %s<br>MIT Off-target score: %.2f<br>Position: %s</small>"
            % (cfdStr, mitScore, posStr)
        )
        if chromDist != None and org != None:
            htmlText2 += "<br><small>Distance from target: %.3f Mbp</small>" % (
                float(chromDist) / 1000000.0
            )
            if org.startswith("mm") or org.startswith("hg") or org.startswith("rn"):
                if chromDist > 20000000:
                    htmlText2 += "<br><small>&gt;20Mbp = unlikely to be in linkage with target</small>"
                else:
                    htmlText2 += "<br><small>&lt;20Mbp= likely to be in linkage with "
                    "target! Even if no linkage: beware of chromosomal rearrangements "
                    "when using this guide!</small>"
                    inLinkage = True

    hasLast12Mm = last12MmCount > 0
    return htmlText1 + htmlText2, hasLast12Mm, inLinkage


def parsePos(text):
    """parse a string of format chr:start-end:strand and return a 4-tuple
    Strand defaults to + and end defaults to start+23
    """
    if text != None and len(text) != 0 and text != "?":
        fields = text.split(":")
        if len(fields) == 2:
            chrom, posRange = fields
            strand = "+"
        else:
            chrom, posRange, strand = fields
        posRange = posRange.replace(",", "")
        if "-" in posRange:
            start, end = posRange.split("-")
            start, end = int(start), int(end)
            if start > end:
                start, end = end, start
                strand = "-"
        else:
            # if the end position is not specified (as by default done by UCSC outlinks), use start+23
            start = int(posRange)
            end = start + 23
    else:
        chrom, start, end, strand = "", 0, 0, "+"
    return chrom, start, end, strand


def annotateOfftargets(org, countDict, guideSeq, pam, inputPos):
    """for a given guide sequence, return a list of tuples that
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
    repCount = 0  # if repCount for a guide is !=0, then the guide should not be used. repCount is then the number
    # of matches for the guide in the genome (not looking at the PAM)

    # to calculate the aggregate CFD score, the CFD of offtargets with up to n mimatches are summed
    # from https://doi.org/10.1016/j.xgen.2026.101190
    aggrThreshold = 1
    aggrCfdScores = []

    # for each edit distance, get the off targets and iterate over them
    foundOneOntarget = False
    isSaCas9 = pamIsSaCas9(pam)
    isCpf1 = pamIsCpf1(pam)

    for editDist in range(0, maxMMs + 1):
        # print countDict,"<p>"
        matches = countDict.get(editDist, [])
        # print otCounts,"<p>"
        last12MmOtCount = 0

        # create html and score for every offtarget
        otCount = 0
        for (
            chrom,
            start,
            end,
            otSeq,
            strand,
            segType,
            geneNameStr,
            totalAlnCount,
            isRep,
        ) in matches:
            # if repCount is > 0, then this means that the guide should not be used and we cannot
            # even get any off-targets
            if (totalAlnCount > MAXOCC) or (totalAlnCount > 1 and isRep):
                repCount = totalAlnCount  # any off-target with this condition will trigger the whole guide to be suppressed

            # skip on-targets
            if segType != "":
                segTypeDesc = segTypeConv[segType]
                geneDesc = segTypeDesc + ":" + geneNameStr
                geneDesc = geneDesc.replace("|", "-")
            else:
                geneDesc = geneNameStr

            # is this not an off-target but the on-target?
            # if we got a genome position, use it. Otherwise use a random off-target with 0MMs
            # as the on-target ("auto-ontarget" mode)
            if (
                editDist == 0
                and repCount == 0
                and (
                    (chrom == inChrom and start >= inStart and end <= inEnd)
                    or (inChrom == "" and foundOneOntarget == False)
                )
            ):
                foundOneOntarget = True
                ontargetDesc = geneDesc
                continue

            otCount += 1
            guideNoPam = guideSeq[: len(guideSeq) - len(pam)]
            otSeqNoPam = otSeq[: len(otSeq) - len(pam)]

            # debug
            # print(len(guideNoPam))
            # print(len(otSeqNoPam))

            if len(otSeqNoPam) == 19:
                otSeqNoPam = (
                    "A" + otSeqNoPam
                )  # should not change the score a lot, weight0 is very low
                guideNoPam = "A" + guideNoPam

            if isCpf1:
                # Cpf1 has no off-target scores yet
                mitScore = 0.0
                cfdScore = 0.0
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
                if pam == "NGG" and otSeq[-2:] != "GG":
                    mitScore = mitScore * 0.2

                # CFD score must include the PAM
                cfdScore = calcCfdScore(guideSeq, otSeq)

            mitOtScores.append(mitScore)
            if cfdScore != -1:
                cfdScores.append(cfdScore)

            if editDist <= aggrThreshold:
                aggrCfdScores.append(cfdScore)

            posStr = "%s:%d-%s:%s" % (chrom, start + 1, end, strand)
            if chrom == inChrom:
                dist = abs(start - inStart)
            else:
                dist = None

            parNum = isInPar(org, chrom, start, end)
            if parNum is not None:
                posStr += " PAR%s" % parNum

            alnHtml, hasLast12Mm, inLinkage = makeAlnStr(
                org, guideSeq, otSeq, pam, mitScore, cfdScore, posStr, dist
            )
            if not hasLast12Mm:
                last12MmOtCount += 1
            posList.append(
                (
                    otSeq,
                    mitScore,
                    cfdScore,
                    editDist,
                    posStr,
                    geneDesc,
                    alnHtml,
                    inLinkage,
                )
            )

        last12MmCounts.append(str(last12MmOtCount))
        # create a list of number of offtargets for this edit dist
        otCounts.append(str(otCount))

    # calculate the guide scores
    if pamIsCpf1(pam):
        guideScore = -1
        guideCfdScore = -1
    else:
        if repCount > 0:
            guideScore = 0
            guideCfdScore = 0
        else:
            guideScore = calcMitGuideScore(sum(mitOtScores))

            if doCfdFix:
                guideCfdScore = calcCfdGuideScore(sum(cfdScores))
            else:
                guideCfdScore = calcMitGuideScore(sum(cfdScores))

            # aggrCfdScore = round(sum(aggrCfdScores), 2)

    # obtain the off-target info: coordinates, descriptions and off-target counts
    if repCount > 0:
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

    return (
        posList,
        otDescStr,
        guideScore,
        guideCfdScore,
        # aggrCfdScore,
        last12DescStr,
        ontargetDesc,
        repCount,
    )


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
    if guideSeq != saGuide:
        sys.path.append("bin/src/pairwise-library-screen")
        import predictSingle

        saGuide = guideSeq
        saScorer = predictSingle.SaCas9Scorer(len(guideSeq))

    # to be compatible with the MIT score, has to be in the range 0-100
    # for the MIT aggregate guide specificity score
    return 100.0 * saScorer.calcScore(guideSeq, otSeq)


# MIT offtarget scoring, "Hsu score"

# aka Matrix "M"
hitScoreM = [
    0,
    0,
    0.014,
    0,
    0,
    0.395,
    0.317,
    0,
    0.389,
    0.079,
    0.445,
    0.508,
    0.613,
    0.851,
    0.732,
    0.828,
    0.615,
    0.804,
    0.685,
    0.583,
]


def calcHitScore(string1, string2):
    "see 'Scores of single hits' on http://crispr.mit.edu/about"
    # The Patrick Hsu weighting scheme
    # S. aureus requires 21bp long guides. We fudge by using only last 20bp
    matrixStart = 0
    maxDist = 19

    assert string1[0].isupper()
    assert len(string1) == len(string2)
    # for nmCas9 and a few others with longer guides, we limit ourselves to 20bp
    if len(string1) > 20:
        string1 = string1[-20:]
        string2 = string2[-20:]
    # for 19bp guides, we fudge a little, but first pos has no weight anyways
    elif len(string1) == 19:
        string1 = "A" + string1
        string2 = "A" + string2
    # for shorter guides, I'm not sure if this score makes sense anymore, we force things
    elif len(string1) < 19:
        matrixStart = 20 - len(string1)
        maxDist = len(string1) - 1

    assert len(string1) == len(string2)

    dists = []  # distances between mismatches, for part 2
    mmCount = 0  # number of mismatches, for part 3
    lastMmPos = None  # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(matrixStart, len(string1)):
        if string1[pos] != string2[pos]:
            mmCount += 1
            if lastMmPos != None:
                dists.append(pos - lastMmPos)
            score1 *= 1 - hitScoreM[pos]
            lastMmPos = pos
    # 2nd part of the score
    if mmCount < 2:  # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists) / len(dists)
        score2 = 1.0 / (((maxDist - avgDist) / float(maxDist)) * 4 + 1)
    # 3rd part of the score
    if mmCount == 0:  # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount**2)

    score = score1 * score2 * score3 * 100
    return score


def calcMitGuideScore(hitSum):
    """Sguide defined on http://crispr.mit.edu/about
    Input is the sum of all off-target hit scores. Returns the specificity of the guide.
    """
    score = 100 / (100 + hitSum)
    score = int(round(score * 100))
    return score


def calcCfdGuideScore(hitSum):
    "suggested by Nicholas Parkinson"
    norm_score = 100.0 * 100.0 / (hitSum)
    return norm_score


# === SOURCE CODE cfd-score-calculator.py provided by John Doench =====
# The CFD score is an improved specificity score


def get_mm_pam_scores():
    """ """
    import pickle

    dataDir = join(dirname(__file__), "CFD_Scoring")
    mm_scores = pickle.load(open(join(dataDir, "mismatch_score.pkl"), "rb"))
    pam_scores = pickle.load(open(join(dataDir, "pam_scores.pkl"), "rb"))
    return (mm_scores, pam_scores)


# Reverse complements a given string
def revcom(s):
    basecomp = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A"}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return "".join(letters)


# Calculates CFD score
def calc_cfd(wt, sg, pam):
    # mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace("T", "U")
    wt = wt.replace("T", "U")
    s_list = list(sg)
    wt_list = list(wt)
    for i, sl in enumerate(s_list):
        if wt_list[i] == sl:
            score *= 1
        else:
            key = "r" + wt_list[i] + ":d" + revcom(sl) + "," + str(i + 1)
            score *= mm_scores[key]
    score *= pam_scores[pam]
    return score


mm_scores, pam_scores = None, None


def calcCfdScore(guideSeq, otSeq):
    """based on source code provided by John Doench
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
        mm_scores, pam_scores = get_mm_pam_scores()
    wt = guideSeq.upper()
    off = otSeq.upper()
    m_wt = re.search("[^ATCG]", wt)
    m_off = re.search("[^ATCG]", off)
    if (m_wt is None) and (m_off is None):
        pam = off[-2:]
        sg = off[:20]
        cfd_score = calc_cfd(wt, sg, pam)
        if doCfdFix:
            cfd_score = cfd_score * 100.0
        return cfd_score

    return -1


# ==== END CFD score source provided by John Doench

# --- END OF SCORING ROUTINES


def getSizeFname(genome):
    "return name of chrom.sizes file"
    genomeDir = genomesDir  # make local
    sizeFname = "%(genomeDir)s/%(genome)s/%(genome)s.sizes" % locals()
    return sizeFname


def parseChromSizes(genome):
    "return chrom sizes as dict chrom -> size"
    sizeFname = getSizeFname(genome)
    ret = {}
    for line in open(sizeFname).read().splitlines():
        fields = line.split()
        chrom, size = fields[:2]
        ret[chrom] = int(size)
    return ret


def extendAndGetSeq(db, chrom, start, end, strand, oldSeq, flank=FLANKLEN, noPerfectMatch=None):
    """extend (start, end) by flank and get sequence for it using twoBitTwoFa.
    Return None if not possible to extend.
    #>>> extendAndGetSeq("hg19", "chr21", 10000000, 10000005, "+", flank=3)
    #'AAGGAATGTAG'
    """
    assert (
        "|" not in chrom
    )  # we are using | to split info in BED files. | is not allowed in the fasta
    chromSizes = parseChromSizes(db)
    maxEnd = chromSizes[chrom] + 1

    start -= flank
    end += flank
    if start < 0 or end > maxEnd:
        return None

    genomeDir = genomesDir
    twoBitFname = "%(genomeDir)s/%(db)s/%(db)s.2bit" % locals()
    progDir = binDir
    genome = db
    cmd = (
        "%(progDir)s/twoBitToFa %(genomeDir)s/%(genome)s/%(genome)s.2bit stdout -seq='%(chrom)s' -start=%(start)s -end=%(end)s"
        % locals()
    )
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, encoding="utf8")
    seqStr = proc.stdout.read()
    proc.wait()
    if proc.returncode != 0:
        errAbort("Could not run '%s'. Return code %s" % (cmd, str(proc.returncode)))
    faFile = StringIO(seqStr)
    seqs = parseFasta(faFile)
    assert len(seqs) == 1
    seq = list(seqs.values())[0].upper()

    if strand == "-":
        seq = revComp(seq)

    genomeSeq = seq[FLANKLEN: (FLANKLEN + len(oldSeq))].upper()

    if oldSeq.upper() != genomeSeq and noPerfectMatch is None:
        logging.warning(
            "Input sequence has SNPs compared to genome, not returning extended seq:"
        )
        logging.warning("- Input sequence:  %s" % oldSeq)
        logging.warning("- Genome sequence: %s" % genomeSeq)
        logging.warning(
            "- Diff String    : %s" % highlightMismatches(oldSeq, genomeSeq, 0)
        )
        return None
    # ? make sure that user annotations, like added Ns, are retained in the long sequence
    # fixedSeq = seq[:100]+oldSeq+seq[-100:]
    # assert(len(fixedSeq)==len(seq))
    return seq


def getExtSeq(
    seq, start, end, strand, extUpstream, extDownstream, extSeq=None, extFlank=FLANKLEN
):
    """extend (start,end) by extUpstream and extDownstream and return the subsequence
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
    assert start >= 0
    assert end <= len(seq)

    # extend
    if strand == "+":
        extStart, extEnd = start - extUpstream, end + extDownstream
    else:
        extStart, extEnd = start - extDownstream, end + extUpstream

    # check for out of bounds and get seq
    if extStart >= 0 and extEnd <= len(seq):
        logging.debug("using input seq, pos %d-%d" % (extStart, extEnd))
        subSeq = seq[extStart:extEnd]
    else:
        if extSeq is None or extSeq == "?":
            return None
        # lift to extSeq coords and get seq
        extStart += extFlank
        extEnd += extFlank
        assert extStart >= 0
        assert extEnd <= len(extSeq)
        subSeq = extSeq[extStart:extEnd]
        logging.debug("using extended seq, pos %d-%d" % (extStart, extEnd))

    if strand == "-":
        logging.debug("revcomp'ing result")
        subSeq = revComp(subSeq)

    # check that the extended sequence really contains the whole input seq
    # e.g. when user has added nucleotides to a otherwise matching sequence
    # if seq.upper() not in subSeq.upper():
    # debug("seq is not in extSeq")
    # subSeq = None

    logging.debug(
        "Got -%d/+%d-extended seq for (%d, %d, %s) = %s. Result: %s."
        % (extUpstream, extDownstream, start, end, strand, seq[start:end], subSeq)
    )
    return subSeq


def pamStartToGuideRange(startPos, strand, pamLen):
    """given a PAM start position and its strand, return the (start,end) of the guide.
    Coords can be negative or exceed the length of the input sequence.
    """
    if not pamIsFirst:
        if strand == "+":
            return (startPos - GUIDELEN, startPos)
        else:  # strand is minus
            return (startPos + pamLen, startPos + pamLen + GUIDELEN)
    else:
        if strand == "+":
            return (startPos + pamLen, startPos + pamLen + GUIDELEN)
        else:  # strand is minus
            return (startPos - GUIDELEN, startPos)


def htmlHelp(text):
    "show help text with tooltip or modal dialog"
    className = "tooltipster"
    if "href" in text:
        className = "tooltipsterInteract"

    print(
        """<img style="padding-bottom: 3px; height:1.1em; width:1.0em" src="%simage/info-small.png" class="help %s" title="%s" />"""
        % (HTMLPREFIX, className, text)
    )


def htmlWarn(text):
    "show help text with tooltip"
    print(
        """<img style="height:0.9em; width:0.8em; padding-bottom: 2px" src="%simage/warning-32.png" class="help tooltipster" title="%s" />"""
        % (HTMLPREFIX, text)
    )


def readRestrEnzymes():
    """parse restrSites.txt and
    return as dict length -> list of (name, suppliers, seq)"""
    fname = "restrSites.txt"
    enzList = {}
    for line in open(join(baseDir, fname)):
        if line.startswith("#"):
            continue
        seq, name, suppliers = line.rstrip("\n").rstrip("\r").split("\t")
        suppliers = tuple(suppliers.split(","))
        enzList.setdefault(len(seq), []).append((name, suppliers, seq))
    return enzList


def patMatch(seq, pat, notDegPos=None):
    """return true if pat matches seq, both have to be same length
    do not match degenerate codes at position notDegPos (0-based)
    """
    assert len(seq) == len(pat)
    for x in range(0, len(pat)):
        patChar = pat[x]
        nuc = seq[x]

        assert patChar in "MKYRACTGNWSDVBH"
        assert nuc in "MKYRACTGNWSDXVBH"

        if notDegPos is not None and x == notDegPos and patChar != nuc:
            return False

        if nuc == "X":
            return False

        if patChar == "N":
            continue

        if patChar == "H" and nuc in "ACT":
            continue
        if patChar == "D" and nuc in "AGT":
            continue
        if patChar == "B" and nuc in "CGT":
            continue
        if patChar == "V" and nuc in "ACG":
            continue

        if patChar == "W" and nuc in "AT":
            continue
        if patChar == "S" and nuc in "GC":
            continue
        if patChar == "M" and nuc in "AC":
            continue
        if patChar == "K" and nuc in "TG":
            continue
        if patChar == "R" and nuc in "AG":
            continue
        if patChar == "Y" and nuc in "CT":
            continue

        if patChar != nuc:
            return False

    return True


def findSite(seq, restrSite):
    """return the positions where restrSite matches seq
    seq can be longer than restrSite
    Do not allow degenerate characters to match at position len(restrSite) in seq
    """
    posList = []
    for i in range(0, len(seq) - len(restrSite) + 1):
        subseq = seq[i : i + len(restrSite)]

        # JP does not want any potential site to be suppressed
        # if i<len(restrSite):
        # isMatch = patMatch(subseq, restrSite, len(restrSite)-i-1)
        # else:
        # isMatch = patMatch(subseq, restrSite)
        isMatch = patMatch(subseq, restrSite)

        if isMatch:
            posList.append((i, i + len(restrSite)))
    return posList


def matchRestrEnz(allEnzymes, guideSeq, pamSeq, pamPlusSeq, pamPat):
    """return list of enzymes that overlap the -3 position in guideSeq
    returns dict (name, pattern, suppliers) -> list of matching positions
    """
    matches = defaultdict(set)

    if pamPlusSeq is None:
        pamPlusSeq = "XXXXX"  # make sure that we never match a restriction site outside the seq boundaries

    fullSeq = concatGuideAndPam(guideSeq, pamSeq, pamPlusSeq)
    # print guideSeq, pamSeq, pamPlusSeq, fullSeq, "<br>"

    for siteLen, sites in allEnzymes.items():
        if pamIsCpf1(pamPat):
            # most modified position: 4nt from the end
            # see http://www.nature.com/nbt/journal/v34/n8/full/nbt.3620.html
            # Figure 1
            startSeq = len(fullSeq) - 4 - pamPlusLen - (siteLen) + 1
        else:
            # most modified position for Cas9: 3bp from the end
            startSeq = len(fullSeq) - len(pamSeq) - 3 - pamPlusLen - (siteLen) + 1

        seq = fullSeq[startSeq:].upper()
        for name, suppliers, restrSite in sites:
            posList = findSite(seq, restrSite)
            if len(posList) != 0:
                liftOffset = startSeq
                posList = [(liftOffset + x, liftOffset + y) for x, y in posList]
                matches.setdefault((name, restrSite, suppliers), set()).update(posList)
    return matches


def calcGlobScore(guideSeq, pamSeq, MitScore, CfdScore, effs, GC, freeE, globEffScore):
    "Calculate a global score for a sgRNA based of Specificity, Efficientcy, GC content and free-energy"

    MitScaled = MitScore / 100
    # CfdScaled = CfdScore/100

    effScore = effs[globEffScore]

    # for Cas12a, get the selected production method, independent of the global effScore to use
    selGlobEffScore = cgiParams.get("globEffScore")

    if globEffScore == "rs3":
        effScore = (effScore + 200) / 400
    else:
        effScore = effScore / 100

    # no specificity scores available for Cas12a : use predicted efficiency only

    if globEffScore == "seqDeepCpf1" or globEffScore == "EVA":
        mainScore = 100 * effScore
    else:
        mainScore = 100 * (0.60 * MitScaled + 0.40 * effScore)
    # coefficients will need to be adjusted

    # penalties

    grafType = crisporEffScores.getGrafType(guideSeq)

    if grafType and globEffScore != "seqDeepCpf1":
        if grafType == "tt" and globEffScore == "rs3":
            mainScore -= 25

    if "TTTT" in guideSeq and (
        globEffScore == "rs3"
        or (globEffScore == "seqDeepCpf1" and selGlobEffScore == "rs3")
    ):
        mainScore -= 25

    if GC < 0.25 or GC > 0.75:
        mainScore -= 25

    if freeE < -3.6 and globEffScore != "EVA":
        if freeE < -6:
            mainScore -= 20
        else:
            mainScore -= 10

    return mainScore


def calcEVAscore(EVAlike, MIT):
    "add the MIT score weight to the EVA-like score"

    if MIT >= 75:
        MIT = 75

    EVAfull = EVAlike + 0.1784 * MIT
    return EVAfull


def calcInsertDistance(
    insertIdx, pamStart, pamSeq, guideStart, guideSeq, strand, kiType, insertSeq, pamPat
):
    """returns the distance between the cut site and the insertion site
    if one of the 15 bases at the extremity of the guide don't overlap with the insertion site,
    returns doRecoding = True (need to recode the donor DNA to prevent its cleavage)
    """

    pamWindowStart = pamStart
    pamWindowEnd = pamStart + len(pamSeq) - 1

    if strand == "+":
        guideWindowStart = guideStart + (len(guideSeq) - 15)
        guideWindowEnd = guideStart + len(guideSeq)
    else:
        guideWindowStart = guideStart
        guideWindowEnd = guideStart + 15
    if pamIsFirst:
        if strand == "+":
            cutPos = guideStart + 19
        else:
            cutPos = pamStart - 18
    else:
        if strand == "+":
            cutPos = pamStart - 3
        else:
            cutPos = guideStart + 4

    if kiType == "deletion" or kiType == "replacement":
        editEnd = insertIdx + len(insertSeq)

        # overlap between guide recoding window / pam and deletion
        # deletion start overlaps guide seq \ guide seq in deletion \ deletion end overlaps guide seq
        cutInGuide = (
            (insertIdx >= guideWindowStart and insertIdx <= guideWindowEnd)
            or (insertIdx <= guideWindowStart and editEnd >= guideWindowEnd)
            or (editEnd >= guideWindowStart and editEnd <= guideWindowEnd)
        )

        cutInPam = (
            (insertIdx >= pamWindowStart and insertIdx <= pamWindowEnd)
            or (insertIdx <= pamWindowStart and editEnd >= pamWindowEnd)
            or (editEnd >= pamWindowStart and editEnd <= pamWindowEnd)
        )
    else:
        cutInGuide = insertIdx >= guideWindowStart and insertIdx <= guideWindowEnd
        cutInPam = insertIdx >= pamWindowStart and insertIdx <= pamWindowEnd

    # check if the substitution alters the pam sequence
    # need to atake replacement into account too
    if kiType == "substitution" and cutInPam:
        # position of the substitution relative to the pam sequence
        editPos = abs(pamWindowStart - insertIdx)
        if strand == "-":
            pamPat = revComp(pamPat)
        pamBase = pamPat[editPos]
        if (
            pamBase == "N"
            or (pamBase == "R" and insertSeq in ["G", "A"])
            or (pamBase == "Y" and insertSeq in ["C", "T"])
            or (pamBase == "V" and insertSeq in ["G", "C", "A"])
        ):
            doRecoding = True
        else:
            doRecoding = False

    elif kiType != "substitution":
        if cutInGuide:
            doRecoding = False
        elif cutInPam:
            if strand == "-":
                pamPat = revComp(pamPat)
            editPos = abs(pamWindowStart - insertIdx)
            pamBase = pamPat[editPos:]
            if insertSeq[0: len(pamBase)] == pamBase:
                doRecoding = True
            else:
                doRecoding = False
        else:
            doRecoding = True
    else:
        doRecoding = True

    # insertDistance should be 0 if the cut site is within the deletion
    # cutUpstream is to determine the template strand for ssODN design
    # onPosupstream/downstream is to display insertDistance even
    # if there is no need for strand selection

    if kiType == "deletion" or kiType == "replacement":
        if cutPos >= insertIdx and cutPos <= editEnd:
            insertDistance = 0
            cutUpstream = "onPos"
        elif cutPos < insertIdx:
            insertDistance = insertIdx - cutPos
            if insertIdx - cutPos <= 10:
                cutUpstream = "onPosupstream"
            else:
                cutUpstream = "upstream"
        elif cutPos > editEnd:
            insertDistance = cutPos - editEnd
            if cutPos - editEnd <= 10:
                cutUpstream = "onPosdownstream"
            else:
                cutUpstream = "downsteam"
    else:
        insertDistance = abs(insertIdx - cutPos)
        if insertIdx - cutPos > 10:
            cutUpstream = "upstream"
        elif insertIdx - cutPos < -10:
            cutUpstream = "downstream"
        elif insertIdx - cutPos < 10 and insertDistance != 0:
            cutUpstream = "onPosupstream"
        elif insertIdx - cutPos > -10 and insertDistance != 0:
            cutUpstream = "onPosdownstream"
        else:
            cutUpstream = "onPos"

    return insertDistance, doRecoding, cutUpstream


def mergeGuideInfo(
    seq,
    startDict,
    pamPat,
    otMatches,
    inputPos,
    effScores,
    sortBy=None,
    org=None,
    exonId=None,
    globEffScore=None,
    pamFullName=None,
    insertIdx=None,
    kiType=None,
    insertSeq=None,
    getSuppInfo=False,
    stopGuides=None,
    allEditData=None,
    beFilter=None
):
    """
    merges guide information from the sequence, the efficiency scores and the off-targets.
    creates rows with too many fields. needs refactoring.

    for each pam in startDict, retrieve the guide sequence next to it and score it
    sortBy can be "main", "effScore", "mhScore", "oofScore" or "pos"
    """
    allEnzymes = readRestrEnzymes()
    guideData = []
    guideScores = {}
    hasNotFound = False
    pamIdToSeq = {}
    if getSuppInfo:
        pamIdToSuppInfo = {}
    if pamFullName:
        pamPat = setupPamInfo(pamFullName)

    pamSeqs = list(
        flankSeqIter(seq.upper(), startDict, len(pamPat), True, exonId, pamFullName)
    )

    # transform into a PAM-based dict
    if allEditData:
        editData = buildEditData(allEditData)

        # get a list of all models in base editor mode to sort by model
        # Base editors have different windows, so models may differ for a given guide
        usedBeModelSet = set()
        for editTpl in editData.values():
            for edits in editTpl:
                for edit in edits:
                    _, _, effs, _ = edits
                    for model, _ in effs:
                        usedBeModelSet.add(model)
        # reassign usedBeModel to a list using the order in modelToEnzyme
        usedBeModels = []
        for model in modelToEnzyme.keys():
            for usedModel in usedBeModelSet:
                if model in usedModel:
                    usedBeModels.append(usedModel)
    else:
        editData = None

    if beFilter:
        orgPamList, bePamIds = beFilter

    for pamId, pamStart, guideStart, strand, guideSeq, pamSeq, pamPlusSeq in pamSeqs:
        # matches in genome
        # one desc in last column per OT seq

        # in KO / STOP mode, filter guides that don't introduce a STOP codon
        if stopGuides and pamId not in stopGuides.keys():
            continue

        # in KI mode, discard BE guides that can't introduce the substitution,
        # while keeping all HDR guides for the selected PAM list
        if beFilter and pamPat not in orgPamList and pamId not in bePamIds:
            continue

        if pamId in otMatches:
            pamMatches = otMatches[pamId]
            guideSeqFull = concatGuideAndPam(guideSeq, pamSeq)
            mutEnzymes = matchRestrEnz(allEnzymes, guideSeq, pamSeq, pamPlusSeq, pamPat)
            (
                posList,
                otDesc,
                guideScore,
                guideCfdScore,
                # aggrCfdScore,
                last12Desc,
                ontargetDesc,
                repCount,
            ) = annotateOfftargets(org, pamMatches, guideSeqFull, pamPat, inputPos)
            if repCount != 0:
                guideScore = 0
                guideCfdScore = 0

            # print(guideCfdScore, aggrCfdScore, "<br>")

        # no off-targets found?
        else:
            # assign a CFD score of 100 (not for substitutions or short deletions / replacements)
            posList, otDesc, guideScore = None, "Not found", -1
            guideCfdScore = -1

            hasNotFound = True
            last12Desc = ""
            mutEnzymes = []
            ontargetDesc = ""
            repCount = 0
            seq34Mer = None

        gcFrac = gcContent(guideSeq)

        effScoring = effScores.get(pamId, {})

        beScoring = {}
        beOutcomes = {}

        # in Ki mode, HDR and BE tables share the same guideData
        if editData and pamId in editData:
            edits = editData[pamId]
            # sorting by outcome is not needed : just add the models for the current guide
            for _, _, _, outcomes in edits:
                for beModel, outcome in outcomes:
                    outcome.sort(key=lambda x: x[1], reverse=True)
                    beOutcomes[beModel] = outcome

            for model in usedBeModels:
                # assign 0 if the model isn't available for the current guide to sort by model
                beScore = -1
                for pos, base, effs, outcomes in edits:
                    for beModel, beEff in effs:
                        if beModel == model:
                            beScore = float(beEff)
                beScoring[model] = beScore
            # print(beScoring, beOutcomes, "<br>")

        # make old jobs compatible
        freeEnergy = effScoring.get("freeEnergy", 0)

        if "EVA" in effScoring:
            EVAlike = effScoring["EVA"]
            if EVAlike == "NA" or EVAlike == 0:
                EVAscore = 0
            else:
                EVAscore = calcEVAscore(EVAlike, guideScore)
            effScoring["EVA"] = EVAscore
        # in knock-in mode, get the distance between cut site and insertion site
        if insertIdx is not None:
            insertDistance, doRecoding, cutUpstream = calcInsertDistance(
                insertIdx,
                pamStart,
                pamSeq,
                guideStart,
                guideSeq,
                strand,
                kiType,
                insertSeq,
                pamPat,
            )
            effScoring["insertDistance"] = insertDistance
        else:
            insertDistance = None
            doRecoding = None
            cutUpstream = None

        if effScoring:
            if pamIsFirst:
                globEffScore = "seqDeepCpf1"
            if globEffScore is None:
                globEffScore = "EVA"
            if globEffScore in effScoring:
                mainScore = calcGlobScore(
                    guideSeq,
                    pamSeq,
                    guideScore,
                    guideCfdScore,
                    effScoring,
                    gcFrac,
                    freeEnergy,
                    globEffScore,
                )
            else:
                mainScore = 0
        else:
            mainScore = 0

        guideRow = [
            guideScore,
            guideCfdScore,
            effScoring,
            pamStart,
            guideStart,
            strand,
            pamId,
            guideSeq,
            pamSeq,
            posList,
            otDesc,
            last12Desc,
            mutEnzymes,
            ontargetDesc,
            repCount,
            gcFrac,
            freeEnergy,
            doRecoding,
            cutUpstream,
            mainScore,
            beScoring,
            beOutcomes
        ]
        guideData.append(guideRow)
        guideScores[pamId] = guideScore
        pamIdToSeq[pamId] = guideSeq
        if getSuppInfo:
            pamIdToSuppInfo[pamId] = (doRecoding, insertDistance)
    # when the function is not called in a loop, sort now
    if exonId is None and pamFullName is None and insertIdx is None:
        sortGuideData(guideData, sortBy)

    # should return pamIdToRecode as None instead ?
    if getSuppInfo:
        return guideData, guideScores, hasNotFound, pamIdToSeq, pamIdToSuppInfo

    return guideData, guideScores, hasNotFound, pamIdToSeq


def sortGuideData(guideData, sortBy, returnGuideData=False, exonSort=False):
    "sorts the guide data according to the value of sortBy"

    if sortBy == "main":
        sortFunc = lambda row: row[19]
        reverse = True
    elif sortBy == "pos":
        sortFunc = lambda row: row[3]
        reverse = False
    elif sortBy == "offCount":
        sortFunc = lambda row: len(row[9])
        reverse = False
    elif sortBy == "cfdSpec":
        sortFunc = operator.itemgetter(1)
        reverse = True
    elif sortBy == "spec" or sortBy is None:
        sortFunc = lambda row: row[0]
        reverse = True
    elif sortBy.startswith("beScore."):
        sortFunc = lambda row: row[20].get(sortBy.lstrip("beScore."), 0)
        reverse = True
    elif sortBy is not None and not sortBy.endswith("pec"):
        sortFunc = lambda row: row[2].get(sortBy, 0)
        if sortBy == "insertDistance":
            reverse = False
        else:
            reverse = True
    else:
        errAbort("Unknown sortBy value. This is a bug. Please contact us.")

    guideData.sort(reverse=reverse, key=sortFunc)

    if exonSort:
        getPrefix = lambda row: int((row[6].split(".")[0]))
        guideData.sort(key=getPrefix)


def sortPairedGuides(pairedGuides, pairSortBy):

    reverse = True
    if pairSortBy == "nickDist":
        sortFunc = lambda pairData: pairData[2]
        reverse = False
    elif pairSortBy == "meanGlob":
        sortFunc = lambda pairData: pairData[3]
    elif "rv" in pairSortBy:
        if pairSortBy[2:] == "CFD":
            sortFunc = lambda pairData: pairData[0][2]
        elif "Main" in pairSortBy:
            sortFunc = lambda pairData: pairData[0][3]
        elif "Eff" in pairSortBy:
            sortFunc = lambda pairData: pairData[0][4].get(pairSortBy[5:])
    elif "fw" in pairSortBy:
        if pairSortBy[2:] == "CFD":
            sortFunc = lambda pairData: pairData[1][2]
        elif "Main" in pairSortBy:
            sortFunc = lambda pairData: pairData[1][3]
        elif "Eff" in pairSortBy:
            sortFunc = lambda pairData: pairData[1][4].get(pairSortBy[5:])
    else:
        errAbort("Unknown pairSortBy value. This is a bug. Please contact us.")

    pairedGuides.sort(key=sortFunc, reverse=reverse)


def printDownloadTableLinks(batchId, addTsv=False, nonClassicMode=None):
    print('<div id="downloads" style="text-align:left">')
    print("Download as Excel tables: ", end=" ")
    print(
        '<a href="crispor.py?batchId=%s&download=guides&format=xls">Guides</a>&nbsp;/&nbsp;'
        % batchId,
        end=" ",
    )
    if not pamIsFirst and not saCas9Mode:
        print(
            '<a href="crispor.py?batchId=%s&showAllScores=1&download=guides&format=xls">Guides, all scores</a>&nbsp;/&nbsp;'
            % batchId,
            end=" ",
        )
    print(
        '<a href="crispor.py?batchId=%s&download=offtargets&format=xls">Off-targets</a>'
        % batchId,
        end=" ",
    )
    if not nonClassicMode:
        print(
            (
                '&nbsp;/&nbsp;<a href="crispor.py?batchId=%s&satMut=1">Saturating mutagenesis assistant</a><br>'
                % batchId
            )
        )
    # print "<small>Plasmid Editor: ",
    # print '<a href="crispor.py?batchId=%s&download=genbank">Guides</a></small>' % batchId,

    if addTsv:
        print("<small>Tab-sep format: ", end=" ")
        print(
            '<a href="crispor.py?batchId=%s&download=guides&format=tsv">Guides</a>&nbsp;/&nbsp;'
            % batchId,
            end=" ",
        )
        print(
            '<a href="crispor.py?batchId=%s&download=offtargets&format=tsv">Off-targets</a></small>'
            % batchId,
            end=" ",
        )

    print("</div>")


def hasGeneModels(org):
    "return true if this organism has gene model information"
    geneFname = join(genomesDir, org, org + ".segments.bed")
    return isfile(geneFname)


def getTableColumnWidths(pam, pamFullName, scoreNames, mutScoreNames, usedBeModels=None):
    """Return a dict of px widths shared by printTableHead() and showGuideTable().
    This is the single source of truth for column widths so the two (intentionally
    separate) tables stay column-aligned across every PAM/score/mode combination."""
    if (
        len(scoreNames) == 2
        or pamIsCpf1(pam)
        or (pamIsSaCas9(pam) and pamFullName is None)
    ):
        effTotal = 150
    else:
        effTotal = 270

    if len(mutScoreNames) <= 1:
        outcomeTotal = 45
    else:
        outcomeTotal = 67

    nEff = max(len(scoreNames), 1)
    nOutcome = max(len(mutScoreNames), 1)

    beEffColWidth = 50
    if usedBeModels:
        nBeModels = max(len(usedBeModels), 1)
        beEffTotal = max(beEffColWidth * nBeModels, 150)
    else:
        beEffTotal = 0
        nBeModels = 1

    return {
        "pos": 80,
        "guide": 235 if pamFullName else 245,
        "distance": 110,
        "globalScore": 150,
        "mitSpec": 80,
        "cfdSpec": 60,
        "effTotal": effTotal,
        "effCol": effTotal // nEff,
        "outcomeTotal": outcomeTotal,
        "outcomeCol": outcomeTotal // nOutcome,
        "beEffTotal": beEffTotal,
        "beEffCol": beEffTotal // nBeModels,
        "beOutcome": 350,
        "offTargets": 117,
        "browser": 500,
    }


def _visualColumns(pam, pamFullName, showColumns, scoreNames, mutScoreNames, editData, usedBeModels, colWidths):
    """Yield (dataColId, widthPx) for every visual column, in body-row order.
    Used both for emitting the shared <colgroup> and for computing total width."""
    yield ("pos", colWidths["pos"])
    yield ("guide", colWidths["guide"])
    if pamFullName and editData is None:
        yield ("distance", colWidths["distance"])
    yield ("global", colWidths["globalScore"])
    if not pamIsCpf1(pam):
        yield ("mit", colWidths["mitSpec"])
    if "cfdGuideScore" in showColumns:
        yield ("cfd", colWidths["cfdSpec"])
    if editData is None:
        for scoreName in scoreNames:
            if scoreName in ("oof", "proxGc"):
                continue
            yield ("eff-" + scoreName, colWidths["effCol"])
    if "proxGc" in scoreNames:
        yield ("proxGc", colWidths["effCol"])
    if not baseEditor and editData is None and not pamFullName:
        for mutScoreName in mutScoreNames:
            yield ("outcome-" + mutScoreName, colWidths["outcomeCol"])
    if editData is not None and usedBeModels is not None:
        for model in usedBeModels:
            yield ("beEff-" + re.sub(r"\s+", "", model), colWidths["beEffCol"])
        yield ("beOutcome", colWidths["beOutcome"])
    yield ("offtargets", colWidths["offTargets"])
    yield ("browser", colWidths["browser"])


def printOtColgroup(
    pam, pamFullName, showColumns, scoreNames, mutScoreNames, editData, colWidths, usedBeModels=None
):
    """Emit the <colgroup> that both the header and body tables share.
    Under table-layout:fixed, <col> widths are authoritative, so an identical
    colgroup in both tables pins their columns to the same pixel widths."""
    print("<colgroup>")
    for colId, width in _visualColumns(
        pam, pamFullName, showColumns, scoreNames, mutScoreNames, editData, usedBeModels, colWidths
    ):
        print('<col data-col-id="%s" style="width:%dpx;">' % (colId, width))
    print("</colgroup>")


def getOtTableTotalWidth(
    pam, pamFullName, showColumns, scoreNames, mutScoreNames, editData, colWidths, usedBeModels=None
):
    "sum of all per-column widths — used to set table width + container min-width consistently"
    return sum(
        w
        for _, w in _visualColumns(
            pam, pamFullName, showColumns, scoreNames, mutScoreNames, editData, usedBeModels, colWidths
        )
    )


def printTableHead(
    pam,
    batchId,
    chrom,
    org,
    varHtmls,
    showColumns,
    geneId,
    pamFullName=None,
    nonClassicMode=None,
    editData=None,
    usedBeModels=None
):
    "print guide score table description and columns"
    # one row per guide sequence

    if not pamIsCpf1(pam) or pamFullName:
        print(
            """<div class='substep'>Ranked by default from highest to lowest Global Score. Click on a column title to rank by a specific score.<br>"""
        )
        print("</div>")
    if pamFullName:
        print(
            """Hover over the PAM in the first column of the table to show information about its corresponding enzyme<br> """
        )
        # print("""<b>Our recommendation:</b> Use Fusi for in-vivo (U6) transcribed guides, Moreno-Mateos for in-vitro (T7) guides injected into Zebrafish/Mouse oocytes.<br>""")
        print(
            """If you use this website, please cite our <a href="https://academic.oup.com/nar/article/46/W1/W242/4995687">paper in NAR 2018</a>."""
        )
        print(
            "Too much information? Look at the <a target=_blank href='manual/'>CRISPOR manual</a>.<p>"
        )
    if not pamFullName:
        print(
            "Highlighted in yellow are the 3 non-overlapping guides with the highest selected score"
        )
        print(
            """<img src="%simage/info-small.png" title="A guide is defined as overlapping with another one if their region 10bp upstream of PAM to 3bp downstream of the guide insersects. This region represents the occupancy of the RNP on the DNA strand (see LINK), so you can use these guides simultaneously." class="tooltipster">"""
            % HTMLPREFIX
        )

    printDownloadTableLinks(batchId, nonClassicMode=nonClassicMode)

    print(
        """
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
        div.style.width="auto";
        div.style.height="auto";
        div.style.border="1px solid black";
        div.style.padding="10px";
        div.style.position="fixed";
        div.style.backgroundColor="white";
        div.style.left=x+"px";
        div.style.top=y+"px";

        var pos = parseInt(this.getAttribute("pos"));
        var exonId = parseInt(this.getAttribute("exonId"));
        var nucl = this.textContent;
        if (nucl.toUpperCase()==="T")
            origNucl = "C";
        else if (nucl.toUpperCase()==="G")
            origNucl = "C"
        else if (nucl.toUpperCase()==="C")
            origNucl = "G"
        else
            origNucl = "G";

        var htmls=[];
        htmls.push("The following guides can mutate "+origNucl+" to "+nucl+" at position "+pos+":<br><small>Note : only the most frequent outcome is shown.<br>Click on the edit to go to the table below and see all predicted outcomes.</small>");

        var guides = editData[exonId][pos][nucl];
        guides.sort( function (a, b) { a[6] - b[6] } ); // sort by the first efficiency score
        for (var i=0; i<guides.length; i++) {
            guide = guides[i];
            pamId = guide[0];
            pamStrand = pamId.slice(-1);
            guideSeq = guide[1];
            pam = guide[2];
            // mutPos is relative to the guide, not the extended guide
            mutPos = guide[3] + 4;
            if (pamStrand === "-") {
                mutPos = 30 - mutPos;
            }

            effs = guide[4];
            outcomes = guide[5];
            // specScore = guide[6];

            htmls.push("<p style='font-weight: bold;'>Guide ID : "+pamId+"</p>")

            // subtable to display the proportion of each edit
            // A graph (similar to ForeCastBE web) may be more readable

            htmls.push("<table class='editTable'>");


            htmls.push("<th>Model</th><th>Predicted efficiency</th><th>Most frequent outcome</th><th>Predicted Frequency</th>");
            let count = 0;
            for (outcome of outcomes) {

                let modelName = outcome[0];
                htmls.push("<tr><td>" + modelName + "</td>");
                let eff = effs[count][1] * 100;
                htmls.push("<td>" + eff.toFixed(2) + "</td>");

                count += 1;

                let modelVals = outcome[1];

                modelVals.sort((a, b) => b[1] - a[1] ); // sort by frequency

                let outCount = 0;

                for (edit of modelVals) {

                    outCount += 1;
                    let seq = edit[0];
                    let freq = edit[1]*100;

                    // position of the guide + edit + PAM

                    let guideSpan = "<span style='background-color: rgba(0, 0, 255, 0.35)'>";
                    let editSpan = "<span style='background-color: rgba(255, 255, 0, 0.35)'>";
                    let pamSpan = "<span style='background-color: rgba(0, 255, 255, 0.35)'>";
                    let spanEnd = "</span>";

                    if (pamStrand === "+") {
                        let ext1 = seq.slice(0, 4);
                        let guide1 = guideSpan + seq.slice(4, mutPos) + spanEnd;
                        let editNucl = editSpan + seq.slice(mutPos, mutPos+1) + spanEnd;
                        let guide2 = guideSpan + seq.slice(mutPos+1, 24) + spanEnd;
                        let pam = pamSpan + seq.slice(24, 27) + spanEnd;
                        let ext2 = seq.slice(27, 30);
                        finalSeq = ext1 + guide1 + editNucl + guide2 + pam + ext2;
                    } else if (pamStrand === "-") {
                        let ext1 = seq.slice(0, 3);
                        let pam = pamSpan + seq.slice(3, 6) + spanEnd;
                        let guide1 = guideSpan + seq.slice(6, mutPos-1) + spanEnd;
                        let editNucl = editSpan + seq.slice(mutPos-1, mutPos) + spanEnd;
                        let guide2 = guideSpan + seq.slice(mutPos, 26) + spanEnd;
                        let ext2 = seq.slice(26, 30);
                        finalSeq = ext1 + pam + guide1 + editNucl + guide2 + ext2;
                    } else {
                        finalSeq = "";
                    }

                    // only show the most frequent outcome
                    if (outCount > 1) {
                        continue;
                    }
                    htmls.push("<td>"+finalSeq+"</td><td>"+freq.toFixed(2)+" %</td></td></tr>");
                }
            }
            htmls.push("</table>");

        }

        htmls.push("</table>");
        $(div).append(htmls.join(""));
        document.body.appendChild(div);
    }

    function emptyBeRows() {
    // hides the rows for which no resuts can be shown based on the selected BE models (see below)

        $(".beModel").closest("tr").each(function () {
            var hasVisible = $(this).find(".beModel").filter(function () {
                return this.style.display !== "none";
                }).length > 0;
                $(this).toggle(hasVisible);
        });

    }

    // base editing model columns the user has unchecked (keyed by the
    var hiddenBeModels = new Set();

    var MIN_BE_GROUP_WIDTH = 150;

    function resizeOtTables() {
        ['otTableHeader', 'otTable'].forEach(function (id) {
            var table = document.getElementById(id);
            if (!table) return;
            var beCols = [];
            table.querySelectorAll('col[data-col-id]').forEach(function (col) {
                var colId = col.getAttribute('data-col-id');
                if (colId.indexOf('beEff-') === 0) {
                    // remember the design width the first time we touch this col
                    if (col.dataset.baseWidth === undefined)
                        col.dataset.baseWidth = parseFloat(col.style.width) || 0;
                    var hidden = hiddenBeModels.has(colId.slice('beEff-'.length));
                    col.style.width = (hidden ? 0 : parseFloat(col.dataset.baseWidth)) + 'px';
                    beCols.push(col);
                }
            });
            // keep the beEff group from collapsing to 0 so its spanning
            // "Predicted editing efficiency" header stays visible; pad the first
            // beEff col (its cells are visibility:hidden, so it just reads as empty)
            if (beCols.length) {
                var beWidth = beCols.reduce(function (s, c) { return s + (parseFloat(c.style.width) || 0); }, 0);
                if (beWidth < MIN_BE_GROUP_WIDTH)
                    beCols[0].style.width = ((parseFloat(beCols[0].style.width) || 0) + MIN_BE_GROUP_WIDTH - beWidth) + 'px';
            }
            // recompute the table width from the (now adjusted) column widths
            var total = 0;
            table.querySelectorAll('col[data-col-id]').forEach(function (col) {
                total += parseFloat(col.style.width) || 0;
            });
            table.style.width = total + 'px';
        });
    }

    function showBeModelResults(checkbox, modelId) {
    // hides / shows base editing outcomes and scores for a given model
        // The beEff column is hidden two ways: (1) collapse its <col> to width 0 in
        // resizeOtTables() to reclaim the space, and (2) set its th/td to
        // visibility:hidden here. We use visibility, NOT display:none: display:none
        // removes a cell from the column flow (shifting every column to its right),
        // whereas visibility:hidden keeps the cell in flow but hides its content and
        // any rotated-label overflow. The outcome <div>s inside the beOutcome cell
        // are block elements, so display:none is fine for them.
        const $models = $('[name="' + modelId + '"]');
        const $cells = $('.beEffCol-' + modelId);

        if (checkbox.checked) {
            $models.show();
            $cells.css('visibility', 'visible');
            hiddenBeModels.delete(modelId);
        } else {
            $models.hide();
            $cells.css('visibility', 'hidden');
            hiddenBeModels.add(modelId);
        }
        emptyBeRows();
        resizeOtTables();
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

    function changeGlobEffScore() {
        var newScore = $('input[name="globEffScore"]:checked').val();
        var searchParams = new URLSearchParams(window.location.search);
        searchParams.set('globEffScore', newScore);
        window.location.search = searchParams.toString();
    }

    function openOtPrimers(batchId, pamId) {
        var url = "crispor.py?batchId=" + batchId + "&pamId=" + pamId + "&otPrimers=1";
        if ($("#onlyExonBox").prop("checked")) {
            url += "&onlyExons=1";
        }
        if ($("#onlySameChromBox").prop("checked")) {
            url += "&onlyChrom=1";
        }
        window.open(url, '_blank');
    }

    </script>

    """
    )

    print('<div id="guideTableScroll" style="overflow-x:auto; max-width:100vw;">')

    colWidths = getTableColumnWidths(pam, pamFullName, scoreNames, mutScoreNames, usedBeModels=usedBeModels)
    tableWidth = max(
        1650,
        getOtTableTotalWidth(
            pam, pamFullName, showColumns, scoreNames, mutScoreNames, editData, colWidths, usedBeModels=usedBeModels
        ),
    )

    print(
        """<div class="otTableWrap" style="width: 100%%; min-width: %dpx; display: table;">"""
        % tableWidth
    )
    print(
        '<table id="otTableHeader" style="background:white; table-layout:fixed; width: %dpx; border-collapse: collapse;">'
        % tableWidth
    )
    printOtColgroup(pam, pamFullName, showColumns, scoreNames, mutScoreNames, editData, colWidths, usedBeModels=usedBeModels)

    print("""<thead style="position:sticky;">""")
    print(
        '<tr style="top: 1px; z-index:2; box-shadow: inset 0 1px black; border-left:5px solid black; border-bottom: none; background-color:#F0F0F0; background-clip: padding-box;">'
    )

    print(
        '<th data-col-id="pos" style="top: 1px; z-index:2; box-shadow: inset -1px 0 black; width:%dpx; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=pos" class="tooltipster" title="Click to sort the table by the position of the PAM site">Position/<br>Strand</a>'
        % (colWidths["pos"], batchId)
    )
    htmlHelp(
        "You can click on the links in this column to highlight the <br>PAM site in the sequence viewer at the top of the page."
    )
    print("</th>")

    print(
        '<th data-col-id="guide" style="top: 1px; z-index:2; box-shadow: inset -1px 0 black; width:%dpx; border-bottom:none">Guide Sequence + <i>PAM</i><br>'
        % colWidths["guide"]
    )

    if not pamFullName:
        print("+ Restriction Enzymes")
        htmlHelp(
            "Restriction enzymes can be very useful for screening mutations induced by the guide RNA using PCR and Restrictrion frament length polymorphism (RFLP).<br>Enzyme sites shown here overlap the main cleavage site 3bp 5' to the PAM.<br>Digestion of the PCR product with these enzymes will not cut the product if the genome was mutated by Cas9. This is a lot easier than screening with the T7 assay, Surveyor or sequencing."
        )
        print("<br>")

    if varHtmls is not None:
        print(" + Variants")
        htmlHelp(
            "Variants that overlap the guide sequence are shown. You can change the variant database with the drop-down box above the sequence viewer at the top of the page."
        )
        print("<br>")

    if not pamFullName:
        print("""<small>""")
        print(
            """<input type="checkbox" class="prefixBox" id="onlyWithGBox" onchange="onlyWith('G')">Only G-"""
        )
        print(
            """<input type="checkbox" class="prefixBox" id="onlyWithGGBox" onchange="onlyWith('GG')">Only GG-"""
        )
        print(
            """<input type="checkbox" class="prefixBox" id="onlyWithABox" onchange="onlyWith('A')">Only A-"""
        )
        htmlHelp(
            "The three checkboxes allow you to show only guides that start with GG-, G- or A-. While we recommend prefixing a 20bp guide with G for U6 expression with spCas9, some protocols recommend using only guides with a G- prefix for U6 and A- for U3."
        )
        print("""</small>""")

    if pamFullName and editData is None:
        print(
            """ <th data-col-id="distance" style="top: 0; z-index:2;  box-shadow: inset -1px 0 black; width:%dpx; border-bottom:none;">
            <a href="crispor.py?batchId=%s&sortBy=insertDistance" class="tooltipster" title="The distance between the cut site (3bp 5' of the PAM on the non-target strand for spCas9 and 18bp 3' of the PAM for Cas12a (Cpf1) and the edit site. Click to sort this table by this distance (default)">Distance between cut site and editing site</a><br>
            </th>"""
            % (colWidths["distance"], batchId)
        )
        isDefaultText = ""
    else:
        isDefaultText = " (default)"
    print(
        """<th data-col-id="global" style="top: 0; z-index:2;  box-shadow: inset -1px 0 black; width:%dpx; border-bottom:none;">
        <a href="crispor.py?batchId=%s&sortBy=main" class="tooltipster" title="Click to sort the table by this score%s. Hover over the (i) bubble on the right to get more information about how this score is calculated.">Global Score</a>"""
        % (colWidths["globalScore"], batchId, isDefaultText)
    )

    if pamFullName:
        globScoreAddStr = "Note : the way this score is calculated differs for Cas12a / Cas9 enzymes. If the search was done with a list of PAMs, the global score isn't comparable between different enzyme types."
    else:
        globScoreAddStr = ""
    htmlHelp(
        """
            The global score is used to rank the guides based on their specificity and efficiency. It ranges from 0 to 100 and is calculated as :<br><br>
            <i>0.6*MIT specificity + 0.4*efficiency - penalties</i><br><br>
Each score is normalized before calculation.<br>
For Cas12a (Cpf1) enzymes, this score is only the predicted efficiency minus penalties, since no specificity score is available for these enzymes.<br>
You can adapt the global score to your delivery method (select below), which changes the efficiency score as well as the penalties.<br><br>

<u>Efficiency scores :</u>
    <ul>
        <li>Rule set 3 is used for guides transcribed <i>in vivo</i>.</li>
        <li>Moreno-Mateos is used for guides transcribed <i>in vitro</i>.</li>
        <li>EVA activity score is used for synthetic guides.</li>
        <li>deepCpf1 is used for Cas12a guides (penalties are still applied relative to the selected production method).</li>
    </ul>
<u>Penalties :</u>
    <ul>
        <li>GC content of lower than 25%% or higher than 75%% : -25 ( <a href='https://doi.org/10.1126/science.1246981'>Wang et al, 2014</a>)).</li>
        <li>Minimum free energy of < -3.6 / -6 kcal/mol : -7.5 / -15 (-15 / -30 for synthetic guides, <a href='https://doi.org/10.1038/s41467-025-59947-0'>Riesenberg et al, 2025</a>).</li>
        <li>Presence of a a stretch of four T that terminates transcription : -25 (only for guides transcribed <i>in vivo</i>).</li>
        <li>Presence of a 'TT' motif : -25 (only for spCas9 guides transcribed <i>in vivo</i>, <a href='https://doi.org/10.1016/j.celrep.2019.01.024'>Graf et al, 2019</a>).</li>
        <li>presence of a 'GCC' motif -40 (only for spCas9 guides, <a href='https://doi.org/10.1016/j.celrep.2019.01.024'>Graf et al, 2019</a>).</li>
    </ul><br>
    %s
    """
        % globScoreAddStr
    )
    print("""<small style="align-text: bottom;"><br> """)
    print("""<small>Select a production method</small><br>""")
    globEffScore = cgiParams.get("globEffScore", "EVA")
    useScores = [
            ("rs3", "cell culture U6", "<i>In vivo</i> transcription from a U6 promoter."),
        ("crisprScan", "T7 transcription", "<i>In vitro</i> transcription from a T7 promoter."),
        ("EVA", "Chemical synthesis", "Chemically synthetized guide RNA."),
    ]
    for scoreVal, scoreLabel, title in useScores:
        checked = "checked" if scoreVal == globEffScore else ""
        print(
            """ <input type=radio name="globEffScore" value="%s" %s onchange="changeGlobEffScore()"/><span class="tooltipsterInteract" title="%s">%s</span><br> """
            % (scoreVal, checked, title, scoreLabel)
        )
    print("</small>")
    print("</th>")
    if not pamIsCpf1(pam):
        print(
            '<th data-col-id="mit" style="top: 0; z-index:2;  box-shadow: inset -1px 0 black; width:%dpx; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=spec" class="tooltipster" title="Click to sort the table by specificity score. Hover over the (i) bubble on the right to get more information about the specificity score.">MIT Specificity Score</a>'
            % (colWidths["mitSpec"], batchId)
        )
        if pamIsSaCas9(pam):
            htmlHelp(
                "The higher the specificity score, the lower are off-target effects in the genome.<br>This specificity score has been adapted for SaCas9 and based on the off-target scores shown on mouse-over. The algorithm was provided by Josh Tycko. Like the MIT score for spCas9, it is aggregated from all off-target scores and ranges 0-100. See <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6063963/'>Tycko et al. Nat Comm 2018</a> for details."
            )
        elif pamFullName == "multipam":
            htmlHelp(
                "The higher the specificity score, the lower are off-target effects in the genome.<br>The specificity score ranges from 0-100 and measures the uniqueness of a guide in the genome. See <a href='http://dx.doi.org/10.1038/nbt.2647'>Hsu et al. Nat Biotech 2013</a>. We recommend values &gt;50, where possible. See <a target=_blank href='manual/#offs'>the CRISPOR manual</a>"
            )
        else:
            htmlHelp(
                "The higher the specificity score, the lower are off-target effects in the genome.<br>The specificity score ranges from 0-100 and measures the uniqueness of a guide in the genome. See <a href='http://dx.doi.org/10.1038/nbt.2647'>Hsu et al. Nat Biotech 2013</a>. We recommend values &gt;50, where possible. See <a target=_blank href='manual/#offs'>the CRISPOR manual</a>"
            )
        print("</th>")

    if "cfdGuideScore" in showColumns:
        print(
            '<th data-col-id="cfd" style="top: 0; z-index:2; box-shadow: inset -1px 0 black; width:%dpx; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=cfdSpec" class="tooltipster" title="Click to sort the table by CFD specificity score">CFD Spec. score</a>'
            % (colWidths["cfdSpec"], batchId)
        )
        htmlHelp(
            "The CFD specificity score, inspired by guidescan.com, behaves like the MIT specificity score, but it is based on the more accurate CFD off-target model, from <a href='http://www.nature.com/nbt/journal/v34/n2/full/nbt.3437.html'>Doench 2016</a>, which is also used by Crispor to rank the off-targets. The CFD specificity score takes into account the identity of mismatches, and correlates better than the MIT score with the total off-target cleavage fraction of a guide, see <a target=_blank href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6731277/'>Tycko et al, Nat Comm 2019</a> and also the <a target=_blank href='/manual/#faq'>CRISPOR manual</a>."
        )
        print("</th>")
    if editData is None:
        if (
            len(scoreNames) == 2
            or pamIsCpf1(pam)
            or pamIsSaCas9(pam)
            and pamFullName is None
        ):
            print(
                '<th style="top: 0; z-index:2;  box-shadow: inset -1px 0 black; width:%dpx; height:100px; border-bottom:none" colspan="%d">Predicted Efficiency'
                % (colWidths["effTotal"], len(scoreNames))
            )
        else:
            if editData is None:
                effColName = "Predicted Efficiency"
            else:
                effColName = "Predicted Nuclease Efficiency"

            print(
                '<th style="top: 0; z-index:2;  box-shadow: inset -1px 0 black; width:%dpx; border-bottom:none" colspan="%d">%s'
                % (colWidths["effTotal"], len(scoreNames), effColName)
            )  # -1 because proxGc is in scoreNames but has no column

        htmlHelp(
            "The higher the efficiency score, the more likely is cleavage at this position. For details on the scores, mouseover their titles below.<br>Note that these predictions are not very accurate, they merely enrich for more efficient guides by a factor of 2-3 so you have to test a few guides to see the effect. <a target=_blank href='manual/#onEff'>Read the CRISPOR manual</a>"
        )

        if not pamIsCpf1(pam) and not pamIsSaCas9(pam) and not geneId:
            if cgiParams.get("showAllScores", "0") == "0":
                print(
                    (
                        """<br><a style="font-size:12px" href="%s" class="tooltipsterInteract" title="By default, only the two most relevant scores are shown, based on our study <a href='http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2'>Haeussler et al. 2016</a>. Click this link to show all efficiency scores.">Show all scores</a>"""
                        % cgiGetSelfUrl({"showAllScores": "1"}, anchor="otTable")
                    )
                )
                scoreDescs["crisprScan"][0] = "Mor.-Mateos"
            else:
                print(
                    (
                        """<br><a style="font-size:12px" href="%s" class="tooltipsterInteract" title="Show only the two main scores">Show main scores</a>"""
                        % cgiGetSelfUrl({"showAllScores": None}, anchor="otTable")
                    )
                )

        print("</th>")

    mhColName = "Outcome"
    if not baseEditor and editData is None and not pamFullName:
        if len(mutScoreNames) <= 1:
            mhColName = ""

        colSpan = len(mutScoreNames)
        print(
            '<th colspan=%d style="top: 0; z-index:2; box-shadow: inset -1px 0 black; width:%dpx; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=oof" class="tooltipster" title="Prediction of the DNA sequence after strand break repair. Click to sort the table by frameshift/out-of-frame scores. Hover over the score names to show information about a particular score. Click a score number to see the predicted indel pattern around the guide.">%s</a>'
            % (colSpan, colWidths["outcomeTotal"], batchId, mhColName)
        )
        # htmlHelp(scoreDescs["oof"][1])
        # print "<small>%s</small>" % oofDesc
        print("</th>")

    if editData and usedBeModels:
        print('<th data-col-id="beEffs" colspan="%d" style="top: 0; z-index:2; box-shadow: inset -1px 0 black; width:%dpx; height: 275px; border-bottom:none">' % (len(usedBeModels), colWidths["beEffTotal"]))
        print('Predicted editing efficiency')
        htmlHelp("""This column shows the predicted base editing efficiencies for several deep-learning models.<br>
                 The scores represents the expected percentages of edited reads after sequencing the target locus (total editing).<br>
                 You can click on each column to sort the table by the corresponding score.""")
        print('</th>')

        print('<th data-col-id="beOutcome" style="top: 0; z-index:2; box-shadow: inset -1px 0 black; width:%dpx; border-bottom:none">' % colWidths["beOutcome"])
        print('Predicted outcome sequences')
        htmlHelp("""This column shows the predicted frequency of each possible edit for this guide sequence.<br>
        The guide and PAM sequences are highlighted in blue and cyan, respectively. Edited bases are highlighted in yellow.<br>
                 Note that outcome sequence prediction vary according to the model, as explained below.<br>
                 <ul>
                     <li>FORECasT-BE shows the frequency of each possible edit, taken individually. The percentage represents the proportion of reads containing this edit, <b>relative to the number of edited reads</b>.</li>
                     <li>DeepBE shows the frequency of each possible combination of edits in the target sequence. Here, the percentage represents the proportion of reads containing this particular sequence, <b>relative to the total number of reads</b>.</li>
                 </ul>""")
        print('<br><small>Show / hide model predictions for :<br>')
        # print buttons to show / hide results for each model

        allBeModels = {
                "A &#8594 G Base Editors": ["ForecastBe - ABE", "DeepBe - DeepNG-BE_17m", "DeepBe - DeepNG-BE_8e"],
                "C &#8594 T Base Editors": ["ForecastBe - CBE", "DeepBe - DeepNG-BE_Ss", "DeepBe - DeepNG-BE_YE1"],
                "C &#8594 G Base Editors": ["DeepBe - DeepNG-BE_mini", "DeepBe - DeepNG-BE_CGBE1", "DeepBe - DeepNG-BE_Bi"]
                }

        beModels = set()
        for editList in editData.values():
            for editTpl in editList:
                effs = editTpl[3]
                for model, eff in effs:
                    beModels.add(model)
        for ezType, modelList in allBeModels.items():
            for i, model in enumerate(modelList):
                if model in beModels:
                    # replace the model name with the corresponding enzyme
                    ezName = modelToEnzyme[model.split(" - ")[1]][0]
                    modelStr = model.split(" - ")[0] + " - " + ezName
                    modelHtml = re.sub(r"\s+", "", model)
                    if i == 0:
                        print("""<p style="font-weight: bold;">%s</p>""" % ezType)
                    print("""<input type="checkbox" id="selectBeModel-%s" checked autocomplete="off" onchange="showBeModelResults(this, '%s')"/>%s<br>""" % (modelHtml, modelHtml, modelStr))

        print("</small>")
        print('</th>')

    print(
        '<th data-col-id="offtargets" style="top: 0; z-index:2; box-shadow: inset -1px 0 black; width:%dpx; border-bottom:none"><a href="crispor.py?batchId=%s&sortBy=offCount" class="tooltipster" title="Click to sort the table by number of off-targets">Off-targets for <br>0-1-2-3-4 mismatches<br></a><span style="color:grey">+ next to PAM </span>'
        % (colWidths["offTargets"], batchId)
    )

    altPamsHelp = [pam]
    if pam in offtargetPams:
        altPamsHelp.extend(offtargetPams[pam])

    htmlHelp(
        "For each number of mismatches, the number of off-targets is indicated.<br>Example: 1-3-20-50-60 means 1 off-target with 0 mismatches, 3 off-targets with 1 mismatch, <br>20 off-targets with 2 mismatches, etc.<br>The CRISPOR website only searches up to four mismatches (use the command line version for 5 or 6). Off-targets are considered if they are flanked by one of these motifs: %s .<br>Shown in grey are the off-targets that have no mismatches in the 12 bp adjacent to the PAM. These are the most likely off-targets."
        % (", ".join(altPamsHelp))
    )

    print("</th>")
    print(
        '<th data-col-id="browser" style="top: 0; z-index:2; box-shadow: inset -1px 0 black; width:%dpx; border-bottom:none">Genome Browser links to matches sorted by CFD off-target score'
        % colWidths["browser"]
    )
    htmlHelp(
        "For each off-target the number of mismatches is indicated and linked to a genome browser. <br>Matches are ranked by CFD off-target score (see Doench 2016 et al) from most to least likely.<br>Matches can be filtered to show only off-targets in exons or on the same chromosome as the input sequence.<br>On most organisms, you can click the links below to open a window with a genome browser at this position."
    )

    print("<br><small>")
    print('<input type="hidden" name="batchId" value="%s">' % batchId)

    if hasGeneModels(org):
        print(
            """<input type="checkbox" id="onlyExonBox" onchange="onlyExons()">in exons only"""
        )
    else:
        print(
            '<small title="When this genome was loaded into CRISPOR, gene models were not available. Contact us if you want to filter for off-targets in exons and think that a gene models are now available for this genome." style="color:grey">No exons.</small>'
        )

    if chrom != "":
        if chrom[0].isdigit():
            chrom = "chrom " + chrom
        print(
            """<input type="checkbox" id="onlySameChromBox" onchange="onlySameChrom()"> on %s only (chromosome of target sequence)"""
            % chrom
        )
    else:
        print('<small style="color:grey">&nbsp;No match, no chrom filter</small>')

    print("</small>")
    print("</th>")
    print("</tr>")

    # subheaders: emit exactly one cell per visual column, in the same order as
    # _visualColumns (the single source of truth shared with the colgroup and body).
    # Columns with a real subheader (eff/proxGc/outcome/beEff) get a rotated <th>;
    # every other column gets an empty offset <th>. Driving the row from
    # _visualColumns guarantees the header stays aligned with the body in every mode,
    # and tagging each cell with its data-col-id means hidden beEff columns collapse
    # together with their <col> automatically.
    print(
        '<tr style="position: sticky; top: 125px; z-index:25; box-shadow: inset -1px 0 black; border-top:none; border-bottom: none; border-left: solid black 5px; background-color:#f0f0f0">'
    )

    # map the whitespace-stripped beEff id (used in data-col-id) back to the full model name
    modelByHtml = {}
    if usedBeModels:
        modelByHtml = {re.sub(r"\s+", "", m): m for m in usedBeModels}

    # empty offset cell for columns with no rotated subheader (pos, guide, global, ...)
    emptyTh = '<th data-col-id="%s" style="position: sticky; top: 125px; z-index:25; box-shadow: inset -1px 0 black; border-top:none"></th>'

    for colId, colWidth in _visualColumns(
        pam, pamFullName, showColumns, scoreNames, mutScoreNames, editData, usedBeModels, colWidths
    ):
        if colId.startswith("eff-") and colId[len("eff-"):] in scoreDescs:
            scoreName = colId[len("eff-"):]
            scoreLabel, scoreDesc = scoreDescs[scoreName]
            print(
                '<th data-col-id="eff-%s" style="position: sticky; top: 125px; z-index:25; box-shadow: inset -1px 0 black; width: 10px; border: none; border-top:none; border-right: none" class="rotate"><div><span><a title="%s" class="tooltipsterInteract" href="crispor.py?batchId=%s&sortBy=%s">%s</a></span></div></th>'
                % (scoreName, scoreDesc, batchId, scoreName, scoreLabel)
            )
        elif colId == "proxGc":
            # the ProxGC score
            print(
                """<th data-col-id="proxGc" style="position: sticky; top: 125px; z-index:25; box-shadow: inset -1px 0 black; border: none; border-top:none; border-right: none; border-left:none" class="rotate">"""
            )
            print("""<div><span style="border-bottom:none">""")
            print(
                """<a title="This column shows two heuristics based on observations rather than computational models: <a href='http://www.cell.com/cell-reports/abstract/S2211-1247%2814%2900827-4'>Ren et al</a> 2014 obtained the highest cleavage in Drosophila when the final 6bp contained &gt;= 4 GCs, based on data from 39 guides. <a href='http://www.genetics.org/content/early/2015/02/18/genetics.115.175166.abstract'>Farboud et al.</a> obtained the highest cleavage in C. elegans for the 10 guides that ended with -GG, out of the 50 guides they tested.<br>The column contains + if the final GC count is &gt;= 4 and GG if the guide ends with GG." href="crispor.py?batchId=%s&sortBy=finalGc6" class="tooltipsterInteract">Prox GC</span></div></th>"""
                % (batchId)
            )
        elif colId.startswith("outcome-"):
            scoreName = colId[len("outcome-"):]
            scoreLabel, scoreDesc = scoreDescs[scoreName]
            print(
                '<th data-col-id="outcome-%s" style="position: sticky; top: 125px; z-index:25; box-shadow: inset -1px 0 black; width: 10px; border-top:none; border-right: none" class="rotate"><div><span><a title="%s" class="tooltipsterInteract" href="crispor.py?batchId=%s&sortBy=%s">%s</a></span></div></th>'
                % (scoreName, scoreDesc, batchId, scoreName, scoreLabel)
            )
        elif colId.startswith("beEff-"):
            modelHtml = colId[len("beEff-"):]
            model = modelByHtml[modelHtml]
            ezName, modelDesc = modelToEnzyme[model.split(" - ")[1]]
            modelStr = model.split(" - ")[0] + " - " + ezName
            print(
                '<th data-col-id="beEff-%s" class="rotate beEffCol-%s" style="position: sticky; top: 125px; z-index:25; box-shadow: inset -1px 0 black; width: %dpx; border: none; border-top:none; border-right: none"><div><span><a title="%s" class="tooltipsterInteract" href="crispor.py?batchId=%s&sortBy=%s">%s</a></span></div></th>'
                % (modelHtml, modelHtml, colWidths["beEffCol"], modelDesc, batchId, "beScore." + model, modelStr)
            )
        else:
            # pos, guide, distance, global, mit, cfd, beOutcome, offtargets, browser
            print(emptyTh % colId)

    print("</tr>")
    print("</thead>")
    print("</table>")

    print("""
    <script>
    // save the states of checkboxes on page reload
    // NOTE: this must run after the WHOLE table (incl. the otTable body rows that
    // hold the .beEffCol-* score cells) has been parsed, otherwise the change
    // handlers find no cells to hide. $(document).ready waits for the full DOM.
    $(document).ready(function() {
        var $checkboxes = $('input[type="checkbox"][id]');
        $checkboxes.each(function() {
            var savedState = localStorage.getItem('checkbox-' + this.id);
            if (savedState !== null) {
                this.checked = savedState === 'true';
                $(this).trigger('change');
            }
        });
        resizeOtTables();

        $checkboxes.on('change', function() {
            localStorage.setItem('checkbox-' + this.id, this.checked);
        });
    });

    </script>


    """)
    print("</div>")


def scoreToColor(guideScore):
    if guideScore is None:
        color = ("#000000", "black")
    elif guideScore > 50:
        color = ("#32cd32", "green")
    elif guideScore > 20:
        color = ("#ffff00", "yellow")
    elif guideScore == -1:
        color = ("#000000", "black")
    else:
        color = ("#aa0114", "red")
    return color


def hexToRgb(hexCode):
    "convert hex color to RGB in UCSC format, https://stackoverflow.com/questions/29643352/converting-hex-to-rgb-value-in-python"
    hexCode = hexCode.lstrip("#")
    return ",".join(tuple(str(int(hexCode[i : i + 2], 16)) for i in (0, 2, 4)))


def makeOtBrowserLinks(otData, chrom, dbInfo, pamId):
    "return a list with the html texts of the offtarget links"
    links = []

    i = 0
    for otSeq, score, cfdScore, editDist, pos, gene, alnHtml, inLinkage in otData:
        cssClasses = ["tooltipster"]
        if not gene.startswith("exon:"):
            cssClasses.append("notExon")
        if pos.split(":")[0] != chrom:
            cssClasses.append("diffChrom")
        if inLinkage:
            cssClasses.append("inLinkage")

        classStr = ""
        if len(cssClasses) != 0:
            classStr = ' class="%s"' % " ".join(cssClasses)

        link = makeBrowserLink(dbInfo, pos, gene, alnHtml, ["tooltipster"])
        editDist = str(editDist)
        links.append("""<div%(classStr)s>%(editDist)s:%(link)s</div>""" % locals())

    return links


def filterOts(otDatas, minScore):
    "remove all offtargets with score < minScore"
    newList = []
    for otData in otDatas:
        score = otData[1]
        if score > minScore:
            newList.append(otData)
    return newList


def findOtCutoff(otData):
    "try cutoffs 0.5, 1.0, 2.0, 3.0 until not more than 20 offtargets left"
    for cutoff in [0.3, 0.5, 1.0, 2.0, 3.0, 10.0, 99.9]:
        otData = filterOts(otData, cutoff)
        if len(otData) <= 30:
            return otData, cutoff

    if len(otData) > 30:
        return otData[:30], None

    return otData, 1000


def findDoubleOts(otCoords1, otCoords2):
    """
    using off-target coordinates of two paired guides (double nicking strategy),
    finds potential off-target double nicking sites
    """

    doubleOts = []

    otTpl1 = [parsePos(posStr) for posStr in otCoords1]
    otTpl2 = [parsePos(posStr) for posStr in otCoords2]

    # temporary solution : not scalable
    # need to make an aglorithm that scales linearly with off-target count

    for chrom1, start1, end1, strand1 in otTpl1:
        for chrom2, start2, end2, strand2 in otTpl2:
            if chrom2 != chrom1:
                continue
            else:
                mean1 = (start1 + end1) // 2
                mean2 = (start2 + end2) // 2

                if abs(mean2 - mean1) < 118:
                    doubleOts.append(
                            ("%s:%s-%s:%s" % (chrom1, start1, end1, strand1),
                             "%s:%s-%s:%s" % (chrom2, start2, end2, strand2))
                            )

    return doubleOts


def printNote(s):
    print(
        '<div style="text-align:left; background-color: aliceblue; padding:5px; border: 1px solid black"><strong>Note:</strong>'
    )
    print(s)
    print("</div>")


def printWarning(s):
    print(
        '<div style="text-align:left; background-color: #FFDDDD; padding:5px; border: 1px solid black"><strong>Warning:</strong>'
    )
    print(s)
    print("</div>")


def printNoEffScoreFoundWarn(effScoresCount, pam):
    if effScoresCount == 0 and not pamIsCpf1(pam):
        note = "No guide could be scored for efficiency. This happens when the input sequence is shorter than 100bp and there is no genome available to extend it or if there is simply not guide socring method. In the first case, please add flanking 50bp on both sides of the input sequence and submit this new, longer sequence. For the second case, you can contact me and suggest an efficiency scoring method, send me the published paper in this case."
        printNote(note)


def showPairedGuidesTable(pairedGuides, annotParams, params, batchId):

    selGeneModel = cgiParams.get("geneModelSelection")
    scriptName = basename(__file__)
    pairSortBy = params.get("pairSortBy", "meanGlob")

    headerCss = """style = "background-color:#F0F0F0;" """
    headerCssSmall = """style = "background-color:#F0F0F0; width: 65px;" """
    headerCssCenter = """style = "background-color:#F0F0F0; width: 65px; text-align: center;" """
    htmlPrefix = HTMLPREFIX

    colspan = 5

    sortPairedGuides(pairedGuides, pairSortBy)

    print("""
        <style>
            .selRow { outline: 4px solid #ff7f04; }
        </style>

        <script>

            function showHdrGuide(pamId) {
                let allTables = document.getElementsByName("guideTablePanel");
                let allButtons = document.getElementsByName("tableSelectButton");

                for (button of allButtons) {
                    if (button.id == 'hdrSelect') {
                        button.className = "assistantButton active tooltipsterInteract";
                    } else {
                        button.className = "assistantButton tooltipsterInteract";
                    }
                }

                // show HDR table and highlight the row corresponding to the selected guide
                for (table of allTables) {
                    if (table.id === 'hdrTable') {
                        table.style.display = "block";


                    } else {
                        table.style.display = "none";
                    }}

                let selected = document.querySelectorAll('.selRow');
                for (let i = 0; i < selected.length; i++) {
                    selected[i].classList.remove('selRow');
                }

                let guideRow = document.getElementById(pamId);
                guideRow.classList.add('selRow');
                guideRow.scrollIntoView();

                }
        </script>
    """)

    print("""<div name="guideTablePanel" id="pairTable">""")

    print("""
    <p>The double nicking strategy relies on using a Cas9 nickase with a pair of guides that flank the edit site, to generate a double strand break. This way, you can use guides that are more distant to position of the edit compared to the single guide, DSB-based method. The guides are in a PAM-out orientation (the guide upstream of the edit is design from the non-target strand, and vice-versa).<br>
    This method is useful if no guides are found close to the edit site, as editing efficiency quickly decrease by this distance.<br>Note that the position of the edit between the two guides doesn't affect efficiency, as long as the distance between nicks is within 40-68bp. For more information, see <a href='https://doi.org/10.1038/s41598-021-98965-y' target='blank'>Schubert et al. 2021</a>.</p>
    """)

    print("<table>")
    print("""
    <thead>
    <tr>
        <th colspan=%(colspan)s %(headerCssCenter)s>Reverse guide</th>

        <th colspan= %(colspan)s %(headerCssCenter)s>Forward guide</th>

        <th rowspan=2 %(headerCssSmall)s><a class="tooltipsterInteract" title="Click to sort the table by the mean of global scores of the pair of guides." href="crispor.py?batchId=%(batchId)s&pairSortBy=meanGlob">Mean of global scores </a></th>

<th rowspan=2 %(headerCss)s><a class="tooltipsterInteract" title="Click to sort the table by this distance." href="crispor.py?batchId=%(batchId)s&pairSortBy=nickDist">Distance between nicks</a> <img src="%(htmlPrefix)simage/info-small.png" title="The distance betweeen nick sites, from 40 to 118bp. Nickase D10A shows optimal editing efficiency for distances between 40 and 68bp (51-68bp for D840A). While pairs of guides with a distance of more than 68bp can be used with D10A, we don't recommend using these with D840A. For more information, see <a href='https://doi.org/10.1038/s41598-021-98965-y' target='blank'>Schubert et al. 2021</a>, supplementary figure 2." class="tooltipster"><br>
<small style="margin-top: 24px; width: 20px;">Click to highlight the corresponding PAMs on the sequence viewer</small>
</th>

        <th rowspan=2 %(headerCss)s>Off-targets double nicking sites <img src="%(htmlPrefix)simage/info-small.png" title="This column shows pairs of off-target sites that are within 118bp of each other, regardless of the orientation. While using a nickase mitagates the risk of off-target effects, such off-target sites may be subject to double strand breaks." class="tooltipster">
</th>

        <th rowspan=2 %(headerCss)s>Design link</th>
    </tr>
    <tr>
        <th %(headerCss)s>Guide Sequence + <i>PAM</i></th>
        <th %(headerCssSmall)s><a class="tooltipsterInteract" title="Click to sort the table by the CFD score of the reverse guide." href="crispor.py?batchId=%(batchId)s&pairSortBy=rvCFD">CFD score</a></th>
        <th %(headerCssSmall)s><a class="tooltipsterInteract" title="Click to sort the table by the global score of the reverse guide." href="crispor.py?batchId=%(batchId)s&pairSortBy=rvEffMain">Global Score</a></th>
        <th %(headerCss)s><a class="tooltipsterInteract" title="Click to sort the table by the Rule Set 3 score of the reverse guide." href="crispor.py?batchId=%(batchId)s&pairSortBy=rvEffrs3">rs3</a></th>
        <th %(headerCss)s><a class="tooltipsterInteract" title="Click to sort the table by the EVA activity score of the reverse guide." href="crispor.py?batchId=%(batchId)s&pairSortBy=rvEffEVA">EVA</a></th>

        <th %(headerCss)s>Guide Sequence + <i>PAM</i></th>
        <th %(headerCssSmall)s><a class="tooltipsterInteract" title="Click to sort the table by the CFD score of the forward guide." href="crispor.py?batchId=%(batchId)s&pairSortBy=fwCFD">CFD score</a></th>
        <th %(headerCssSmall)s><a class="tooltipsterInteract" title="Click to sort the table by the global score of the forward guide." href="crispor.py?batchId=%(batchId)s&pairSortBy=fwEffMain">Global Score</a></th>
        <th %(headerCss)s><a class="tooltipsterInteract" title="Click to sort the table by the Rule Set 3 score of the forward guide." href="crispor.py?batchId=%(batchId)s&pairSortBy=fwEffrs3">rs3</a></th>
        <th %(headerCss)s><a <a class="tooltipsterInteract" title="Click to sort the table by the EVA activity score of the forward guide." href="crispor.py?batchId=%(batchId)s&pairSortBy=fwEffEVA">EVA</a></th>
    </tr>
    </thead>
          """ % locals())

    for i, (guideInfo1, guideInfo2, nickDist, meanScore) in enumerate(pairedGuides):

        # pamId, MIT, CFD, globalScore, effScores, offTargets, guideSeq, pamSeq, cutUpstream, doRecoding
        pamId1, mit1, cfd1, glob1, effScores1, otData1, guideSeq1, pamSeq1, doRecoding1, cutUpstream1 = guideInfo1
        pamId2, mit2, cfd2, glob2, effScores2, otData2, guideSeq2, pamSeq2, doRecoding2, cutUpstream2 = guideInfo2

        otCoords1 = [otTpl[4] for otTpl in otData1]
        otCoords2 = [otTpl[4] for otTpl in otData2]

        doubleOts = findDoubleOts(otCoords1, otCoords2)

        if len(doubleOts) == 0:
            otText = "no double off-targets"
        else:
            otText = doubleOts

        donorParams = {
            "batchId": batchId,
            "pamId": pamId1,
            "pam": "NGG",
            "doRecoding": doRecoding1,
            "cutUpstream": cutUpstream1,
            "doubleNicking": True,
            "insertDistance": effScores1.get("insertDistance"),
            "revPamId": pamId1,
            "revDoRecoding": doRecoding1,
            "revCfd": cfd1,
            "fwPamId": pamId2,
            "fwDoRecoding": doRecoding2,
            "fwCfd": cfd2
        }

        if selGeneModel and selGeneModel not in ("None", "noGenes"):
            donorParams["geneModelSelection"] = selGeneModel

        donorParams.update(annotParams)

        donorLink = '&nbsp;<a href="%s?%s" target="blank"><strong>Design Donor DNA</strong></a>' % (scriptName, urllib.parse.urlencode(donorParams))

        print("""<tr>""")
        print("""<td>
              <span style="font-family: Source Code Pro;" >%s <i>%s</i></span><br>
              <a onclick="showHdrGuide('%s')">show on main table</a>
              </td>""" % (guideSeq1, pamSeq1, pamId1))
        print("""<td>%d</td>""" % cfd1)
        print("""<td>%d</td>""" % glob1)
        print("""<td>%d</td>""" % effScores1.get("rs3"))
        print("""<td>%d</td>""" % effScores1.get("EVA"))

        print("""<td>
              <span style="font-family: Source Code Pro;" >%s <i>%s</i></span><br>
              <a onclick="showHdrGuide('%s')">show on main table</a>
              </td>""" % (guideSeq1, pamSeq1, pamId2))
        print("""<td>%d</td>""" % cfd2)
        print("""<td>%d</td>""" % glob2)
        print("""<td>%d</td>""" % effScores2.get("rs3"))
        print("""<td>%d</td>""" % effScores2.get("EVA"))

        print("""<td>%s</td>""" % round(meanScore, 2))
        print("""<td><a style="font-weight: bold;" class="guidePair" data-id1="%s" data-id2="%s">%s bp</a></td>""" % (pamId1, pamId2, nickDist))
        print("""<td>%s</td>""" % otText)
        print("""<td>%s</td>""" % donorLink)
        print("<tr>")
    print("</table></div>")


def showGuideTable(
    guideData,
    pam,
    otMatches,
    dbInfo,
    batchId,
    org,
    chrom,
    varHtmls,
    geneId=None,
    pamFullName=None,
    koMethod=None,
    pamWindow=None,
    exonSelect=None,
    annotParams=None,
    editData=None
):
    "shows table of all PAM motif matches"
    if pamFullName:
        batchInfo = readBatchAsDict(batchId)
        multipam = batchInfo["multipam"]
        annotParams = annotParams or {}

    if koMethod is not None:
        if koMethod == "excision":
            print(
                "<br><div class='title'>Guide sequences for upstream / downstream regions of %s with PAM %s</div>"
                % (geneId, pam)
            )
        elif koMethod == "promoter":
            print(
                "<br><div class='title'>Guide sequences for upstream / downstream regions of %s promoter with PAM %s</div>"
                % (geneId, pam)
            )
        elif koMethod == "splicing":
            print(
                "<br><div class='title'>Guide sequences for %s exons junctions with PAM %s</div>"
                % (geneId, pam)
            )

        else:
            print(
                "<br><div class='title'>Guide sequences for %s exons with PAM %s</div>"
                % (geneId, pam)
            )

    elif pamFullName:
        if editData is None:
            print(
                """ <br> <div class="title">Guide sequences for PAMs %s bp around the edit site</div>"""
                % pamWindow
            )
        else:
            print(
                """ <br> <div class="title">Guide sequences for base editing</div>"""
            )

    else:
        print("<br><div class='title'>Predicted guide sequences for PAMs</div>")

    if pamFullName or koMethod:
        nonClassicMode = True
    else:
        nonClassicMode = False

    usedBeModels = None
    if editData:
        # dict to display the name of the enzyme instead of the model

        # get a list of all models in base editor mode to get the number of columns to display
        # Base editors have different windows, so models may differ for a given guide
        usedBeModelSet = set()
        for editTpl in editData.values():
            for edits in editTpl:
                for edit in edits:
                    _, _, effs, _ = edits
                    for model, _ in effs:
                        usedBeModelSet.add(model)
        # reassign usedBeModel to a list using the order in modelToEnzyme
        usedBeModels = []
        for model in modelToEnzyme.keys():
            for usedModel in usedBeModelSet:
                if model in usedModel:
                    usedBeModels.append(usedModel)

    global scoreNames
    if geneId and not pamFullName:
        if pamIsCpf1(pam):
            scoreNames = cpf1ScoreNames
        elif pamIsSaCas9(pam):
            scoreNames = saCas9ScoreNames
        else:
            scoreNames = koCas9ScoreNames
    elif pamFullName:
        if multipam == "20bp-NGG":
            scoreNames = ["rs3", "EVA", "crisprScan"]
        else:
            scoreNames = ["rs3", "EVA", "crisprScan", "seqDeepCpf1", "najm"]
    elif cgiParams.get("showAllScores", "0") == "1":
        scoreNames = allScoreNames
    showColumns = set()
    # show the CFD guide score?
    if pamIsSpCas9(pam):
        showColumns.add("cfdGuideScore")

    if koMethod is not None and editData is None:
        showPamWarning(pam)
    showNoGenomeWarning(dbInfo)
    printTableHead(
        pam,
        batchId,
        chrom,
        org,
        varHtmls,
        showColumns,
        geneId,
        pamFullName=pamFullName,
        nonClassicMode=nonClassicMode,
        editData=editData,
        usedBeModels=usedBeModels
    )

    count = 0
    effScoresCount = 0
    showProxGcCol = "proxGc" in scoreNames

    # single source of truth shared with printTableHead so the two tables stay aligned
    colWidths = getTableColumnWidths(pam, pamFullName, scoreNames, mutScoreNames, usedBeModels=usedBeModels)
    effTotalWidth = colWidths["effTotal"]
    effColWidth = colWidths["effCol"]
    beEffTotalWidth = colWidths["beEffTotal"]
    beEffColWidth = colWidths["beEffCol"]
    outcomeColWidth = colWidths["outcomeCol"]
    tableWidth = max(
        1650,
        getOtTableTotalWidth(
            pam, pamFullName, showColumns, scoreNames, mutScoreNames, editData, colWidths, usedBeModels=usedBeModels
        ),
    )

    highlightedGuidesIds = []
    highlightedGuidesPos = []

    print(
        '<div id="guideRowsScroll" style="overflow-y:auto; overflow-x: hidden; max-height: 75vh; min-width: %dpx;">'
        % tableWidth
    )
    print(
        """<div class="otTableWrap" style="width: 100%%; min-width: %dpx;">"""
        % tableWidth
    )
    print(
        '<table id="otTable" style="table-layout: fixed; width: %dpx; border-collapse: collapse;">'
        % tableWidth
    )
    printOtColgroup(pam, pamFullName, showColumns, scoreNames, mutScoreNames, editData, colWidths, usedBeModels=usedBeModels)
    print("<tbody>")

    for guideIdx, guideRow in enumerate(guideData):
        (
            guideScore,
            guideCfdScore,
            effScores,
            pamStart,
            guideStart,
            strand,
            pamId,
            guideSeq,
            pamSeq,
            otData,
            otDesc,
            last12Desc,
            mutEnzymes,
            ontargetDesc,
            repCount,
            gcFrac,
            freeEnergy,
            doRecoding,
            cutUpstream,
            mainScore,
            beScoring,
            beOutcomes
        ) = guideRow

        guidelen = len(guideSeq)
        pamlen = len(pamSeq)

        # don't show the row outside the selected edit / PAM distance in KI mode
        if pamWindow and pamFullName and abs(effScores["insertDistance"] > pamWindow) and editData is None:
            continue

        # in KI / HDR mode, don't show guides with alternative PAMs used for Base Editing
        if pamFullName:
            if editData is None:
                orgPamList = multiPamDict[multipam][0]
                currentPam = pamId.split('.')[0]
                if currentPam not in orgPamList:
                    continue

        # don't show the rows that don't correspond to an edit (base editing mode)
        if editData and pamId not in editData:
            continue

        # flag the 3 best non-overlapping guides
        if koMethod:
            inExon = int(pamId.split(".")[0])
        else:
            inExon = None
        pamStart = guideRow[3]

        if len(highlightedGuidesIds) >= 3 or koMethod is None or editData:
            highlight = False
        else:
            # In KO mode, define Cas9 occupancy region to highlight the best 3 non-overlapping guides
            overlapStart = guideStart - 3 if strand == "+" else pamStart - 10
            overlapEnd = (
                pamStart + pamlen + 10 if strand == "+" else guideStart + guidelen + 3
            )

            for hgOverlapStart, hgOverlapEnd, exon in highlightedGuidesPos:
                if (
                    (overlapStart > hgOverlapStart and overlapStart < hgOverlapEnd)
                    or (overlapEnd > hgOverlapStart and overlapEnd < hgOverlapEnd)
                ) and (inExon is None or exon == inExon):
                    overlap = True
                    highlight = False
                    break
            else:
                overlap = False
            if not overlap or guideIdx == 0:
                highlightedGuidesIds.append(pamId)
                highlightedGuidesPos.append((overlapStart, overlapEnd, inExon))
                highlight = True

        if highlight and pamFullName is None:
            backgroundColor = "#ffffb3"
        else:
            backgroundColor = "#ffffff"

        color = scoreToColor(guideScore)[0]

        classStr = cssClassesFromSeq(guideSeq)
        if geneId is not None and koMethod is not None:
            exonId = int(pamId.split(".")[0])
            # /!\ the exon Ids to reference the rows (for filtering) are 0-based
            # the Ids to show as text are 1-based
            if koMethod == "splicing":
                if exonId % 2 == 0:
                    originalExon = (exonId + 1) // 2
                else:
                    originalExon = exonId // 2
                classStr += " exonRow exon-" + str(originalExon)
            else:
                classStr += " exonRow exon-" + str(exonId)
        print(
            '<tr id="%s" class="%s" style="border-left: 5px solid %s; background-clip: padding-box;">'
            % (pamId, classStr, color)
        )

        # position and strand
        # print '<td id="%s">' % pamId
        print(
            """<td style="width:%dpx; background-color:%s;">"""
            % (colWidths["pos"], backgroundColor)
        )
        print('<a href="#list%s">' % (pamId))
        print(str(pamStart + 1) + " /")
        if strand == "+":
            print("fw")
        else:
            print("rev")

        # in multiseq / multipam mode, the exon number or PAM name are prepended to pamId
        if pamId[0] != "s":
            print("<br>")
            pamPrefix = pamId.split(".")[0]
            if pamPrefix.isdigit():
                exonId = int(pamPrefix)
                if koMethod in ["excision", "promoter"]:
                    if exonId == 0:
                        print("upstream")
                    else:
                        print("downstream")

                elif koMethod == "splicing":
                    if exonId % 2 == 0:
                        if exonSelect and exonSelect.isnumeric():
                            # the text corresponding to the exon differs from exonId (index of multiseq)
                            originalExonText = int(exonSelect)
                        else:
                            originalExonText = (exonId + 1) // 2
                        # make 1-based only before printing
                        print(
                            "exon %d<br>splicing acceptor site" % (originalExonText + 1)
                        )
                    else:
                        if exonSelect and exonSelect.isnumeric():
                            originalExonText = int(exonSelect)
                        else:
                            originalExonText = exonId // 2
                        print("exon %d<br>splicing donor site" % (originalExonText + 1))
                else:
                    print("in exon %d" % (exonId + 1))
            else:
                for desc in pamDesc:
                    if desc[0] == pamPrefix:
                        pamText = desc[1]
                        break
                print("</a>")
                print(
                    """<div class="tooltipster" title="%s">with PAM %s</div>"""
                    % (pamText, pamPrefix)
                )
        else:
            pamPrefix = None

        if pamPrefix is None or pamPrefix.isdigit():
            print("</a>")

        print("</td>")

        # sequence with variants and PCR primer link

        print(
            """<td style="width:%dpx; background-color:%s;">"""
            % (colWidths["guide"], backgroundColor)
        )
        print("<small>")

        # guide sequence + PAM sequence
        if pamIsFirst:
            fullGuideHtml = "<tt><i>" + pamSeq + "</i> " + guideSeq + "</tt>"
            spacePos = len(pamSeq)
        else:
            fullGuideHtml = "<tt>" + guideSeq + " <i>" + pamSeq + "</i></tt>"
            spacePos = len(guideSeq)
        print(fullGuideHtml)
        print("<br>")

        # variant-string
        if varHtmls is not None:
            varFound = False
            varStrs = []
            guideHtmlStart = min(guideStart, pamStart)
            guideHtmls = varHtmls[
                guideHtmlStart: guideHtmlStart + len(guideSeq) + len(pamSeq)
            ]
            if strand == "-":
                guideHtmls = list(reversed(guideHtmls))

            for i in range(len(guideHtmls)):
                html = guideHtmls[i]
                if html != ".":
                    varFound = True
                if i == spacePos:
                    varStrs.append("&nbsp;")
                varStrs.append(html)
            print(("<tt style='color:#888888'>%s</tt><br>" % ("".join(varStrs))))

        if pamFullName and doRecoding and editData is None:
            text = "If you use this guide, the donor DNA needs to be recoded, as the 15bp at the end of the guide don't overlap with the insertion site. Without recoding of the donor this guide sequence will hybridize to it, resulting in its cleavage."
            htmlWarn(text)
            print("Donor needs recoding")
            print("<br>")
        if "TTTT" in guideSeq.upper():
            text = "This guide contains the sequence TTTT. It cannot be transcribed with a U6 or U3 promoter, as TTTT terminates the transcription."
            htmlWarn(text)
            print(" Not with U6/U3")
            print("<br>")

        if pam == "NGG":
            grafType = crisporEffScores.getGrafType(guideSeq)
            if grafType:
                if grafType == "tt":
                    grafText = "The guide ends with TTC or TTT or contains only T and C in the last four nucleotides and more than 2 Ts or at least one TT and one T or C ('TT-motif'). These guides should be avoided in polymerase III (Pol III)-based gene editing experiments requiring high sgRNA expression levels."
                elif grafType == "gcc":
                    grafText = "The guide ends with [AGT]GCC or GCCT ('GCC motif'). These sgRNAs appear to be inefficient irrespective of the delivery method and should thus be generally avoided."

                text = (
                    "This guide contains one of the motifs described by <a target=_blank href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6352712/'>Graf et al, Cell Reports 2019</a>. %s "
                    % grafText
                )
                htmlWarn(text)
                print(" Inefficient")
                print("<br>")

        if gcFrac > 0.75:
            text = "This sequence has a GC content higher than 75%.<br>In the data of Tsai et al Nat Biotech 2015, the two guide sequences with a high GC content had almost as many off-targets as all other sequences combined. We do not recommend using guide sequences with such a high GC content."
            htmlWarn(text)
            print(" High GC content")
            print("<br>")

        if gcFrac < 0.25:
            text = "This sequence has a GC content lower than 25%.<br>In the data of Wang/Sabatini/Lander Science 2014, guides with a very low GC content had low cleavage efficiency."
            htmlWarn(text)
            print(" Low GC content<br>")
            print("<br>")

        if len(mutEnzymes) != 0 and not pamFullName:
            print("<div style='margin-top: 3px'>Enzymes: <i>", end=" ")
            print(", ".join([x.split("/")[0] for x, y, z in list(mutEnzymes.keys())]))
            print("</i></div>")
        else:
            print("<br>")

        scriptName = basename(__file__)
        if pamFullName and batchInfo.get("posStr") != "?" and editData is None:
            donorParams = {
                "batchId": batchId,
                "pamId": pamId,
                "pam": pam,
                "doRecoding": doRecoding,
                "cutUpstream": cutUpstream,
                "insertDistance": effScores.get("insertDistance"),
            }
            # geneModelSelection is omitted from annotParams for the manual case;
            # carry it from cgiParams so donorDesignPage sees useManualAnnotation
            selGeneModel = cgiParams.get("geneModelSelection")
            if selGeneModel and selGeneModel not in ("None", "noGenes"):
                donorParams["geneModelSelection"] = selGeneModel
            donorParams.update(annotParams)
            print(
                '&nbsp;<a href="%s?%s" target="_blank"><strong>Design Donor DNA</strong></a>'
                % (scriptName, urllib.parse.urlencode(donorParams))
            )
            print("<br>")
        if otData is not None and repCount == 0:
            print(
                (
                    '&nbsp;<a href="%s?batchId=%s&pamId=%s&pam=%s" target="_blank"><strong>Cloning / PCR primers</strong></a><br>'
                    % (scriptName, batchId, urllib.parse.quote(str(pamId)), pam)
                )
            )

        # Display RNAfold secondary structure plot
        if freeEnergy < -3.6:
            text = (
                "This sequence can form one or several secondary structures (minimum free energy = %(freeEnergy)s kcal/mol). This is of importance for synthetic gRNAs. For more information, read <a target=blank href='https://doi.org/10.1038/s41467-025-59947-0'>Riesenberg et al, Nat. Commun. 2025</a>."
                % locals()
            )
            htmlWarn(text)
        # else:
        # print("""<nobr><small style="border-radius: 50%; height:8px; width: 8px; background-color: #00ff00; color: #00ff00;">----</small>""")
        print(
            (
                '&nbsp;<a href="%s?batchId=%s&pamId=%s&guideSeq=%s&freeEnergy=%s" target="_blank"><strong>Predicted secondary structure</strong></a>'
                % (
                    scriptName,
                    batchId,
                    urllib.parse.quote(str(pamId)),
                    guideSeq,
                    freeEnergy,
                )
            )
        )

        print("</small>")
        print("</td>")

        # in knock-in mode, show the distance between the cut site and insertion site
        if pamFullName and editData is None:
            insertDistance = effScores.get("insertDistance")

            print('<td style="width:%dpx;">' % colWidths["distance"])
            print("%d bp" % insertDistance)
            print("</td>")

        # Global Score
        print(
            """ <td style="width:%dpx; background-color:%s;"> """
            % (colWidths["globalScore"], backgroundColor)
        )
        if guideScore == None or mainScore == "NA":
            print("No matches")
        else:
            print("%d" % mainScore)
        print("</td>")

        # off-target score, aka specificity score aka MIT score
        if not pamIsCpf1(pam):
            print(
                """ <td style="width:%dpx; background-color:%s;"> """
                % (colWidths["mitSpec"], backgroundColor)
            )
            if guideScore == None:
                print("No matches")
            else:
                print("%d" % guideScore)
            print("</td>")

        # guide score based on CFD scores, aka guidescan score
        if "cfdGuideScore" in showColumns:
            print(
                """ <td style="width:%dpx; background-color:%s;"> """
                % (colWidths["cfdSpec"], backgroundColor)
            )
            if guideCfdScore == None:
                print("No matches")
            else:
                print("%d" % guideCfdScore)
            print("</td>")

        # eff scores
        if effScores == None:
            print(
                '<td colspan="%d" style="width:%dpx; background-color:%s;">Too close to end</td>'
                % (len(scoreNames), effTotalWidth, backgroundColor)
            )
            htmlHelp(
                "The efficiency scores require some flanking sequence<br>This guide does not have enough flanking sequence in your input sequence and could not be extended as it was not found in the genome.<br>"
            )
        elif editData is None:
            for scoreName in scoreNames:
                # out-of-frame and prox. gc need special treatment
                if scoreName in ["oof", "proxGc"]:
                    continue
                score = effScores.get(scoreName, None)
                if score != None:
                    effScoresCount += 1
                # in multipam mode, all scores for all enzymes are displayed, and non available scores are 0.
                if score == None or (score == 0 and pamFullName is not None):
                    print(
                        """<td style="width:%dpx; background-color:%s;">--</td>"""
                        % (effColWidth, backgroundColor)
                    )
                elif scoreName == "ssc":
                    # save some space
                    numStr = "%.1f" % (float(score))
                    print(
                        """<td style="width:%dpx; font-size:small; background-color:%s;">%s</td>"""
                        % (effColWidth, backgroundColor, numStr)
                    )
                elif scoreDigits.get(scoreName, 0) == 0:
                    print(
                        """<td style="width:%dpx; background-color:%s;">%d</td>"""
                        % (effColWidth, backgroundColor, int(score))
                    )
                else:
                    print(
                        """<td style="width:%dpx; background-color:%s;">%0.1f</td>"""
                        % (effColWidth, backgroundColor, float(score))
                    )
            # print "<!-- %s -->" % seq30Mer

        if showProxGcCol:
            print(
                """<td style="width:%dpx; background-color:%s;">"""
                % (effColWidth, backgroundColor)
            )
            # close GC > 4
            finalGc = int(effScores.get("finalGc6", -1))
            if finalGc == 1:
                print("+")
            elif finalGc == 0:
                print("-")
            else:
                print("--")

            # main motif is "NGG" and last nucleotides are GGNGG
            if int(effScores.get("finalGg", 0)) == 1:
                print("<br>")
                print("<small>-GG</small>")
            print("</td>")

        if not baseEditor and editData is None and not pamFullName:
            for mutScoreName in mutScoreNames:
                print(
                    """<td style="width:%dpx; background-color:%s;">"""
                    % (outcomeColWidth, backgroundColor)
                )
                oofScore = str(effScores.get(mutScoreName, None))

                if mutScoreName == "oof":
                    scoreDesc = "out-of-frame deletions"
                else:
                    scoreDesc = "frameshift mutations"

                if oofScore == None or oofScore == "None":
                    print("--")
                else:
                    print(
                        """<a href="%s?batchId=%s&pamId=%s&showMh=%s" target=_blank class="tooltipster" title="This score indicates how likely %s are. Click to show the induced deletions based on the micro-homology around the cleavage site.">%s</a>"""
                        % (
                            myName,
                            batchId,
                            urllib.parse.quote(pamId),
                            mutScoreName,
                            scoreDesc,
                            oofScore,
                        )
                    )
                    # print """<br><br><small><a href="%s?batchId=%s&pamId=%s&showMh=1" target=_blank class="tooltipster">Micro-homology</a></small>""" % (myName, batchId, pamId)
                print("</td>")

        # outcome sequences (for base editing)
        # print sub-columns for each Base editor model (of there is no model for the current guide, shows "-")

        # new method : BE scores are stored in guideData
        if editData:
            for model, beEff in beScoring.items():
                modelHtml = re.sub(r"\s+", "", model)
                if beEff == -1:
                    beEff = '-'
                else:
                    beEff = str(round(beEff * 100, 2))
                print("""<td class="beEffCol-%s" style="width:10px; background-color: %s;">%s</td>""" % (modelHtml, backgroundColor, beEff))
                print("</td>")

            print("""<td style="width:%dpx; background-color: %s;">""" % (colWidths["beOutcome"], backgroundColor))
            # or display a barplot ?
            # make a window for each model, can be toggled in printTableHead by clicking on the appropriate model

            guideSpan = "<span style='background-color: rgba(0, 0, 255, 0.35)'>"
            editSpan = "<span style='background-color: rgba(255, 255, 0, 0.35)'>"
            pamSpan = "<span style='background-color: rgba(0, 255, 255, 0.35)'>"
            spanEnd = "</span>"

            for model, outcomes in beOutcomes.items():

                mainModel, subModel = model.split(" - ")
                ezStr = mainModel + " - " + modelToEnzyme[subModel][0]

                # remove spaces in model names to use it as an id
                modelHtml = re.sub(r"\s+", "", model)
                outcomes.sort(key=lambda x: x[1], reverse=True)
                print("""<div class="beModel" name="%s" style="margin-top: 8px;">""" % modelHtml)
                print(ezStr, "<br>")
                for i, outcome in enumerate(outcomes):

                    # show the top 5 outcomes by default
                    if i == 3 and len(outcomes) > 3:
                        cssPamId = (pamId.replace("-", "minus").replace("+", "plus").replace(".", "_"))
                        print("""<div id="%s-%s" class="otMore" style="display: none;">""" % (cssPamId, modelHtml))

                    outSeq, prop = outcome
                    propStr = str(round(100 * prop, 2)) + " %"

                    if strand == "+":
                        pamRange = range(24, 27)
                        guideRange = range(4, 24)
                    else:
                        pamRange = range(3, 6)
                        guideRange = range(6, 26)

                    spanList = []
                    # or only show the target codon and the edits at this position
                    for j, base in enumerate(outSeq):
                        if base.isupper():
                            spanList.append(editSpan + base + spanEnd)
                        elif j in pamRange:
                            spanList.append(pamSpan + base + spanEnd)
                        elif j in guideRange:
                            spanList.append(guideSpan + base + spanEnd)
                        else:
                            spanList.append(base)
                    outcomeHtml = ''.join(spanList)
                    print("%s - %s<br>" % (outcomeHtml, propStr))
                else:
                    if i > 3:
                        print("</div>")
                        print(
                            """<a id="%s-%sMoreLink" class="otMoreLink" onclick="showAllOts('%s-%s')">"""
                            % (cssPamId, modelHtml, cssPamId, modelHtml)
                        )
                        print("show all outcomes...</a>")
                        print(
                            """<a id="%s-%sLessLink" class="otLessLink" style="display:none" onclick="showLessOts('%s-%s')">"""
                            % (cssPamId, modelHtml, cssPamId, modelHtml)
                        )
                        print("show less outcomes...</a>")
                print("</div>")

        print("</td>")

        # mismatch description
        print(
            """<td style="width:%dpx; background-color:%s;">"""
            % (colWidths["offTargets"], backgroundColor)
        )
        # otCount = sum([int(x) for x in otDesc.split("/")])
        if otData == None:
            # no genome match
            print(otDesc)
            htmlHelp(
                "This exact sequence was not found in the genome.<br>If you have pasted a cDNA multi-exon sequence, note that sequences that overlap a splice site cannot be used as guide sequences. If you only have a cDNA sequence, please BLAST or BLAT your sequence first against the genome, then use the resulting exon from the genome for CRISPOR.<br>This warning also appears if you have selected the wrong or no genome."
            )
        elif repCount > 0:
            print("Repeat")
            htmlHelp(
                "At <= 4 mismatches, %d alignments were found in the genome for this sequence, without looking at the PAM sequence around these alignments.<br>This guide is a repeated region, it is too unspecific.<br>Usually, CRISPR cannot be used to target repeats. Also, note that sequences that include long repeats will make the CRISPOR website slow. You can mask repeats with Ns to speed up the search."
                % repCount
            )
        else:
            print(otDesc)
            print("<br>")

            # mismatch description, last 12 bp
            print('<small style="color:grey">' + last12Desc + "</small><br>")
            otCount = len(otData)
            print("<br><small>%d off-targets</small>" % otCount)
        print("</td>")

        # links to offtargets
        print(
            """<td style=" width: %dpx; background-color:%s;"><small>"""
            % (colWidths["browser"], backgroundColor)
        )
        if otData != None:
            if len(otData) > 500 and len(guideData) > 1:
                otData, cutoff = findOtCutoff(otData)
                if cutoff == None:
                    print(
                        "More than 1000 off-targets, showing only top "
                        + str(len(otData))
                    )
                else:
                    print(
                        "More than 500 off-targets, showing %d with score &gt;%0.1f "
                        % (len(otData), cutoff)
                    )

                htmlHelp(
                    "This guide sequence has a high number of off-targets, its use is discouraged.<br>To show all off-targets, paste only the guide sequence into the input sequence box."
                )

            otLinks = makeOtBrowserLinks(otData, chrom, dbInfo, pamId)

            print("\n".join(otLinks[:3]))
            if len(otLinks) > 3:
                cssPamId = (
                    pamId.replace("-", "minus").replace("+", "plus").replace(".", "_")
                )
                cssPamId = cssPamId + "More"
                print(
                    '<div id="%s" class="otMore" style="display:none; width:100%%;">'
                    % cssPamId
                )
                print("\n".join(otLinks[3:]))

                print(
                    """<a style="float:right;text-decoration:underline" href="#" onclick="openOtPrimers('%s', '%s'); return false;" id="%s">"""
                    % (batchId, urllib.parse.quote(pamId), cssPamId)
                )
                print("<strong>Off-target primers</strong></a>")

                print("</div>")

                print(
                    """<a id="%sMoreLink" class="otMoreLink" onclick="showAllOts('%s')">"""
                    % (cssPamId, cssPamId)
                )
                print("show all...</a>")

                print(
                    """<a id="%sLessLink" class="otLessLink" style="display:none" onclick="showLessOts('%s')">"""
                    % (cssPamId, cssPamId)
                )
                print("show less...</a>")

        print("</small></td>")

        print("</tr>")
        count = count + 1
    print("</tbody>")
    print("</table>")
    print("</div>")  # body otTableWrap
    print("</div>")  # guideRowsScroll
    print("</div>")  # guideTableScroll

    # Column widths drift when the two tables reflow independently (colspans in
    # the header vs one-cell-per-column in the body). This script copies each
    # rendered header column width onto the matching body <col> on load/resize.
    print(
        """
    <script type="text/javascript">
    (function() {
        function syncOtTableColumns() {
            var header = document.getElementById('otTableHeader');
            var body = document.getElementById('otTable');
            if (!header || !body) return;
            var bodyCols = body.querySelectorAll('col[data-col-id]');
            if (!bodyCols.length) return;
            for (var i = 0; i < bodyCols.length; i++) {
                var colId = bodyCols[i].getAttribute('data-col-id');
                var headerCell = header.querySelector('[data-col-id="' + colId + '"]');
                if (!headerCell) continue;
                var rect = headerCell.getBoundingClientRect();
                if (rect.width > 0) {
                    bodyCols[i].style.width = rect.width + 'px';
                }
            }
            // re-collapse any columns the user has hidden (and reset table width)
            // so a window resize / reflow doesn't re-expand them
            if (typeof resizeOtTables === 'function') resizeOtTables();
        }
        function schedule() {
            if (window.requestAnimationFrame) {
                window.requestAnimationFrame(syncOtTableColumns);
            } else {
                setTimeout(syncOtTableColumns, 0);
            }
        }
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', schedule);
        } else {
            schedule();
        }
        window.addEventListener('load', schedule);
        window.addEventListener('resize', schedule);
    })();
    </script>
    """
    )

    printDownloadTableLinks(batchId, addTsv=True, nonClassicMode=nonClassicMode)

    if editData is None:
        printNoEffScoreFoundWarn(effScoresCount, pam)


def linkLocalFiles(listFname):
    """write a <link> statement for each filename in listFname. Version them via mtime
    (-> browser cache)
    """
    for fname in open(listFname).read().splitlines():
        fname = fname.strip()
        if not isfile(fname):
            fname = join(HTMLDIR, fname)
            if not isfile(fname):
                print("missing: %s<br>" % fname)
                continue
        mTime = str(os.path.getmtime(fname)).split(".")[0]  # seconds is enough
        if fname.endswith(".css"):
            # url = fname.replace("/var/www/", "http://tefor.net/")
            print(
                "<link rel='stylesheet' media='screen' type='text/css' href='%s%s?%s'/>"
                % (HTMLPREFIX, fname, mTime)
            )


def printHeader(batchId, title):
    "print the html header"

    print(
        """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">"""
    )
    print("<html><head>")

    if title == None:
        if batchName != "":
            print("""<title>CRISPOR - %s</title>""" % batchName)
        else:
            print("""<title>CRISPOR</title>""")
    else:
        print("""<title>%s</title>""" % title)
    print(
        """
<meta name='description' content='Design CRISPR guides with off-target and efficiency predictions, for more than 100 genomes.'/>
<meta http-equiv='Content-Type' content='text/html; charset=utf-8' />
<meta property='fb:admins' content='692090743' />
<meta name="google-site-verification" content="OV5GRHyp-xVaCc76rbCuFj-CIizy2Es0K3nN9FbIBig" />
<meta property='og:type' content='website' />
<meta property='og:url' content='http://crispor.gi.ucsc.edu/' />
<meta property='og:image' content='http://crispor.gi.ucsc.edu/image/CRISPOR.png' />
<script src="https://cdn.jsdelivr.net/npm/clipboard@2/dist/clipboard.min.js"></script>

"""
    )

    # load jquery and select2 from local copy, not from CDN, for offline use
    print(
        """<script src='%sjs/jquery.min.js'></script>
<script src='%sjs/jquery-ui.min.js'></script>
"""
        % (HTMLPREFIX, HTMLPREFIX)
    )

    print('<script src="%sjs/select2.min.js"></script>' % HTMLPREFIX)
    print(
        '<link rel="stylesheet" type="text/css" href="%sstyle/select2.min.css" />'
        % HTMLPREFIX
    )

    # print('<link rel="stylesheet" href="//fonts.googleapis.com/css?family=Roboto:300,300italic,700,700italic" />')
    # print('<link rel="stylesheet" type="text/css" href="https://cdnjs.cloudflare.com/ajax/libs/normalize/5.0.0/normalize.min.css" />')
    # print('<link rel="stylesheet" href="//cdn.rawgit.com/milligram/milligram/master/dist/milligram.min.css">')

    linkLocalFiles("includes.txt")

    print(
        '<link ref="stylesheet" type="text/css" href="%sstyle/exon-box.css" />'
        % HTMLPREFIX
    )
    print(
        '<link rel="stylesheet" type="text/css" href="%sstyle/tooltipster.css" />'
        % HTMLPREFIX
    )
    print(
        '<link rel="stylesheet" type="text/css" href="%sstyle/tooltipster-shadow.css" />'
        % HTMLPREFIX
    )
    print(
        '<link rel="stylesheet"  href="https://cdnjs.cloudflare.com/ajax/libs/chosen/1.6.2/chosen.css" />'
    )
    print(
        '<link rel="stylesheet"  href="https://cdn.jsdelivr.net/npm/source-code-pro@2.38.0/source-code-pro.css" />'
    )

    # the UFD combobox, https://code.google.com/p/ufd/wiki/Usage
    # patched to allow mouse wheel
    # https://code.google.com/p/ufd/issues/detail?id=86&q=mouse%20wheel
    print(
        '<script type="text/javascript" src="%sjs/jquery.ui.ufd.js"></script>'
        % HTMLPREFIX
    )
    # print '<link rel="stylesheet" type="text/css" href="%sstyle/ufd-base.css" />' % HTMLPREFIX
    print(
        '<link rel="stylesheet" type="text/css" href="%sstyle/plain.css" />'
        % HTMLPREFIX
    )
    print(
        '<link rel="stylesheet" type="text/css"  href="%sstyle/jquery-ui.css" />'
        % HTMLPREFIX
    )
    print(
        '<link rel="stylesheet" type="text/css"  href="%sstyle/assistant.css" />'
        % HTMLPREFIX
    )

    print('<script type="text/javascript" src="js/jquery.tooltipster.min.js"></script>')
    print(
        '<script src="https://cdnjs.cloudflare.com/ajax/libs/chosen/1.6.2/chosen.jquery.min.js"></script>'
    )

    # override the main TEFOR css
    print(
        """
<style>

select { font-size: 80%; }

body {
   text-align: left;
   /* float: left; */
   min-width: 1200px;
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
"""
    )

    # activate tooltipster
    # theme: 'tooltipster-shadow',
    # activate jqueryUI tooltips
    print(
        """
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
    </script>"""
    )

    # style of Jquery UI tooltips, default style is div.ui-tooltip
    print(
        """<style>
        .alignStyle {
            background-color: #FFFFFF;
            width: 350px;
            max-width: 400px;
            height: 110px;
            position : absolute;
            text-align: left;
            border:1px solid #cccccc;
        }
            </style>"""
    )

    # style from https://css-tricks.com/rotated-table-column-headers/ to rotate table headers
    print(
        """<style>
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
      """
    )

    # if we're showing all scores, we have very little space, so turn by
    # 90degrees. otherwise, we can afford 45 degrees, which is easier to read.
    # if cgiParams.get("showAllScores", "0")=="1":
    print(
        """
    -webkit-transform: rotate(-90);
    -moz-transform: rotate(270deg);
    -ms-transform: rotate(270deg);
    -o-transform: rotate(270deg);
    transform: rotate(270deg);
    width: 25px;"""
    )

    # else:
    # print("""
    # -webkit-transform: rotate(-45deg);
    # -moz-transform: rotate(315deg);
    # -ms-transform: rotate(315deg);
    # -o-transform: rotate(315deg);
    # transform: rotate(315deg);
    # width: 25px;""")

    print(
        """
       }

       th.rotate > div > span {
         /* border-bottom: 1px solid #ccc; */
         padding: 0px 3px;
         white-space: nowrap;
       }
    </style>"""
    )

    print("</head>")

    print('<body id="wrapper">')


def firstFreeLine(lineMasks, y, start, end):
    "recursively search for first free line to place a feature (start, end)"
    # print "first free line called with y", y, "<br>"
    if y >= len(lineMasks):
        return None
    lineMask = lineMasks[y]
    for x in range(start, end):
        # print "checking pos", x, "<br>"
        if lineMask[x] != 0:
            return firstFreeLine(lineMasks, y + 1, start, end)
    return y
    # return None


def layoutPamLines(pamLines, seqLen, stopGuides=None):
    """place pam motifs on lines so they don't overlap"""
    # max number of lines in y direction to draw
    MAXLINES = 60
    # amount of free space around each feature
    SLOP = 2

    # bitmask, one per line, 1 = we have a feature here, 0 = no feature here
    lineMasks = []
    for i in range(0, MAXLINES):
        lineMasks.append([0] * (seqLen + 10))

    # dict with lineCount (0...MAXLINES) -> list of (start, strand) tuples
    ftsByLine = defaultdict(list)
    maxY = 0

    # sort by start position to pack them nicely
    pamLines.sort(key=lambda x: x[0])

    for startFt, endFt, label, strand, pamId in pamLines:

        if stopGuides is not None and pamId not in stopGuides:
            continue

        y = firstFreeLine(lineMasks, 0, startFt, endFt)
        if y is None:
            # errAbort("not enough space to plot features")
            continue

        # fill the current mask
        mask = lineMasks[y]
        maskStart = max(startFt - SLOP, 0)
        maskEnd = min(endFt + SLOP, seqLen)
        for i in range(maskStart, maskEnd):
            mask[i] = 1

        maxY = max(y, maxY)
        ft = (startFt, endFt, label, strand, pamId)
        ftsByLine[y].append(ft)
    return ftsByLine, maxY


def distrOnLines(seq, startDict, featLen, pam, exonId=None, stopGuides=None):
    """given a dict with start -> (start,end,name,strand) and a motif len, create lines of annotations such that
    the motifs don't overlap on the lines
    """
    pamLines = getPamLines(seq, startDict, featLen, pam, exonId)
    return layoutPamLines(pamLines, len(seq), stopGuides=stopGuides)


def getPamLines(seq, startDict, featLen, pam, exonId=None, pamFullName=None):
    """generate motif tuples for a single PAM type"""
    pamLines = []
    # Cannot use Unicode here: these symbols are not part of the
    # monospace font on some platforms and therefore their width
    # is not the same as the other characters
    arrNE = "/"
    arrSE = "\\"

    for start in sorted(startDict):
        end = start + featLen
        strand = startDict[start]

        ftSeq = seq[start:end]
        if strand == "+":
            if pamIsFirst:
                if pamIsCas12max(
                    pam
                ):  # modify cleavage sites to 14-16 (target strand) and 24nt (non-target strand) / hfCas12Max
                    label = "%s" % (ftSeq) + ".............%s%s%s.......%s" % (
                        arrNE,
                        arrNE,
                        arrNE,
                        arrSE,
                    )
                else:
                    label = "%s" % (ftSeq) + ".................%s....%s" % (
                        arrNE,
                        arrSE,
                    )
                startFt = start
                endFt = start + len(label)
            else:
                # label = '%s..%s'%(seq[start-3].lower(), ftSeq)
                # label = '---%s'%(ftSeq)
                # label = '&#45;&#45;&#45;%s'%(ftSeq)
                label = "&#8722;&#8722;&#8722;%s" % (ftSeq)
                startFt = start - 3
                endFt = end
        else:
            if pamIsFirst:
                if pamIsCas12max(
                    pam
                ):  # modify cleavage sites to 14-16 (target strand) and 24nt (non-target strand) / hfCas12Max
                    spc1 = "......."
                    spc2 = "............."
                    labelPrefix = "%s%s%s%s%s%s" % (
                        arrSE,
                        spc1,
                        arrNE,
                        arrNE,
                        arrNE,
                        spc2,
                    )
                    label = labelPrefix + ftSeq
                else:
                    spc1 = "...."
                    spc2 = "................."
                    labelPrefix = "%s%s%s%s" % (arrSE, spc1, arrNE, spc2)
                    label = labelPrefix + ftSeq
                startFt = start - len(labelPrefix)
                endFt = startFt + len(label)
            else:
                # label = '%s..%s'%(ftSeq, seq[end+2].lower())
                label = "%s&#45;&#45;&#45;" % (ftSeq)
                startFt = start
                endFt = end + 3

        if exonId is not None:
            pamId = "%d.s%d%s" % (exonId, start, strand)
        elif pamFullName is not None:
            pamId = "%s.s%d%s" % (pamFullName, start, strand)
        else:
            pamId = "s%d%s" % (start, strand)
        pamLines.append((startFt, endFt, label, strand, pamId))
    return pamLines


def writePamFlank(seq, startDict, pam, faFname, pamFullName=None, beFilter=None):
    "write pam flanking sequences to fasta file, optionally with versions where each nucl is removed"
    # print "writing pams to %s<br>" % faFname

    logging.info("writing guides for a single sequence")

    if beFilter:
        orgPamList, bePamIds = beFilter

    faFh = open(faFname, "w")
    for (
        pamId,
        pamStart,
        guideStart,
        strand,
        flankSeq,
        pamSeq,
        pamPlusSeq,
    ) in flankSeqIter(
        seq, startDict, len(pam), True, exonId=None, pamFullName=pamFullName
    ):
        if beFilter and pam not in orgPamList and pamId not in bePamIds:
            continue

        logging.info("wrote PAM %s to fasta" % pamId)

        faFh.write(">%s\n%s\n" % (pamId, flankSeq))
    faFh.close()


def appendMultiPamFlank(seq, startDict, pam, faFh, exonId):
    """appends pam flanking sequences to fasta file, to be used in a loop processing multiple sequences
    adds a unique id and the exonId (= exon number) to the fasta header"""

    # no longer needed, need to fuse with writePamFlank()

    for (
        pamId,
        pamStart,
        guideStart,
        strand,
        flankSeq,
        pamSeq,
        pamPlusSeq,
    ) in flankSeqIter(seq, startDict, len(pam), True, exonId):
        faFh.write(">%s\n%s\n" % (pamId, flankSeq))


def writeMultiFasta(multiseq, faFname, pam):
    """in KO mode, writes all guides corresponding to multiple sequences in a fasta file
    the fasta header is >PAM.seqId.exonId.
    """
    faFh = open(faFname, "w")
    for exonId, seq in multiseq:
        uppSeq = seq.upper()
        startDict, endSet = findAllPams(uppSeq, pam, exonId)
        appendMultiPamFlank(seq, startDict, pam, faFh, exonId)
    faFh.close()


def runCmd(cmd, ignoreExitCode=False, useShell=True):
    "run shell command, check ret code, replaces BIN and SCRIPTS special variables"
    if useShell:
        cmd = cmd.replace("$BIN", binDir)
        cmd = cmd.replace("$PYTHON", sys.executable)
        cmd = cmd.replace("$SCRIPT", scriptDir)
        cmd = "set -o pipefail; " + cmd
        executable = "/bin/bash"
    else:
        cmd = [
            x.replace("$BIN", binDir)
            .replace("$PYTHON", sys.executable)
            .replace("$SCRIPT", scriptDir)
            for x in cmd
        ]
        executable = None

    debug("Running %s" % cmd)
    ret = subprocess.call(cmd, shell=useShell, executable=executable)
    if ret != 0 and not ignoreExitCode:
        if not useShell:
            cmd = " ".join(cmd)
        if commandLineMode:
            logging.error("Error: could not run command %s." % cmd)
            sys.exit(1)
        else:
            print("Server error: could not run command %s, error %d.<p>" % (cmd, ret))
            print(
                "please send us an email, we will fix this error as quickly as possible. %s "
                % contactEmail
            )
            raise
            sys.exit(0)


def isAltChrom(chrom):
    """return true is chrom name looks like it's not on the primary assembly. This is mostly relevant for hg38.
    examples: chr6_*_alt (hg38)
    """
    return chrom.endswith("_alt")


def parseOfftargets(db, batchId, onTargetChrom=""):
    """parse a bed file with annotataed off target matches from overlapSelect,
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
    bedFname = batchBase + ".bed.gz"
    # example input:
    # chrIV 9864393 9864410 s41-|-|5|ACTTGACTG|0    chrIV   9864303 9864408 ex:K07F5.16
    # chrIV   9864393 9864410 s41-|-|5|ACTGTAGCTAGCT|9999    chrIV   9864408 9864470 in:K07F5.16
    debug("reading offtargets from %s" % bedFname)

    # first sort into dict (pamId,chrom,start,end,editDist,strand)
    # -> (segType, segName)
    pamData = {}

    # ifh = open(bedFname) # switched to gzip compression in Dec 2018, converted old files with bash script
    try:
        ifh = gzip.open(bedFname, "rt")
    except FileNotFoundError:
        print(
            "Off-target results from this link were temporarily removed to save space. "
        )
        linkUrl = "crispor.py?batchId=" + batchId
        print(
            "Please go back to <a href='%s'>your job page</a>, which will rerun the job, then try the link again."
            % linkUrl
        )
        exit(0)

    maxOtLines = 500000
    count = 0
    for line in ifh:
        fields = line.rstrip("\n").split("\t")
        count += 1
        if count > maxOtLines:
            print(
                "Error: More than %d off-targets. CRISPR has trouble with handling extremely unspecific inputs. Please email us at "
                "%s and discuss. Your input sequence most likely includes a recent L1HS, SVA or similar repetitive elements. "
                "It is hard to design guides for these, you can try to re-run your input sequence, but without the repetitive element. "
                % (maxOtLines, contactEmail)
            )
            sys.exit(1)
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
        # print pamId, strand, editDist, seq, chrom, start, end, name, segment, "<br>"

        if targetIsAlt:
            if (
                editDist == "0"
                and onTargetChrom.split("_")[0] == chrom.split("_")[0]
                and not pamId in skippedPams
            ):
                logging.debug(
                    "altChrom edge case: target is on alt-chrom, skipping a single 0-mismatch off-target on primary chrom"
                )
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
        if len(nameFields) > 5:
            totalAlnCount = int(nameFields[4])
            if len(nameFields) > 6:
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
        if parNum is not None and chrom == "chrX":
            logging.debug("off-target on PAR region, skipping")
            continue

        # if a offtarget overlaps an intron/exon or ig/exon boundary it will
        # appear twice; in this case, we only keep the exon offtarget
        if otKey in pamData and segType != "ex":
            logging.debug("skipping off-target: ex/ig boundary edge case")
            continue
        pamData[otKey] = (segType, segName)

    # index by pamId and edit distance
    indexedOts = defaultdict(dict)
    for otKey, otVal in pamData.items():
        pamId, chrom, start, end, editDist, seq, strand, totalAlnCount, isRep = otKey
        segType, segName = otVal
        otTuple = (
            chrom,
            start,
            end,
            seq,
            strand,
            segType,
            segName,
            totalAlnCount,
            isRep,
        )
        indexedOts[pamId].setdefault(editDist, []).append(otTuple)

    return indexedOts


class ConsQueue:
    """a pseudo job queue that does nothing but report progress to the console"""

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

        if start > 1000000:
            startStr = "%.2f Mbp" % (float(start) / 1000000)
        else:
            startStr = "%.2f Kbp" % (float(start) / 1000)
        desc = "%s %s" % (chrom, startStr)

        ofh.write(line.rstrip("\n"))
        ofh.write("\t")
        ofh.write(desc)
        ofh.write("\n")
    ofh.close()


def findAllGuides(seq, pam, exonId=None, pamFullName=None):
    startDict, endSet = findAllPams(seq, pam, exonId)
    pamInfo = list(flankSeqIter(seq, startDict, len(pam), False, exonId, pamFullName))
    return pamInfo


def extractMutScores(scoreDict, pamIds):
    "make a list of the guide-related outcome scores in the order of pamIds"
    res = []
    for pamId in pamIds:
        res.append(scoreDict[pamId][0])
    return res


def calcSaveEffScores(
    batchId, seq, extSeq, pam, queue, writeHeader, seqNumber=None, exonId=None, stopGuides=None
):
    """given a sequence and an extended sequence, get all potential guides
    with pam, extend them to 100mers and score them with various eff. scores.
    Return a
    list of rows [headers, (guideSeq, 100mer, score1, score2, score3,...), ... ]

    Also write the results to a database so they can be retrieved later.

    extSeq can be None, if we were unable to extend the sequence
    """
    seq = seq.upper()
    if extSeq:
        extSeq = extSeq.upper()

    pamInfo = findAllGuides(seq, pam, exonId)

    pamIds = []
    guides = []
    longSeqs = []

    for pamId, startPos, guideStart, strand, guideSeq, pamSeq, pamPlusSeq in pamInfo:

        # in KO / stop mode, don't calculate the scores for guides that can't introduce a STOP codon
        if stopGuides and pamId not in stopGuides:
            continue

        logging.debug("PAM ID: %s - guideSeq %s" % (pamId, guideSeq))
        gStart, gEnd = pamStartToGuideRange(startPos, strand, len(pam))
        longSeq = getExtSeq(
            seq, gStart, gEnd, strand, 50 - GUIDELEN, 50, extSeq
        )  # +-50 bp from the end of the guide
        if longSeq != None:
            longSeqs.append(longSeq)
            pamIds.append(pamId)
            guides.append(guideSeq + pamSeq)
    logging.info(longSeqs)
    if len(longSeqs) > 0 and doEffScoring:
        global scoreNames
        enz = None

        if pamIsCpf1(pam) and not pam == "NGTN":
            enz = "cpf1"
            scoreNames = cpf1ScoreNames
        elif pamIsSaCas9(pam):
            enz = "sacas9"
            scoreNames = saCas9ScoreNames
        # for spcas9, we use the extended list for the calculation
        if enz is None:
            if exonId is None:
                scoreNames = allScoreNames
            else:
                scoreNames = koCas9ScoreNames

        effScores = crisporEffScores.calcAllScores(
            longSeqs, enzyme=enz, scoreNames=scoreNames
        )

        if not baseEditor:

            # these are slow algorithms, so store the results for later
            queue.startStep(batchId, "outcome", "Calculating editing outcomes")
            mutScores = crisporEffScores.calcMutSeqs(
                pamIds, longSeqs, enz, scoreNames=mutScoreNames
            )
            saveOutcomeData(batchId, mutScores, seqNumber)

            # for output and sorting, it's easier to treat the outcome-derived scores like an efficiency score
            for mutScoreName in mutScoreNames:
                if mutScoreName in mutScores:
                    effScores[mutScoreName] = extractMutScores(
                        mutScores[mutScoreName], pamIds
                    )

        # make sure the "N bug" reported by Alberto does never happen again:
        # we must get back as many scores as we have sequences
        for scoreName, scores in effScores.items():
            if len(scores) != len(longSeqs):
                print("Internal error when calculating score %s" % scoreName)
                assert False
    else:
        effScores = {}

    activeScoreNames = list(effScores.keys())

    # in knock out mode, don't write the header if no guides were found for the current exon
    if not effScores:
        writeHeader = False

    # reformat to rows, write all scores to file
    assert len(pamIds) == len(guides) == len(longSeqs)
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

    if writeHeader is True:
        headerRow = ["guideId", "guide", "longSeq"]
        headerRow.extend(activeScoreNames)
        rows.insert(0, headerRow)

    return rows


def calcMultiSaveEffScores(batchId, seq, extSeq, pam, queue, pamFullName, iterIdx, beFilter=None):
    """given a sequence and an extended sequence, get all potential guides
    with pam, extend them to 100mers and score them with various eff. scores.
    Return a
    list of rows [headers, (guideSeq, 100mer, score1, score2, score3,...), ... ]

    Also write the results to a database so they can be retrieved later.

    extSeq can be None, if we were unable to extend the sequence
    """

    # need to reduce possible scores
    # need to create the same rows for all PAMs (if the score isn't calculated for this pam, assign the value 0 / NA)
    # move the score calculation to a shared function, and row creation to a specialized one ?

    if beFilter:
        orgPamList, bePamIds = beFilter

    # logging.info("ORGPAMLIST : %s ; BEPAMIDS : %s" % (orgPamList, bePamIds))
    seq = seq.upper()
    if extSeq:
        extSeq = extSeq.upper()
    pam = setupPamInfo(pamFullName)

    pamInfo = findAllGuides(seq, pam, exonId=None, pamFullName=pamFullName)

    pamIds = []
    guides = []
    longSeqs = []

    for pamId, startPos, guideStart, strand, guideSeq, pamSeq, pamPlusSeq in pamInfo:

        # in KI mode, the list of PAMs was extended for base editing guides
        # filter out the guides that can't introduce the substitution,
        # but only for PAMs that are only used for base editing
        if beFilter and pam not in orgPamList and pamId not in bePamIds:
            continue

        logging.info("PAM ID: %s - guideSeq %s" % (pamId, guideSeq))
        gStart, gEnd = pamStartToGuideRange(startPos, strand, len(pam))

        longSeq = getExtSeq(
            seq, gStart, gEnd, strand, 50 - GUIDELEN, 50, extSeq
        )  # +-50 bp from the end of the guide
        # Always append, even if longSeq is None (to keep track of all guides)
        longSeqs.append(longSeq)
        pamIds.append(pamId)
        guides.append(guideSeq + pamSeq)

    # Define consistent output columns for all PAMs

    outputScoreNames = list(allScoreNames) + ["seqDeepCpf1", "najm", "oof", "lindel"]

    if baseEditor and outputScoreNames:
        outputScoreNames.remove("oof")
        outputScoreNames.remove("lindel")

    if len(longSeqs) > 0 and doEffScoring:
        enz = None
        if pamIsCpf1(pam) and not pam == "NGTN":
            enz = "cpf1"
        elif pamIsSaCas9(pam):
            enz = "sacas9"

        global scoreNames

        if enz == "cpf1":
            scoreNames = cpf1ScoreNames
        elif enz == "sacas9":
            scoreNames = saCas9ScoreNames
        else:
            scoreNames = allScoreNames

        global mutScoreNames
        mutScoreNames = ["oof"]  # may remove this ?

        # Filter valid sequences for scoring (cannot score None)
        validIndices = [i for i, s in enumerate(longSeqs) if s is not None]
        validLongSeqs = [longSeqs[i] for i in validIndices]
        validPamIds = [pamIds[i] for i in validIndices]

        # Initialize effScores with 0 for all guides (default value)
        effScores = {}
        for score in outputScoreNames:
            effScores[score] = [0] * len(longSeqs)

        if len(validLongSeqs) > 0:  # doesn't execute for saCas9 PAMS!
            validEffScores = crisporEffScores.calcAllScores(
                validLongSeqs, enzyme=enz, scoreNames=scoreNames
            )
            logging.info("the valid effscores are")
            logging.info(validEffScores)
            # these are slow algorithms, so store the results for later
            queue.startStep(batchId, "outcome", "Calculating editing outcomes")
            mutScores = crisporEffScores.calcMutSeqs(
                validPamIds, validLongSeqs, enz, scoreNames=mutScoreNames
            )
            saveOutcomeData(batchId, mutScores)

            # for output and sorting, it's easier to treat the outcome-derived scores like an efficiency score
            if "oof" in mutScores:
                validEffScores["oof"] = extractMutScores(mutScores["oof"], validPamIds)

            # Map valid scores back to the full list
            for scoreName in validEffScores:
                if scoreName in effScores:
                    for i, validIdx in enumerate(validIndices):
                        effScores[scoreName][validIdx] = validEffScores[scoreName][i]

    else:
        effScores = {}
        if len(longSeqs) > 0:
            for score in outputScoreNames:
                effScores[score] = [0 for i in range(len(longSeqs))]

    # reformat to rows, write all scores to file
    assert len(pamIds) == len(guides) == len(longSeqs)
    rows = []
    for i, (guideId, guide, longSeq) in enumerate(zip(pamIds, guides, longSeqs)):
        row = [guideId, guide, str(longSeq)]  # longSeq can be None
        for scoreName in outputScoreNames:
            scoreList = effScores[scoreName]
            if len(scoreList) > 0:
                row.append(scoreList[i])
            else:
                row.append("noScore?")
        rows.append(row)

    # avoid creating multiple headers for each pam
    if iterIdx == 0:
        headerRow = ["guideId", "guide", "longSeq"]
        headerRow.extend(outputScoreNames)
        rows.insert(0, headerRow)
    return rows


def writeRow(ofh, row):
    "write list to file as tab-sep row"
    row = [str(x) for x in row]
    ofh.write("\t".join(row))
    ofh.write("\n")


def createBatchEffScoreTable(
    batchId,
    queue,
    outFname,
    outFh=None,
    seqInMultiseq=None,
    extSeqInMultiSeq=None,
    seqNumber=None,
    exonId=None,
    stopGuides=None
):
    """annotate all potential guides with efficiency scores and write to file.
    tab-sep file for easier debugging, no pickling
    """

    # Todo: why don't we get these from the caller as arguments instead of reading the batch?
    batchInfo = readBatchAsDict(batchId)
    if seqInMultiseq:
        seq = seqInMultiseq
        extSeq = extSeqInMultiSeq
    else:
        seq = batchInfo["seq"]
        extSeq = batchInfo.get(
            "extSeq"
        )  # cannot always extend a sequence, e.g. when no perfect match

    pam = batchInfo["pam"]
    pam = setupPamInfo(pam)
    seq = seq.upper()
    if extSeq:
        extSeq = extSeq.upper()

    if seqInMultiseq:
        # the function is called in a loop : the file is appended to
        guideFh = outFh
        # if the there are no guides for the first exon, eff scores headers are missing
        # write the headers only if there are scores
        if (guideFh.tell() == 0 and exonId > 0) or exonId == 0:
            writeHeader = True
        else:
            writeHeader = False
    else:
        guideFh = open(outFname, "w")
        logging.info("wrote scores to file")
        writeHeader = True

    guideRows = calcSaveEffScores(
        batchId, seq, extSeq, pam, queue, writeHeader, seqNumber, exonId, stopGuides=stopGuides
    )

    for row in guideRows:
        writeRow(guideFh, row)

    if seqInMultiseq:
        logging.info(
            "Wrote %s rows to eff scores to %s for exon %s"
            % (len(guideRows), guideFh.name, exonId)
        )
    else:
        guideFh.close()


def createMultiBatchEffScoreTable(
    batchId, Fname, queue, pam, pamFullName, extSeq, iterIdx, beFilter=None
):

    pam = setupPamInfo(pamFullName)
    batchInfo = readBatchAsDict(batchId)
    seq = batchInfo["seq"]
    extSeq = batchInfo.get("extSeq")

    guideFh = open(Fname, "w")

    logging.info("writing eff scores for PAM %s" % pam)
    seq = seq.upper() if seq is not None else seq
    guideRows = calcMultiSaveEffScores(
        batchId, seq, extSeq, pam, queue, pamFullName, iterIdx, beFilter=beFilter
    )

    for row in guideRows:
        writeRow(guideFh, row)


def readEffScores(batchId, multipam=None):
    "parse eff scores from tab sep file and return as dict pamId -> dict of scoreName -> value"

    effScoreFname = join(batchDir, batchId) + ".effScores.tab"
    seqToScores = {}
    if isfile(effScoreFname) and os.path.getsize(effScoreFname) > 0:

        for row in lineFileNext(open(effScoreFname)):
            scoreDict = {}
            rowDict = row._asdict()
            # the first three fields are the pamId, shortSeq, longSeq, they are not scores
            allScoreNames = row._fields[3:]
            for scoreName in allScoreNames:
                score = rowDict[scoreName]
                if score == "None" or score == 0:
                    score = 0
                elif "." in score or "e" in score:
                    score = float(score)
                else:
                    score = int(score)
                scoreDict[scoreName] = score
            seqToScores[row.guideId] = scoreDict

    return seqToScores


def findOfftargetsBwa(
    queue, batchId, batchBase, faFname, genome, pamDesc, bedFname, pamFullName=None
):
    "align faFname to genome and create matchedBedFname"
    matchesBedFname = batchBase + ".matches.bed"
    saFname = batchBase + ".sa"

    if pamFullName:
        pam = setupPamInfo(pamFullName)
    else:
        pam = setupPamInfo(pamDesc)

    pamLen = len(pam)
    genomeDir = genomesDir  # make var local, see below

    open(matchesBedFname, "w")  # truncate to 0 size

    # increase MAXOCC if there is only a single query, but only in CGI mode
    # if len(parseFasta(open(faFname)))==1 and not commandLineMode:
    # global MAXOCC
    # global maxMMs
    # MAXOCC=max(HIGH_MAXOCC, MAXOCC)
    # maxMMs=HIGH_maxMMs

    maxDiff = maxMMs
    queue.startStep(
        batchId, "bwa", "Alignment of potential guides, mismatches <= %d" % maxDiff
    )
    convertMsg = "Converting alignments"
    seqLen = GUIDELEN

    bwaM = MFAC * MAXOCC  # -m is queue size in bwa
    cmd = (
        "$BIN/bwa aln -o 0 -m %(bwaM)s -n %(maxDiff)d -k %(maxDiff)d -N -l %(seqLen)d %(genomeDir)s/%(genome)s/%(genome)s.fa %(faFname)s > %(saFname)s"
        % locals()
    )
    runCmd(cmd)

    queue.startStep(batchId, "saiToBed", convertMsg)
    maxOcc = MAXOCC  # create local var from global
    # EXTRACTION OF POSITIONS + CONVERSION + SORT/CLIP
    # the sorting should improve the twoBitToFa runtime
    python = sys.executable
    cmd = (
        "$BIN/bwa samse -n %(maxOcc)d %(genomeDir)s/%(genome)s/%(genome)s.fa %(saFname)s %(faFname)s | $SCRIPT/xa2multi.pl | %(python)s $SCRIPT/samToBed %(pam)s %(seqLen)d | sort -k1,1 -k2,2n | $BIN/bedClip stdin %(genomeDir)s/%(genome)s/%(genome)s.sizes stdout >> %(matchesBedFname)s "
        % locals()
    )
    runCmd(cmd)

    filtMatchesBedFname = batchBase + ".filtMatches.bed"
    queue.startStep(batchId, "filter", "Removing matches without a PAM motif")
    altPats = ",".join(offtargetPams.get(pam, ["na"]))
    bedFnameTmp = bedFname + ".tmp"
    altPamMinScore = str(ALTPAMMINSCORE)
    shmFaFname = join("/dev/shm", genome + ".fa")

    # EXTRACTION OF SEQUENCES + ANNOTATION - big headache!!
    # twoBitToFa was 15x slower than python's twobitreader, after markd's fix it is better
    # but bedtools uses an fa.idx file and also mmap, so is a LOT fasterr
    # arguments: guideSeq, mainPat, altPats, altScore, passTotalAlnCount
    if isfile(shmFaFname):
        logging.info("Using bedtools and genome fasta on ramdisk, %s" % shmFaFname)
        cmd = (
            "time bedtools getfasta -s -name -fi %(shmFaFname)s -bed %(matchesBedFname)s -fo /dev/stdout | $SCRIPT/filterFaToBed %(faFname)s %(pam)s %(altPats)s %(altPamMinScore)s > %(filtMatchesBedFname)s"
            % locals()
        )
    else:
        cmd = (
            "time $BIN/twoBitToFa %(genomeDir)s/%(genome)s/%(genome)s.2bit stdout -bed=%(matchesBedFname)s | %(python)s $SCRIPT/filterFaToBed %(faFname)s %(pam)s %(altPats)s %(altPamMinScore)s > %(filtMatchesBedFname)s"
            % locals()
        )
    # cmd = "$SCRIPT/twoBitToFaPython %(genomeDir)s/%(genome)s/%(genome)s.2bit %(matchesBedFname)s | $SCRIPT/filterFaToBed %(faFname)s %(pam)s %(altPats)s %(altPamMinScore)s %(maxOcc)d > %(filtMatchesBedFname)s" % locals()
    runCmd(cmd)

    segFname = "%(genomeDir)s/%(genome)s/%(genome)s.segments.bed" % locals()

    # if we have gene model segments, annotate them, otherwise just use the chrom position
    if isfile(segFname):
        queue.startStep(batchId, "genes", "Annotating matches with genes")
        cmd = (
            "cat %(filtMatchesBedFname)s | $BIN/overlapSelect %(segFname)s stdin stdout -mergeOutput -selectFmt=bed -inFmt=bed | cut -f1,2,3,4,8 | gzip > %(bedFnameTmp)s "
            % locals()
        )
        runCmd(cmd)
    else:
        queue.startStep(
            batchId, "chromPos", "Annotating matches with chromosome position"
        )
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
    "generate all possible variants of sequence at 1bp-distance"
    seqs = []
    for i in range(0, len(seq)):
        for l in "ACTG":
            if l == seq[i]:
                continue
            newSeq = seq[:i] + l + seq[i + 1 :]
            seqs.append((i, seq[i], l, newSeq))
    return seqs


def expandIupac(seq):
    """expand all IUPAC characters to nucleotides, returns list.
    >>> expandIupac("NY")
    ['GC', 'GT', 'AC', 'AT', 'TC', 'TT', 'CC', 'CT']
    """
    # http://stackoverflow.com/questions/27551921/how-to-extend-ambiguous-dna-sequence
    d = {
        "A": "A",
        "C": "C",
        "B": "CGT",
        "D": "AGT",
        "G": "G",
        "H": "ACT",
        "K": "GT",
        "M": "AC",
        "N": "GATC",
        "S": "CG",
        "R": "AG",
        "T": "T",
        "W": "AT",
        "V": "ACG",
        "Y": "CT",
        "X": "GATC",
    }
    seqs = []
    for i in product(*[d[j] for j in seq]):
        seqs.append("".join(i))
    return seqs


def writeBowtieSequences(inFaFname, outFname, pamPat):
    """write the sequence and one-bp-distant-sequences + all possible PAM sequences to outFname
    Return dict querySeqId -> querySeq and a list of all
    possible PAMs, as nucleotide sequences (not IUPAC-patterns)
    """
    ofh = open(outFname, "w")
    outCount = 0
    inCount = 0
    guideSeqs = {}  # 20mer guide sequences
    qSeqs = (
        {}
    )  # 23mer query sequences for bowtie, produced by expanding guide sequences
    allPamSeqs = expandIupac(pamPat)
    for seqId, seq in parseFastaAsList(open(inFaFname)):
        inCount += 1
        guideSeqs[seqId] = seq
        for pamSeq in allPamSeqs:
            # the input sequence + the PAM
            newSeqId = "%s.%s" % (seqId, pamSeq)
            newFullSeq = seq + pamSeq
            ofh.write(">%s\n%s\n" % (newSeqId, newFullSeq))
            qSeqs[newSeqId] = newFullSeq

            # all one-bp mutations of the input sequence + the PAM
            for nPos, fromNucl, toNucl, newSeq in makeVariants(seq):
                newSeqId = "%s.%s.%d:%s>%s" % (seqId, pamSeq, nPos, fromNucl, toNucl)
                newFullSeq = newSeq + pamSeq
                ofh.write(">%s\n%s\n" % (newSeqId, newFullSeq))
                qSeqs[newSeqId] = newFullSeq
                outCount += 1
    ofh.close()
    logging.debug(
        "Wrote %d variants+expandedPam of %d sequences to %s"
        % (outCount, inCount, outFname)
    )
    return guideSeqs, qSeqs, allPamSeqs


def applyModifStr(seq, modifStrs, strand):
    """bowtie: given a list of pos:toNucl>fromNucl and a seq, return the original seq.
    position is 0-based
    >>> applyModifStr("ACAATAAGACATAAACATATCGG", "14:T>A,21:A>G,22:C>G".split(","), "+")
    'ACAATAAGACATAATCATATCAC'
    """
    seq = list(seq)
    for modifStr in modifStrs:
        # logging.debug( modifStr)
        pos, toFromNucl = modifStr.split(":")
        fromNucl, toNucl = toFromNucl.split(">")
        pos = int(pos)
        if strand == "-":
            fromNucl = revComp(fromNucl)
        seq[pos] = fromNucl
    return "".join(seq)


def parseRefout(tmpDir, guideSeqs, pamLen):
    """parse all .map file in tmpDir and return as list of chrom,start,end,strand,guideSeq,tSeq"""
    fnames = glob.glob(join(tmpDir, "*.map"))

    # while parsing, make sure we keep only the hit with the lowest number of mismatches
    # to the guide. Saves time when parsing.
    posToHit = {}
    hitBestMismCount = {}
    for fname in fnames:
        for line in open(fname):
            # s20+.17:A>G     -       chr8    26869044        CCAGCACGTGCAAGGCCGGCTTC IIIIIIIIIIIIIIIIIIIIIII 7       4:C>G,13:T>G,15:C>G
            (
                guideIdWithMod,
                strand,
                chrom,
                start,
                tSeq,
                weird,
                someScore,
                alnModifStr,
            ) = line.rstrip("\n").split("\t")

            guideId = guideIdWithMod.split(".")[0]
            modifParts = alnModifStr.split(",")
            if modifParts == [""]:
                modifParts = []
            mismCount = len(modifParts)
            hitId = (guideId, chrom, start, strand)
            oldMismCount = hitBestMismCount.get(hitId, 9999)
            if mismCount < oldMismCount:
                hit = (
                    mismCount,
                    guideIdWithMod,
                    strand,
                    chrom,
                    start,
                    tSeq,
                    modifParts,
                )
                posToHit[hitId] = hit
                hitBestMismCount[hitId] = mismCount  # thanks to github user mbsimonovic

    ret = []
    for guideId, hit in posToHit.items():
        mismCount, guideIdWithMod, strand, chrom, start, tSeq, modifParts = hit
        if strand == "-":
            tSeq = revComp(tSeq)
        guideId = guideIdWithMod.split(".")[0]
        guideSeq = guideSeqs[guideId]
        genomeSeq = applyModifStr(tSeq, modifParts, strand)
        start = int(start)
        bedRow = (
            guideId,
            chrom,
            start,
            start + GUIDELEN + pamLen,
            strand,
            guideSeq,
            genomeSeq,
        )
        ret.append(bedRow)

    return ret


def getEditDist(str1, str2):
    """return edit distance between two strings of equal length
    >>> getEditDist("HIHI", "HAHA")
    2
    """
    assert len(str1) == len(str2)
    str1 = str1.upper()
    str2 = str2.upper()

    editDist = 0
    for c1, c2 in zip(str1, str2):
        if c1 != c2:
            editDist += 1
    return editDist


def findOfftargetsBowtie(queue, batchId, batchBase, faFname, genome, pamPat, bedFname):
    "align guides with pam in faFname to genome and write off-targets to bedFname"
    tmpDir = batchBase + ".bowtie.tmp"
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

    genomePath = abspath(join(genomesDir, genome, genome))
    oldCwd = os.getcwd()

    # run bowtie
    queue.startStep(batchId, "bowtie", "aligning with bowtie")
    os.chdir(tmpDir)  # bowtie writes to hardcoded output filenames with --refout
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
    maxOcc = MAXOCC  # meaning in BWA: includes any PAM, in bowtie we have the PAM in the input sequence
    cmd = (
        "$BIN/bowtie -e 1000 %(genomePath)s -f %(bwFaFname)s  -v 3 -y -t -k %(maxOcc)d -m %(maxOcc)d dummy --max tooManyHits.txt --mm --refout --maxbts=2000 -p 4"
        % locals()
    )
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
        logging.debug("PAM seq: %s of %s" % (genomePamSeq, tSeq))
        if genomePamSeq in altPamSeqs:
            minScore = ALTPAMMINSCORE
        elif genomePamSeq in allPamSeqs:
            minScore = MINSCORE
        else:
            logging.debug(
                "Skipping off-target for %s: %s:%d-%d" % (guideId, chrom, start, end)
            )
            continue

        logging.debug("off-target minScore = %f" % minScore)

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
        guideId = guideId.split(".")[0]  # full guide ID looks like s33+.0:A>T
        name = (
            guideId
            + "|"
            + strand
            + "|"
            + str(editDist)
            + "|"
            + tSeq
            + "|"
            + str(guideHitCount)
            + "|"
            + str(otScore)
        )
        row = [chrom, str(start), str(end), name]
        # this way of collecting the features will remove the duplicates
        otKey = (chrom, start, end, strand, guideId)
        logging.debug("off-target key is %s" % str(otKey))
        offTargets[otKey] = row

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
    genomeDir = genomesDir  # make local var
    segFname = "%(genomeDir)s/%(genome)s/%(genome)s.segments.bed" % locals()

    # annotate with genome locus names
    cmd = (
        "$BIN/overlapSelect %(segFname)s %(tempBedPath)s stdout -mergeOutput -selectFmt=bed -inFmt=bed | cut -f1,2,3,4,8 > %(tmpAnnotOffsPath)s"
        % locals()
    )
    runCmd(cmd)

    shutil.move(tmpAnnotOffsPath, bedFname)
    queue.startStep(batchId, "done", "Job completed")

    if DEBUG:
        logging.info("debug mode: Not deleting %s" % tmpDir)
    else:
        shutil.rmtree(tmpDir)


def processSubmission(faFname, genome, pamDesc, bedFname, batchBase, batchId, queue):
    """search fasta file against genome, filter for pam matches and write to bedFName
    optionally write status updates to work queue. Remove faFname.
    """
    batchInfo = readBatchAsDict(batchId)

    if genome == "noGenome":
        posStr = "?"
    elif (
        "batchName" in batchInfo and batchInfo["batchName"].count(":") == 2
    ):  # chrom:start-end:strand
        posStr = batchInfo["batchName"]
    else:
        queue.startStep(
            batchId,
            "bwasw",
            "Searching genome for one 100% identical match to input sequence",
        )
        posStr = findPerfectMatch(batchId)

    batchInfo["posStr"] = posStr

    if posStr != "?":
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
        effScoresFname = join(batchDir, batchId + ".effScores.tab")
        queue.startStep(batchId, "effScores", "Calculating guide efficiency scores")
        createBatchEffScoreTable(batchId, queue, effScoresFname)

    if genome == "noGenome":
        # skip the off-target search entirely
        open(bedFname, "w")  # create a 0-byte file to signal job completion
        queue.startStep(batchId, "done", "Job completed")
        return

    if useBowtie:
        findOfftargetsBowtie(
            queue, batchId, batchBase, faFname, genome, pamDesc, bedFname
        )
    else:
        findOfftargetsBwa(queue, batchId, batchBase, faFname, genome, pamDesc, bedFname)

    if not DEBUG:
        os.remove(faFname)

    return bedFname


def getStopEditData(genome, seq, pam, batchId, koMethod, koGeneId, exonId, exonPosStr, stopGuides):
    """
    To be used in KO / STOP mode, called in a loop for each exon :
    Searches for guides that can introduce STOP codons in the current exon.
    Searches for all potential edits, then gets the annotation of the coding sequence and filters out STOP guides.
    Guides are scored after filtering.
    """

    extendPos = GUIDELEN + 6
    geneModels = getGeneModels(genome)
    if geneModels:
        for model, modelStr in geneModels:
            exonInfo, maxTransIdLen = getExonInfo(
                genome, model, exonPosStr, extendPos=extendPos
            )
            for transId, sym in list(exonInfo.keys()):
                if transId in koGeneId or koGeneId in transId or sym == koGeneId.rstrip("~SYM"):
                    # get the first transcript in common exons mode
                    selTransId = transId
                    break

        startDict, endSet = findAllPams(seq, pam, exonId)
        pamSeqs = list(flankSeqIter(seq, startDict, len(pam), True, exonId=exonId))

        stopEnzymes = ["CBE", "CGBE"]
        wins = [model["win"] for ez in stopEnzymes for model in allBeModels[ez]]
        beWinStart = min(win[0] for win in wins)
        beWinEnd = max(win[1] for win in wins)

        # get a list of potential edits
        _, editData = makeEditLines(
            seq, pamSeqs, beWinStart, beWinEnd, None, exonId
        )
        # logging.info("Edit Data : %s" % editData)
        editInfo = (pam, editData)

        # get guides that can introduce a STOP codon
        _, newStopGuides = makeExonLines(exonInfo, seq, selTransId, koMethod, editInfo=editInfo)

        if len(newStopGuides) > 0:
            # calculate the scores
            editLines, newEditData = makeEditLines(
                seq,
                pamSeqs,
                beWinStart,
                beWinEnd,
                None,
                exonId,
                stopGuides=newStopGuides,
                batchId=batchId,
                pam=pam
            )
        else:
            newEditData = {}

        return newEditData, newStopGuides


def processMultiSeqSubmission(
    multiseq, genome, pam, batchBase, batchId, queue, koMethod
):
    """In Ko mode, writes the off-target and efficiency scores for multiple sequences :
    - writes efficiency scores to batchId.effScores.tab in a loop
    - writes off-targets from a fasta file containing the PAMs of all sequences
    """

    pam = setupPamInfo(pam)

    # global maxMMs
    # maxMMs=4

    # if otBedFname and effScores exists -> use these files

    batchInfo = readBatchAsDict(batchId)

    logging.info("%s exons will be processed" % len(multiseq))

    bedFname = batchBase + ".bed.gz"
    effScoresFname = batchBase + ".effScores.tab"
    bedFnameTmp = bedFname + ".tmp"
    effScoresFnameTmp = effScoresFname + ".tmp"
    faFname = batchBase + ".fa"

    # use the effscores and offtarget files if they already exist
    if isfile(effScoresFname) and isfile(bedFname):
        return bedFname, effScoresFname

    else:

        batchInfo["koMethod"] = koMethod
        batchInfo["exonSeqs"] = []
        batchInfo["extSeqList"] = []
        batchInfo["exonPosStr"] = []

        koGeneId = batchInfo["koGeneId"]

        if koMethod == "stop":
            stopGuides, allEditData = {}, {}

            # only increase sequence length for restrictive PAMs
            if pam not in ["NGN", "NRN"]:
                global MAXSEQLEN
                MAXSEQLEN = 1e4
        else:
            stopGuides = None

        guideFh = open(effScoresFnameTmp, "w")
        # write temp eff scores per sequence
        for seqNumber, (exonId, exonPosStr) in enumerate(multiseq):

            chrom, start, end, strand = parsePos(exonPosStr)

            # extend the exon sequence to find PAMs that result in a cut max 6bp of a splice site
            if not pamIsFirst:
                # in stop mode, extend further to include PAMs that con edit a splice site from the - strand
                if koMethod == "stop":
                    flank = GUIDELEN + 6
                else:
                    flank = GUIDELEN - 6
                extStart = start - flank
                extEnd = end + flank
                extPosStr = "%s:%s-%s:%s" % (chrom, extStart, extEnd, strand)

                extSeq = getSeq(genome, extPosStr)
                # make intron seq lowercase
                seq = (
                    extSeq[0:flank].lower()
                    + extSeq[flank:-flank].upper()
                    + extSeq[-flank:].lower()
                )
                extSeq = extendAndGetSeq(genome, chrom, extStart, extEnd, strand, seq)

            else:
                # get the sequence of the current exon
                seq = getSeq(genome, exonPosStr)
                seq = seq.upper()
                # get a 100bp-extended version of the input seq
                extSeq = extendAndGetSeq(genome, chrom, start, end, strand, seq)

            # in stop mode, discard the guides that can't introduce a STOP codon / splice site mutation before calculating the scores
            if koMethod == "stop":
                newEditData, newStopGuides = getStopEditData(genome, seq, pam, batchId, koMethod, koGeneId,
                                                             exonId, exonPosStr, stopGuides)
                if len(newStopGuides) > 0:
                    allEditData.update(newEditData)
                    stopGuides.update(newStopGuides)

            batchInfo["exonSeqs"].append((exonId, seq))
            logging.info("Exon %s is %s bp long" % (exonId, len(seq)))

            if extSeq is None:
                batchInfo["exonPosStr"] = "?"
                # this can only happen if there is a 100%-M match but small SNPs in it compared to the input sequence
                # so the extension of the input fails.
                # in this case, we also invalidate the position, as there was no perfect match and the user
                # has to do something to fix it
            else:
                logging.debug(
                    "100pb-extended seq (len: %d) is: %s" % (len(extSeq), extSeq)
                )

                # if exonPos could be extended, get the exon seq from extSeq
                seq = extSeq[FLANKLEN:-FLANKLEN]
                batchInfo["extSeqList"].append(extSeq)
                batchInfo["exonPosStr"].append(exonPosStr)

            queue.startStep(batchId, "effScores", "Calculating guide efficiency scores")
            createBatchEffScoreTable(
                batchId, queue, None, guideFh, seq, extSeq, seqNumber, exonId, stopGuides=stopGuides
            )

        guideFh.close()

        # if no STOP codons can be introduced with SpCas9, switch to SpRY
        if koMethod == "stop" and pam != "NRN" and len(stopGuides) == 0:
            logging.info("no guides found for PAM %s: search with SpRY (pam NRN)" % pam)
            pam = "NRN"
            pam = setupPamInfo(pam)
            batchInfo["pam"] = pam

            # rewrite eff scores and editData
            stopGuides, allEditData = {}, {}
            os.remove(effScoresFnameTmp)

            guideFh = open(effScoresFnameTmp, "w")
            for seqNumber, (exonId, exonPosStr) in enumerate(multiseq):
                extSeq = batchInfo["extSeqList"][exonId]
                seq = batchInfo["exonSeqs"][exonId][1]

                newEditData, newStopGuides = getStopEditData(genome, seq, pam, batchId, koMethod, koGeneId,
                                                             exonId, exonPosStr, stopGuides)
                if len(newStopGuides) > 0:
                    allEditData.update(newEditData)
                    stopGuides.update(newStopGuides)

                queue.startStep(batchId, "effScores", "Calculating guide efficiency scores")
                createBatchEffScoreTable(
                    batchId, queue, None, guideFh, seq, extSeq, seqNumber, exonId, stopGuides=stopGuides
                )
            guideFh.close()

        if stopGuides:
            batchInfo["stopGuides"] = stopGuides
            # persist the accumulated, stop-filtered edit data for all exons,
            # not the per-exon `editData` left over from the loop above
            editData = allEditData
        else:
            editData = None

        writeBatchAsDict(batchInfo, batchId, editData=editData)
        writeMultiFasta(batchInfo["exonSeqs"], faFname, pam)

        # find offtargets and append them to the main file
        if useBowtie:
            findOfftargetsBowtie(
                queue, batchId, batchBase, faFname, genome, pam, bedFnameTmp
            )
        else:
            findOfftargetsBwa(
                queue, batchId, batchBase, faFname, genome, pam, bedFnameTmp
            )

        if not DEBUG:
            if isfile(faFname):
                os.remove(faFname)

        # create the final file to end the job
        shutil.move(bedFnameTmp, bedFname)
        shutil.move(effScoresFnameTmp, effScoresFname)

        return bedFname, effScoresFname


def processMultiPamSubmission(genome, seq, posStr, multipam, batchBase, batchId, queue):
    """In KI mode :
    For each PAM in multiPamDesc, creates a fasta file containing guides for each sequence in multiseq.
    Then, search these files against genome, filter for pam matches and append to bedFName.
    optionally write status updates to work queue. Remove faFname.
    """

    # allow up to 3 mismatches for offtarget search
    global maxMMs
    maxMMs = 4
    batchInfo = readBatchAsDict(batchId)

    insertIdx = batchInfo["insertIdx"]
    kiType = batchInfo["kiType"]
    insertSeq = batchInfo["insertSeq"]
    noPerfectMatch = batchInfo.get("noPerfectMatch")

    # kiType = batchInfo["kiType"]

    if seq is None and posStr is None:
        return None, None

    pamList = set(multiPamDict[multipam][0])

    logging.info("PAMs are : %s" % ", ".join(pamList))

    bedFname = batchBase + ".bed.gz"
    effScoresFname = batchBase + ".effScores.tab"

    if isfile(bedFname) and isfile(effScoresFname):
        return bedFname, effScoresFname

    bedFnameTmp = bedFname + ".tmp"
    effScoresFnameTmp = effScoresFname + ".tmp"

    # Clear output offtargets and effscores file
    open(bedFnameTmp, "wb").close()
    open(effScoresFnameTmp, "w").close()

    # initialize sequence data
    uppSeq = seq.upper()
    chrom, start, end, strand = parsePos(posStr)

    if posStr == "?":
        extSeq = None
        batchInfo["extSeq"] = "?"
    elif seq == "null":
        seq = None
        batchInfo["seq"] = None
    else:
        extSeq = extendAndGetSeq(genome, chrom, start, end, strand, seq, noPerfectMatch=noPerfectMatch)
        batchInfo["extSeq"] = extSeq
        if extSeq is None:
            batchInfo["extSeq"] = "?"

    # Base editing can be used for these substitutions
    useBaseEditor = False
    if kiType == "substitution":
        global baseEditor
        fromToNucl = (seq[insertIdx].upper(), insertSeq.upper())

        substInfo = (insertIdx, insertSeq)
        enzyme = None

        for ez, editList in possibleEdits.items():
            if fromToNucl in editList:
                fw, rev = editList
                if fromToNucl == fw or fromToNucl == rev:
                    enzyme = ez

                # baseEditor is reinitialized in setupPamInfo
                useBaseEditor = True

    if useBaseEditor is True:
        # search for alternative PAMs to design guides for base editing
        for pamVariant in pamVariantModels.keys():
            pamList.add(pamVariant)

    writeBatchAsDict(batchInfo, batchId)
    # do eff scoring and off target search for each pam
    pamSeqs = []

    # first, check PAMs that can be used with BE
    for i, pamFullName in enumerate(pamList):

        pam = setupPamInfo(pamFullName)

        if useBaseEditor:
            # reset the global here
            baseEditor = True

            startDict, endSet = findAllPams(seq.upper(), pam)
            singlePamSeqs = list(
                flankSeqIter(
                    seq,
                    startDict,
                    len(pam),
                    True,
                    exonId=None,
                    pamFullName=pamFullName,
                )
            )
            pamSeqs.extend(singlePamSeqs)

    if useBaseEditor:
        # beWinStart, beWinEnd = getBeWin(DEFAULTBEWIN)
        # get the largest editing window for this enzyme type (results will be filtered later)
        enzList = allBeModels[enzyme]
        beWinStart = min([enzDict["win"][0] for enzDict in enzList])
        beWinEnd = max([enzDict["win"][1] for enzDict in enzList])

        _, editData = makeEditLines(
            seq, pamSeqs, beWinStart, beWinEnd, None, substInfo=substInfo, enzyme=enzyme, extSeq=extSeq
        )
        writeBatchAsDict({}, batchId, editData=editData)

        bePamIds = buildEditData(editData).keys()
    else:
        bePamIds = None

    # PAM list selected from the menu (not extended to include BE pams)
    # score and search off-targets for all of these PAMs
    # for PAMs to be used with BE, only score guides that can be used to introduce the substitution
    orgPamList = multiPamDict[multipam][0]

    for i, pamFullName in enumerate(pamList):

        # set globals for this PAM
        pam = setupPamInfo(pamFullName)

        iterBatchBase = "%s.%d" % (batchBase, i)
        faFname = iterBatchBase + ".fa"
        tempBedFname = iterBatchBase + ".bed"
        tempEffScoresFname = iterBatchBase + ".effScores.tab"

        # calculate eff scores and append them to the main file
        if doEffScoring:
            queue.startStep(batchId, "effScores", "Calculating guide efficiency scores")
            createMultiBatchEffScoreTable(
                batchId, tempEffScoresFname, queue, pam, pamFullName, extSeq, i, beFilter=(orgPamList, bePamIds)
            )  # simplify the scores table

        if isfile(tempEffScoresFname):
            with open(effScoresFnameTmp, "a") as mainScores:
                with open(tempEffScoresFname, "r") as tempScores:
                    shutil.copyfileobj(tempScores, mainScores)
            if not DEBUG:
                os.remove(tempEffScoresFname)

        logging.info("searching for off-targets with PAM %s" % pamFullName)
        startDict, endSet = findAllPams(uppSeq, pam)
        writePamFlank(seq, startDict, pam, faFname, pamFullName, beFilter=(orgPamList, bePamIds))

        # find offtargets and append them to the main file
        # if useBowtie:
        #    findOfftargetsBowtie(queue, batchId, iterBatchBase, faFname, genome, pam, tempBedFname)
        # else:
        findOfftargetsBwa(
            queue,
            batchId,
            iterBatchBase,
            faFname,
            genome,
            pam,
            tempBedFname,
            pamFullName,
        )

        if isfile(tempBedFname):
            with open(bedFnameTmp, "ab") as mainBed:
                with open(tempBedFname, "rb") as tempBed:
                    shutil.copyfileobj(tempBed, mainBed)

            if not DEBUG:
                os.remove(tempBedFname)

        if not DEBUG:
            if isfile(faFname):
                os.remove(faFname)


    # create the final file to end the job
    shutil.move(bedFnameTmp, bedFname)
    shutil.move(effScoresFnameTmp, effScoresFname)

    return bedFname, effScoresFname


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
    Record = namedtuple("tsvRec", headers)

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
            # raise Exception("wrong field count in line %s" % line)
            continue
        # convert fields to correct data type
        yield rec


allGenomes = None


def readGenomes():
    "return list of all genomes supported"
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
                addStr = "UCSC "
            else:
                addStr = ""
            genomes[row.name] = (
                row.scientificName
                + " - "
                + row.genome
                + " - "
                + addStr
                + row.description
            )

    genomes = list(genomes.items())
    genomes.sort(key=operator.itemgetter(1))
    allGenomes = genomes
    return allGenomes


def readAnnGenomes():
    "returns a list of all genomes with a genePred file"

    myDir = dirname(__file__)
    genomesDir = join(myDir, "genomes")

    annGenomes = []
    for subDir in os.listdir(genomesDir):
        subPath = join(genomesDir, subDir)
        if os.path.isdir(subPath):
            gpFiles = [f for f in os.listdir(subPath) if f.endswith(".gp")]
            if gpFiles:
                orgName = os.path.basename(subDir)
                annGenomes.append(orgName)

    return annGenomes


def printOrgDropDown(lastorg, genomes):
    "print the organism drop down box."
    print(
        '<select id="genomeDropDown" class style="width: 85%; max-width:600px;" name="org" tabindex="2">'
    )
    print("<option ")
    if lastorg == "noGenome":
        print("selected ")
    print(
        'value="noGenome">-- No Genome: no specificity, only cleavage efficiency scores (max. len 25kbp)</option>'
    )

    for db, desc in genomes:
        print("<option ")
        if db == lastorg:
            print("selected ")
        print('value="%s">%s</option>' % (db, desc))

    print("</select>")
    # print ('''
    # <script type="text/javascript">
    # $("#genomeDropDown").ufd({maxWidth:350, listWidthFixed:false});
    # </script>''')
    print("""<br>""")

    print(
        """<script>
    $("#genomeDropDown").chosen();
    $(".chosen-choices li").css("background","red");
    </script>
    """
    )


def dbsearchGene(params, onlySymbol=False, commonExons=False):
    """outputs gene IDs corresponding to the selected organism in json format
    optionally returns the selected gene symbol and its corresponding transcript.
    When commonExon is true, add "common exons" as an option for each gene symbol,
    with the gene symbol as value
    """
    print("Content-type: application/json\n")

    org = cgiGetStr(params, "org", "")
    term = cgiGetStr(params, "term", "").lower()

    if not org or not term:
        print(json.dumps({"results": []}))
        sys.exit(0)

    genomePath = join(genomesDir, org)

    gpFiles = [f for f in os.listdir(genomePath) if f.endswith(".gp")]
    if not gpFiles:
        print(
            json.dumps(
                {"results": [{"id": "0", "text": "No gene files for this organism"}]}
            )
        )
        sys.exit(0)

    matches = []

    # GenePred format : 0: gene name ; 1: chr ; 2: strand ; 3: TSS ; 4: TES ; 5: CDS start ; 6: CDS end ; 7: nb. exons ; 8: exons start ; 9: exons end ; 10: Score ; 11: alt Name.
    for gpFile in gpFiles:
        if "Select" in gpFile:
            isMane = True
        else:
            isMane = False
        with open(join(genomePath, gpFile), "r") as genePred:
            foundGenes = 0
            # stop the search at 100 matches for non specific terms (e.g LOC)
            for line in genePred:
                if foundGenes > 100:
                    break
                cols = line.strip().split("\t")
                mainId = cols[0]
                cdsStart = int(cols[5])
                cdsEnd = int(cols[6])
                if cdsEnd - cdsStart > 3:
                    coding = True
                else:
                    coding = False
                isAltName = len(cols) > 11
                if isAltName:
                    altId = cols[11]
                else:
                    altId = ""
                if len(cols) >= 14 and coding:
                    exFrames = cols[14]
                else:
                    exFrames = ""
                if isAltName and term in altId.lower() and onlySymbol is True:
                    foundGenes += 1
                    exonCount = cols[7]
                    matches.append((altId, mainId, exonCount, exFrames, coding, isMane))

    # select2 data formats : https://select2.org/data-sources/formats

    symDicts = []
    seen = set()
    # merge transcripts belonging to the same gene

    for match in matches:
        sym, transcript, exonCount, exFrames, coding, isMane = match
        if coding:
            codingStr = ""
        else:
            codingStr = ", non coding"
            exFrames = ""
        if isMane:
            matchText = "Mane select transcript for %s (%s - %s exons%s)" % (
                sym,
                transcript,
                exonCount,
                codingStr,
            )
        else:
            matchText = "%s (%s exons%s)" % (transcript, exonCount, codingStr)

        maneCount = 0
        if sym in seen:
            for symDict in symDicts:
                selSym = symDict["text"]
                transList = symDict["children"]
                if selSym == sym:
                    transDict = {
                        "id": transcript,
                        "text": matchText,
                        "exonCount": exonCount,
                        "exFrames": exFrames,
                    }
                    if commonExons and (len(transList) - maneCount) == 2:
                        transList.insert(
                            0,
                            {
                                "id": "%s~SYM" % sym,
                                "text": "Search exons common to all transcripts in %s"
                                % sym,
                                "exonCount": 0,
                                "exFrames": "",
                            },
                        )
                    # show the MANE select transcript at the top of the options
                    if isMane:
                        maneCount += 1
                        transList.insert(0, transDict)
                    else:
                        transList.append(transDict)

        else:
            symTransList = [
                {
                    "id": transcript,
                    "text": matchText,
                    "exonCount": exonCount,
                    "exFrames": exFrames,
                }
            ]

            symDicts.append({"text": sym, "children": symTransList})

            seen.add(sym)

    print(json.dumps({"results": symDicts[:30]}))
    sys.exit(0)


def printGeneSelection(paramName, onlySymbol=False, commonExons=False):
    """ "
    prints the searchable dropdown menu for gene/transcript selection.
    """

    if commonExons:
        ajaxStr = "geneSearchCommon"
    else:
        ajaxStr = "geneSearch"

    scriptName = basename(__file__)
    print(
        """
    <select class="js-select-gene" name="%s" id="geneSelection" style="width: 80%%;"></select>
    <input type="hidden" name="exonCount" id="exonCountVal">
    <br>
    <script>
    $(document).ready(function() {
        var gene_select = $('.js-select-gene');
        gene_select.select2({
            ajax: {
                url: '%s',
                dataType: 'json',
                delay: 250,
                data: function (params) {
                    return {
                        ajax: '%s',
                        term: params.term, // search term
                        org: $('#genomeDropDown').val() // get currently selected organism
                    };
                },
                processResults: function (data) {
                    return {
                        results: data.results
                    };
                },
                cache: true
            },
            placeholder: 'Type a gene symbol and select a corresponding transcript',
            minimumInputLength: 2,
        });

        $('#genomeDropDown').on('change', function() {
            // clear the selection when organism changes
            gene_select.val(null).trigger('change');
            $('#exonCountVal').val('');
        });

        gene_select.on('select2:select', function (e) {
            var data = e.params.data;
            if (data.exonCount) {
                $('#exonCountVal').val(data.exonCount);
            };
            if (data.exFrames) {
                $('#exFrames').val(data.exFrames);
            };
        });

    });
    </script>

    """
        % (paramName, scriptName, ajaxStr)
    )


def printPamDropDown(lastpam, name=None):

    if name is None:
        name = "pam"

    print(
        '<select class="js-example-basic-single" id="pamDropDown" style="float:left; width: 85%%;" name="%s" tabindex="3">'
        % name
    )
    for key, value in pamDesc:
        print("<option ")
        if key == lastpam:
            print("selected ")
        print('value="%s">%s</option>' % (key, value))
    print("</select>")
    print(
        """
    <script>
        $(document).ready(function() {
        $('.js-example-basic-single').select2();
        });
    </script>"""
    )


def printForm(params):
    "print html input form for classic mode"
    scriptName = basename(__file__)

    genomes = readGenomes()
    annGenomes = readAnnGenomes()

    haveHuman = False
    for g in genomes:
        if g[0] == "hg19":
            haveHuman = True

    # The returned cookie is available in the os.environ dictionary
    cookies = http.cookies.SimpleCookie(os.environ.get("HTTP_COOKIE"))
    if "lastorg" in cookies and "lastseq" in cookies and "lastpam" in cookies:
        lastorg = cookies["lastorg"].value
        lastseq = cookies["lastseq"].value
        lastpam = cookies["lastpam"].value
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

    if "exonCountVal" in params:
        exonCountVal = params["exonCountVal"]
    seqName = ""
    if "seqName" in params:
        seqName = params["seqName"]

    # print('''March 6 2023: Sorry, no CRISPOR on the new UCSC-based-server (with RS3 scores) today. Too many performance problems on the new server. We were able to renew the old server. Please use the <a href="http://37.187.154.234">old server</a> temporarily.''')
    # sys.exit(0)

    print(
        """
<form id="main-form" method="post" action="%s">

<div style="display:grid; clear:both; width: %%; min-width: 1500px; grid-template-columns: 42%% 58%%; grid-template-rows: auto auto; place-self:center; justify-self:center; space:20px; padding:12px;">

<div class="windowstep subpanel" style="width:100%%; grid-column:2; grid-row:-1/1;">
    <div class="substep">
        <div class="title">
            Step 3
        </div>


        Sequence name (optional): <input type="text" name="name" size="20" value="%s"><br>

        Enter a single genomic sequence, &lt; %d bp, typically an exon
        <img src="%simage/info-small.png" title="CRISPOR conserves the lowercase and uppercase format of your sequence, allowing to highlight sequence features of interest such as ATG or STOP codons.<br>Avoid using cDNA sequences as input, CRISPR guides that straddle splice sites are unlikely to work.<br>You can paste a single >23bp sequence and even multiple sequences, separated by N characters." class="tooltipster">
    <br>
    <small><a href="javascript:clearInput()">Clear Box</a> - </small>
    <small><a href="javascript:resetToExample()">Reset to default</a></small>
    </div>

    <textarea tabindex="1" style="width:98%%" name="seq" rows="12" autocorrect="off" spellcheck="false"
              placeholder="Paste here the genomic - not a cDNA - sequence of the exon you want to target. The sequence has to include the PAM site for your enzyme of interest, e.g. NGG. Maximum size %d bp. If you only have a cDNA, please BLAST or BLAT the cDNA first to find the right exon sequence for CRISPOR.">%s</textarea>
      <small>Text case is preserved, e.g. you can mark ATGs with lowercase.<br>Instead of a sequence, you can paste a chromosome range, e.g. chr1:11,130,540-11,130,751</small>

<details id="geneSelection" style = "margin-top:12px;">
    <summary>Click here to enter a gene ID and select a target exon instead</summary>
          """
        % (scriptName, seqName, MAXSEQLEN, HTMLPREFIX, MAXSEQLEN, lastseq)
    )

    printGeneSelection("koGeneId")

    print("""<div style="margin-top:12px; margin-bottom:12px">""")
    print(
        """<select class="js-select-hidden" name="exonSelect" id="exonSelect" style="width:40%%"> """
    )

    print(
        """
    </select>
<img src="%simage/info-small.png" title="Due to technical constraints, small exons (<23 bp) will be omitted, and large exons (>2300 bp) splitted into several smaller ones." class="tooltipsterInteract"><br>
    </div>
"""
        % HTMLPREFIX
    )
    print(
        """
    <small>Currently, %d out of %d genomes are annotated with genes. If your genome isn't included, paste a sequence above.</small>
    </details>
<div style="width:100%%; margin-top: 25px; margin-left:50px; text-align:center; display:block">
    <input type="submit" name="submit" value="SUBMIT" tabindex="4" style="height:40px; width:100px;"/>
</div>

</div>

<div class="windowstep subpanel" style="width:90%%; grid-column:1; grid-row:1;">
    <details id="classic1" open>
    <summary><small>Show / Hide step 1</small></summary>
    <div class="substep" style="margin-bottom: 1px">
        <div class="title" style="cursor:pointer;" onclick="$('#helpstep2').toggle('fast')">
            Step 1
        </div>
        Select a genome
    </div>

    """
        % (len(annGenomes), len(genomes))
    )

    printOrgDropDown(lastorg, genomes)

    print(
        """
    <div id="trackHubNote" style="margin-bottom:5px">
    <small>Note: pre-calculated exonic guides for this species are on the <a id='hgTracksLink' target=_blank href="">UCSC Genome Browser</a>.</small>
    </div>
    """
    )
    print(
        '<small style="float:left">We have %d genomes, but not yours? Search <a href="https://www.ncbi.nlm.nih.gov/assembly">NCBI assembly</a> and send a GCF_/GCA_ ID to <a href="mailto:%s">CRISPOR support</a>.</small>'
        % (len(genomes), contactEmail)
    )
    print(
        """
    </div>
    </details>
    <div class="windowstep subpanel" style="width:90%%; grid-column:1; grid-row:2;">
    <details id="classic2" open>
    <summary><small>Show / Hide step 2</small></summary>

    <div class="substep">
        <div class="title" style="cursor:pointer;" onclick="$('#helpstep3').toggle('fast')">
            Step 2
            <img src="%simage/info-small.png" title="The most common system uses the NGG PAM recognized by Cas9 from S. <i>pyogenes</i>. The VRER and VQR mutants were described by <a href='http://www.nature.com/nature/journal/vaop/ncurrent/abs/nature14592.html' target='_blank'>Kleinstiver et al</a>, Cas9-HF1 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4851738/'>Kleinstiver 2016</a>, eSpCas1.1 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4714946/'>Slaymaker 2016</a>, Cpf1 by <a href='http://www.cell.com/abstract/S0092-8674(15)01200-3'>Zetsche 2015</a>, SaCas9 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/25830891/'>Ran 2015</a> and KKH-SaCas9 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/26524662/'>Kleinstiver 2015</a>, modified As-Cpf1s by <a href='http://biorxiv.org/content/early/2016/12/04/091611'>Gao et al. 2017</a>." class="tooltipsterInteract">
        </div>
        Select a Protospacer Adjacent Motif (PAM)
    </div>
    """
        % HTMLPREFIX
    )

    printPamDropDown(lastpam)

    print(
        """<br>See <a target=_blank href="manual/manual.html#enzymes">notes on enzymes</a> in the manual.<br>
    <details id="customPAM" style = "margin-top:12px;">
        <input name = "customPAM" placeholder="PAM (at least 1 non-N, 3-8 nt)" style="height:22px; width:180px;" onkeydown="handleEnter(event)"></input>
        <select name = "customType" class="js-example-basic-single" style="width:25%%">
            <option value="">Select enzyme type</option>
            <option value="Cas9">Cas9</option>
            <option value="Cas12a">Cas12a (Cpf1)</option>
            <!--<option value="CBE">Cytosine base editor</option>-->
            <!--<option value="ABE">Adenine Base Editor (ABE)</option>-->
        </select>
        <input type="range" id="customGUIDELEN" name="customGUIDELEN" value="20" min="16" max="30" style="vertical-align:middle; width:15%%;" oninput="this.nextElementSibling.value = this.value">
        <output>20</output> nt guide
        </select>
        <summary>Using a Custom enzyme ? Click here to enter its specifications <img src="%simage/info-small.png" title="We cannot guarrantee the accuracy of efficiency scores for custom enzymes. If you want to add a new score or enzyme, please send an email at %s" class="tooltipsterInteract"> </summary>
    </details>

    </div>
    </details>
    </div>
          """
        % (HTMLPREFIX, contactEmail)
    )

    print(
        """
<script>
/* set the dropbox to hg19 and paste the example sequence into the input box. */
function resetToExample() {
    $("textarea[name='seq']").val("%s");
    $("#genomeDropDown").val("%s").trigger("chosen:updated");
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

<script>
    $(document).ready(function() {
    $('.js-example-basic-single').select2();
});
</script>

<script>
function handleEnter(event) {
    if (event.keyCode === 13) {
        event.preventDefault();
        event.target.blur();
    }
}
</script>

<script>
// get the exon number for the selected geneID
$(document).ready(function() {
    $('.js-select-gene').on('select2:select', function (e) {
        var data = e.params.data;
        var exonSelect = $('#exonSelect');
        exonSelect.empty();

        if (data.exonCount) {
            for (var i = 0; i < data.exonCount; i++) {
                j = i+1;

                var exonText = 'find guides for exon ' + j;
                var option = new Option(exonText, i, false, false);
                exonSelect.append(option);
            }
            exonSelect.trigger('change');
        }
    });
});
</script>

<script>
// select2 with hidden search box
$(".js-select-hidden").select2({
     minimumResultsForSearch: -1,
     placeholder: 'Select the exon to target'})
</script>

<script>
// save the states of detail elements on page reload
(function() {
    var $details = $('details[id]');
    $details.each(function() {
        var savedState = localStorage.getItem('details-' + this.id);
        if (savedState !== null) {
            this.open = savedState === 'true';
        }
    });

    $details.on('toggle', function() {
        localStorage.setItem('details-' + this.id, this.open);
    });
})();
</script>

</form>
    """
        % (DEFAULTSEQ, DEFAULTORG)
    )


def readBatchAsDict(batchId):
    "return contents of batch as a dictionary or None"
    batchBase = join(batchDir, batchId)
    jsonFname = batchBase + ".json"
    if isfile(jsonFname):
        params = json.load(open(jsonFname))
    else:
        # fix for the test version were the job archive is empty
        if not isfile(batchArchive) or os.path.getsize(batchArchive) == 0:
            return None
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


def writeBatchAsDict(batchInfo, batchId, editData=None):
    batchBase = join(batchDir, batchId)
    tmpFname = batchBase + ".json.tmp"

    if len(batchInfo) > 0:
        ofh = open(tmpFname, "w")
        json.dump(batchInfo, ofh)
        ofh.close()

        jsonFname = batchBase + ".json"
        os.rename(tmpFname, jsonFname)
        logging.debug("Wrote batch info to %s: %s" % (jsonFname, batchInfo))

    if editData is not None:

        editTmpFname = batchBase + ".editData.json.tmp"
        efh = open(editTmpFname, "w")
        json.dump(editData, efh)
        efh.close()

        editFname = batchBase + ".editData.json"
        os.rename(editTmpFname, editFname)
        logging.debug("Wrote edit data info to %s" % (editFname))


def readBatchParams(batchId):
    """given a batchId, return the genome, the pam, the input sequence and the
    chrom pos and extSeq, a 100bp-extended version of the input sequence.
    Returns None for pos if not found."""
    params = readBatchAsDict(batchId)
    if params != None:
        return (
            params.get("seq"),
            params["org"],
            params.get("pam"),
            params.get("posStr"),
            params.get("extSeq"),
            params.get("multiseq"),
            params.get("koMethod"),
            params.get("geneModel"),
            params.get("koGeneId"),
            params.get("multipam"),
        )

    # FROM HERE UP TO END OF FUNCTION: legacy cold for old batches pre-end-2016 (no json files back then)
    # remove in 2017
    batchBase = join(batchDir, batchId)
    inputFaFname = batchBase + ".input.fa"
    if not isfile(inputFaFname):
        errAbort(
            "Could not find the batch %s. We cannot keep Crispor runs for more than "
            "a few months. Please resubmit your input sequence via"
            ' <a href="crispor.py">the query input form</a>' % batchId
        )

    ifh = open(inputFaFname, encoding="utf8")
    ifhFields = ifh.readline().replace(">", "").strip().split()
    if len(ifhFields) == 2:
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

    return inSeq, genome, pamSeq, position, extSeq, None, None, None, None, None


def gzipStr(s):
    "compress a string with gzip and return"
    out = StringIO()
    with gzip.GzipFile(fileobj=out, mode="w") as f:
        f.write(s)
    return out.getvalue()


def gunzipStr(s):
    "uncompress a string with gzip and return"
    print(len(s), type(s), dir(s))
    f = gzip.GzipFile(StringIO(s))
    result = f.read()
    f.close()
    return result


def openDbm(dbFname, mode):
    "some distributions don't include the dbm module anymore"
    # import dbm.ndbm
    # dbMod = dbm
    # import dbm.gnu
    # dbMod = gdbm
    # import semidbm
    # lmdbm is faster than everything else: https://pypi.org/project/lmdbm/
    # though semidbm is not bad either
    # Also see leveldb Wiki page
    # import lmdb
    # import dbm.dumb
    from lmdbm import Lmdb

    db = Lmdb.open(dbFname + ".lmdb", mode)
    if mode == "c":
        # hack, I don't know how to set the permissions on the open call
        cmd = ["chmod", "-R", "a+rw", dbFname + ".lmdb"]
        runCmd(cmd, useShell=False)
    return db


def saveOutcomeData(batchId, data, seqNumber=None):
    """save outcome data of batch. data is a dictionary with key = score name
    if the function is called in a loop, merge files instead of overwriting
     seqNumber = index in loop"""
    batchBase = join(batchDir, batchId)
    dbFname = batchBase
    db = openDbm(dbFname, "c")

    # conn = sqlite3.connect(dbFname, "w")
    # c = conn.cursor()
    # c.execute('''CREATE TABLE outcomes (id text PRIMARY KEY, data blob))''' % scoreName)
    # c.commit()

    for scoreName, data in data.items():
        if scoreName in db:
            # if run in a loop, merge to the existing data
            existing_compressed = db[scoreName]
            existing_data = json.loads(
                zlib.decompress(existing_compressed).decode("utf8")
            )
            existing_data.update(data)
            db[scoreName] = zlib.compress(json.dumps(existing_data).encode("utf8"))
        else:
            # No existing data, just store the new data
            db[scoreName] = zlib.compress(json.dumps(data).encode("utf8"))

    db.close()

    # c.commit()


def readOutcomeData(batchId, scoreName):
    """open outcome data of batch, key is score name"""
    batchBase = join(batchDir, batchId)
    # conn = sqlite3.connect(dbFname, "r")
    # c = conn.cursor()
    # binData = c.execute("SELECT data FROM outcomes where id=?", scoreName)
    # try:
    #    import dbm.ndbm
    #    db = dbm.ndbm.open(batchBase, "r") # dbm always adds .db to the file name
    # except:
    #    # old batches on crispor.org are still using gdbm
    # dbFname = batchBase+".dbm"
    #    import dbm.gnu
    #    db = dbm.gnu.open(dbFname, "r")
    # import dbm.dumb
    # db = dbm.dumb.open(dbFname, "r") # dbm always adds .db to the file name
    db = openDbm(batchBase, "r")
    dbObj = db[scoreName]
    jsonStr = zlib.decompress(dbObj)
    data = json.loads(jsonStr)
    db.close()
    return data


def findAllPams(seq, pam, exonId=None):
    """find all matches for PAM and return as dict startPos -> strand and a set
    of end positions. The start positions for the negative strand are for the
    rev-complemented PAM
    """
    seq = seq.upper()
    startDict, endSet = findPams(seq, pam, "+", {}, set(), exonId)
    startDict, endSet = findPams(seq, revComp(pam), "-", startDict, endSet, exonId)

    if pam in multiPams:
        for pam2 in multiPams[pam]:
            startDict, endSet = findPams(seq, pam2, "+", startDict, endSet, exonId)
            startDict, endSet = findPams(
                seq, revComp(pam2), "-", startDict, endSet, exonId
            )

    return startDict, endSet


def newBatch(batchName, seq, org, pam):
    """obtain a batch ID and write seq/org/pam to their files.
    Return batchId.
    """
    batchId = makeTempBase(seq, org, pam, batchName)
    batchBase = join(batchDir, batchId)
    jsonFname = batchBase + ".json"

    if isfile(jsonFname):
        return batchId

    batchData = {}
    batchData["org"] = org
    batchData["pam"] = pam
    batchData["batchName"] = batchName
    batchData["seq"] = seq
    batchData["posStr"] = ""

    writeBatchAsDict(batchData, batchId)
    return batchId


def newMultiSeqBatch(
    batchName,
    multiseq,
    org,
    pam,
    koMethod=None,
    geneModel=None,
    koGeneId=None,
    assist=None,
    exonSelect=None,
):
    """obtain a batch ID and write essential params to their files.
    Return batchId.
    """

    allSeq = "".join([seq[1] for seq in multiseq])
    # add exonSelect and the selected geneId to get a unique batchId
    if exonSelect is not None:
        allSeq += exonSelect
    allSeq += koGeneId
    batchId = makeTempBase(allSeq, org, pam, batchName)
    batchBase = join(batchDir, batchId)
    jsonFname = batchBase + ".json"

    if isfile(jsonFname):
        return batchId

    else:
        batchData = {}
        batchData["org"] = org
        batchData["pam"] = pam
        batchData["batchName"] = batchName
        # batchData["seq"] = seq
        batchData["posStr"] = ""
        batchData["multiseq"] = multiseq
        batchData["koMethod"] = koMethod
        batchData["geneModel"] = geneModel
        batchData["koGeneId"] = koGeneId
        batchData["assist"] = assist
        batchData["exonSelect"] = exonSelect

        writeBatchAsDict(batchData, batchId)

        return batchId


def newMultiPamBatch(
    batchName,
    seq,
    posStr,
    org,
    multipam,
    insertSeq,
    insertIdx,
    assist=None,
    koGeneId=None,
    insertPos=None,
    kiType=None,
    tagNames=None,
    geneModel=None,
    nonCoding=None,
    clippedSeq=None,
    noPerfectMatch=None
):

    if seq:
        concatSeq = seq + insertSeq + str(insertIdx)
    else:
        concatSeq = posStr + insertSeq + str(insertIdx)

    batchId = makeTempBase(concatSeq, org, multipam, batchName)
    batchBase = join(batchDir, batchId)
    jsonFname = batchBase + ".json"

    if isfile(jsonFname):
        return batchId
    else:
        batchData = {}
        batchData["seq"] = seq
        batchData["posStr"] = posStr
        batchData["org"] = org
        batchData["multipam"] = multipam
        batchData["insertSeq"] = insertSeq
        batchData["insertIdx"] = insertIdx
        batchData["batchName"] = batchName
        batchData["assist"] = assist
        if nonCoding:
            batchData["nonCoding"] = nonCoding
        if koGeneId:
            batchData["koGeneId"] = koGeneId
        if insertPos:
            batchData["insertpos"] = insertPos
        if kiType:
            batchData["kiType"] = kiType
        if tagNames:
            batchData["tagNames"] = tagNames
        if geneModel:
            batchData["geneModel"] = geneModel
        if clippedSeq:
            batchData["clippedSeq"] = clippedSeq
        if noPerfectMatch:
            batchData["noPerfectMatch"] = noPerfectMatch

        writeBatchAsDict(batchData, batchId)

        return batchId


def readDbInfo(org):
    "return a dbInfo object with the columsn in the genomeInfo.tab file"
    myDir = dirname(__file__)
    genomesDir = join(myDir, "genomes")
    infoFname = join(genomesDir, org, "genomeInfo.tab")
    if not isfile(infoFname):
        return None
    dbInfo = next(lineFileNext(open(infoFname)))
    return dbInfo


def printQueryNotFoundNote(dbInfo):
    print(
        "<div class='title'>Query sequence, not found in the selected genome, %s (%s)</div>"
        % (dbInfo.scientificName, dbInfo.name)
    )
    print(
        "<div class='substep' style='border: 1px black solid; padding:5px; background-color: aliceblue'>"
    )
    print(
        "<strong>Warning:</strong> The query sequence was not found in the selected genome."
    )
    print("This can be a valid query, e.g. a GFP sequence.<br>")
    print(
        "If not, you might want to check if you selected the right genome for your query sequence.<br>"
    )
    print(
        "Use a tool like <a target=_blank href='http://genome.ucsc.edu/cgi-bin/hgBlat'>BLAT</a> to check if the "
        "sequence really has a 100% identical match in the target genome.<p>"
    )
    print(
        "When reading the list of guide sequences and off-targets below, bear in mind that in case that the input sequence is really in the genome and just has a few differences, the software will use the first found match as the on-target as it cannot distinguish 0-mismatch off-targets from 0-mismatch on-targets. In this case, the specificity scores of guide sequences are too low. In other words, some guides may be fine, the problem may just be that the on-target is shown as an off-target. <br>"
    )
    print(
        "Because there is no flanking sequence available, the guides in your sequence that are within 50bp of the ends will have no efficiency scores. The efficiency scores will instead be shown as '--'. Include more flanking sequence > 50bp to obtain these scores."
    )
    print("</div>")


def submitMultiSearch(batchId, org, pamDesc, mode):
    "sends a job to the worker to process multiple PAMs or seqs"

    if mode == "multiseq" or mode == "multipam":
        batchBase = join(batchDir, batchId)
        otBedFname = batchBase + ".bed.gz"
        effScoresFname = batchBase + ".effScores.tab"

        q = JobQueue()
        q.openSqlite()
        ip = os.environ.get("REMOTE_ADDR", "noIp")

        wasOk = q.addJob(mode, batchId, "ip=%s,org=%s,pam=%s" % (ip, org, pamDesc))

        if not wasOk:
            print("CRISPOR job %s failed-running..." % batchId)
            pass
        q.close()

        if isfile(otBedFname) and isfile(effScoresFname):
            return otBedFname, effScoresFname


def getOfftargets(seq, org, pamDesc, batchId, startDict, queue):
    """write guides to fasta and run bwa or use cached results.
    Return name of the BED file with the matches or None if not yet available.
    Write progress status updates to queue object.
    """
    pam = setupPamInfo(pamDesc)
    assert "-" not in pam

    batchBase = join(batchDir, batchId)
    otBedFname = batchBase + ".bed.gz"

    batchInfo = readBatchAsDict(batchId)

    flagFile = batchBase + ".running"
    if isfile(flagFile):
        errAbort(
            "This sequence is still being processed. Please wait for ~20 seconds "
            "and try again, e.g. by reloading this page. If you see this message for "
            "more than 2-3 minutes, please send an email to %s. Thanks!" % contactEmail
        )

    if (
        not batchInfo
        or not isfile(otBedFname)
        or commandLineMode
        or not "posStr" in batchInfo
        or (batchInfo["posStr"] == "" and not batchInfo["org"] == "noGenome")
    ):  # pre-4.8 batches don't have a posStr at all
        # write potential PAM sites to file
        faFname = batchBase + ".fa"
        writePamFlank(seq, startDict, pam, faFname)
        if commandLineMode:
            processSubmission(
                faFname, org, pamDesc, otBedFname, batchBase, batchId, queue
            )
        else:
            q = JobQueue()
            q.openSqlite()
            ip = os.environ.get("REMOTE_ADDR", "noIp")

            if ip == "195.176.112.240":
                errAbort("IP address blocked.")

            wasOk = q.addJob(
                "search", batchId, "ip=%s,org=%s,pam=%s" % (ip, org, pamDesc)
            )
            if not wasOk:
                print("CRISPOR job %s failed-running..." % batchId)
                pass
            q.close()
            return None

    return otBedFname


def showPamWarning(pam):
    if pamIsCpf1(pam):
        print(
            '<div style="text-align:left; border: 1px solid; background-color: aliceblue; padding: 3px">'
        )
        print("<strong>Note:</strong> You are using the Cpf1 enzyme or related enzyme.")
        print(
            "While there is an efficiency score specificially for Cpf1, there is no off-target ranking algorithm available in the literature, to our knowledge. We use Hsu and CFD scores below for off-target ranking, but they were developed for spCas9. There is not enough data yet to support their usefulness for Cpf1. Contact us for more info if you need to rank Cpf1 off-targets for validation or if you have a dataset that could elucidate this question. We are showing out-of-frame scores, but they are based on micro-homology that assumes a spCas9 cut site, so most likely the out-of-frame scores are not accurate for the staggered cut of Cpf1 either."
        )
        print("</div>")
    # elif pamIsSaCas9(pam):
    # print '<div style="text-align:left; border: 1px solid; background-color: aliceblue; padding: 3px">'
    # print "<strong>Note:</strong> Your query is using a Cas9 from S. aureus.<br>"
    # print "Please note that while the efficiency scoring was built for saCas9, the off-target ranking below and specificity scores are based on CFD/Hsu models, which were developed for spCas9. The ranking of off-targets could be very inaccurate. If you have a saCas9 off-target dataset, you can contact us for further info, we are only aware of the BLESS dataset by <a href='https://www.nature.com/articles/nature14299' target=_blank>Ran et al. 2015</a>.<br>As for out-of-frame and micro-homology, this model is also based on spCas9, but <a target=_blank href='https://www.nature.com/articles/nature14299'>Ran et al 2015</a> showed that the saCas9 cleavage pattern looks identical to spCas9's, so the OOF micro-homology model should work with saCas9."
    # print '</div>'
    elif not pamIsSpCas9(pam) and not pamIsSaCas9(pam):
        print(
            '<div style="text-align:left; border: 1px solid; background-color: aliceblue; padding: 3px">'
        )
        print(
            "<strong>Warning:</strong> Your query involves a Cas9 that is not from S. Pyogenes and is also not Cpf1 nor saCas9."
        )
        print(
            "Please bear in mind that specificity and efficiency scores were designed using data with S. Pyogenes Cas9 and will very likely not be applicable to this particular Cas9. There is nothing we can do about this, we are unaware of a published dataset for this enzyme. If you know one, please contact us. Also contact us if you think another one of the existing scoring model would be more appropriate for this enzyme.<br>"
        )
        print("</div>")

    if pam == "NNNNACA":
        printNote(
            "You selected the old version of the CjCas9 PAM. You may want to select the more recent "
            + "PAMs from the menu on the first page, based on the study by "
            + "<a target=_blank href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5473640/'>Kim et al 2017</a>."
        )

    if pam == "NGN":
        printNote(
            "You have selected the NGN pam for xCas9. While this PAM is documented to work, "
            + "if you read the paper in detail, you will notice that the editing efficiency is much lower. "
            + "For optimal efficiency, consider going back and switching to the 'high-efficiency' xCas9 PAM."
        )

    if pam == "NGK":
        printNote(
            "You have selected the most efficient PAM for xCas9. You can also select the more general/flexible"
            + " NGN PAM from the menu when you submit your job. If you read the xCas9 paper in detail, you will find "
            + "that NGN is not as efficient though."
        )

    if pam == "TTTN":
        printWarning(
            "You selected TTTN as the PAM for Cpf1. "
            + "This is not the best PAM. The actual PAM is TTTV, as shown in Fig. 2a of "
            + "<a href='https://www.ncbi.nlm.nih.gov/pubmed/27992409'>Kim HK et al. Nat Meth 2017</a>.<br>"
        )


def showNoGenomeWarning(dbInfo):
    if dbInfo == None:
        printNote(
            'As there is no genome that can be used to get flanking sequence for your sequence, efficiency scores 50bp from the start or the end of your sequence cannot be calculated and are shown as "--". If needed, extend the input sequence and retry.'
        )


def getSeq(db, posStr, maxlen=True, minlen=True):
    """
    given a database name and a string with the position as chrom:start-end, return the sequence as
    a string.
    """
    chrom, start, end, strand = parsePos(posStr)

    if end - start > MAXSEQLEN and db != "noGenome" and maxlen is True:
        errAbort(
            "Input sequence range too long. Please retry with a sequence range shorter than %d bp."
            % MAXSEQLEN
        )
    genomeDir = genomesDir  # pull in global var
    twoBitFname = getTwoBitFname(db)
    binPath = join(binDir, "twoBitToFa")

    chromSizes = parseChromSizes(db)
    if chrom not in chromSizes:
        errAbort(
            "Sorry, the chromosome '%s' is not valid in the genome %s. Check upper/lowercase, e.g. for most mammalian genomes, "
            "it is chrX not chrx, and chr1, not Chr1." % (html.escape(chrom), db)
        )
    if start < 0 or end < 0 or start > chromSizes[chrom] or end > chromSizes[chrom]:
        errAbort(
            "Sorry, the coordinates '%d-%d' are not valid in the genome %s. Coordinates must not be outside chromosome boundaries or less than 0."
            % (start, end, db)
        )

    cmd = [
        binPath,
        twoBitFname,
        "-seq=" + chrom,
        "-start=" + str(start),
        "-end=" + str(end),
        "stdout",
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    seqStr = proc.stdout.read()

    retCode = proc.wait()
    if retCode != 0:
        errAbort(
            "Error on sequence retrieval. This looks like a bug. Please contact us and tell us the input and genome, we will fix this error."
        )

    # remove fasta header line
    lines = seqStr.decode("utf8").splitlines()
    if len(lines) > 0:
        lines.pop(0)
    seq = "".join(lines)
    if minlen and len(seq) < 23:
        errAbort(
            "Sorry, the sequence range %s on genome %s is not longer than 23bp. To find a valid CRISPR/Cas9 site, one needs at least a 23bp long sequence."
            % (db, posStr)
        )

    if strand == "-":
        seq = revComp(seq)

    return seq


def printStatus(batchId, msg):
    "print status, not using any Ajax"
    q = JobQueue()
    q.openSqlite()
    status = q.getStatus(batchId)
    q.close()

    errorState = False
    print(status)  # temporary, for debugging
    if "Traceback" in status:
        print("<!--")
        print(status)
        print("-->")
        status = (
            "An error occured during the processing.<br> Please send an email to %s and tell us that the failing batchId was %s.<br>We can usually fix this quickly. Thanks! <br>If you submit the same sequence/genome/name again, it will not be re-run, Crispor will pickup the old error. We will have to reset it before you can resubmit this particular sequence, so you will have to contact us or change the sequence to get a new job into the system."
            % (contactEmail, batchId)
        )
        errorState = True
    else:
        print('<meta http-equiv="refresh" content="10" >')
        """
        if len(msg) != 0:
            print((msg + "<p>"))
        """
        print("CRISPOR job has been submitted.<p>")

    if status == None:
        status = "Batch completed. Refresh page to show results."

    print(("Job Status: <tt>%s</tt><p>" % status))

    if not errorState:
        print("<p><small>This page will refresh every 10 seconds</small><br>")
        print(
            (
                "<p><small>If you see this message for longer than 5 minutes, please <a href='mailto:%s'>contact us</a>."
                % contactEmail
            )
        )


def readVarDbs(db):
    """find all possible variant VCFs and return as list of (shortLabel, fname, label, hasAF)
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
            if len(fields) == 4:
                shortLabel, fname, desc, hasAF = fields
            else:
                errAbort("not four fields in vcfDescs.txt: %s" % fields)

            fpath = join(genomesDir, db, fname)
            if not isfile(fpath):
                print("Error: Cannot find VCF file %s" % fpath)
                continue
            hasAF = hasAF == "1"
            ret.append((shortLabel, fname, desc, hasAF))
    return ret


def parseVcfInfo(info):
    "parse a VCF info string and return as a dict key=val"
    fs = info.split(";")
    ret = {}
    for field in fs:
        parts = field.split("=", maxsplit=1)
        if len(parts) == 2:
            key, val = parts
        elif len(parts) == 1:
            key = parts[0]
            val = True
        ret[key] = val
    return ret


def findVariantsInRange(vcfFname, chrom, start, end, strand, minFreq):
    """find variants that overlap the position.
    varDb is a tuple of label, vcfFname
    return as a dict relative position -> (chrom, pos, refAllele, altAllele, list of info-dicts)
    special position is "label" which is the label of the variant db
    """
    minFreq = float(minFreq)
    seqLen = end - start
    if not isfile(vcfFname):
        errAbort("%s not found" % vcfFname)
    tb = tabix.open(vcfFname)
    chrom = chrom.replace("chr", "")
    try:
        records = tb.query(chrom, start + 1, end)  # VCF is 1-based
    except tabix.TabixError:
        sys.stderr.write(
            "Chromosome in query does not exist in VCF file? chrom: %s, VCF file: %s\n"
            % (chrom, vcfFname)
        )
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
            # afList = infoDict["AF"].split(",")
            # altAllList = altAll.split(",")
            # newAltAllList = []
            # newAfList = []
            # for af, altAll in zip(afList, altAllList):
            # afNum = float(af)
            # if afNum < minFreq:
            # continue
            # newAltAllList.append(altAll)
            # newAfList.append(af)
            ##if len(newAltAllList)==0:
            # continue
            # altAll = ",".join(newAltAllList)
            # infoDict["AF"] = ",".join(newAfList)
            if allFreq != None:
                attribs["freq"] = allFreq
            relPos = int(varPos) - 1 - start
            if strand == "-":
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
    "show a little dropdown menu so user can get annotated sequence in genbank format"
    print("""<div style="padding-top:4px"><small>Download for: """)

    htmls = []

    baseUrl = "crispor.py?batchId=%s" % batchId

    myUrl = baseUrl + "&download=serialcloner"
    html = (
        "<a href='%s'>SerialCloner</a> (<a target=_blank href='http://serialbasics.free.fr/Serial_Cloner-Download.html'>free</a>)"
        % myUrl
    )
    htmls.append(html)

    myUrl = baseUrl + "&download=ape"
    html = (
        '<a href="%s">ApE</a> (<a target=_blank href="http://biologylabs.utah.edu/jorgensen/wayned/ape/">free</a>)'
        % myUrl
    )
    htmls.append(html)

    myUrl = (
        "http://crispor.tefor.net/crispor.py?batchId=%s&download=genomecompiler"
        % batchId
    )
    # backUrl = "https://designer.genomecompiler.com/plasmid_iframe?file_url=%s#/plasmid" % urllib.quote(myUrl)
    backUrl = (
        "https://designer.genomecompiler.com/plasmid_iframe?file_url=%s#/plasmid"
        % urllib.parse.quote(myUrl)
    )
    html = "<a target=_blank href='%s'>GenomeCompiler</a>" % backUrl
    htmls.append(html)

    myUrl = baseUrl + "&download=benchling"
    html = "<a href='%s'>Benchling</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl + "&download=snapgene"
    html = "<a href='%s'>SnapGene</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl + "&download=geneious"
    html = "<a href='%s'>Geneious</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl + "&download=vnti"
    html = "<a href='%s'>Vector NTI</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl + "&download=lasergene"
    html = "<a href='%s'>LaserGene</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl + "&download=genbank"
    html = "<a href='%s'>Genbank</a>" % myUrl
    htmls.append(html)

    myUrl = baseUrl + "&download=fasta"
    html = "<a href='%s'>FASTA</a>" % myUrl
    htmls.append(html)

    # html = "<div id='copyLink' data-clipboard-target='#seqAsText'>Copy sequence to clipboard</div>"
    # htmls.append(html)

    print(" - ".join(htmls))

    print("</small></div>")


def mapToGenome(seqStart, seqStrand, pamStart, guideStart, guideStrand):
    if pamIsFirst:
        # thick part = PAM comes first
        chromStart = seqStart + pamStart
        thickStart = seqStart + guideStart
        thickEnd = thickStart + GUIDELEN
        chromEnd = thickEnd
    else:
        chromStart = seqStart + guideStart
        thickStart = chromStart
        thickEnd = chromStart + GUIDELEN
        chromEnd = thickEnd

    strands = seqStrand + guideStrand
    chromStrand = "+"
    if strands == "+-" or strands == "-+":
        chromStrand = "-"

    return chromStart, chromEnd, thickStart, thickEnd, chromStrand


def makeCustomTrack(
    org, chrom, seqStart, seqEnd, seqStrand, guideData, batchId, batchName
):
    "create a custom track file for a given batch and return the filename"
    ctDir = join(batchDir, "customTracks")
    if not isdir(ctDir):
        os.makedirs(ctDir)

    baseUrl = join(ctDir, batchId)
    ctFname = join(ctDir, batchId + ".bed")  # temporary bed file
    bbFname = join(ctDir, batchId + ".bb")  # bigBed file
    ctFname = join(ctDir, batchId + ".txt")  # custom track settings
    if isfile(ctFname):
        return bbFname

    seqStart = int(seqStart)
    seqEnd = int(seqEnd)

    rows = []
    for guideRow in guideData:
        (
            guideScore,
            guideCfdScore,
            effScores,
            pamStart,
            guideStart,
            guideStrand,
            pamId,
            guideSeq,
            pamSeq,
            otData,
            otDesc,
            last12Desc,
            mutEnzymes,
            ontargetDesc,
            repCount,
            gcFrac,
            freeEnergy,
            doRecoding,
            cutUpstream,
            mainScore,
            beScoring,
            beOutcomes
        ) = guideRow

        rgb = hexToRgb(scoreToColor(guideScore)[0])
        chromStart, chromEnd, thickStart, thickEnd, chromStrand = mapToGenome(
            seqStart, seqStrand, pamStart, guideStart, guideStrand
        )

        mitScore = str(guideScore)
        fusiScore = str(effScores.get("fusi", -1))
        crisprScanScore = str(effScores.get("crisprScan", -1))
        oofScore = str(effScores.get("oof", -1))

        name = pamId
        bed = [
            chrom,
            chromStart,
            chromEnd,
            name,
            mitScore,
            chromStrand,
            thickStart,
            thickEnd,
            rgb,
            guideSeq,
            pamSeq,
            mitScore,
            fusiScore,
            crisprScanScore,
            oofScore,
            batchId,
        ]
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
    cmd = [
        "$BIN/bedToBigBed",
        "-type=bed9+",
        "-tab",
        "-as=" + asFname,
        ctFname,
        sizeFname,
        bbFname,
    ]
    runCmd(cmd, useShell=False)

    bbUrl = baseUrl + "/%s.bb" % batchId

    ofh = open(ctFname, "w")
    if batchName == "":
        batchName = "Results"
    ofh.write("browser position %s:%d-%d\n" % (chrom, seqStart, seqEnd))
    ofh.write(
        'track type=bigBed name="CRISPOR %(batchName)s" description="CRISPOR Results %(batchName)s %(batchId)s" bigDataUrl=%(bbUrl)s itemRgb=On visibility=pack\n'
        % locals()
    )
    ofh.close()

    hubFname = join(ctDir, batchId+".txt")
    ofh = open(hubFname, "w")
    ofh.write("hub CRISPOR\n")
    ofh.write("shortLabel CRISPOR %s\n" % org)
    ofh.write("longLabel CRISPOR batch %s %s\n" % (org, batchId))
    ofh.write("genomesFile %s\n" % hubFname)
    ofh.write("email crispor@tefor.net\n")
    ofh.write("descriptionUrl http://crispor.org\n")
    ofh.write("genome %s\n" % org)
    ofh.write("trackDb %s\n" % hubFname)
    ofh.write("track crispor%s\n" % batchId)
    # ofh.write("shortLabel M-CAP
    # longLabel M-CAP
    # group genes
    # visibility dense
    # type bigBed 9 +
    ##os.remove(ctFname)

    # return bbFname
    ctUrl = baseUrl + "/%s.txt" % batchId
    return ctUrl, bbUrl


def iterBbLines(bbPath, chrom, start, end, strand):
    "yield bigGenePred rows from bigBed that overlap pos"
    binPath = join(binDir, "bigBedToBed")
    cmd = [
        binPath,
        bbPath,
        "stdout",
        "-chrom=" + chrom,
        "-start=" + str(start),
        "-end=" + str(end),
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, encoding="latin1")
    for line in proc.stdout:
        yield line.split("\t")


def trimExonAndFlip(exStart, exEnd, exStrand, seqLen, seqStrand, exFrame=None):
    """Put the exon into the current sequence window:
    - trim exon to the window (0, seqLen), return None if completely outside the view.
    - reverse the exon coordinates if seqStrand=="-"
    """
    # Flip first
    if seqStrand == "-":
        exStart, exEnd = seqLen - exEnd, seqLen - exStart
        exStrand = "+" if exStrand == "-" else "-"

    # save the unadjusted exFrame, to check if removing the exon destroys the reading frame (even after trimming)
    if exFrame is not None:
        oldExFrame = exFrame
    else:
        oldExFrame = None

    # Trim start
    if exStart < 0:
        if exEnd < 0:
            # the whole exon is outside the view on the left side
            return None, None, None, None, None
        else:
            # if the sequence start inside en exon, the frame can be truncated
            if exFrame is not None and exFrame != -1:
                # take into account the phase shift due to the missing part of the exon
                # use a positive value for the missing length!
                exOffset = abs(exStart) % 3
                exFrame = (exFrame + exOffset) % 3
            exStart = 0

    # Trim end
    if exEnd > seqLen:
        if exStart > seqLen:
            # the whole exon is outside the view on the right side
            return None, None, None, None, None
        else:
            # truncate the end
            exEnd = seqLen

    return exStart, exEnd, exFrame, oldExFrame, exStrand


def getExonInfo(org, geneName, position, extendPos=None):
    """retrieve exon info between position, return format transId -> (exNumber, start, end, exFrame, nextExonFrame, strand)
    - start and end are relative to position!
    - exNumber starts at 0
    - nextExonFrame is relative to strand: for a transcript on the - strand, it is for the more 5' exon
    """
    # bigGenePred format:
    # string chrom;       "Reference sequence chromosome or sca
    # uint   chromStart;  "Start position in chromosome"
    # uint   chromEnd;    "End position in chromosome"
    # string name;        "Name or ID of item, ideally both hum
    # uint score;         "Score (0-1000)"
    # char[1] strand;     "+ or - for strand"
    # uint thickStart;    "Start of where display should be thi
    # uint thickEnd;      "End of where display should be thick
    # uint reserved;       "RGB value (use R,G,B string in inpu
    # int blockCount;     "Number of blocks"
    # int[blockCount] blockSizes; "Comma separated list of bloc
    # int[blockCount] chromStarts; "Start positions relative to
    # string name2;       "Alternative/human readable name"
    # string cdsStartStat; "enum('none','unk','incmpl','cmpl')"
    # string cdsEndStat;   "enum('none','unk','incmpl','cmpl')"
    # int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame
    # string type;        "Transcript type"
    # string geneName;    "Primary identifier for gene"
    # string geneName2;   "Alternative/human readable gene name
    # string geneType;    "Gene type"

    ret = defaultdict(list)
    seqChrom, seqStart, seqEnd, seqStrand = parsePos(position)

    # in KO mode, position is extended by extendPos bases 
    if extendPos is not None and not pamIsFirst:
        seqStart -= extendPos
        seqEnd += extendPos

    seqLen = seqEnd - seqStart
    fname = join(genomesDir, org, geneName + ".bb")

    maxIdLen = 0

    for row in iterBbLines(fname, seqChrom, seqStart, seqEnd, seqStrand):
        (
            chrom,
            chromStart,
            chromEnd,
            name,
            score,
            strand,
            thickStart,
            thickEnd,
            reserved,
            blockCount,
            blockSizes,
            blockStarts,
            name2,
            cdsStartStat,
            cdsEndStat,
            exonFrames,
            tType,
            geneName,
            geneName2,
            geneType,
        ) = row

        chromStart = int(chromStart)
        chromEnd = int(chromEnd)
        thickStart = int(thickStart)
        thickEnd = int(thickEnd)

        blockSizes = [int(x) for x in blockSizes.split(",") if x != ""]
        blockStarts = [int(x) for x in blockStarts.split(",") if x != ""]
        exonFrames = [int(x) for x in exonFrames.split(",") if x != ""]
        assert len(blockSizes) == len(blockStarts) == len(exonFrames)

        symbol = ""
        if geneName2 != "":
            symbol = geneName2

        blocks = list(zip(blockSizes, blockStarts, exonFrames))

        for exIdx, (blockSize, blockStart, exonFrame) in enumerate(blocks):
            exChromStart = chromStart + blockStart
            exChromEnd = exChromStart + blockSize
            exStrand = strand
            nextFrame = None
            if strand == "+":

                if exIdx + 1 < len(blocks):
                    nextFrame = blocks[exIdx + 1][-1]
            else:
                if exIdx > 0:
                    nextFrame = blocks[exIdx - 1][-1]
                # for the - strand, count the exons backwards
                exIdx = (len(blocks) - 1) - exIdx

            # print("next frame", nextFrame, "<br>")

            # figure out exon start/end: special case for UTRs: trim down to CDS start/end
            if exChromStart < thickStart < exChromEnd:
                exStart = thickStart - seqStart
                if strand == "+":
                    exonFrame = 0  # Segments starting at translation start always have frame 0 (only for the + strand)
                # add the UTR as a special exon
                utrStart, utrEnd, _, _, utrStrand = trimExonAndFlip(
                    exChromStart - seqStart, exStart, exStrand, seqLen, seqStrand
                )
                if utrStart is not None:
                    ret[(name, symbol)].append(
                        (-1, utrStart, utrEnd, -1, -1, nextFrame, utrStrand)
                    )
            else:
                exStart = chromStart + blockStart - seqStart
            if exChromStart < thickEnd < exChromEnd:
                if strand == "-":
                    exonFrame = 0

                exEnd = thickEnd - seqStart
                # add the UTR as a special exon
                utrStart, utrEnd, _, _, utrStrand = trimExonAndFlip(
                    exEnd, exChromEnd - seqStart, exStrand, seqLen, seqStrand
                )
                if utrStart is not None:
                    ret[(name, symbol)].append(
                        (-1, utrStart, utrEnd, -1, -1, nextFrame, utrStrand)
                    )
            else:
                exEnd = chromStart + blockStart + blockSize - seqStart
            # print chromStart, chromStart+blockSize, seqStart, "<br>"
            exStart, exEnd, exonFrame, oldExonFrame, exStrand = trimExonAndFlip(
                exStart, exEnd, exStrand, seqLen, seqStrand, exFrame=exonFrame
            )
            if exStart is None:
                # whole exon is outside of view
                continue
            symbol = ""
            if geneName2 != "":
                symbol = geneName2
            ret[(name, symbol)].append(
                (exIdx, exStart, exEnd, exonFrame, oldExonFrame, nextFrame, exStrand)
            )
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

    if pamIsXCas9(pam) and org == "noGenome":
        errAbort(
            "You selected no genome, so only efficiency scoring is active. "
            "You also selected the enzyme xCas9. "
            "This does not work, since no efficiency score has been published yet for xCas9. "
            + "Please contact us if you think there is a xCas9 scoring model available somewhere. Thanks!"
        )

    return minFreq, varDb


def crisprSearch(params):
    "do crispr off target search and eff. scoring"
    if "org" in params:
        db = params["org"]
        twoBitFname = getTwoBitFname(db)
        if not isfile(twoBitFname) and db != "noGenome":
            errAbort(
                "Sorry, a genome assembly called %s is not on Crispor "
                "yet or not anymore. "
                "Please send us an email if you want us to add it." % db
            )

    # retrieve sequence if not provided
    if "pos" in params and not "seq" in params:
        params["seq"] = getSeq(params["org"], params["pos"])

    if "batchId" in params and params["batchId"] is not None:
        # if we're getting only the batchId, extract the parameters from the batch
        # this allows a stable link to a batch that is done
        batchId = params["batchId"]
        (
            seq,
            org,
            pamDesc,
            posStr,
            _,
            multiseq,
            koMethod,
            geneModel,
            koGeneId,
            multipam,
        ) = readBatchParams(batchId)
        if multiseq:
            params["multiseq"] = multiseq
        if multipam:
            params["multipam"] = multipam

        # pamDesc can include additional options, like guidelen and base editor
        # added after the pam, e.g. "NGG-BE1". setupPamInfo(pam) will set the globals
        # based on it
        if seq:
            seq, warnMsg = cleanSeq(seq, org)
        else:
            warnMsg = ""
    else:
        if (
            "customPAM" in params
            and params["customPAM"].strip()
            and "customType" in params
            and "customGUIDELEN" in params
        ):
            pamSeq = params["customPAM"].upper()
            ezType = params["customType"]
            guideLen = params["customGUIDELEN"]

            if len(pamSeq) < 3 or len(pamSeq) > 8:
                print("Please use a PAM sequence between 3nt and 8nt long")
                return

            legalChars = ["N", "A", "T", "G", "C"]  # add Y / R / V ?
            ncount = 0
            for char in pamSeq:
                if char not in legalChars:
                    print("Please use only ATGCN nucleotides in the PAM sequence")
                    return

                if char != "N":
                    ncount += 1
            if ncount < 1:
                print("Please use at least one non-N nucleotide in the PAM sequence")
                return

            pamDesc = "%(pamSeq)s.%(ezType)s.%(guideLen)s.custom" % locals()
        else:
            pamSeq = None
        # this is a new sequence: create a new batch (TODO: and add it to the queue?)

        org = params["org"]
        newBatchName = params.get("name", "")
        multipam = params.get("multipam")
        donorSeq = params.get("donorSeq")

        koGeneId = params.get("koGeneId")
        geneModel = params.get("geneModel")
        # if a new search has been launched via the KI results page, geneModel is in json format
        if isinstance(geneModel, str):
            geneModel = json.loads(geneModel)

        koMethod = params.get("koMethod")
        multiseq = params.get("multiseq")
        exonSelect = params.get("exonSelect")

        if multipam or multiseq:
            if multipam:
                assist = params["assist"]
                seq = params["seq"]
                posStr = params.get("pos")
                insertPos = params.get("insertpos")
                insertSeq, warnMsg = cleanSeq(params["insertSeq"], org)
                insertIdx = int(params["insertIdx"])
                kiType = params.get("kiType")
                tagNames = params.get("tagNames")
                nonCoding = params.get("nonCoding")
                clippedSeq = params.get("clippedSeq")

                # allow to use the coordinates of the best match as posStr
                noPerfectMatch = params.get("noPerfectMatch")

                if isinstance(tagNames, str):
                    tagNames = json.loads(tagNames)
                batchId = newMultiPamBatch(
                    newBatchName,
                    seq,
                    posStr,
                    org,
                    multipam,
                    insertSeq,
                    insertIdx,
                    assist=assist,
                    koGeneId=koGeneId,
                    insertPos=insertPos,
                    kiType=kiType,
                    tagNames=tagNames,
                    geneModel=geneModel,
                    nonCoding=nonCoding,
                    clippedSeq=clippedSeq,
                    noPerfectMatch=noPerfectMatch
                )
            else:
                assist = params["assist"]
                pamDesc = params["pam"]
                batchId = newMultiSeqBatch(
                    newBatchName,
                    multiseq,
                    org,
                    pamDesc,
                    koMethod,
                    geneModel,
                    koGeneId,
                    assist,
                    exonSelect,
                )

        else:
            if not pamSeq:
                pamDesc = params["pam"]
            seq = params["seq"]

            # the "seq" parameter can contain a chrom:start-end position instead of the sequence.
            if re.match(" *[a-zA-Z0-9_-]+: *[0-9, ]+ *- *[0-9,]+(:[+-])? *", seq):
                seq = getSeq(params["org"], seq.replace(" ", "").replace(",", ""))

            seq, warnMsg = cleanSeq(seq, org)

            if len(seq) > MAXSEQLEN2 and (isSlowPam(pamDesc)):
                errAbort(
                    "Sorry, but xCas9, SCanis and enCas12a have so many PAM sites that we are restricting "
                    " the input sequence length to %d bp at the moment to keep the "
                    "web site fast enough. We will revisit this in a few months. Let us know if "
                    "you think this is too short." % MAXSEQLEN2,
                    isWarn=True,
                )

            if len(seq) > MAXSEQLEN3 and (pamDesc in verySlowPams):
                errAbort(
                    "Sorry, but SpRY has so many PAM sites that we are restricting "
                    " the input sequence length to %d bp at the moment to keep the "
                    "web site usable. We will revisit this in a few months. Please let us know if "
                    "you think this is too short or have other ideas how to handle the issue."
                    % MAXSEQLEN3,
                    isWarn=True,
                )

            else:
                batchId = newBatch(newBatchName, seq, org, pamDesc)

        print("<script>")
        print(
            (
                """history.replaceState('crispor.py', document.title, '?batchId=%s');"""
                % (batchId)
            )
        )
        print("</script>")

    multiseq = params.get("multiseq")
    multipam = params.get("multipam")

    if multiseq and "seq" not in params and not multipam:
        mode = "multiseq"
        multiSearchDone = submitMultiSearch(batchId, org, pamDesc, mode)
        if multiSearchDone is None:
            printStatus(batchId, "")
            return
        KoResultsPage(params, batchId, koGeneId)
    elif "multipam" in params:
        mode = "multipam"
        multiSearchDone = submitMultiSearch(batchId, org, multipam, mode)
        if multiSearchDone is None:
            printStatus(batchId, "")
            return
        KiResultsPage(params, batchId)
    else:

        pam = setupPamInfo(pamDesc)

        assert "-" not in pam
        minFreq, varDb = checkOtherArgs(params)

        if len(warnMsg) != 0:
            print(warnMsg + "<p>")

        batchBase = join(batchDir, batchId)

        # read genome info tab file into memory
        dbInfo = readDbInfo(org)
        uppSeq = seq.upper()
        startDict, endSet = findAllPams(uppSeq, pam)
        otDone = getOfftargets(uppSeq, org, pamDesc, batchId, startDict, None)

        if otDone is None:
            # Job has been added to the queue or is not done yet.
            printStatus(batchId, warnMsg)
            return

        classicResultsPage(
            params,
            batchId,
            batchDir,
            warnMsg,
            pamDesc,
            pam,
            seq,
            org,
            dbInfo,
            uppSeq,
            startDict,
            endSet,
            minFreq,
            varDb,
        )


def classicResultsPage(
    params,
    batchId,
    batchDir,
    warnMsg,
    pamDesc,
    pam,
    seq,
    org,
    dbInfo,
    uppSeq,
    startDict,
    endSet,
    minFreq,
    varDb,
):
    """
    parse eff scores and offtargets and prints the results
    """

    # if we reach this, the batch has been processed
    batchInfo = readBatchAsDict(batchId)
    position = batchInfo.get("posStr")  # if there was no match, the posStr key is "?"

    if dbInfo == None:
        print(
            "<div class='title'>No Genome selected, specificity scoring is deactivated</div>"
        )
        print(
            '<div style="text-align:left;"><strong>Note:</strong> There are no predicted off-targets below and all specificity scores are shown in red as their score is 0. <br></div>'
        )
        chrom = ""

    elif position == "?":
        printQueryNotFoundNote(dbInfo)
        chrom = ""
    else:
        genomePosStr = ":".join(position.split(":")[:2])
        chrom, start, end, strand = parsePos(position)
        start = str(int(start) + 1)
        chrom = applyChromAlias(org, chrom)
        oneBasedPosition = "%s:%s-%s" % (chrom, start, end)

        print("<div class='title'><em>")
        if batchName != "":
            print(batchName + ":")

        ctUrl = None
        # if org in ["hg19", "mm10"]:
        # ctUrl = ctBaseUrl+"/%s.txt" % batchId

        print("%s (%s)</em>, " % (dbInfo.scientificName, dbInfo.name))
        print('<span style="text-decoration:underline">')
        # mouseOver = "link to UCSC,Ensembl or Gbrowse Genome Browser"
        mouseOver = None
        if dbInfo.server == "manual":
            mouseOver = "no genome browser link available for this organism"
        print(
            makeBrowserLink(
                dbInfo,
                genomePosStr,
                oneBasedPosition,
                mouseOver,
                ["tooltipster"],
                ctUrl=ctUrl,
            )
            + "</span>, "
        )
        if strand == "+":
            print(" forward genomic strand")
        else:
            print(" reverse genomic strand")
        print("</div>")
        # print " (link to Genome Browser)</div>"

    otMatches = parseOfftargets(org, batchId, chrom)
    effScores = readEffScores(batchId)
    sortBy = params.get("sortBy", "main")
    globEffScore = params.get("globEffScore", "EVA")
    guideData, guideScores, hasNotFound, pamIdToSeq = mergeGuideInfo(
        uppSeq,
        startDict,
        pam,
        otMatches,
        position,
        effScores,
        sortBy,
        org=org,
        globEffScore=globEffScore,
    )
    if len(guideScores) == 0:
        print(
            "Found no possible guide sequence. Make sure that your input sequence is long enough and contains at least one match to the PAM motif %s."
            % pam
        )
        print('<br><a class="neutral" href="crispor.py">')
        print(
            '<div class="button" style="margin-left:auto;margin-right:auto;width:150px;">New Query</div></a>'
        )
        return
    if hasNotFound and not position == "?":
        printNote(
            "At least one of the possible guide sequences was not found "
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
            "identifier (GCA_xxx or GCF_xxx), to crispor@tefor.net"
        )

    chrom, start, end, strand = parsePos(position)

    parNum = isInPar(org, chrom, start, end)
    if parNum is not None:
        print(
            (
                "<div style='text-align:left; background-color: aliceblue; padding:5px; border: 1px solid black'><strong>Note</strong>: The target sequence is in the PAR%s region. The off-targets on chrY's PAR copy have been removed from the off-target search. We treat the PAR regions as a single region, as all guides are assumed to modify both copies.</div>"
                % parNum
            )
        )

    varHtmls, varDbs, varDb = getVariants(
        seq, org, varDb, position, chrom, start, end, strand, minFreq
    )

    showSeqAndPams(
        org,
        seq,
        startDict,
        pam,
        guideScores,
        varHtmls,
        varDbs,
        varDb,
        minFreq,
        position,
        pamIdToSeq,
    )

    showSeqDownloadMenu(batchId)

    showGuideTable(guideData, pam, otMatches, dbInfo, batchId, org, chrom, varHtmls)

    print('<br><a class="neutral" href="crispor.py">')
    print(
        '<div class="button" style="margin-left:auto;margin-right:auto;width:150px;">New Query</div></a>'
    )

    makeCustomTrack(org, chrom, start, end, strand, guideData, batchId, batchName)

    print(
        """
    <script>
        // save the states of detail elements on page reload
        (function() {
            var $details = $('details[id]');
            $details.each(function() {
                var savedState = localStorage.getItem('details-' + this.id);
                if (savedState !== null) {
                    this.open = savedState === 'true';
                }
            });

            $details.on('toggle', function() {
                localStorage.setItem('details-' + this.id, this.open);
            });
        })();
    </script>
    """
    )


def KiResultsPage(params, batchId, download=False):
    """
    Parses and prints the results from Knock-in jobs.
    Optionnally returns the data formatted for downloadFile()

    Note : print calls should only occur when download is False
    """

    batchInfo = readBatchAsDict(batchId)
    org = batchInfo["org"]
    seq = batchInfo["seq"]
    posStr = batchInfo["posStr"]
    uppSeq = seq.upper()
    dbInfo = readDbInfo(org)
    multipam = batchInfo["multipam"]
    insertIdx = batchInfo["insertIdx"]
    insertPos = batchInfo["insertpos"]
    geneId = batchInfo.get("koGeneId")
    geneModel = batchInfo.get("geneModel")
    tagNames = batchInfo.get("tagNames")
    nonCoding = batchInfo.get("nonCoding")
    clippedSeq = batchInfo.get("clippedSeq")
    kiType = batchInfo.get("kiType")
    extSeq = batchInfo.get("extSeq")

    noPerfectMatch = batchInfo.get("noPerfectMatch")

    if kiType:
        insertSeq = batchInfo["insertSeq"]
    else:
        insertSeq = None

    sortBy = params.get("sortBy", "insertDistance")

    otMatches = parseOfftargets(org, batchId)
    effScores = readEffScores(batchId, multipam=multipam)
    globEffScore = params.get("globEffScore", "rs3")

    # increase the default window for substitutions to display PAMs for base editing
    maxPamWindow = max(len(seq[0:insertIdx]), len(seq[insertIdx: len(seq)]))
    if params.get("pairSortBy") is not None:
        defaultWindow = maxPamWindow
    elif kiType == "substitution":
        defaultWindow = 20
    else:
        defaultWindow = 10

    pamWindow = int(params.get("pamWindow", defaultWindow))
    otherPam = params.get("otherPam")

    pamList = set(multiPamDict[multipam][0])

    chrom, start, end, strand = parsePos(posStr)

    minFreq, varDb = checkOtherArgs(params)

    # these are read at function scope (in the mergeGuideInfo loop below), so they
    # must be initialized in both modes; the download path skips the block below
    useBaseEditor = False
    editData = None

    if not download:

        print("""
        <style>
            .pamPinned { background-color: #ffe680; outline: 2px solid #ff7f04; }
        </style>
        <script>
            $(document).ready(function() {
                $('a.guidePair').on('mouseenter', function() {
                    pinPairHighlight(this);
                }).on('click', function() {
                    pinPairHighlight(this);
                    var pamEls = pairPamEls(this);
                    var target = pamEls.length ? pamEls[0]
                                               : document.getElementById('seqStart');
                    if (target) {
                        target.scrollIntoView(
                            {behavior: 'smooth', block: 'center', inline: 'center'});
                    }
                });

                // switching to any guide table other than the pairs table clears
                // the pair PAM highlights, which only make sense in the pairs view.
                $('button[name="tableSelectButton"]').on('click', function() {
                    if (this.id !== 'pairSelect') {
                        clearPinnedHighlights();
                    }
                });
            });

            // collect the PAM elements of a guide pair in the sequence viewer.
            // pamIds contain +/-/. so we look them up via getElementById, not a #selector.
            function pairPamEls(pairEl) {
                return ['data-id1', 'data-id2'].map(function(attr) {
                    return document.getElementById('list' + pairEl.getAttribute(attr));
                }).filter(Boolean);
            }

            // remove every pinned PAM highlight currently on the page.
            function clearPinnedHighlights() {
                var pinned = document.querySelectorAll('.pamPinned');
                for (var i = 0; i < pinned.length; i++) {
                    pinned[i].classList.remove('pamPinned');
                }
            };

            // highlight both PAMs of a guide pair on click, clearing any previous pin.
            function pinPairHighlight(pairEl) {
                clearPinnedHighlights();
                pairPamEls(pairEl).forEach(function(el) {
                    el.classList.add('pamPinned');
                });
            };

        </script>
        """)

        print("""
        <script>
            function showTable(tableId, parentButton, setDist) {
                let targetTable = document.getElementById(tableId);
                let allTables = document.getElementsByName("guideTablePanel");
                let allButtons = document.getElementsByName("tableSelectButton");

                if (setDist === 1 && tableId != 'hdrTable') {
                    // for base editing, prime editing and double nicking, show all PAMs

                    let distButton = document.getElementById("distButton");
                    let distRange = document.getElementById("distRange");

                    if (distRange.value < %d) {
                        distRange.value = %d;
                        distButton.click();
                    }
                }

                for (button of allButtons) {
                    if (button.id === parentButton.id) {
                        // button.style.color = "#ff7f04";
                        button.className = "assistantButton active tooltipsterInteract"
                    } else {
                        // button.style.color = "#8A8278";
                        button.className = "assistantButton tooltipsterInteract";
                    }
                }

                for (table of allTables) {
                    if (table.id === tableId) {
                        table.style.display = "block";
                    } else {
                        table.style.display = "none";
                    }}
                }
        </script>
        """ % (maxPamWindow, maxPamWindow))

        genomePosStr = ":".join(posStr.split(":")[:2])
        chrom, start, end, strand = parsePos(posStr)
        start = str(int(start) + 1)
        chrom = applyChromAlias(org, chrom)
        oneBasedPosition = "%s:%s-%s:%s" % (chrom, start, end, strand)

        # mouseOver = "link to UCSC,Ensembl or Gbrowse Genome Browser"
        mouseOver = None
        if dbInfo.server == "manual":
            mouseOver = "no genome browser link available for this organism"
        if strand == "+":
            strandStr = " forward genomic strand"
        else:
            strandStr = " reverse genomic strand"

        browserLink = makeBrowserLink(
            dbInfo,
            genomePosStr,
            oneBasedPosition,
            mouseOver,
            ["tooltipster"],
            ctUrl=None,
        )

        # should build the link here instead to avoid calling the function twice
        browserUrl = makeBrowserLink(
            dbInfo,
            genomePosStr,
            oneBasedPosition,
            mouseOver,
            ["tooltipster"],
            ctUrl=None,
            returnUrl=True,
        )

        browserLinkHtml = """<div>
        <em>%s</em> (%s),
        %s  %s
        </div>""" % (
            dbInfo.scientificName,
            dbInfo.name,
            browserLink,
            strandStr,
        )

        varHtmls, varDbs, varDb = getVariants(
            seq, org, varDb, posStr, chrom, int(start), int(end), strand, minFreq
        )

        # Base editing can be used for these substitutions
        useBaseEditor = False
        editData = None
        beFilter = None
        if kiType == "substitution":

            seqMsg = "%s -> %s substitution" % (seq[insertIdx].upper(), insertSeq.upper())

            fromToNucl = (seq[insertIdx].upper(), insertSeq)
            for editList in possibleEdits.values():
                if fromToNucl in editList:
                    # load edit data
                    batchBase = join(batchDir, batchId)
                    editFname = batchBase + ".editData.json"
                    if isfile(editFname):
                        editData = json.load(open(editFname))
                        bePamIds = buildEditData(editData).keys()
                    else:
                        # avoid error message if no editData was written (subserver crash)
                        editData = {}
                    if len(editData) > 0:
                        # baseEditor in reinitialized in setupPamInfo
                        useBaseEditor = True
                        orgPamList = multiPamDict[multipam][0]

                        # data to filter out BE guides that can't introduce the substitution,
                        # while keeping all HDR guides for the selected PAM list
                        beFilter = (orgPamList, bePamIds)

                        for pamVariant in pamVariantModels.keys():
                            pamList.add(pamVariant)

        elif kiType == "deletion":
            seqMsg = "%dbp deletion" % (len(batchInfo["insertSeq"]))
        elif kiType == "replacement":
            seqMsg = "replacement of a %s bp sequence" % len(batchInfo["insertSeq"])
        else:
            seqMsg = "knock-in of a %s bp sequence" % len(batchInfo["insertSeq"])

        if geneId:
            if "ENST" in geneId:
                transcriptUrl = (
                    """in %s of <a href="https://www.ensembl.org/Multi/Search/Results?q=%s;site=ensembl;page=1" target="blank">%s</a> """
                    % (insertPos, geneId.split("_")[0], geneId.split("_")[0])
                )
            else:
                transcriptUrl = (
                    """in %s of <a href="https://www.ncbi.nlm.nih.gov/nuccore/%s/" target="blank">%s</a> """
                    % (insertPos, geneId, geneId)
                )
        else:
            chromPos = int(start) + int(insertIdx)
            transcriptUrl = "at position <a target='blank' href='%s'>%d</a> in %s" % (
                browserUrl,
                chromPos,
                chrom,
            )

        print(
            """<div class="title" style="text-align:left; margin-bottom:12px;margin-top:12px;"><i>%s (%s)</i> : <span style="text-decoration:underline">%s %s</span></div><br> """
            % (dbInfo.scientificName, org, seqMsg, transcriptUrl)
        )

        printKiSteps(batchId, step=1, align="left", useBaseEditor=useBaseEditor, insertSeq=insertSeq)

        if nonCoding is not None:
            htmlWarn("Could not get START or STOP codons ")
            print(
                """ No START or STOP codons could be found. This gene is either non-coding, or something is wrong with the annotation of its coding sequence. The insertion will likely not be in-frame, so we suggest manually copying the target sequence by clicking on "Enter a sequence" from the Knock-in menu page (this mode supports the editing of coding or non-coding regions)."""
            )

        if geneModel:
            exonSeqsPlaceholder = []
            print(
                "<p>Show below is the gene model with the edits. Hover on these to get their full name and length</p>"
            )
            printGeneModel(
                geneModel,
                exonSeqsPlaceholder,
                koMethod=None,
                insertSeq=insertSeq,
                insertPos=insertPos,
                kiType=kiType,
                tagNames=tagNames,
            )

    allGuideData = []
    allGuideScores = {}
    allPamIdToSeq = {}
    allPamIdToSuppInfo = {}

    for pamFullName in pamList:
        pam = setupPamInfo(pamFullName)
        startDict, endSet = findAllPams(uppSeq, pam)

        guideData, guideScores, hasNotFound, pamIdToSeq, pamIdToSuppInfo = (
            mergeGuideInfo(
                uppSeq,
                startDict,
                pam,
                otMatches,
                None,
                effScores,
                sortBy,
                org=org,
                exonId=None,
                globEffScore=globEffScore,
                pamFullName=pamFullName,
                insertIdx=insertIdx,
                kiType=kiType,
                insertSeq=insertSeq,
                getSuppInfo=True,
                allEditData=editData,
                beFilter=beFilter
            )
        )

        allGuideData.extend(guideData)
        allGuideScores.update(guideScores.copy())
        allPamIdToSeq.update(pamIdToSeq.copy())
        allPamIdToSuppInfo.update(pamIdToSuppInfo.copy())

    sortGuideData(allGuideData, sortBy=sortBy)

    if download:
        return seq, org, pam, posStr, allGuideData
    else:

        pairedGuides = None
        if len(insertSeq) < 10:
            # select paired guides in PAM-out orientation for double nicking strategy width D10A
            # only documented for small edits
            pairedGuides = getNickPairs(seq, pamList, insertIdx, allGuideData)

        showSeqAndPams(
            org,
            seq,
            None,
            None,
            allGuideScores,
            varHtmls,
            varDbs,
            varDb,
            minFreq,
            posStr,
            allPamIdToSeq,
            browserLinkHtml=browserLinkHtml,
            multiPamInfo=(multipam, pamList, insertIdx, kiType, insertSeq),
            pamIdToSuppInfo=allPamIdToSuppInfo,
            pamWindow=pamWindow,
            otherPam=otherPam,
            strand=strand,
            clippedSeq=clippedSeq,
            geneId=geneId,
            useBaseEditor=useBaseEditor,
            extSeq=extSeq,
            editData=editData,
            batchId=batchId,
            noPerfectMatch=noPerfectMatch,
        )

        showSeqDownloadMenu(batchId)

        # Sequence viewer filtering / display menu
        print("""<details id="results1" open>""")
        print(
            """<summary style="font-weight: bold; font-size: 20px; margin-top: 24px; margin-bottom: 12px;">Sequence viewer display and filtering</summary>"""
        )
        print(
            """
        <div style="width: 100%; display: flex; flex-direction: row; gap: 12px; margin-top: 24px; overflow-x: scroll;">
        """
        )

        print(
            (
                """<form style="display:inline" id="paramForm" action="%s" method="GET">"""
                % basename(__file__)
            )
        )

        print("""<input type="hidden" name="batchId" value="%s"/>""" % batchId)

        print(
            """
        <div class="windowstep subpanel" style="align-items: center; padding: 12px; margin-right: 32px; min-width: 450px; box-shadow: 0 0 15px 4px rgba(47, 129, 203, 0.14);">

        <div style="display: flex; flex-direction: column;">
            <div style="display: flex; flex-direction: row; align-items: top; gap: 5%%; margin-bottom: 12px;">
                <div style="font-weight: bold; font-size: 18px;">Options to modify the display of PAMs on the sequence viewer</div>
            </div>
            <div style="display: flex">
                <div>
                <div style="display: flex; flex-direction: row; align-items: center; gap: 6px;">
                Show PAMs <input id="distRange" type=range value=%s min=10 max=%d name="pamWindow" oninput="this.nextElementSibling.value = this.value"/>
                <output>%s</output> bp from the edit site
                </div>
        """
            % (
                pamWindow,
                maxPamWindow,
                pamWindow,
            )
        )

        print(
            """
        <div style="display: flex; flex-direction: row; align-items: center; gap: 12px;">
        <p>Add PAMs</p>
              """
        )

        print(""" <select style="width:80%; height: 25px;" name="otherPam"> """)
        print("""<option value="">only %s</option>""" % multiPamDict[multipam][1])
        for pamKey in multiPamDict:
            if pamKey == multipam:
                continue
            if pamKey == otherPam:
                selectStr = 'selected="selected"'
            else:
                selectStr = ""
            pamText = multiPamDict[pamKey][1]
            print(
                """ <option %s value="%s">+ %s</option> """
                % (selectStr, pamKey, pamText)
            )
        print(
            """
              </select>
              </div>"""
        )

        print(
            """
                </div>
            </div>
            <div style="margin-top: 24px; display: flex; flex-direction: row; gap: 10%; align-self: center;">
                <button id="distButton" style="width: 80px; height: 32px; display: flex; align-items: center; justify-contents: center;" name="submit" type="submit" value="Update">Update</button>
        """
        )

        if otherPam is not None:
            # get the url to submit a new search with the selected PAM list
            newParams = {
                param: val
                for param, val in batchInfo.items()
                if param not in ["batchId", "multipam"]
            }
            # gene model is a list of lists : convert it to json format
            if "tagNames" in newParams:
                newParams["tagNames"] = json.dumps(newParams["tagNames"])
            if "geneModel" in newParams:
                newParams["geneModel"] = json.dumps(newParams["geneModel"])
            newParams.update({"newSearch": 1, "submit": "submit"})
            newParams.update({"multipam": otherPam})
            actionStr = "?" + urllib.parse.urlencode(newParams)

            print(
                """
                <a href="%s" style="width: 320px; display: flex; align-items: center;">
                    <div class="button" style="display: flex; align-items: center; height: 32px;">New search with the selected PAMs</div>
                </a>
            """
                % actionStr
            )
        print(
            """
            </div>

        </div>
        </div>
        </form>
        """
        )
        print(
            """
        <p style="width: 50%;">
        Note : by default, only guides that result in a DSB less than 10bp away from the edit site are displayed.<br>
        If there are no guides found or you want to examine more possible guides, either increase this theshold, or look for other PAMs in the sequence (in gray on the sequence viewer). Note that guides with other PAMs don't have any specificity or efficiency scores yet. You can submit a new search with the selected PAMs to get the corresponding guides and scores.
        </p>
        </div>
        """
        )

        print("</details>")

        if posStr == "?":
            print(
                """
            <p>Since the input sequence was not found in the genome, we can't show to page to design the donor DNA. To design a donor DNA for this experiment, either do it manually or submit a new search with a genomic sequence corresponding to the genome you selected. If the sequence you entered is already edited, please use the wild-type version. If you want to add a new genome, contact us at %s </p>
            """
                % contactEmail
            )

        # showSeqAndPams above has already validated and written annotation
        # params back into cgiParams, so resolveAnnotationParams returns a
        # clean dict consistent with the currently selected gene model
        annotParams = resolveAnnotationParams(org, seq, posStr)

        print("""
        <div class="assistantMenu" style="margin-bottom: 24px; margin-left: 18px; margin-top: 12px;">
            <button class="assistantButton active tooltipsterInteract" style="font-size: 24px;" name="tableSelectButton" id="hdrSelect" onclick="showTable('hdrTable', this, setDist=0)">guides for HDR-based editing </button>
        """)

        if pairedGuides and len(pairedGuides) > 0:
            print("""
            <button class="assistantButton tooltipsterInteract" style="font-size: 24px;" name="tableSelectButton" id="pairSelect" onclick="showTable('pairTable', this, setDist=1)">pairs of guides for HDR-based editing using a double-nicking strategy</button>
            """)

        if useBaseEditor:
            print("""
            <button class="assistantButton tooltipsterInteract" style="font-size: 24px;" name="tableSelectButton" id="beSelect" onclick="showTable('beTable', this, setDist=1)">guides for base editing</button>
            """)

        print("</div>")

        print("""<div name="guideTablePanel" id="hdrTable" >""")
        showGuideTable(
            allGuideData,
            pam,
            otMatches,
            dbInfo,
            batchId,
            org,
            chrom,
            None,
            geneId=None,
            pamFullName="multipam",
            pamWindow=pamWindow,
            annotParams=annotParams,
        )
        print("</div>")

        if pairedGuides:

            showPairedGuidesTable(pairedGuides, annotParams, params, batchId)

        if editData:
            print("""<div name="guideTablePanel" id="beTable" >""")
            tableEditData = buildEditData(editData, targetPos=insertIdx)
            showGuideTable(
                allGuideData,
                pam,
                otMatches,
                dbInfo,
                batchId,
                org,
                chrom,
                None,
                geneId=None,
                pamFullName="multipam",
                pamWindow=pamWindow,
                annotParams=annotParams,
                editData=tableEditData
            )
            print("</div>")

        print('<br><a class="neutral" href="crispor.py?expType=ki">')
        print(
            '<div class="button" style="margin-left:auto;margin-right:auto;width:150px;">New Query</div></a>'
        )

    if not download:
        print(
            """
        <script>
        // save the states of detail elements on page reload
        (function() {
            var $details = $('details[id]');
            $details.each(function() {
                var savedState = localStorage.getItem('details-' + this.id);
                if (savedState !== null) {
                    this.open = savedState === 'true';
                }
            });

            $details.on('toggle', function() {
                localStorage.setItem('details-' + this.id, this.open);
            });
        })();

        // restore which guide-table panel was selected on page reload.
        $(document).ready(function() {
            var $buttons = $('button[name="tableSelectButton"]');

            const savedId = localStorage.getItem('kiActiveTable');
            const savedButton = savedId ? document.getElementById(savedId) : null;
            if (savedButton && savedButton.getAttribute('name') === 'tableSelectButton') {
                let buttonId = savedButton.getAttribute('id');
                if (buttonId === "pairSelect") {
                   showTable('pairTable', savedButton, setDist=0);
                } else if (buttonId === "beSelect") {
                    showTable('beTable', savedButton, setDist=0)
                } else {
                    showTable('hdrTable', savedButton, setDist=0)
                }
                // will add PE table here
            } else {
                // default view: show the HDR table, hide the rest
                const hdrButton = document.getElementById('hdrSelect');
                if (hdrButton && typeof showTable === 'function') {
                    showTable('hdrTable', hdrButton, setDist=0);
                }
            }

            // remember the active table whenever a button is clicked
            $buttons.on('click', function() {
                localStorage.setItem('kiActiveTable', this.id);
            });
        });

        </script>
              """
        )


def printKiSteps(batchId: str, step=1, backParams=None, align="center", useBaseEditor=False, insertSeq=None):
    """
    prints an interactive recap of the workflow for KI experiments
    """

    if step not in range(1, 4):
        return

    # make the current step bold
    else:
        stepStyles = [
            "opacity: 1; font-weight: 1000;" if i == step else "opacity: 0.5;"
            for i in range(1, 4)
        ]

    backUrl = basename(__file__) + "?" + "batchId=" + batchId

    # first step

    if useBaseEditor:
        guideSelectText = """
        Several CRISPR-based methods can be used to edit this sequence.
        """
        guideHtml = (
            """ <a href="%s" class="tooltipsterInteract" title="%s" style="%s; font-size: 1.25em;">Select an editing stategy</a> """
            % (backUrl, guideSelectText, stepStyles[0])
        )

        print("""
        <script>
            function showBeTable(tableId) {
                let allTables = document.getElementsByName("guideTablePanel");
                let allButtons = document.getElementsByName("tableSelectButton");

                for (button of allButtons) {
                    if (tableId == 'beTable' && button.id === 'beSelect') {
                        button.className = "assistantButton active tooltipsterInteract";
                    } else if (tableId == 'hdrTable' && button.id === 'hdrSelect') {
                        button.className = "assistantButton active tooltipsterInteract";
                    } else {
                        // button.style.color = "#8A8278";
                        button.className = "assistantButton tooltipsterInteract";
                    }
                }

                for (table of allTables) {
                    if (table.id === tableId) {
                        table.style.display = "block";
                        table.scrollIntoView();
                    } else {
                        table.style.display = "none";
                    }}
                }
        </script>
        """)

    else:
        guideSelectText = """
        First, select guides that introduce a DSB as close as the possible to the editing site.<br>
        By default, the table is sorted by this distance and only guides that result in a cut at less than 10bp from the editing site are shown.<br>
        If no guides satisfies this condition, go to the 'Options to modify the display of PAMs on the sequence viewer' box.<br>
        You can either :
        <ul>
            <li>Display guides that cut further away from the editing site (using the cursor).</li>
            <li>See if another enzyme may be more appropriate by showing other PAM patterns.</li>
        </ul>
        Using guides distant to the editing site may result in a low efficiency.<br>
        In this case, using nickases with two guides that flank the editing site (double nicking strategy) may yield better results (see Schubert et al. 2021 - fig. 2).<br>
        """
        guideHtml = (
            """ <a href="%s" class="tooltipsterInteract" title="%s" style="%s; font-size: 1.25em;">Select guide sequences</a> """
            % (backUrl, guideSelectText, stepStyles[0])
        )

    if useBaseEditor and insertSeq:
        editorText = """
        One or more guide sequences can be used with a base editor to introduce this substitution. Hover or click on the 'Edits to %s' field on the sequence viewer to get a summary of the available guides.<br>
        For each type of editor, a prediciton of its efficiency and outcome sequences is shown, allowing you to select the most appropriate deaminase domain.<br>
        For a more detailed explanation of base editor selection rules, read <a href='https://doi.org/10.1038/s41587-023-01792-x' target='blank'>Kim et al. 2023</a>. Depending on the aim of your experiment, you may consider bystander editing (i.e base conversions that occur outside of the intended substitution) in addition to the predicted editing efficiency.
        """ % insertSeq
        editorHtml = """<a class="tooltipsterInteract" title="%s" style="%s; font-size: 1.25em; margin-left: 12px;" onclick=showBeTable('beTable')> Choose a base editor</a>""" % (editorText, stepStyles[1])

    donorDesignText = """
    Once you have selected a suitable guide sequence, click on 'design donor DNA' under the 'guide sequence + PAM' column of the table.<br>
    Note that the target region (guide + PAM) may be identical between the genome and the donor DNA.<br>
    <ul>
        <li>This can result in the cutting of the donor, or re-cutting of the inserted sequence after a successful edit.</li>
        <li>To prevent this, you can choose to recode the donor DNA sequence by introducing blocking mutations in the PAM or the seed region of the spacer.</li>
        <li>In this case, the design will likely be specific to its corresponding guide.</li>
        <li>Otherwise, all guides that don't have the 'donor needs recoding' flag (highlighted in blue on the sequence below) can be used with the same, non recoded donor DNA design.</li>
    <ul>
    """

    # second step
    if step == 3:
        # make the "design donor DNA" step clickable only as a return link in the last step, to avoid confusion
        donorBackUrl = printBackLink(toDonorPage=True, returnUrl=True)

        # if present, include manual annotation params in the url
        if backParams:
            donorBackUrl = donorBackUrl + "&" + urllib.parse.urlencode(backParams)

        donorDesignHtml = (
            """<a href="%s" class="tooltipsterInteract" title="%s" style="%s color: #ff6000; font-size: 1.25em;">&nbsp Design donor DNA</a> """
            % (donorBackUrl, donorDesignText, stepStyles[0])
        )

    else:
        donorDesignHtml = (
            """<div class="tooltipsterInteract" title="%s" style="%s color: #ff6000; font-size: 1.25em;">&nbsp Design donor DNA</div> """
            % (donorDesignText, stepStyles[1])
        )

    # third step
    donorDisplayText = """
    After designing the donor DNA, its sequence will be displayed with some annotations.<br>
    Then, you may trim the sequence to facilitate its synthesis (for example by removing repeated regions), or generate alternative designs with a barcode to check for homozygous editing (<i>not implemented yet</i>).<br>
    """
    donorDisplayHtml = (
        """ <div class="tooltipsterInteract" title="%s" style="%s color: #ff6000; font-size: 1.25em;">&nbsp Visualize and download guide + donor DNA</div> """
        % (donorDisplayText, stepStyles[2])
    )

    # print all steps
    if useBaseEditor:
        print(
            """<div style="text-align: %(align)s;">
                <p>Workflow overview <small>(hover on each step to get details, or click on a previous step to go back)</small></p>
                <div style="display: flex; flex-direction: row; justify-content: %(align)s; gap: 12px; margin-top: 12px;">
                    <div style="margin-top: 14px;">
                        %(guideHtml)s
                    </div>
                    <div style="display: flex; flex-direction: column; gap: 12px">
                        <div style="font-weight: 1000; font-size: 1.25em;">&nbsp &#8599</div>
                        <div style="font-weight: 1000; font-size: 1.25em;">&nbsp &#8600</div>
                    </div>
                    <div style="display: flex; flex-direction: column; gap: 32px; margin-top: 8px;">
                        <div style="margin-top: -12px; display: flex; flex-direction: row;">
                            %(donorDesignHtml)s
                            <div style="font-weigth: bold; font-size: 1.25em;">&nbsp &#8594</div>
                            %(donorDisplayHtml)s
                        </div>
                        <div style="margin-top: -12px; display: flex; flex-direction: row;">
                            %(editorHtml)s
                        </div>
                    </div>
                    </div>
                </div>
            </div>
           """
            % locals()
        )

    else:
        print(
            """<div style="text-align: %s; margin-bottom: 48px;">
                <p>Workflow overview <small>(hover on each step to get details, or click on a previous step to go back)</small></p>
                <div style="display: flex; flex-direction: row; justify-content: %s;">
                    %s
                    <div style="font-weight: 1000; font-size: 1.25em;">&nbsp &#8594</div>
                    %s
                    <div style="font-weigth: bold; font-size: 1.25em;">&nbsp &#8594</div>
                    %s
                </div>
            </div>
           """
            % (align, align, guideHtml, donorDesignHtml, donorDisplayHtml)
        )


def deserialize(inStr, inType="list"):
    """
    deserializes the repr of a tuple or a list into a list
    JSON should be used, but params saved as hidden html inputs can't be loaded as json
    """

    if inType == "tuple":
        lbrd, rbrd = "(", ")"
    else:
        lbrd, rbrd = "[", "]"

    return inStr.strip(lbrd).strip(rbrd).replace("'", "").replace(" ", "").split(",")


def showDonor(
    HA5,
    HA3,
    newInsertSeq,
    recodedArmSeq,
    mutEvents,
    noModel,
    recodeArm,
    HA5repeats,
    HA3repeats,
    params,
):
    """
    Dispays the donor DNA sequence
    needs a lot of reactoring, contains a lot more features than anticipated
    """

    donorType = params["donorType"]
    batchId = params["batchId"]
    batchInfo = readBatchAsDict(batchId)
    donorName = params.get("donorName")
    pamId = params["pamId"]
    guideSeq = params["guideSeq"]
    guideInfo = params["guideInfo"]

    doubleNicking = params.get("doubleNicking")

    if doubleNicking:

        revPamId = params["revPamId"]
        revGuideSeq, revPamSeq, revGuideStrand, revGuideStart = deserialize(params["revGuideInfo"], inType="tuple")
        revDoRecoding = params["revDoRecoding"]
        revCfd = params["revCfd"]

        fwPamId = params["fwPamId"]
        fwGuideSeq, fwPamSeq, fwGuideStrand, fwGuideStart = deserialize(params["fwGuideInfo"], inType="tuple")
        fwDoRecoding = params["fwDoRecoding"]
        fwCfd = params["fwCfd"]

        revGuideSeq, fwGuideSeq = revGuideSeq.upper(), fwGuideSeq.upper()

        doubleNickingParams = {
            "doubleNicking": doubleNicking,
            "revPamId": revPamId,
            "revDoRecoding": revDoRecoding,
            "revCfd": revCfd,
            "fwPamId": fwPamId,
            "fwDoRecoding": fwDoRecoding,
            "fwCfd": fwCfd
            }

    geneId = batchInfo.get("koGeneId")
    kiType = batchInfo.get("kiType")
    org = batchInfo["org"]
    seq = batchInfo["seq"]
    # newInsertSeq = sequence to use in the donor DNA ("" for deletions)
    # insertSeq = deleted sequence in deletion mode
    # outside of deletions, insertSeq = newInsertSeq
    insertSeq = batchInfo["insertSeq"]
    posStr = batchInfo["posStr"]
    insertIdx = int(batchInfo["insertIdx"])
    tagNames = batchInfo.get("tagNames")
    insertPos = batchInfo["insertpos"]
    seq = batchInfo["seq"]

    # save params to include it in the return link in printKiSteps()
    backParams = {}
    for param in ["manualExStart", "manualExEnd", "manualExFrame"]:
        val = params.get(param)
        if val:
            backParams[param] = val

    if doubleNicking:
        backParams.update(doubleNickingParams)

    if ("tagseq" in params and "linkerseq") in params or (
        "markerseq" in params and "expressionSeq" in params and "qTag" in params
    ):
        tagNames, newInsertSeq = getInsertSeq(
            params.get("linkerseq"),
            params.get("tagseq"),
            params.get("markerseq"),
            params.get("expressionSeq"),
            params.get("qTag"),
            batchInfo["insertpos"],
        )

    # change the insert sequence on this page
    if "replaceInsertSeq" in params:
        replaceInsertSeq = params["replaceInsertSeq"]
        replaceInsertSeq = re.sub(r"[^ATGCNatgcn]", "", replaceInsertSeq)
        if not (len(replaceInsertSeq) != len(insertSeq) and kiType == "replacement"):
            newInsertSeq = replaceInsertSeq

    dbInfo = readDbInfo(org)
    chrom, seqStart, seqEnd, strand = parsePos(posStr)

    # this is not good
    # but can't save guideInfo as json as it comes from a hidden html input
    pamSeq, guideStart, pamStrand = deserialize(guideInfo, inType="tuple")
    pamStart = int(pamId.split(".")[1].strip("[s+-]"))

    if donorType == "ss":
        polarity = params["polarity"]
        if polarity == "positive" and strand == "-":
            templateStrand = "+"
        elif polarity == "negative" and strand == "+":
            templateStrand = "-"
        else:
            templateStrand = strand
    else:
        templateStrand = strand

    if recodedArmSeq:
        orgDonor = HA5 + newInsertSeq + HA3
        if recodeArm == "HA3":
            donorSeq = HA5 + newInsertSeq + recodedArmSeq
        else:
            donorSeq = recodedArmSeq + newInsertSeq + HA3
    else:
        donorSeq = HA5 + newInsertSeq + HA3
        orgDonor = donorSeq

    if strand != templateStrand:
        isRev = True
        donorSeq = revComp(donorSeq)
        orgDonor = revComp(orgDonor)
    else:
        isRev = False

    insertCoord = "%s:%s" % (chrom, seqStart + insertIdx)
    if kiType == "substitution":
        kiTypeStr = "%s-%ssubst" % (seq[insertIdx], newInsertSeq)
    elif kiType == "deletion":
        kiTypeStr = "del%sbp" % (len(insertSeq))
    elif kiType == "replacement":
        kiTypeStr = "replace%sbp" % (len(newInsertSeq))
    else:
        kiTypeStr = "insert%sbp" % (len(newInsertSeq))

    # move this to a new function (annotateDonor)

    donorBins, homopolymers = processDonor(donorSeq)
    highlights = []

    if donorType == "ss" and templateStrand != strand:
        editStart = len(HA3)
    else:
        editStart = len(HA5)

    # highlight the substitution
    if kiType == "substitution":
        editEnd = editStart + 1
    # highlight insert
    elif kiType != "deletion":
        editEnd = editStart + len(newInsertSeq)
    else:
        editEnd = editStart + len(insertSeq)

    # coordinates of the PAM + guide
    # + avoid out of bounds coordinates
    pamInDel = False
    if pamStart < insertIdx:
        pamStartCoord = len(HA5) - (insertIdx - pamStart)
        # invert the coordinates for ssODN with reverse polarity
        if donorType == "ss" and templateStrand != strand:
            pamStartCoord = len(donorSeq) - pamStartCoord - len(pamSeq)
        pamEndCoord = pamStartCoord + len(pamSeq)

    # PAM inside the deletion : don't highlight it
    # for deletions, adjust corrdinates with insertSeq
    # for other edits, adjust with insertSeq (replaced sequence), "" for deletions
    elif (
        kiType == "deletion"
        and pamStart >= insertIdx
        and pamStart + len(pamSeq) <= insertIdx + len(insertSeq)
    ):
        pamInDel = True
        pamStartCoord, pamEndCoord = editStart, editStart
        # distance between the PAM and deletion start
        pamDistFromEdit = pamStart - insertIdx
    else:
        if kiType == "deletion":
            pamStartCoord = len(HA5) + (pamStart - insertIdx - len(insertSeq))
        elif kiType in ["substitution", "replacement"]:
            pamStartCoord = len(HA5) + (pamStart - insertIdx)
        else:
            pamStartCoord = len(HA5) + len(newInsertSeq) + (pamStart - insertIdx)
        if donorType == "ss" and templateStrand != strand:
            pamStartCoord = len(donorSeq) - pamStartCoord - len(pamSeq)
        pamEndCoord = pamStartCoord + len(pamSeq)

    if (
        pamStartCoord < len(HA5)
        and pamEndCoord > len(HA5)
        and kiType not in ["substitution", "deletion", "replacement"]
    ):
        pamCoordList = [
            (pamStartCoord, len(HA5)),
            (len(HA5) + len(newInsertSeq), pamEndCoord + len(newInsertSeq)),
        ]
    else:
        pamCoordList = [(pamStartCoord, pamEndCoord)]

    # the logic is getting out of hand
    # 3 cases : --ggPAM-- + strand, - strand with Cas12a or - strand with inverse polarity ssODN
    guideFirst = (
        (pamStrand == "+" and not pamIsFirst and templateStrand == strand)
        or (pamIsFirst and pamStrand == "-" and templateStrand == strand)
        or (pamStrand == "-" and (pamIsFirst or templateStrand != strand))
    )

    # PAM in 3'
    if guideFirst:
        guideStartCoord = pamStartCoord - len(guideSeq)
        guideEndCoord = pamStartCoord

        # insert sequence inside the guide
        if kiType == "insertion" and guideStartCoord < editEnd < guideEndCoord:
            guideStartCoord = guideStartCoord - len(newInsertSeq)

        # OK
        if kiType == "deletion":

            # the guide start inside the deletion : clip its length to the deletion end
            # --xxxxxxxx------
            # ----gggggggPAM---
            if guideStartCoord > editStart and guideEndCoord < editEnd:
                guideStartCoord = editEnd

            #        <--> pamDistFromEdit
            # -------xxxxxxxx--
            # ----gggggggPAM---
            elif guideStartCoord < editStart and pamInDel:
                guideStartCoord = editStart - (len(guideSeq) - pamDistFromEdit)
                guideEndCoord = editStart

            # the guide sequence overlaps the deletion : reduce its length
            # -----xxxxx------
            # ----gggggggPAM---
            elif guideStartCoord < editStart and editEnd < pamStartCoord:
                guideStartCoord += len(insertSeq)

    # PAM in 5'
    else:
        guideStartCoord = pamEndCoord
        guideEndCoord = guideStartCoord + len(guideSeq)

        if kiType == "insertion" and guideStartCoord < editStart < guideEndCoord:
            guideEndCoord = guideEndCoord + len(newInsertSeq)

        # OK
        if kiType == "deletion":
            # ------xxxxxxxxx--
            # ---PAMgggggg-----
            if guideEndCoord > editStart and guideEndCoord < editEnd:
                guideEndCoord = editStart

            # -xxxxxxxxx-------
            # -----PAMggggg----
            elif guideEndCoord > editEnd and pamInDel:
                guideStartCoord = editStart

                # with ssODN in reverse polarity, pamDistFromEdit is calculated as if the PAM is in 3'
                if templateStrand != strand:
                    guideEndCoord = editStart + (len(guideSeq) - pamDistFromEdit)
                else:
                    guideEndCoord = (
                        editStart
                        + len(guideSeq)
                        - (len(insertSeq) - (pamDistFromEdit + len(pamSeq)))
                    )

            # ------xxxx--------
            # ---PAMgggggg-----
            elif guideEndCoord > editEnd and pamEndCoord < editStart:
                guideEndCoord -= len(insertSeq)

    guideCoordList = [(guideStartCoord, guideEndCoord)]

    # calculate the CFD score of the guide + PAM region in the donor DNA
    guideOutsideDonor = False
    if (pamStrand == "+" and templateStrand == strand) or (pamStrand == "-" and templateStrand != strand):

        if guideStartCoord < 0 or pamEndCoord > len(donorSeq):
            guideOutsideDonor = True
        else:
            guideTarget = guideSeq + pamSeq
            guideDonor = donorSeq[guideStartCoord:pamEndCoord]
    else:
        if pamStartCoord < 0 or guideEndCoord > len(donorSeq):
            guideOutsideDonor = True
        else:
            if templateStrand != strand:
                guideTarget = guideSeq + pamSeq
            else:
                guideTarget = guideSeq + pamSeq

            guideDonor = donorSeq[pamStartCoord:guideEndCoord]
            guideDonor = revComp(guideDonor)

    if not guideOutsideDonor:
        donorCfd = calcCfdScore(guideTarget.upper(), guideDonor.upper())
        # donorCfd = 1
    else:
        donorCfd = -1

    if homopolymers and donorType != "ss":
        for start, end in homopolymers:
            highlights.append(
                (
                    start,
                    end,
                    {"background-color": "rgba(0, 0, 0, 0)"},
                    "homopolymerCoord",
                )
            )
    if donorBins and donorType != "ss":
        for start, end, GCfrac in donorBins:
            if GCfrac > 0.8:
                highlights.append(
                    (
                        start,
                        end,
                        {"background-color": "rgba(0, 0, 0, 0)"},
                        "GcRichCoord",
                    )
                )

    # coordinates of repeated regions
    for repeatStart, repeatEnd in HA5repeats:
        highlights.append(
            (
                repeatStart,
                repeatEnd,
                {"background-color": "rgba(0, 0, 0, 0)"},
                "repeatCoord",
            )
        )

    if len(HA3repeats) > 0:
        if kiType in ["substitution", "deletion", "replacement"]:
            offset = len(HA3)
        else:
            offset = len(HA3) + len(newInsertSeq)

    for repeatStart, repeatEnd in HA3repeats:
        highlights.append(
            (
                repeatStart + offset,
                repeatEnd + offset,
                {"background-color": "rgba(0, 0, 0, 0)"},
                "repeatCoord",
            )
        )

    if pamCoordList and guideCoordList:
        for start, end in guideCoordList:
            highlights.append(
                (start, end, {"background-color": "rgba(0, 0, 255, 0.5)"}, "guideCoord")
            )
        for start, end in pamCoordList:
            highlights.append(
                (start, end, {"background-color": "rgba(0, 255, 255, 0.5)"}, "pamCoord")
            )

    # priority to insert sequence
    if kiType == "deletion":
        # highlight deletion start / end in this case
        editStart, editEnd = (editStart - 1, editStart + 1)
    highlights.append(
        (
            editStart,
            editEnd,
            {"background-color": "rgba(255, 255, 0, 0.5)"},
            "insertCoord",
        )
    )

    if not donorName:
        donorName = "Don_%s_%s_%s" % (org, insertCoord, kiTypeStr)
    donorTypeText = "ssODN" if donorType == "ss" else "double stranded donor DNA"

    # if geneId and insertPos == "Cter":
    #     stop = ''.join([base for base in seq if base.isupper()])

    if kiType == "substitution":
        insertText = "substitution"
    elif kiType == "replacement":
        insertText = "replacement sequence"
    else:
        insertText = "insert sequence"

    print(
        """
    <script>
    // from https://www.w3schools.com/howto/howto_js_copy_clipboard.asp
    // blocked if the website is accessed over http
          function copyDonor() {

          var donorSeq = document.getElementById("donorSeq");

          navigator.clipboard.writeText(donorSeq.value).then(function() {

            // alert("Copied the donor DNA sequence");

          }).catch(function(err) {

          console.error("failed to copy: ", err);
          document.execCommand("copy");

          })
        }
    </script>
          """
    )

    if kiType == "substitution":
        seqMsg = "%s -> %s substitution" % (seq[insertIdx], newInsertSeq)
    elif kiType == "deletion":
        seqMsg = "%dbp deletion" % len(insertSeq)
    elif kiType == "replacement":
        seqMsg = "%dbp replacement" % len(insertSeq)
    else:
        seqMsg = "knock-in of a %s bp sequence" % len(newInsertSeq)

    if geneId:
        if "ENST" in geneId:
            transcriptUrl = (
                """in %s of <a href="https://www.ensembl.org/Multi/Search/Results?q=%s;site=ensembl;page=1" target="blank">%s</a> """
                % (insertPos, geneId.split("_")[0], geneId.split("_")[0])
            )
        else:
            transcriptUrl = (
                """in %s of <a href="https://www.ncbi.nlm.nih.gov/nuccore/%s/" target="blank">%s</a> """
                % (insertPos, geneId, geneId)
            )
    else:
        transcriptUrl = "at position %s in %s" % (seqStart + insertIdx, chrom)

    print(
        """<div class="title" style="text-align:center; margin-bottom=50px;margin-top=50px;">%s : %s %s </div><br> """
        % (dbInfo.scientificName, seqMsg, transcriptUrl)
    )

    # print("<small>Black lines = homology arms. Sequence = insert sequence</small>")

    printKiSteps(batchId, step=3, backParams=backParams)

    print("""<input type="hidden" id="donorSeq" value="%s"></input>""" % donorSeq)

    print(
        """<div style='width: 80%; margin-top: 54px; margin-left:20%; margin-right:50%; text-align:left;'>"""
    )

    print("""<h2>Guide sequence : %s</h2> """ % guideSeq.upper())
    print(
        """
          <div style="margin-top: 32px; margin-bottom: 12px; display: flex; align-items: flex-end;">
              <h4 style="align-self: flex-end; margin: 0;">Sequence of the %s : </h4>&nbsp&nbsp
              <button style="align-self: flex-end; margin: 0;" onclick="copyDonor()"><small>Copy sequence to clipboard</small></button><br> """
        % donorTypeText
    )

    print("<form>")

    baseParams = {
            "batchId": batchId,
            "donorSeq": donorSeq,
            "pamId": pamId,
            "guideSeq": guideSeq,
            "donorName": donorName
            }

    if doubleNicking:
        baseParams.update({
            "doubleNicking": doubleNicking,
            "revPamId": revPamId,
            "revGuideSeq": revGuideSeq,
            "fwPamId": fwPamId,
            "fwGuideSeq": fwGuideSeq
            })
    if recodedArmSeq:
        baseParams.update({"recodedDonorSeq": donorSeq})

    printHiddenFields(params, baseParams)

    if doubleNicking:
        pairStr = " pair"
    else:
        pairStr = ""
    print(
        """
              <button name="downloadDonor" value="download" style="align-self: flex-end; margin: 0; margin-left: 12px;"><small>Download fasta (guide%s + donor DNA)</small></button>
          </form>
          </div> """ % pairStr
    )

    # barcode step, to add later
    '''
    if kiType != "substitution":
        print("""
        <div style="margin-bottom: 24px;">
            Output a second DNA donor with a barcode to check for homozygous editing <i>(not implemented yet)</i>:&nbsp
            <input type="checkbox" style="-webkit-box-shadow:  0px 0px 0px 2px #ff6000;
            -moz-box-shadow: 0px 0px 0px 2px #ff6000;
            box-shadow: 0px 0px 0px 2px #ff6000;" id="doBarcode" autocomplete="off"/>
            <img src=" %s image/info-small.png" title="With this option, a second donor DNA with a few mutations will be displayed, so that homozygous editing can be detected by PCR or NGS." class="tooltipsterInteract"><br>
        </div>
        """ % HTMLPREFIX)
    '''

    if kiType == "qTag":
        tagNamesStr = ", ".join([name for name in tagNames if "lox" not in name])
        print(
            """<p>You selected the qTag system (<a target="blank" href="https://doi.org/10.1038/s44318-024-00337-5">Philip et al. 2025</a>) with the following elements : <br>
            <ul>
            <li>%s.</li>
            </ul>
          Note that the donor DNA sequence below is designed for synthesis and only contains the sequence of these elements.<br>
          The plasmid construct is available on <a target="blank" href="https://www.addgene.org/browse/article/28238680/">Addgene</a> and may be more convenient.</p>"""
            % tagNamesStr
        )

    if kiType == "substitution":
        editSpanText = "Substitued base"
    elif kiType == "deletion":
        editSpanText = "Deletion start / end"
    elif kiType == "replacement":
        editSpanText = "Replacement"
    else:
        editSpanText = "Insert sequence"

    print(
        """
    <script>
    function showHighlight(coord, color) {
        // let spans = document.querySelectorAll(`[name="${coord}"], #${coord}, .${coord}`);
        let spans = document.getElementsByName(coord);
        let check = event.target;

        if (check.checked) {
            for (let span of spans) {
                span.style.backgroundColor = color;
            }
        } else {
            for (let span of spans) {
                span.style.backgroundColor = 'rgba(0, 0, 0, 0)';
            }
        }
    }
    </script>
    """
    )

    print(
        """
    <div style="display:flex; flex-direction: row; gap: 24px;">
        <p><span style="background-color: rgba(0, 255, 255, 0.5)">PAM</span></p>
        <p><span style="background-color: rgba(0, 0, 255, 0.5)">Spacer</span></p>
        <p><span style="background-color: rgba(255, 255, 0, 0.5)"><u>%s</u></span></p><br>
    </div>
    Click on the checkboxes below to show the following sequence features:
    <div style="display:flex; flex-direction: row; gap: 2px;">
        <input type="checkbox" autocomplete="off" onchange="showHighlight('repeatCoord', 'rgba(102, 255, 51, 0.5)')"/><p><span class="tooltipsterInteract" title="UCSC genomes are maksed using repeatMakser (for interspersed repeats and low complexity DNA sequences) and TRF (for tandem repeats). Other genomes might have been annotated with different programs." style="background-color: rgba(102, 255, 51, 0.5)">Repeats (TRF + repeatMasker)</span></p>
    """
        % editSpanText
    )

    print(
        """
    <input type="checkbox" autocomplete="off" onchange="showHighlight('homopolymerCoord', 'rgba(255, 0, 0, 0.5)')"/><p><span style="background-color: rgba(255, 0, 0, 0.5)">Homopolymers (10+ A/T or 6+ G/C)</span></p>
    <input type="checkbox" autocomplete="off" onchange="showHighlight('GcRichCoord', 'rgba(153, 51, 255, 0.5)')"/><p><span style="background-color: rgba(153, 51, 255, 0.5)">GC rich (> 80% over 20+ nt)</span></p>
    """
    )
    print(
        """
    </div>
    """
    )

    print(
        """<div style="font-family: Source Code Pro; justify-self:center; margin-bottom:24px;">"""
    )
    fastaWidth = 80
    if len(donorSeq) > fastaWidth:
        maxStrCount = len(str(len(donorSeq)))
        last = 0
        print("".join([" " for i in range(maxStrCount)]) + ">%s" % donorName)
        print("<br>")
        for i in range(fastaWidth, len(donorSeq), fastaWidth):
            print(getHighlightedRow(donorSeq, last, i, highlights))
            print("<br>")
            last = i

        if last < len(donorSeq):
            print(getHighlightedRow(donorSeq, last, len(donorSeq), highlights))
            print("<br>")
    else:
        print(">%s" % donorName)
        print("<br>")
        print(getHighlightedRow(donorSeq, 0, len(donorSeq), highlights))

    print("""</div>""")

    print("""<div style="margin-right: 15%;">""")

    if donorCfd != -1:
        print("<h4>CFD score of the guide against donor DNA : %s</h4>" % round(donorCfd, 2))
    else:
        print("<p>The guide sequence is outside the coordinates of the donor DNA and will unlikely re-cleave it after insertion.</p>")

    if noModel is True:
        print(
            """
            <p>The donor DNA sequence could not be recoded because no gene annotation file could be located for this organism. If you want to add one, contact us %s
            """
            % contactEmail
        )
    if mutEvents:
        print("<h3>Silent mutations introduced to prevent re-cut</h3>")
        print("<h4>Notes on recoding</h4>")
        print(
            """<p>
        Codons in the regions to recode are replaced by the synonymous codon with the highest frequency (codon frequency is calculated based on the longest transcript for each protein-coding gene). The recoded bases are show in uppercase.
        </p>
        """
        )
        recNc = False
        for wt, freq, pos in mutEvents:
            mut, mutFreq, motif = mutEvents[(wt, freq, pos)]
            # successful recoding in non-coding region
            if len(wt) == 1 and len(mut) == 1:
                recNc = True
        if recNc:
            print(
                """<p>
            Some non-coding regions were recoded. In this case, a transition is introduced every 3nt.
            Motifs 5bp upstream and downstream of exon boundaries are avoided
            to minimize the disruption of splicing. The kozak consensus sequence (6bp upstream of the
            START codon) is also avoided.<br>
            </p>"""
            )

        print(
            """<p>
            However, we cannot guarantee that recoding will not affect gene expression, as some
            transcription factor binding sites or RNA methylation motifs may be modified.<br>
            If you want to manually recode the donor, you can find its unmodified sequence in the fasta file (click on the button at the top of the page to download it).
            </p>"""
        )

        printMutEventsTable(mutEvents, HA3, insertSeq, HA5, recodeArm, isRev)

    # form to change the insert sequence
    if kiType != "deletion":

        print("<form>")

        # to verify
        if recodedArmSeq:
            printHiddenFields(
                params,
                {
                    "batchId": batchId,
                    "recodedDonorSeq": donorSeq,
                    "donorSeq": orgDonor,
                    "pamId": pamId,
                    "guideSeq": guideSeq,
                    "donorName": donorName,
                    "donorType": donorType,
                    "replaceInsertSeq": None,
                    "newInsertSeq": None,
                    "tagseq": None,
                    "markerseq": None,
                    "expressionSeq": None,
                    "qTag": None,
                },
            )
        else:
            printHiddenFields(
                params,
                {
                    "batchId": batchId,
                    "donorSeq": donorSeq,
                    "pamId": pamId,
                    "guideSeq": guideSeq,
                    "donorName": donorName,
                    "donorType": donorType,
                    "replaceInsertSeq": None,
                    "newInsertSeq": None,
                    "tagseq": None,
                    "markerseq": None,
                    "expressionSeq": None,
                    "qTag": None,
                },
            )

        print(
            """<div style="display: flex; flex-direction: column; align-items: left; gap: 8px; margin: 0 auto; margin-right: 30%; min-width: 800px;"> """
        )
        if kiType == "tagging":
            print(
                """<p>If you want to select other markers, use the dropdown menu below and click on update button. Updating with nothing selected will restore the original sequence.</p>"""
            )
            printTagsAndLinkers(qTAG=False, tagNames=tagNames)

        elif kiType == "qTag":
            print(
                """<p>If you want to select other elements, use the dropdown menu below and click on the update button. Updating with nothing selected will restore the original sequence.</p>"""
            )
            printTagsAndLinkers(tag=False, tagNames=tagNames)

        elif kiType == "substitution":
            print(
                """<p>If you want to change the substitution, select another base from the dropdown menu below and click on the update button.</p>"""
            )
            print("""<select name="replaceInsertSeq" style="width: 48px;">""")
            for base in [
                base
                for base in ["A", "T", "G", "C"]
                if base not in [insertSeq.upper(), seq[insertIdx]]
            ]:
                print("""<option value="%s">%s</option>""" % (base, base))
            print("</select>")
        elif kiType == "replacement":
            replaceLen = len(insertSeq)
            print(
                """
                <p>If you want to change the replaced sequence, paste a new one here and click on the update button. The sequence should be of the same length as the current one (%(replaceLen)d bp). Updating with an empty box will restore the original sequence.</p>
                <input name="replaceInsertSeq" style="width: 45%%;" maxlength=%(replaceLen)d minlength=%(replaceLen)d placeholder="paste a new sequence replacement here (%(replaceLen)d bp)"/>
            """
                % locals()
            )
        else:
            print(
                """
                <p>If you want to change the insert sequence, paste a new one here and click on the update button. Updating with an empty box will restore the original sequence.</p>
                <textarea name="replaceInsertSeq" maxlength=5000 cols=100 rows=6 placeholder="paste a new insert sequence here (max 5kb)"></textarea>
            """
            )

        # print("""<input type="hidden" name=%s value=%s />""" % ("insertSeq", insertSeq))
        # print("""<input type="hidden" name=%s value=%s />""" % ("donorType", donorType))

        # trimming options
        # to add : minimun length of homology arms
        '''
        if donorType == "ds":
            print("""
                <div id="trimOptions">
                    <strong>Trim the donor to facilitate its synthesis</strong> <i>(not implemented yet)</i><br>
                    <input type="checkbox" name="trimHomopolymers" form="main" value="True" autocomplete="off"/>Trim homopolymers<br>
                    <input type="checkbox" name="trimGC" form="main" value="True" autocomplete="off"/>Trim regions with high GC content<br>
                    <input type="checkbox" name="trimRepeat" form="main" value="True" autocomplete="off"/>Trim repeated sequences<br>
                </div>
                """)
        '''

        print(
            """<button style="align-self: center; width: 125px;" type="submit" value="update">update</button>"""
        )
        print("</form>")

    if donorType == "ss":
        print("<details>")
        print(
            "<summary>Show the free energy and secondary structure of the ssODN</summary>"
        )
        showSecondaryStructure(params, donorSeq=donorSeq)
        print("</details>")

    print("</div>")
    print("</div>")


def getHighlightedRow(seq, rowStart, rowEnd, highlights):
    """
    Highlights regions of a sequence string with HTML spans.

    Args:
        seq: The full sequence string.
        rowtart: The start index of the sequence chunk to process.
        rowEnd: The end index of the sequence chunk to process.
        highlights: A list of tuples, where each tuple is
                    (hlStart, hlEnd, styleDict).
                    style_dict is a dictionary of CSS properties,
                    e.g. {'background-color': 'yellow'}.
    """
    if not highlights:
        return seq[rowStart:rowEnd]
    styles = [{} for _ in range(len(seq))]
    names = ["" for _ in range(len(seq))]
    for hlStart, hlEnd, style_dict, name in highlights:
        for i in range(hlStart, hlEnd):
            if 0 <= i < len(styles):
                styles[i].update(style_dict)
                if name:
                    names[i] = name

    htmlParts = []
    currentStyleStr = ""
    currentName = ""

    # process within the bounds of the sequence
    processingEnd = min(rowEnd, len(seq))

    for i in range(rowStart, processingEnd):
        styleDict = styles[i]
        # Create a canonical string representation for comparison.
        styleStr = "; ".join(sorted(["%s: %s" % (k, v) for k, v in styleDict.items()]))
        name = names[i]

        if styleStr != currentStyleStr or name != currentName:
            if currentStyleStr != "":
                if currentStyleStr == "background-color: rgba(255, 255, 0, 0.5)":
                    htmlParts.append("</u></span>")
                else:
                    htmlParts.append("</span>")
            if styleStr != "":
                if styleStr == "background-color: rgba(255, 255, 0, 0.5)":
                    htmlParts.append("""<span style="%s"><u>""" % styleStr)
                elif name:
                    htmlParts.append(
                        """<span name="%s" style="%s">""" % (name, styleStr)
                    )
                else:
                    htmlParts.append("""<span style="%s">""" % styleStr)

            currentStyleStr = styleStr
            currentName = name

        htmlParts.append(seq[i])

    if currentStyleStr != "":
        htmlParts.append("</u></span>")

    return "".join(htmlParts)


def buildEditData(jsonData, stopGuides=None, targetPos=None):
    """
    converts edit data from a dict based on edit position to a dict based on pamIds
    returns {pamId: [scores]}, to be used in showGuideTable()

    showGuideTable() only renders the first edit of each guide, so the edit at the
    guide's target position is moved to the front: otherwise a lower-position
    bystander edit (possibly with a restricted set of models) would be shown
    instead of the edit that actually introduces the desired change.
    In KO/stop mode the target is the per-guide stop position (stopGuides[pamId][0]);
    in KI/substitution mode it is the single insertion position (targetPos)."""

    editData = {}

    for exonId, exonIdDict in jsonData.items():
        for pos, posDict in exonIdDict.items():
            for base, pamIdList in posDict.items():
                for editTpl in pamIdList:
                    pamId, _, _, _, effs, outcomes, _ = editTpl
                    if len(outcomes) == 0:
                        continue
                    if pamId in editData:
                        editData[pamId].append((pos, base, effs, outcomes))
                    else:
                        editData[pamId] = [(pos, base, effs, outcomes)]

    # put the target-position edit first for each guide
    if stopGuides or targetPos is not None:
        for pamId, edits in editData.items():
            if stopGuides:
                stopInfo = stopGuides.get(pamId)
                if stopInfo is None:
                    continue
                guideTargetPos = stopInfo[0]
            else:
                guideTargetPos = targetPos
            guideTargetPos = int(guideTargetPos)
            edits.sort(key=lambda edit: int(edit[0]) != guideTargetPos)

    return editData


def getNickPairs(seq, pamList, insertIdx, guideData):
    """
    select pairs of guides in PAM-out orientation to be used for HDR-based editing with D10A
    Nick sites should be 40-68bp apart
    The edit position may be positioned anywhere inside the nick sites
    See https://doi.org/10.1038/s41598-021-98965-y
    """

    pairedGuides = []
    minDist, maxDist = 40, 118
    # pamId, MIT, CFD, globalScore, effScores, offTargets, guideSeq, pamSeq, cutUpstream, doRecoding
    scoreById = {row[6]: (row[6], row[0], row[1], row[19], row[2], row[9], row[7], row[8], row[17], row[18]) for row in guideData}

    for pam in pamList:
        if pam != "NGG":
            continue
        pam = setupPamInfo(pam)
        startDict, endSet = findAllPams(seq, pam)

    leftGuides = []
    rightGuides = []

    for pamStart, pamStrand in startDict.items():
        pamId = "NGG.s%d%s" % (pamStart, pamStrand)
        if pamStrand == "+" and pamStart - 3 > insertIdx:
            nickPos = pamStart - 3
            rightGuides.append((pamId, nickPos))
        elif pamStrand == "-" and pamStart + 6 < insertIdx:
            nickPos = pamStart + 6
            leftGuides.append((pamId, nickPos))

    # print(leftGuides, "<br>", rightGuides)

    # select pairs with adequate distance between nick sites
    for leftPamId, leftNickPos in leftGuides:

        leftPamInfo = scoreById[leftPamId]

        for rightPamId, rightNickPos in rightGuides:
            nickDist = rightNickPos - leftNickPos

            if minDist <= nickDist <= maxDist:

                rightPamInfo = scoreById[rightPamId]
                meanScore = 0.5 * (leftPamInfo[3] + rightPamInfo[3])
                pairedGuides.append((leftPamInfo, rightPamInfo, nickDist, meanScore))

    # sort by the mean of global scores
    pairedGuides.sort(key=lambda x: x[3], reverse=True)

    return pairedGuides


def printMutEventsTable(mutEvents, HA3, insertSeq, HA5, recodeArm, isRev):
    """displays the silent mutations introduced by recoding"""

    codonTable = buildCodonTable()

    # sort the mutation by position
    mutEvents = dict(sorted(mutEvents.items(), key=lambda x: x[0][2]))
    coding = False
    for wt, freq, pos in mutEvents:
        if len(wt) == 3:
            coding = True
            break

    # display rules for recoding in non-coding regions

    print("<table>")
    if coding:
        print(
            """
        <tr>
            <th>Amino acid</th>
            <th>Wild type</th>
            <th>Original codon frequency</th>
            <th>Recoded</th>
            <th>Recoded codon frequency</th>
            <th>Position</th>
            <th>Region</th>
        </tr>
        """
        )
    else:
        print(
            """
        <tr>
            <th>Region</th>
            <th>Wild type</th>
            <th>Mutated</th>
            <th>Position</th>
            <th>Region</th>
        </tr>
        """
        )

    for wt, freq, pos in mutEvents:
        mut, muFreq, motif = mutEvents[(wt, freq, pos)]
        if len(wt) == 3:
            aaStr = "<strong>%s</strong>" % codonTable[wt]
        else:
            if mut == "No (in 5'UTR)" and coding is False:
                aaStr = "5'UTR"
            elif mut == "No (in a splicing site)" and coding is False:
                aaStr = "splicing site"
            else:
                aaStr = "non-coding"

        if motif == "pam":
            regionTxt = "PAM motif"
        elif motif == "seed":
            regionTxt = "PAM-proximal end of the guide"
        elif motif == "gap":
            regionTxt = "between the cut site and edit site"

        if recodeArm == "HA3" and pos is not None:
            if isRev:
                fullPos = len(HA3) - pos - 3
            else:
                fullPos = len(HA5) + len(insertSeq) + pos
        elif pos is not None:
            if isRev:
                fullPos = (len(HA5) + len(insertSeq) + len(HA3)) - pos
            else:
                fullPos = pos
        else:
            fullPos = pos

        mut, mutFreq, motif = mutEvents[(wt, freq, pos)]
        print("<tr>")
        print(""" <td>%s</td> """ % aaStr)
        print(""" <td>%s</td> """ % wt)
        if coding:
            print(""" <td>%s</td> """ % freq)
        print(""" <td>%s</td> """ % mut)
        if coding:
            print(""" <td>%s</td> """ % mutFreq)
        print(""" <td>%s</td> """ % fullPos)
        print(""" <td>%s</td> """ % regionTxt)
        print("</tr>")
    print("</table>")

    print("""</div>""")


def KoResultsPage(params, batchId, koGeneId, download=False):
    """read offtargets and effcores files generated by the multiseq job,
    aggregates the data in a loop and display it. Optionally returns guide data to prepare for downloadFile()
    """

    batchInfo = readBatchAsDict(batchId)
    org = batchInfo["org"]
    dbInfo = readDbInfo(org)
    pam = batchInfo["pam"]
    pam = setupPamInfo(pam)

    koMethod = batchInfo["koMethod"]
    if koMethod == "stop":
        # reset the global here in case the PAM isn't correctly assigned
        global baseEditor
        baseEditor = True
        # beWinStart, beWinEnd = getBeWin(cgiParams.get("beWin", DEFAULTBEWIN))
        batchBase = join(batchDir, batchId)
        editFname = batchBase + ".editData.json"
        if isfile(editFname):
            allEditData = json.load(open(editFname))
        else:
            allEditData = {}
    else:
        allEditData = None

    stopGuides = batchInfo.get("stopGuides")
    geneModel = batchInfo["geneModel"]
    koGeneId = batchInfo["koGeneId"]
    if "SYM" in koGeneId:
        koGeneId = koGeneId.split("~")[0]
        commonExons = True
    else:
        commonExons = False
    exonSeqs = batchInfo["exonSeqs"]
    exonPosStr = batchInfo["exonPosStr"]
    exonSelect = batchInfo.get("exonSelect")

    sortBy = params.get("sortBy", "main")
    globEffScore = params.get("globEffScore", "rs3")
    selGeneModel = params.get("geneModelSelection")

    otMatches = parseOfftargets(org, batchId)
    effScores = readEffScores(batchId)

    minFreq, varDb = checkOtherArgs(params)

    allGuideData = []
    allGuideScores = {}
    allPamIdToSeq = {}
    if not download:

        print("""<div class="title" style="text-align:left; margin-bottom: 50px;">""")
        print("<i>%s (%s)</i></em> : " % (dbInfo.scientificName, dbInfo.name))
        if koMethod == "frameshift":
            titleText = "introducing a frameshift mutation"
        elif koMethod == "stop":
            titleText = "introducing premature STOP codons or disrupting a splice site with base editing"
        elif koMethod == "excision":
            titleText = "excision of the gene locus"
        elif koMethod == "promoter":
            titleText = "excision of the promoter"
        elif koMethod == "splicing":
            titleText = "disrupting splicing"
        else:
            titleText = ""

        if "ENST" in koGeneId:
            transcriptUrl = (
                "https://www.ensembl.org/Multi/Search/Results?q=%s;site=ensembl;page=1"
                % koGeneId.split("_")[0]
            )
        else:
            transcriptUrl = "https://www.ncbi.nlm.nih.gov/nuccore/%s/" % koGeneId

        if commonExons and koMethod not in ["excision", "promoter"]:
            commonExonStr = " (common exons)"
        else:
            commonExonStr = ""
        print('<span style="text-decoration:underline">')
        print(
            """Knock-out of <a href="%s" target="_blank">%s</a>%s </span> by %s"""
            % (transcriptUrl, koGeneId, commonExonStr, titleText)
        )
        print("""</div>""")

        print("""<div style="margin-bottom: 24px;">""")

        if koMethod in ["frameshift", "stop"] or (
            koMethod == "splicing" and len(exonPosStr) > 2
        ):
            printGeneModel(geneModel, exonSeqs, koMethod, commonExons=commonExons)
        print(
            """<p>Below are the target and PAM sequences. Lowercase bases corresponds to an extension of the target exon sequence (to identify guides that bind exon boundaries and induce a DSB in the exon).<p>"""
        )
        if koMethod == "frameshift":
            print(
                """
            In-frame methionine codons are highlighted in green, to avoid selecting guides that could result in a DSB upstream of an alternative START codon.<br>
                  """
            )

        print(
            """Colors <span style="color:#32cd32; text-shadow: 1px 1px 1px #bbb">green</span>, <span style="color:#ffff00; text-shadow: 1px 1px 1px #888">yellow</span> and <span style="text-shadow: 1px 1px 1px #f01; color:#aa0014">red</span> indicate high, medium and low specificity of the PAM's guide sequence in the genome.<br>"""
        )
        print(
            "Click on a match for the PAM %s below to show its %d bp-long guide sequence.<br>"
            % (pam, GUIDELEN)
        )

        print("</div>")

        geneModels, selGeneModel, selTransId = getSelGeneModel(org, noGenes=False)
        if geneModels:
            if koMethod == "stop":
                extendPos = GUIDELEN + 6
            else:
                extendPos = GUIDELEN - 6
            for model, modelStr in geneModels:
                exonInfo, maxTransIdLen = getExonInfo(
                    org, model, exonPosStr[0], extendPos=extendPos
                )
                for transId, sym in list(exonInfo.keys()):
                    # for common exons, koGeneId is a gene symbol : select the first transcript by default
                    if transId in koGeneId or koGeneId in transId or sym == koGeneId:
                        selGeneModel = model
                        if commonExons:
                            selTransId = "allTrans"
                        else:
                            selTransId = transId
                        break

        if koMethod in ["excision", "promoter"]:
            # for experiments with a pair of guides, show two results pages

            print(
                "<p>This knock-out method is based on a large deletion, resulting from two DSBs introduced by a pair of guides. Click on the buttons below to show guides for the regions upstream or downstream of the deletion.</p>"
            )

            print(
                """
            <script>
                function showResults(region, parentButton) {
                    const displayUpstream = document.getElementById('displayUpstream');
                    const displayDownstream = document.getElementById('displayDownstream');
                    const allButtons = document.getElementsByName('pairButton');

                    for (button of allButtons) {
                        if (button.id === parentButton.id) {
                            button.className = "assistantButton active tooltipsterInteract";
                        } else {
                            button.className = "assistantButton tooltipsterInteract";
                        };
                    }

                    // TODO : save button states on page reload
                    const upstreamClicked = document.querySelector("#showUpStream");
                    const downstreamClicked = document.querySelector("#showDownStream");

                    if (region === 'up') {
                        displayUpstream.style.display = 'block';
                        displayDownstream.style.display = 'none';
                    } else {
                        displayUpstream.style.display = 'none';
                        displayDownstream.style.display = 'block';
                    }
                };
            </script>
           """
            )

            print(
                """
                <div class="assistantMenu" style="margin-bottom: 24px; margin-left: 18px; margin-top: 12px;">
                        <button class="assistantButton tooltipsterInteract active" id="showUpstream" name="pairButton" value="up" style="font-size: 24px;" onclick="showResults(this.value, this)">Show results for the upstream region</button>
                        <button  class="assistantButton tooltipsterInteract" id="showDownstream" name="pairButton" value="down" style="font-size: 24px;" onclick="showResults(this.value, this)">Show results for the downstream region</button>
                </div>
                  """
            )

    for exonSeqInfo, posStr in zip(exonSeqs, exonPosStr):

        exonId, seq = exonSeqInfo
        uppSeq = seq.upper()
        startDict, endSet = findAllPams(uppSeq, pam, exonId)
        chrom, start, end, strand = parsePos(posStr)

        if not download:
            start = str(int(start) + 1)
            chrom = applyChromAlias(org, chrom)
            oneBasedPosition = "%s:%s-%s:%s" % (chrom, start, end, strand)

            browserLink = makeBrowserLink(
                dbInfo, posStr, oneBasedPosition, None, ["tooltipster"]
            )

            varHtmls, varDbs, varDb = getVariants(
                seq, org, varDb, posStr, chrom, int(start), int(end), strand, minFreq
            )

            if (
                exonId == 0
                and koMethod not in ["excision", "promoter"]
                and (baseEditor or varDb)
            ):
                print("<form>")
                print("""<input type="hidden" value=%s name="batchId" />""" % batchId)
                if baseEditor:
                    print("""<details id="results4" open autocomplete="off">""")
                    print(
                        """<summary style="font-weight: bold; font-size: 20px; margin-top: 24px; margin-bottom: 12px;">Base editing information</summary>"""
                    )
                    print("<i>Note : this feature is still in early development</i>")
                    print(
                        "<p>Show below the sequence are the possible edits, using this base editor with the selected modification window.<br>"
                    )
                    if koMethod == "stop":
                        print(
                            """
                              <ul>
                                    <li>Edits in red result in the introduction of a premature STOP codon.</li>
                                    <li>Edits in grey corresponds to "bystander" edits (i.e, additional edits that can occur when using the same guide).</li>
                              </ul>
                              """
                        )
                    print(
                        """Hover on an edit to show the corresponding guides and their predicted efficiencies and outcome sequences.<br>
                             Clicking on the edit will redirect to the row of the guide with the highest predicted efficiency.<br>"""
                    )

                    '''
                    print("Base Editor modification window:")
                    selBeWin = "%s-%s" % (beWinStart, beWinEnd)
                    print(
                        (
                            """<input type="text" name="beWin" size="10" value="%s">"""
                            % selBeWin
                        )
                    )
                    print(
                        """<input style="height:18px;margin:0px;font-size:10px;line-height:normal" type="submit" name="submit" value="Update">"""
                    )
                    '''

                    print("</p>")
                    print("</details>")

                if varDb is not None:
                    print("Variant database:")
                    varDbList = [
                        (b, c) for a, b, c, d in varDbs
                    ]  # only keep fname+label
                    printDropDown("varDb", varDbList, varDb)

                    if minFreq == 0.0:
                        minFreq = "0.0"
                    else:
                        minFreq = str(minFreq)

                    # pull out the hasAF field for this varDb
                    varDbHasAF = False
                    for shortLabel, fname, desc, hasAF in varDbs:
                        if fname == varDb:
                            varDbHasAF = hasAF
                            break

                    if varDbHasAF:
                        print("""&nbsp; Min. frequency: """)
                        print(
                            (
                                """<input type="text" name="minFreq" size="8" value="%s">"""
                                % minFreq
                            )
                        )
                    print(
                        """<input style="height:18px;margin:0px;font-size:10px;line-height:normal" type="submit" name="submit" value="Update">"""
                    )
                    print(
                        (
                            "<small style='margin-left:30px'><a href='mailto:%s'>Missing a variant database? We can add it.</a></small>"
                            % contactEmail
                        )
                    )
                print("</form>")

        if (
            koMethod == "splicing"
            and exonId % 2 == 0
            and not exonSelect.isnumeric()
            and not download
        ):
            originalExon = (exonId + 1) // 2
            print(
                """<h3 class="exonGroup%s" name="exonDisplay">Guides for exon %d</h3>"""
                % (originalExon, (originalExon + 1))
            )

        if koMethod in ["excision", "promoter"] and not download:
            if exonId == 0:
                print("""<div id="displayUpstream" style="display: block;"> """)
                print("<h2>Guides in the upstream region</h2>")
            elif exonId == 1:
                print("""<div id="displayDownstream" style="display: none;"> """)
                print("<h2>Guides in the downstream region</h2>")

        guideData, guideScores, hasNotFound, pamIdToSeq = mergeGuideInfo(
            uppSeq,
            startDict,
            pam,
            otMatches,
            posStr,
            effScores,
            sortBy,
            org=org,
            exonId=exonId,
            globEffScore=globEffScore,
            stopGuides=stopGuides,
            allEditData=allEditData
        )
        if koMethod in ["excision", "promoter"]:
            sortGuideData(guideData, sortBy)

        # return a single guideData object when downloading data from pairs of guides (otherwise, guideData are kept separated)
        if (
            len(guideData) > 0
            and koMethod not in ["excision", "promoter"]
            or (download and koMethod in ["excision", "promoter"])
        ):
            allGuideData.extend(guideData)
            allGuideScores.update(guideScores.copy())
            allPamIdToSeq.update(pamIdToSeq.copy())

        if not download:

            showExonAndPams(
                batchId,
                org,
                seq,
                startDict,
                pam,
                guideScores,
                varHtmls,
                varDbs,
                varDb,
                minFreq,
                posStr,
                pamIdToSeq,
                exonId,
                koMethod,
                browserLink,
                selGeneModel=selGeneModel,
                selTransId=selTransId,
                exonSelect=exonSelect,
                stopGuides=stopGuides,
            )

            # for methods that require a pair of guides, two tables are shown
            if koMethod in ["excision", "promoter"]:
                showGuideTable(
                    guideData,
                    pam,
                    otMatches,
                    dbInfo,
                    batchId,
                    org,
                    chrom,
                    None,
                    koGeneId,
                    koMethod=koMethod,
                )
                print("</div>")

    if not download:
        showSeqDownloadMenu(batchId)
        if baseEditor:
            printJson("editData", allEditData)
            editData = buildEditData(allEditData, stopGuides=stopGuides)
        else:
            editData = None

    # for experiements using a pair of guides, sort the table by each target sequence
    # if koMethod in ["excision", "promoter"]:
    #    sortGuideData(allGuideData, sortBy, exonSort=True)
    # else:

    # handle sorting the guide data for pairs of guide in download mode

    if download is False and koMethod not in ["excision", "promoter"]:
        sortGuideData(allGuideData, sortBy)
        showGuideTable(
            allGuideData,
            pam,
            otMatches,
            dbInfo,
            batchId,
            org,
            chrom,
            None,
            koGeneId,
            koMethod=koMethod,
            exonSelect=exonSelect,
            editData=editData
        )

    if download is False:
        print('<br><a class="neutral" href="crispor.py?expType=ko">')
        print(
            '<div class="button" style="margin-left:auto;margin-right:auto;width:150px;">New Query</div></a>'
        )

    else:
        return exonSeqs, org, pam, exonPosStr, allGuideData


def getVariants(seq, org, varDb, position, chrom, start, end, strand, minFreq):
    "returns the variant information to be displayed in showSeqAndPams()"

    # get list of variant databases
    varLabel = None
    varDbs = readVarDbs(org)

    if len(varDbs) > 0 and not position == "?":
        if varDb is None:
            varDb = varDbs[0][1]

        # pull out label of the variant database
        varLabel = None
        for shortLabel, varKey, lab, hasAF in varDbs:
            if varKey == varDb:
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

    return varHtmls, varDbs, varDb


def printGeneModel(
    geneModel,
    exonSeqs,
    koMethod=None,
    insertSeq=None,
    insertPos=None,
    kiType=None,
    tagNames=None,
    commonExons=False,
):
    """
    Displays the gene model, from CDS start to CDS end
    Optionally make target exons as buttons"""

    if koMethod in ["frameshift", "stop"] and exonSeqs:
        thirdLen = 0
        for feature in geneModel:
            if feature[0] == "exon":
                length = feature[2]
                thirdLen += length
        if koMethod == "stop":
            # thirdLen = 2 * math.ceil(thirdLen / 3)
            pass
        else:
            thirdLen = math.ceil(thirdLen / 3)
            thirdLen = max(thirdLen, 250)

        lastseq = str(exonSeqs[-1][1])

        # remove flanking sequences
        lastseq = "".join(base for base in lastseq if base.isupper())
        lastLen = len(lastseq)

    exonSeqs = dict(exonSeqs)
    singleExon = len(exonSeqs) == 1
    if singleExon:
        defaultHeight = "7vw"
    else:
        defaultHeight = "5vw"

    print(
        """
        <script>
function toggleExonSeq(selectedValue) {
    const allExonDisplays = document.querySelectorAll('[name="exonDisplay"]');
    const exonHeights = document.getElementsByName('exonPamSeq');
    const exonFilterMsg = document.querySelectorAll('[name="exonFilterMsg"]');

    if (selectedValue === 'all') {
        for (const exonSeq of allExonDisplays) {
            exonSeq.style.display = 'block';
        }
        for (const exonHeight of exonHeights) {
            exonHeight.style.height = '%s';
        }
        for (const exonMsg of exonFilterMsg) {
            exonMsg.style.display = 'none';
        }
        $('#otTable tr.guideRow').show();

    } else {
        // hide all exons
        for (const displayElement of allExonDisplays) {
            displayElement.style.display = 'none';
        }
        for (const exonHeight of exonHeights) {
            exonHeight.style.height = '7vw';
        }
        for (const exonMsg of exonFilterMsg) {
            exonMsg.style.display = 'none';
        }

        // show all matching exons (e.g. both donor and acceptor for splicing)
        var groupNum = selectedValue.replace('exon', '');
        var targetElements = document.querySelectorAll('.exonGroup' + groupNum);
        for (const el of targetElements) {
            el.style.display = 'block';
        }

        // show the filter message
        var targetMsgs = document.querySelectorAll('.exonFilterMsg' + groupNum);
        for (const msg of targetMsgs) {
            msg.style.display = 'block';
        }

        // Filter guide table
        var exonClass = selectedValue.replace('exon', 'exon-');
        $('#otTable tr.guideRow').hide();
        $('#otTable tr.' + exonClass).show();
    }
}
        </script>
            """
        % defaultHeight
    )

    # scroll to the stop codon for insertions in C-terminal
    if insertPos == "Cter":
        print(
            """
            <script>
            window.addEventListener('DOMContentLoaded', function() {
              const div = document.getElementById('geneModel');
                div.scrollLeft = div.scrollWidth - div.clientWidth;
                });
            </script>
            """
        )

    if exonSeqs:
        if commonExons:
            commonStr = "common"
        else:
            commonStr = ""
        if not singleExon:
            showAllExonsButton = (
                """, or
            <button name="exonSelect" value="all" onclick=toggleExonSeq(this.value)
            style="width:110spx; height:25px"><small>show all %s exons targeted</small></button>
            """
                % commonStr
            )
        # for genes with a single exon, don't display le button "show all exons"
        else:
            showAllExonsButton = "."
        print(
            """<div style="margin-top:8px; margin-bottom:8px"> Below is the gene model. Click on an exon to show the corresponding guides%s</div>"""
            % showAllExonsButton
        )

        shownExonMsgs = set()
        for exonId in exonSeqs:
            if koMethod == "splicing":
                if exonId % 2 == 0:
                    originalExon = (exonId + 2) // 2
                else:
                    originalExon = (exonId + 1) // 2
                groupId = originalExon - 1
                exonLabel = originalExon
            else:
                groupId = exonId
                exonLabel = exonId + 1

            if groupId not in shownExonMsgs and not singleExon:
                shownExonMsgs.add(groupId)
                print(
                    """<div name="exonFilterMsg" class="exonFilterMsg%d" style="display:none;"> the results are currently filtered for Exon %s </div>"""
                    % (groupId, exonLabel)
                )

    print(
        """ <div id="geneModel" style="
          width: 1650px;
          overflow-x:scroll;
          display:flex;
          align-items:center;
          padding:5px;
          white-space: nowrap;
          height: 100px;
          border: 0.5px solid lightgray;
          border-radius: 8px;"> """
    )
    currentLen = 0

    if kiType and insertPos:
        tagSeqList = []
        tagBoxes = []
        if tagNames:
            for tagName in tagNames:
                tagSeq = taggingSeqs[tagName]
                tagSeqList.append((tagName, tagSeq))
            for tagName, tagSeq in tagSeqList:
                color = tagToColor[tagName]
                tagMouseOver = 'class="tooltipsterInteract" title="%s (%dbp)"' % (
                    tagName,
                    len(tagSeq),
                )
                if len(tagSeq) < 25:
                    tagName = tagName[0] + "."
                elif len(tagName) > 0 and len(tagSeq) / len(tagName) < 10:
                    tagName = tagName[0:3].strip(" ") + "."

                tagBox = """<div
                    %s
                    style="
                    width: %dpx;
                    height: 25px;
                    background-color: %s;
                    box-shadow: 0 0 0 0.5px gray;
                    clip-path: polygon(0 0, calc(100%% - 10px) 0, 100%% 50%%, calc(100%% - 10px) 100%%, 0 100%%);
                    text-align:center;
                    display: flex;
                    align-items: center;
                    justify-content: center;
                    flex-shrink: 0;
                    padding-right: 10px;
                    ">
                    %s
                    </div> """ % (
                    tagMouseOver,
                    len(tagSeq),
                    color,
                    tagName,
                )
                tagBoxes.append(tagBox)
        else:
            insertSeqText = "Insert sequence (%s bp)" % len(insertSeq)
            tagMouseOver = (
                'class="tooltipsterInteract" title="Custom insert sequence of %d bp"'
                % len(insertSeq)
            )

            if len(insertSeq) < 25:
                insertSeqText = "s."
            elif len(insertSeq) < 50:
                insertSeqText = "%sbp" % len(insertSeq)
            elif len(insertSeq) < 100:
                insertSeqText = "Insert (%s bp)" % len(insertSeq)
            elif len(insertSeq) < 150:
                insertSeqText = "Insert seq. (%s bp)" % len(insertSeq)

            tagBox = """<div
                    %s
                    style="
                    width:%dpx;
                    height: 25px;
                    background-color: #ffff66;
                    box-shadow: 0 0 0 0.5px gray;
                    clip-path: polygon(0 0, calc(100%% - 10px) 0, 100%% 50%%, calc(100%% - 10px) 100%%, 0 100%%);
                    text-align:center;
                    display: flex;
                    align-items: center;
                    justify-content: center;
                    flex-shrink: 0;
                    padding-right: 10px;
                    ">
                    %s
                    </div> """ % (
                tagMouseOver,
                len(insertSeq),
                insertSeqText,
            )
            tagBoxes.append(tagBox)

    for i, feature in enumerate(geneModel):

        featureType = feature[0]
        featureId = int(feature[1])
        featureIdOneBased = featureId + 1
        length = int(feature[2])

        if i == 0:
            print("<small>CDS start&nbsp&nbsp</small>")

        if i == 0 and kiType and insertPos == "Nter":
            for tagBox in tagBoxes:
                print(tagBox)
        if featureType == "exon":
            currentLen += length
            if i == 0:
                borderRadius = "0 8px 8px 0"
            elif i == len(geneModel) - 1:
                borderRadius = "8px 0 0 8px"
            else:
                borderRadius = "8px"
            if commonExons:
                exonStr = "common coding exon"
            else:
                exonStr = "coding exon"

            exonStrMed = "cod. exon"
            exonStrSmall = "exon"
            exonStrMin = "ex."

            fullExonTitle = "%s %d (%d bp)" % (exonStr, featureIdOneBased, length)
            if length >= 150 and not commonExons or length >= 200 and commonExons:
                exonText = "<small>%s </small>" % fullExonTitle
            if length >= 95 and not commonExons or length >= 135 and commonExons:
                exonText = "<small>%s %d</small>" % (exonStr, featureIdOneBased)
            elif length >= 75:
                exonText = "<small>%s %d</small>" % (exonStrMed, featureIdOneBased)
            elif length >= 45:
                exonText = "<small>%s%d</small>" % (exonStrSmall, featureIdOneBased)
            elif length >= 20:
                exonText = "<small>%s%d</small>" % (exonStrMin, featureIdOneBased)
            else:
                exonText = ""

            exonMouseOver = 'class="tooltipsterInteract" title="%s"' % fullExonTitle

            if exonSeqs:
                if koMethod in ["frameshift"]:
                    isSplittedExon = featureId + 1 == len(exonSeqs) and lastLen < length
                    isTargetExon = currentLen <= thirdLen and length >= GUIDELEN
                else:
                    isTargetExon = True
                    isSplittedExon = False

            # Only print the gene model, without highlighting target exons
            else:
                isSplittedExon = None
                isTargetExon = None

            if (isSplittedExon or isTargetExon) and featureId in exonSeqs:

                # don't make a button for exons in which no guide sequences were found
                if isSplittedExon:
                    if lastLen <= 50:
                        exonText = "<small>ex.%d</small>" % featureIdOneBased
                    print(
                        """<button name="exonSelect" value="exon%d" onclick=toggleExonSeq(this.value)
                        %s
                        style="
                        padding:0;
                        margin:0;
                        box-sizing: border-box;
                        width:%dpx;
                        height: 27px;
                        border: 0.5px solid gray;
                        border-right: 0px;
                        border-radius: 8px 0 0 8px;
                        text-align:center;
                        display: flex;
                        align-items: center;
                        justify-content: center;
                        flex-shrink: 0;
                        ">
                        %s
                        </button> """
                        % (featureId, exonMouseOver, lastLen, exonText)
                    )
                    print(
                        """<div
                        %s
                        style="
                        width:%dpx;
                        height: 25px;
                        border: 0.5px solid gray;
                        border-left: 0px;
                        border-radius: 0 8px 8px 0;
                        background-color: #eeeeee;
                        text-align:center;
                        display: flex;
                        align-items: center;
                        justify-content: center;
                        flex-shrink: 0;
                        "></div> """
                        % (exonMouseOver, (length - lastLen))
                    )

                elif isTargetExon:

                    print(
                        """<button name="exonSelect" value="exon%d" onclick=toggleExonSeq(this.value)
                        %s
                        style="
                        padding:0;
                        margin:0;
                        box-sizing: border-box;
                        width:%dpx;
                        height: 27px;
                        border: 0.5px solid gray;
                        border-radius: %s;
                        text-align:center;
                        display: flex;
                        align-items: center;
                        justify-content: center;
                        flex-shrink: 0;
                        ">
                        %s
                        </button> """
                        % (featureId, exonMouseOver, length, borderRadius, exonText)
                    )
            else:
                print(
                    """<div
                    %s
                    style="
                    width:%dpx;
                    height: 25px;
                    border: 0.5px solid gray;
                    border-radius: %s;
                    background-color: #eeeeee;
                    text-align:center;
                    display: flex;
                    align-items: center;
                    justify-content: center;
                    flex-shrink: 0;
                    ">
                    %s
                    </div> """
                    % (exonMouseOver, length, borderRadius, exonText)
                )

        elif featureType == "intron":

            print(
                """ <div class="block_1" style="
                  height: 1px;
                  border: 1px solid #C4C4C4;
                  width: 25px;
                  text-align:center;
                  flex-shrink: 0;
                  "></div> """
            )
            print("<small>/ %d bp /</small>" % length)
            print(
                """ <div class="block_1" style="
                  height: 1px;
                  border: 1px solid #C4C4C4;
                  width: 25px;
                  text-align:center;
                  flex-shrink: 0;
                  "></div> """
            )
        if i == len(geneModel) - 1:
            if kiType and insertPos == "Cter":
                for tagBox in tagBoxes:
                    print(tagBox)
            print("<small>&nbsp&nbspCDS end</small>")

    print("</div><br>")


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
    # print("""<a href='crispor.py'><img style='width:70px' src='%simage/2021-Logo-Do-3.jpg' alt='UCSC Logo'></a>""" % (HTMLPREFIX))

    print(
        """<div style=" width: 100%; display: flex; flex-direction: row; justify-content: space-between;">"""
    )
    print(
        """<div style='margin-top:10px;'><a href='crispor.py'>&nbsp;&larr; Back to CRISPOR homepage</a></div>"""
    )
    print("""<div style="display: flex; flex-direction: row;">""")
    print(
        """<a href='crispor.py'><img style='width:150px; margin-left:25px' src='%simage/2021-Logo-Do-3.jpg' alt='UCSC Logo'></a>"""
        % (HTMLPREFIX)
    )
    print(
        """
        <div style="margin-top: 6px;" >
        <a class="tooltipsterInteract" title="CELPHEDIA (The National Infrastructure for model organisms in health and biomedical research) is a national operational research infrastructure distributed over the French territory.<br>Its mission is to support academic and industrial scientific community to accelerate discoveries in biology and improve biomedical research. To this end, CELPHEDIA operates in 3 main activities with respect of ethical principles and animal welfare.<br>
        <ul>
            <li>Standardized service offers, in the areas of creation, functional exploration, archiving and distribution of animal models, necessary for fundamental research and preclinical approaches: rodents with the mouse as the leader, non-human primates and non-mammals including aquatic vertebrates.</li>
            <li>Research and development activity for new technological offers.</li>
            <li>Training courses adapted to users needs either for the use of animals in research with respect to institutional regulations or to develop specific technological skills.</li>
        </ul>" href='https://celphedia.eu/en/' target="_blank"><img style='width:150px; margin-left:25px;' src='%simage/logo_Celphedia.jpg' alt='Celphedia'></a>
        </div>
        """
        % (HTMLPREFIX)
    )
    print("</div>")
    print("</div>")

    print('<div id="bd">')
    print('<div class="centralpanel" style="margin-left:0px">')
    print('<div class="subpanel" style="background:transparent;box-shadow:none;">')
    print(
        '<div class="contentcentral" style="margin-left:0px; width:100%; background:none">'
    )


def printTeforBodyStart():
    print(
        """<div style="display: flex; flex-direction: row; gap:2%; align-items: center; padding: 12px; margin-left: 24px;">"""
    )
    # print """<a href='http://genome.ucsc.edu'><img style='vertical-align: top; height: 40px' src='%s/image/ucscBioinf.jpg' alt=''></a>""" % (HTMLPREFIX)

    '''
    print(
        """<a href='crispor.py'><img style='width:150px; margin-left:25px' src='%simage/2021-Logo-Do-3.jpg' alt='UCSC Logo'></a>"""
        % (HTMLPREFIX)
    )
    print(
        """
        <a class="tooltipsterInteract" title="CELPHEDIA (The National Infrastructure for model organisms in health and biomedical research) is a national operational research infrastructure distributed over the French territory.<br>Its mission is to support academic and industrial scientific community to accelerate discoveries in biology and improve biomedical research. To this end, CELPHEDIA operates in 3 main activities with respect of ethical principles and animal welfare.<br>
        <ul>
            <li>Standardized service offers, in the areas of creation, functional exploration, archiving and distribution of animal models, necessary for fundamental research and preclinical approaches: rodents with the mouse as the leader, non-human primates and non-mammals including aquatic vertebrates.</li>
            <li>Research and development activity for new technological offers.</li>
            <li>Training courses adapted to users needs either for the use of animals in research with respect to institutional regulations or to develop specific technological skills.</li>
        </ul>" href='https://celphedia.eu/en/' target="_blank"><img style='width:150px; margin-left:25px' src='%simage/logo_Celphedia.jpg' alt='Celphedia'></a>
        """ % (HTMLPREFIX)
    )

    print(
        """<a href='crispor.py'><img style='width:70px; margin-left:25px' src='%simage/logo_tefor.png' alt=''></a>"""
        % (HTMLPREFIX)
    )
    print("""<div style="margin-left: auto; margin-right: 5%%; display: flex; flex-direction: row; gap: 12px; min-width: 250px;">
                <a style="border-right: solid 1px lightgrey; padding: 0 8px 0 0;" target=_blank href="/manual/">Documentation</a>
                <a style="border-right: solid 1px lightgrey; padding: 0 8px 0 0;" href="https://academic.oup.com/nar/article/46/W1/W242/4995687">Citation</a>
                <a href='mailto:%s'>Contact us</a>
            </div>
        """ % contactEmail)

    print("</div>")
    '''

    print('<div id="bd">')
    print('<div class="centralpanel" style="margin-left:0px">')
    print('<div class="subpanel" style="background:transparent;box-shadow:none;">')
    print(
        '<div class="contentcentral" style="margin-left:0px; width:100%; background:none">'
    )


def printReleaseNote():
    print(
        '<div style="clear:both; text-align:center; margin-top: 12px; margin-bottom: 12px;">'
    )
    print(
        """<i>%s See <a href="doc/changes.html">Full list of changes</a></i><br>"""
        % releaseNote
    )
    print("</div>")


def printTeforBodyEnd():

    print(
        '<div style="clear:both; text-align:center; margin-left: auto; margin-right: 5%%; display: flex; flex-direction: row; gap: 12px; min-width: 250px;">'
    )
    print(
        """<div style="margin-left: auto; margin-right: auto;">
                Version %s -
                <a style="border-right: solid 1px lightgrey; padding: 0 8px 0 0; margin-right: 4px;" target=_blank href="/manual/">Documentation</a>
                <a style="border-right: solid 1px lightgrey; padding: 0 8px 0 0; margin-right: 4px;" href="https://academic.oup.com/nar/article/46/W1/W242/4995687">Citation</a>
                <a style="border-right: solid 1px lightgrey; padding: 0 8px 0 0; margin-right: 4px;" href="downloads/">Downloads / local installation</a>
                <a style="border-right: solid 1px lightgrey; padding: 0 8px 0 0; margin-right: 4px;" href="https://github.com/maximilianh/crisporWebsite/blob/master/LICENSE.txt">License</a>
                <a href='mailto:%s'>Contact us</a>
        """
        % (versionStr, contactEmail)
    )

    print("</div>")
    print(
        """
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

    </script>


<script>

// keep scroll position on reload
// from https://stackoverflow.com/questions/17642872

    document.addEventListener("DOMContentLoaded", function(event) {
        var scrollpos = localStorage.getItem('scrollpos');
        if (scrollpos) {
            window.scrollTo({top: parseInt(scrollpos),
                             left: 0,
                             behavior: 'auto'
                             });
            localStorage.removeItem('scrollpos');
        };
    });

    window.onbeforeunload = function(e) {
        localStorage.setItem('scrollpos', window.scrollY);
    };

    document.addEventListener('submit', function(e) {
        localStorage.setItem('scrollpos', window.scrollY);
    }, true);

</script>
"""
    )

    print(
        """
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');
  ga('create', 'UA-38389239-1', 'auto');
  ga('send', 'pageview');
</script>
"""
    )


pamIdRe = re.compile(r"s([0-9]+)([+-])g?([0-9]*)")


def intToExtPamId(pamId, multiseq=None, multipam=None, koMethod=None):
    "convert the internal pam Id like s20+ to the external one, like 21Forw. Handles multiseq/multipam prefixes."
    if multiseq or multipam:
        pamInfo = pamId.split(".")
        pamPrefix = pamInfo[0]
        pamIdCore = pamInfo[1]
    else:
        pamIdCore = pamId

    match = pamIdRe.match(pamIdCore)
    if not match:
        return pamId  # should not happen with current pamId formats

    pamPos, strand, rest = match.groups()
    if strand == "+":
        strDesc = "forw"
    else:
        strDesc = "rev"
    guideDesc = str(int(pamPos) + 1) + strDesc

    if multiseq:
        pamPrefix = int(pamPrefix)
        if koMethod == "frameshift":
            return "exon%d_%s" % (pamPrefix + 1, guideDesc)
        elif koMethod in ["excision", "promoter"]:
            if pamPrefix == 0:
                rowStr = "upstream"
            else:
                rowStr = "downstream"
            return "%s_%s" % (rowStr, guideDesc)
        elif koMethod == "splicing":
            if pamPrefix % 2 == 0:
                originalExon = (pamPrefix + 2) // 2
                return "exon%d_donorSite_%s" % (originalExon, guideDesc)
            else:
                originalExon = (pamPrefix + 1) // 2
                return "exon%d_acceptorSite_%s" % (originalExon, guideDesc)

    elif multipam:
        return "%s_%s" % (pamPrefix, guideDesc)
    else:
        return guideDesc


def concatGuideAndPam(guideSeq, pamSeq, pamPlusSeq=""):
    "return guide+pam or pam+guide, depending on pamIsFirst"
    if pamIsFirst:
        return pamPlusSeq + pamSeq + guideSeq
    else:
        return guideSeq + pamSeq + pamPlusSeq


def makeGuideHeaders(multipam=None):
    "return list of the headers of the guide output file"
    headers = list(tuple(guideHeaders))  # make a copy of the list

    logging.debug("active scoreNames: %s" % scoreNames)
    tableScoreNames = list(tuple(scoreNames))
    if not pamIsFirst:
        tableScoreNames.extend(mutScoreNames)

    headers.append("global-Score")

    if multipam:
        tableScoreNames.insert(0, "insertDistance")
        headers.append("insertDistance")

    for scoreName in tableScoreNames:
        if scoreName != "insertDistance":
            headers.append(scoreDescs[scoreName][0] + "-Score")

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


def iterGuideRows(
    guideData,
    addHeaders=False,
    seqId=None,
    satMutOpt=None,
    minSpec=None,
    minFusi=None,
    multipam=None,
    multiseq=None,
    koMethod=None,
):
    "yield rows from guide data. Need to know if for Cpf1 or not"
    headers, tableScoreNames = makeGuideHeaders(multipam=multipam)

    if satMutOpt:
        headers.append("Oligonucleotide")
        headers.append("AdapterHandle+PrimerFw")
        headers.append("AdapterHandle+PrimerRev")
        (
            oligoPrefix,
            oligoSuffix,
            primerFwPrefix,
            primerRevPrefix,
            batchId,
            genome,
            position,
            ampLen,
            tm,
        ) = satMutOpt

        otMatches = parseOfftargets(genome, batchId)

        guideData.sort(
            key=operator.itemgetter(3)
        )  # sort by position, makes more sense here

    if seqId != None:
        headers.insert(0, "#seqId")
    else:
        headers[0] = "#" + headers[0]

    headers.append("grafType")

    if addHeaders:
        yield headers

    for guideRow in guideData:
        (
            guideScore,
            guideCfdScore,
            effScores,
            pamStart,
            guideStart,
            strand,
            pamId,
            guideSeq,
            pamSeq,
            otData,
            otDesc,
            last12Desc,
            mutEnzymes,
            ontargetDesc,
            repCount,
            gcFrac,
            freeEnergy,
            doRecoding,
            cutUpstream,
            mainScore,
            beScoring,
            beOutcomes
        ) = guideRow
        if minSpec and guideScore < minSpec:
            continue

        if not effScorePass(effScores, minFusi):
            continue

        otCount = 0
        if otData != None:
            otCount = len(otData)

        guideDesc = intToExtPamId(
            pamId, multipam=multipam, multiseq=multiseq, koMethod=koMethod
        )

        fullSeq = concatGuideAndPam(guideSeq, pamSeq)
        row = [guideDesc, fullSeq, guideScore, guideCfdScore, otCount, ontargetDesc]
        row.append(mainScore)

        for scoreName in tableScoreNames:
            row.append(effScores.get(scoreName, "NotEnoughFlankSeq"))

        grafType = crisporEffScores.getGrafType(guideSeq)
        if grafType is None:
            grafType = "GrafOK"
        row.append(grafType)

        if satMutOpt:
            oligo = oligoPrefix + guideSeq + oligoSuffix
            row.append(oligo)

            chrom, start, end, strand, gene, isUnique = findOntargetPos(
                otMatches, pamId, position
            )
            lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq, ampRange, flankSeq, addTags = (
                designPrimer(genome, chrom, start, end, strand, 0, batchId, ampLen, tm)
            )

            fwPrimer = lSeq
            if fwPrimer != None:
                fwPrimer = primerFwPrefix + fwPrimer

            revPrimer = rSeq
            if revPrimer != None:
                revPrimer = primerRevPrefix + revPrimer

            row.append(fwPrimer)
            row.append(revPrimer)

        row = [str(x) for x in row]
        if seqId != None:
            row.insert(0, seqId)
        yield row


def iterOfftargetRows(guideData, addHeaders=False, skipRepetitive=True, seqId=None):
    "yield bulk offtarget rows for the tab-sep download file"
    otRows = []

    headers = list(offtargetHeaders)  # clone list
    if seqId:
        headers.insert(0, "seqId")

    if addHeaders:
        otRows.append(headers)

    skipCount = 0

    for guideRow in guideData:
        (
            guideScore,
            guideCfdScore,
            effScores,
            pamStart,
            guideStart,
            strand,
            pamId,
            guideSeq,
            pamSeq,
            otData,
            otDesc,
            last12Desc,
            mutEnzymes,
            ontargetDesc,
            repCount,
            gcFrac,
            freeEnergy,
            doRecoding,
            cutUpstream,
            mainScore,
            beScoring,
            beOutcomes
        ) = guideRow

        if otData != None:
            otCount = len(otData)
            if otCount > 5000 and skipRepetitive:
                skipCount += otCount
                continue

            for (
                otSeq,
                mitScore,
                cfdScore,
                editDist,
                pos,
                gene,
                alnHtml,
                inLinkage,
            ) in otData:
                gene = gene.replace(",", "_").replace(";", "-")
                chrom, start, end, strand = parsePos(pos)
                guideDesc = intToExtPamId(pamId)
                mismStr = highlightMismatches(guideSeq, otSeq, len(pamSeq))
                fullSeq = concatGuideAndPam(guideSeq, pamSeq)
                row = [
                    guideDesc,
                    fullSeq,
                    otSeq,
                    mismStr,
                    editDist,
                    mitScore,
                    cfdScore,
                    chrom,
                    start,
                    end,
                    strand,
                    gene,
                ]
                if seqId:
                    row.insert(0, seqId)
                row = [str(x) for x in row]
                otRows.append(row)

    if skipCount != 0:
        newRow = [""] * len(headers)
        newRow[0] = (
            "# %d off-targets are not shown: guides with more than 5000 off-targets were considered too repetitive"
            % skipCount
        )
        otRows.insert(0, newRow)

    return otRows


def writeHtmlTable(rows, outFile):
    "write list of rows to outFile in html format"
    outFile.write(
        '<link rel="stylesheet" href="https://unpkg.com/purecss@1.0.0/build/pure-min.css" integrity="sha384-nn4HPE8lTHyVtfCBi5yW9d20FjT8BJwUXyWZT9InLYax14RDjBj46LmSztkmNP9w" crossorigin="anonymous">\n'
    )
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
    "write table to ofh, currently writes xls files as tsv"
    if fileFormat == "tsv" or fileFormat == "xls" or fileFormat == "csv":
        sep = "\t"
        if fileFormat == "csv":
            sep = ","
        for row in rows:
            ofh.write(sep.join(row))
            ofh.write("\n")
        ofh.close()
    else:
        writeHtmlTable(rows, ofh)


def xlsWrite(
    rows,
    title,
    outFile,
    colWidths,
    fileFormat,
    seq,
    org,
    pam,
    position,
    batchId,
    optFields=None,
    multipam=None,
):
    """given rows, writes a XLS binary stream to outFile, if xlwt is available
    Otherwise writes a tab-sep file.
    colWidths is a list of widths of columns, in Arial characters.
    """
    if xlwtLoaded and fileFormat in ["xls"]:
        seqStyle = xlwt.easyxf("font: name Courier New")
        charSize = 269  # see http://reliablybroken.com/b/2011/10/widths-heights-with-xlwt-python/
        wb = xlwt.Workbook()
        ws = wb.add_sheet(title)

        ws.write(0, 0, "# Name")
        ws.write(0, 1, batchName)
        ws.write(1, 0, "# Sequence")
        ws.write(1, 1, seq)
        if multipam:
            ws.write(3, 0, "#PAMS")
            # add BE pams for substitutions
            pamList = multiPamDict[multipam][0]

            ws.write(3, 1, ", ".join([pam for pam in pamList]))
        else:
            ws.write(3, 0, "# PAM")
            ws.write(3, 1, pam)
        ws.write(2, 0, "# Genome")
        ws.write(2, 1, org)
        ws.write(4, 0, "# Position")
        ws.write(4, 1, position)

        ws.write(5, 0, "# Version")
        # http://stackoverflow.com/questions/4530069/python-how-to-get-a-value-of-datetime-today-that-is-timezone-aware
        FORMAT = "%Y-%m-%dT%H:%M:%S%Z"
        dateStr = time.strftime(FORMAT, time.localtime())
        ws.write(5, 1, "CRISPOR %s, %s" % (versionStr, dateStr))

        ws.write(6, 0, "# Results")
        url = "http://crispor.gi.ucsc.edu/crispor.py?batchId=%s" % batchId
        # ws.write(6, 1, xlwt.Formula('HYPERLINK("%s";"Link")' % (url)))
        ws.write(6, 1, url)

        startRow = 7
        curRow = startRow
        if optFields is not None:
            for key, val in optFields.items():
                ws.write(curRow, 0, "# %s" % key)
                ws.write(curRow, 1, val)
                curRow += 1

        skipRows = curRow + 1

        seqCols = [1, 7, 8, 9]  # columns with sequences -> fixed width font

        for rowCount, row in enumerate(rows):
            if rowCount == 65534 - startRow:
                ws.write(
                    rowCount + skipRows,
                    0,
                    "WARNING: cannot write more than 65535 rows to an Excel file. Switch to .tsv format to get all off-targets.",
                )
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
                    ws.write(rowCount + skipRows, colCount, col, seqStyle)
                else:
                    ws.write(rowCount + skipRows, colCount, col)

        # set sizes in characters per column
        for colId, colWidth in enumerate(colWidths):
            ws.col(colId).width = charSize * colWidth

        try:
            outFile.flush()  # flush out the Content-type header
            wb.save(outFile.buffer)
        except:
            print("error")

    elif fileFormat == "html":
        writeHtmlTable(rows, outFile)
    else:
        # raw ASCII tsv output mode
        sep = "\t"
        if fileFormat == "csv":
            sep = ","
        for row in rows:
            outFile.write(sep.join(row))
            outFile.write("\n")

    outFile.flush()


def seqToGenbankLines(seq):
    """chunk sequence string into lines each with six parts of 10bp, return as a list
    >>> seqToGenbankLines("aacacacatggtacacactgactagctagctacgatccagtacgatcgacgtagctatcgatcgatcgatcgactagcta")
    ['aacacacatg gtacacactg actagctagc tacgatccag tacgatcgac gtagctatcg', 'atcgatcgat cgactagcta']
    """
    # first chunk into 10bp parts
    parts = [seq[i : i + 10] for i in range(0, len(seq), 10)]

    # put into lines of 6*10 bp
    lines = []
    for i in range(0, len(parts), 6):
        lines.append(" ".join(parts[i : i + 6]))
    return lines


def writeLn(fh, line, indent=None, doWrap=True):
    "write line to file, using \r\n"
    if indent == None:
        fh.write(line)
        fh.write("\r\n")
    else:
        if doWrap:
            lineSize = 80 - indent
        else:
            lineSize = 10000
        parts = [line[i : i + lineSize] for i in range(0, len(line), lineSize)]
        spacer = "".join(([" "] * indent))
        for p in parts:
            fh.write(spacer + p)
            fh.write("\r\n")


def genbankWrite(batchId, fileFormat, desc, seq, org, position, pam, guideData, ofh):
    "write a description of the current job in genbank format to ofh"
    if fileFormat == "serialcloner":
        # a bug in serial cloner means that we cannot use the linear format
        writeLn(
            ofh,
            """LOCUS       %s    %d bp      DNA     circular   1/1/17"""
            % (desc, len(seq)),
        )
    elif fileFormat == "snapgene":
        writeLn(
            ofh,
            "LOCUS       Exported                 239 bp ds-DNA     linear   SYN 22-MAR-2017",
        )
    else:
        writeLn(
            ofh,
            """LOCUS       %s    %d bp      DNA     linear   1/1/17"""
            % (desc, len(seq)),
        )

    batchUrl = "http://crispor.gi.ucsc.edu/crispor.py?batchId=%s" % batchId
    seqDesc1 = """Sequence exported from CRISPOR.org V%s""" % versionStr
    seqDesc2 = "Genome %s, position %s. View full CRISPOR results at %s" "" % (
        org,
        position,
        batchUrl,
    )

    if fileFormat in ["ape"]:
        seqDesc2 += " Features indicate PAMs. Click the little triangles next to the features to show scores and guide sequences."

    if fileFormat == "genomecompiler":
        # genomecompiler plasmid viewer can't show more than ~30 characters as the seq definition line
        writeLn(ofh, """DEFINITION  %s""" % desc)
    else:
        writeLn(ofh, """DEFINITION  %s""" % seqDesc1)
        writeLn(ofh, seqDesc2, indent=12)
        writeLn(ofh, """             Export for: %s""" % (fileFormat))

    # writeLn(ofh, """ACCESSION""")
    # writeLn(ofh, """VERSION""")
    writeLn(ofh, """SOURCE      %s""" % org)
    writeLn(ofh, """  ORGANISM  %s""" % org)

    if fileFormat in ["serialcloner"]:
        writeLn(ofh, """COMMENT     Serial Cloner Genbank Format""")
        writeLn(ofh, """COMMENT     SerialCloner_Type=DNA""")
        writeLn(
            ofh, """COMMENT     SerialCloner_Comments=%s""" % seqDesc1 + " " + seqDesc2
        )
        writeLn(ofh, """COMMENT     SerialCloner_Ends=0,0,,0,""")
    elif fileFormat in ["ape", "genomecompiler", "vnti"]:
        # ape does not show the definition line, only the comment block
        writeLn(ofh, """COMMENT     %s""" % seqDesc1)
        writeLn(ofh, seqDesc2, indent=12)

    if fileFormat == "snapgene":
        # this tells snapgene that the file is really in snapgene format
        writeLn(ofh, """KEYWORDS    snapgene3""")
        writeLn(ofh, """REFERENCE   1  (bases 1 to 239)""")
        writeLn(ofh, """  AUTHORS   .""")
        writeLn(ofh, """  TITLE     Direct Submission""")
        writeLn(
            ofh,
            """  JOURNAL   Exported Wednesday, Mar 22, 2017 from SnapGene Viewer 3.3.3""",
        )
        writeLn(ofh, """          http://www.snapgene.com""")

    writeLn(ofh, """FEATURES             Location/Qualifiers""")

    i = 1
    guideData.sort(key=operator.itemgetter(3))  # sort by position

    for guideRow in guideData:
        (
            guideScore,
            guideCfdScore,
            effScores,
            pamStart,
            guideStart,
            strand,
            pamId,
            guideSeq,
            pamSeq,
            otData,
            otDesc,
            last12Desc,
            mutEnzymes,
            ontargetDesc,
            repCount,
            gcFrac,
            freeEnergy,
            doRecoding,
            cutUpstream,
            mainScore,
            beScoring,
            beOutcomes
        ) = guideRow

        guideUrl = batchUrl + "&pamId=" + urllib.parse.quote(pamId) + "&pam=" + pam

        otCount = 0
        if otData != None:
            otCount = len(otData)

        fullSeq = concatGuideAndPam(guideSeq, pamSeq)

        if fileFormat in ["geneious", "snapgene"]:
            # this code annotates only the guide sequence, suggested by Alyce Chen
            start = guideStart + 1
            end = start + len(guideSeq) - 1
        else:
            # most viewers don't handle overlaps well. We highlight only the PAM in these cases
            if strand == "+":
                start = pamStart + 1
                end = start + len(pamSeq)
            else:
                start = pamStart + 1
                end = start + len(pamSeq)

        colorHex, colorName = scoreToColor(guideScore)
        guideName = intToExtPamId(pamId)
        urlLink = "full details/primers at %s" % guideUrl

        if pamIsSpCas9(pam):
            fusiScore = str(effScores.get("fusi", -1))
            crisprScanScore = str(effScores.get("crisprScan", -1))
            descStr = "%s: Spec %s, Eff %s/%s" % (
                guideName,
                guideScore,
                fusiScore,
                crisprScanScore,
            )
            guideSeqDescSeq = (
                "Guide %s MIT-Spec %s, Eff Doench2016 %s, Eff Mor.-Mat. %s"
                % (guideSeq, guideScore, fusiScore, crisprScanScore)
            )
            longDesc = (
                "MIT-Specificity score: %s, Efficiency Doench2016 = %s, Efficiency Moreno-Mateos = %s, guide sequence: %s, "
                "full details/primers at %s"
                % (guideScore, fusiScore, crisprScanScore, guideSeq, guideUrl)
            )
        else:
            mainEffName = scoreNames[0]
            descStr = "%s: Spec %s, Eff %s" % (
                guideName,
                guideScore,
                effScores.get(mainEffName, -1),
            )

            shortStrList = []
            longStrList = []
            for scoreName in scoreNames:
                scoreVal = str(effScores.get(scoreName, -1))
                shortStrList.append("Eff %s %s" % (scoreName, scoreVal))
                longStrList.append(
                    "Efficiency %s = %s" % (scoreDescs[scoreName][0], scoreVal)
                )

            guideSeqDescSeq = "Guide %s MIT-Spec %s, %s" % (
                guideSeq,
                guideScore,
                ", ".join(shortStrList),
            )
            longDesc = "MIT-Specificity score: %s, %s, guide sequence: %s, %s" % (
                guideScore,
                ", ".join(longStrList),
                guideSeq,
                urlLink,
            )

        featType = "misc_feature"

        if fileFormat in ["genomecompiler"]:
            # genomecompiler does not store colors in the genbank file.
            # we use the feature type to simulate colors
            if colorName == "green":
                featType = "promoter"
            elif colorName == "red":
                featType = "terminator"
            elif colorName == "yellow":
                featType = "misc_feature"

        if strand == "+":
            writeLn(ofh, """     %s    %d..%d""" % (featType, start, end))
        else:
            writeLn(ofh, """     %s    complement(%d..%d)""" % (featType, start, end))

        if fileFormat in ["serialcloner"]:
            # serialcloner has no label
            writeLn(ofh, '''                     /note="%s"''' % descStr)
            writeLn(
                ofh,
                """                     /SerialCloner_Color=&h%s"""
                % colorHex.replace("#", ""),
            )
            writeLn(ofh, """                     /SerialCloner_Show=True""")
            writeLn(ofh, """                     /SerialCloner_Protect=True""")
            writeLn(ofh, """                     /SerialCloner_Arrow=True""")
            writeLn(ofh, """                     /SerialCloner_Desc=%s""" % longDesc)
        elif fileFormat in ["ape"]:
            # ape can use multiple note lines
            writeLn(ofh, '''                     /locus_tag="%s"''' % descStr)
            writeLn(
                ofh, """                     /note=MIT Specificity: %s""" % guideScore
            )
            writeLn(
                ofh,
                """                     /note=Efficiency: Doench2016 %s Mor-Mat. %s"""
                % (str(effScores["fusi"]), str(effScores["crisprScan"])),
            )
            writeLn(ofh, """                     /ApEinfo_fwdcolor=%s""" % colorHex)
            writeLn(ofh, """                     /ApEinfo_revcolor=%s""" % colorHex)
            writeLn(
                ofh,
                """                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0""",
            )
            writeLn(ofh, """                     /note=Guide %s""" % guideSeq)
        elif fileFormat == "benchling":
            writeLn(ofh, '''                     /note="%s"''' % guideSeqDescSeq)
            writeLn(ofh, """                     /ApEinfo_fwdcolor=%s""" % colorHex)
            writeLn(ofh, """                     /ApEinfo_revcolor=%s""" % colorHex)
            writeLn(
                ofh,
                """                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0""",
            )
        elif fileFormat in ["genomecompiler"]:
            writeLn(ofh, '''                     /label="%s"''' % guideName)
            writeLn(ofh, '''                     /note="%s"''' % guideSeqDescSeq)
        elif fileFormat in ["snapgene"]:
            writeLn(ofh, '''                     /label="%s"''' % guideName)
            writeLn(ofh, '''                     /note="%s"''' % longDesc)
            if strand == "+":
                writeLn(
                    ofh,
                    '''                     /note="color: %s; direction: RIGHT"'''
                    % colorHex,
                )
            else:
                writeLn(
                    ofh,
                    '''                     /note="color: %s; direction: LEFT"'''
                    % colorHex,
                )
        # vector NTI treats attributes as a key-val list
        # lasergene shows all attributes, vector NTI is the most complete way to show all data
        elif fileFormat in ["vnti", "lasergene"]:
            writeLn(ofh, '''/label="%s"''' % guideName, indent=21)
            writeLn(ofh, '''/note="%s"''' % descStr, indent=21)
            writeLn(ofh, '''/MITSpecScore="%s"''' % str(guideScore), indent=21)
            # longDesc = "MIT-Specificity score: %s, Efficiency Doench2016 = %s, Efficiency Moreno-Mateos = %s, guide sequence: %s, full details/primers at %s" % (guideScore, str(effScores["fusi"]), str(effScores["crisprScan"]), guideSeq, guideUrl)
            writeLn(ofh, '''/Doench2016Eff="%s"''' % str(effScores["fusi"]), indent=21)
            writeLn(
                ofh, '''/Mor-MateosEff="%s"''' % str(effScores["crisprScan"]), indent=21
            )
            writeLn(ofh, '''/guide_sequence="%s"''' % guideSeq, indent=21)
            writeLn(ofh, '''/url="%s"''' % guideUrl, indent=21)
        elif fileFormat in ["geneious"]:
            writeLn(
                ofh, '''/label="%s"''' % guideName, indent=21, doWrap=False
            )  # geneious translates \n to spaces, breaks link
            writeLn(
                ofh, '''/MIT-spec_score="%s"''' % guideScore, indent=21, doWrap=False
            )
            writeLn(ofh, '''/guide_sequence="%s"''' % guideSeq, indent=21, doWrap=False)
            writeLn(ofh, '''/note="%s"''' % urlLink, indent=21, doWrap=False)
            if "crisprScan" in effScores:
                writeLn(
                    ofh,
                    '''/Efficiency="%s"''' % effScores["crisprScan"],
                    indent=21,
                    doWrap=False,
                )
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
    "write the http header for attachments"
    if doDownload:
        print("Content-Type: application/octet-stream")
        print('Content-Disposition: attachment; filename="%s"' % (fname))
    else:
        print("Content-type: text/html")
    print("")  # = end of http headers
    sys.stdout.flush()


def buildPoolOptions(barcodeId, custPrefix="", custSuffix=""):
    "return a list of pool settings and a dictionary with pool options"
    barcodeDict = dict(satMutBarcodes)
    barcodeLabel = barcodeDict[barcodeId]

    if int(barcodeId) == 0:
        barcodePre, barcodePost = "", ""
    else:
        barcodePre, barcodePost = barcodeLabel.split()[-1].split("/")

    optFields = OrderedDict()
    optFields["Subpool Barcode"] = barcodeLabel

    optFields["Subpool Prefix"] = barcodePre
    optFields["Subpool Suffix"] = barcodePost

    if custPrefix != "" and custSuffix != "":
        optFields["Custom Prefix"] = custPrefix
        optFields["Custom Suffix"] = custSuffix
        optFields["Oligonucl. structure"] = (
            "Subpool Prefix + Custom Prefix + sgRNA + Custom Suffix + Subpool Suffix"
        )

        fullPrefix = barcodePre + custPrefix
        fullSuffix = barcodePost + custSuffix
    else:

        prePrimer = "GGAAAGG"
        pre = "ACGAAACACCG"
        post = "GTTTTAGAGCTAGAAATA"
        postPrimer = "GCAAGTTAAAATAAGGC"
        fullPrefix = barcodePre + prePrimer + pre
        fullSuffix = post + postPrimer + barcodePost

        optFields["pLentiGuidePre primer"] = prePrimer
        optFields["pLentiGuidePre"] = pre

        optFields["pLentiGuidePost"] = post
        optFields["pLentiGuidePost primer"] = postPrimer

        optFields["Oligonucl. structure"] = (
            "Subpool Prefix + pLentiGuidePre Primer + pLentiGuidePre + sgRNA + pLentiGuide Post + pLentiGuidePost Primer + Subpool Suffix"
        )

    primerFwPrefix = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
    primerRevPrefix = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

    satMutOpt = [fullPrefix, fullSuffix, primerFwPrefix, primerRevPrefix]

    return satMutOpt, optFields


def writeSatMutFile(
    barcodeId, ampLen, tm, batchId, minSpec, minFusi, fileFormat, outFh
):
    "write saturating mutagenesis table of all guides, scores, primers and amplicons around them to outFh"
    seq, org, pam, position, guideData = readBatchAndGuides(batchId)

    satMutOpt, optFields = buildPoolOptions(barcodeId)
    satMutOpt.extend([batchId, org, position, ampLen, tm])

    guideRows = iterGuideRows(
        guideData,
        addHeaders=True,
        satMutOpt=satMutOpt,
        minSpec=minSpec,
        minFusi=minFusi,
    )
    xlsWrite(
        guideRows,
        "guides",
        outFh,
        [20, 28, 10, 10, 10, 10, 10, 60, 21, 21],
        fileFormat,
        seq,
        org,
        pam,
        position,
        batchId,
        optFields=optFields,
    )


def readBatchAndGuides(batchId, globEffScore=None):
    "parse the input file, the batchId-json file and the offtargets and link everything together"
    (
        seq,
        org,
        pam,
        position,
        extSeq,
        multiseq,
        koMethod,
        geneModel,
        koGeneId,
        multipam,
    ) = readBatchParams(batchId)
    chrom, _, _, _ = parsePos(position)
    pam = setupPamInfo(pam)
    uppSeq = seq.upper()

    startDict, endSet = findAllPams(uppSeq, pam)

    # effScoreFname = join(batchDir, batchId + ".effScores.tab")

    otMatches = parseOfftargets(org, batchId, chrom)
    effScores = readEffScores(batchId)
    guideData, guideScores, hasNotFound, pamIdToSeq = mergeGuideInfo(
        uppSeq,
        startDict,
        pam,
        otMatches,
        position,
        effScores,
        org=org,
        globEffScore=globEffScore,
    )
    return seq, org, pam, position, guideData


def writeOntargetAmpliconFile(
    outType, batchId, ampLen, tm, ofh, fileFormat="tsv", minSpec=0, minFusi=0
):
    """design primers with approx ampLen and tm around each guide's target.
    outType can be "primers" or "amplicons"
    """
    (
        inSeq,
        db,
        pamPat,
        position,
        extSeq,
        multiseq,
        koMethod,
        geneModel,
        koGeneId,
        multipam,
    ) = readBatchParams(batchId)
    chrom, _, _, _ = parsePos(position)
    otMatches = parseOfftargets(db, batchId, chrom)

    startDict, endSet = findAllPams(inSeq, pamPat)
    pamSeqs = list(flankSeqIter(inSeq, startDict, len(pamPat), True))

    allEffScores = readEffScores(batchId)
    guideData, guideScores, hasNotFound, pamIdToSeq = mergeGuideInfo(
        inSeq,
        startDict,
        pamPat,
        otMatches,
        position,
        allEffScores,
        sortBy="pos",
        org=db,
    )

    if outType == "primers":
        headers = [
            "#guideId",
            "forwardPrimer",
            "leftPrimerTm",
            "revPrimer",
            "revPrimerTm",
            "ampliconSequence",
            "guideSequence",
        ]
    else:
        headers = ["#guideId", "ampliconSequence", "guideSequence"]

    rows = []
    rows.append(headers)
    # ofh.write("\t".join(headers))
    # ofh.write("\n")

    # for pamId, pamStart, guideStart, strand, guideSeq, pamSeq, pamPlusSeq in pamSeqs:
    for (
        guideScore,
        guideCfdScore,
        effScores,
        pamStart,
        guideStart,
        guideStrand,
        pamId,
        guideSeq,
        pamSeq,
        otData,
        otDesc,
        last12Desc,
        mutEnzymes,
        ontargetDesc,
        repCount,
        gcFrac,
        freeEnergy,
        doRecoding,
        cutUpstream,
        mainScore,
        beScoring,
        beOutcomes
    ) in guideData:

        if guideScore < minSpec:
            continue
        if not effScorePass(effScores, minFusi):
            continue

        chrom, start, end, strand, gene, isUnique = findOntargetPos(
            otMatches, pamId, position
        )
        # effScores = allEffScores.get(pamId, None)

        note = ""
        if not isUnique:
            note = "warning: guide has no unique match in genome"

        lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq, ampRange, flankSeq, addTags = (
            designPrimer(db, chrom, start, end, strand, 0, batchId, ampLen, tm)
        )

        pamName = intToExtPamId(pamId)
        if outType == "primers":
            row = [pamName, lSeq, lTm, rSeq, rTm, targetSeq, guideSeq]
        else:
            row = [pamName, targetSeq, guideSeq]

        row = [str(x) for x in row]
        rows.append(row)

    writeTable(fileFormat, rows, ofh)


def writeTargetSeqs(guideData, ofh, minSpec=None, minFusi=None):
    "write the guide sequences and their pam to ofh"
    for guideRow in guideData:
        (
            guideScore,
            guideCfdScore,
            effScores,
            pamStart,
            guideStart,
            strand,
            pamId,
            guideSeq,
            pamSeq,
            otData,
            otDesc,
            last12Desc,
            mutEnzymes,
            ontargetDesc,
            repCount,
            gcFrac,
            freeEnergy,
            doRecoding,
            cutUpstream,
            mainScore,
            beScoring,
            beOutcomes
        ) = guideRow
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
    "write the guide locations and their sequences to ofh, in the format for crisprSurf"
    seqChrom, seqStart, seqEnd, seqStrand = parsePos(position)

    rows = []
    header = ["Chr", "Start", "Stop", "sgRNA_Sequence", "Strand", "sgRNA_Type"]
    rows.append(header)

    for guideRow in guideData:
        (
            guideScore,
            guideCfdScore,
            effScores,
            pamStart,
            guideStart,
            guideStrand,
            pamId,
            guideSeq,
            pamSeq,
            otData,
            otDesc,
            last12Desc,
            mutEnzymes,
            ontargetDesc,
            repCount,
            gcFrac,
            freeEnergy,
            doRecoding,
            cutUpstream,
            mainScore,
            beScoring,
            beOutcomes
        ) = guideRow
        if minSpec and guideScore < minSpec:
            continue
        if not effScorePass(effScores, minFusi):
            continue

        chromStart, chromEnd, _, _, chromStrand = mapToGenome(
            seqStart, seqStrand, pamStart, guideStart, guideStrand
        )

        row = [seqChrom, chromStart, chromEnd, chromStrand, guideSeq, "observation"]
        row = [str(x) for x in row]
        rows.append(row)

    # crisprsurf takes only csv files
    if fileFormat == "tsv":
        fileFormat = "csv"
    writeTable(fileFormat, rows, ofh)


def fastaWrite(seqId, seq, fh, width=80):
    """output fasta seq to file object, break to 80 char width"""
    fh.write(">" + seqId + "\n")
    if len(seq) > width:
        last = 0
        for l in range(width, len(seq), width):
            fh.write(seq[last:l])
            fh.write("\n")
            last = l
        fh.write(seq[last : len(seq)])
    else:
        fh.write(seq)
    fh.write("\n")


def downloadDonor(params):

    donorSeq = params["donorSeq"]
    recodedDonorSeq = params.get("recodedDonorSeq")
    guideSeq = params["guideSeq"]
    pamId = params["pamId"]
    donorName = params["donorName"]

    doubleNicking = params.get("doubleNicking")
    if doubleNicking:
        revPamId = params["revPamId"]
        revGuideSeq = params["revGuideSeq"]
        fwPamId = params["fwPamId"]
        fwGuideSeq = params["fwGuideSeq"]

    writeHttpAttachmentHeader(donorName + ".fa")
    if doubleNicking:
        fastaWrite(revPamId, revGuideSeq, sys.stdout)
        print(" ")
        fastaWrite(fwPamId, fwGuideSeq, sys.stdout)
    else:
        fastaWrite(pamId, guideSeq.upper(), sys.stdout)
    print(" ")
    if recodedDonorSeq:
        fastaWrite(donorName + "_Original", donorSeq, sys.stdout)
        print(" ")
        fastaWrite(donorName + "_Recoded", recodedDonorSeq, sys.stdout)
    else:
        fastaWrite(donorName, donorSeq, sys.stdout)


def downloadFile(params):
    " "
    global scoreNames

    batchId = params["batchId"]
    batchInfo = readBatchAsDict(batchId)
    if batchInfo:
        multiseq = batchInfo.get("multiseq")
        multipam = batchInfo.get("multipam")
        koMethod = batchInfo.get("koMethod")
    else:
        multiseq = None
        multipam = None
        koMethod = None

    if multiseq:
        koGeneId = batchInfo.get("koGeneId")
        exonSeqs, org, pam, exonPosStr, guideData = KoResultsPage(
            params, batchId, koGeneId, download=True
        )
        if koMethod in ["excision", "promoter"]:
            label = {0: "upstream", 1: "downstream"}
            seq = ", ".join(["%s: %s" % (label[int(s[0])], s[1]) for s in exonSeqs])
        elif koMethod == "splicing":
            seqList = []
            for s in exonSeqs:
                exonId = s[0]
                exonSeq = s[1]
                if exonId % 2 == 0:
                    originalExon = (exonId + 2) // 2
                    seqList.append(
                        "exon %s (splicing donor site): %s" % (originalExon, exonSeq)
                    )
                else:
                    originalExon = (exonId + 1) // 2
                    seqList.append(
                        "exon %s (splicing acceptor site): %s" % (originalExon, exonSeq)
                    )
            seq = ", ".join(seqList)
        else:
            seq = ", ".join(["exon %s: %s" % (int(s[0]) + 1, s[1]) for s in exonSeqs])

        position = koGeneId if koGeneId else ""
    elif multipam:
        seq, org, pam, position, guideData = KiResultsPage(
            params, batchId, download=True
        )

    else:
        globEffScore = params.get("globEffScore")
        seq, org, pam, position, guideData = readBatchAndGuides(
            batchId, globEffScore=globEffScore
        )

    if batchName != "":
        queryDesc = batchName + "_"
    else:
        queryDesc = ""

    if position is None or position == "?":
        queryDesc += org + "-unknownLoc"
    else:
        queryDesc += org + "-" + position.strip(":+-").replace(":", "-")
        # print org, position, queryDesc

    fileType = params["download"]
    # the assistant cannot set the argument as it uses multi-submit buttons
    if fileType == "useGet":
        for key in params:
            if key.startswith("get-"):
                fileType = key.split("-")[1]
                break

    fileFormat = params.get("format", "tsv")
    if not fileFormat in ["tsv", "xls", "html"]:
        errAbort("Invalid fileFormat argument")

    doDownload = True
    if fileFormat == "html":
        doDownload = False

    if fileType == "guides":
        if cgiParams.get("showAllScores", "0") == "1":
            if not pamIsCpf1(pam) and not pamIsSaCas9(pam):
                scoreNames = allScoreNames
        writeHttpAttachmentHeader("guides_%s.%s" % (queryDesc, fileFormat), doDownload)
        xlsWrite(
            iterGuideRows(
                guideData,
                addHeaders=True,
                multipam=multipam,
                multiseq=multiseq,
                koMethod=koMethod,
            ),
            "guides",
            sys.stdout,
            [9, 28, 10, 10],
            fileFormat,
            seq,
            org,
            pam,
            position,
            batchId,
            multipam=multipam,
        )

    elif fileType == "offtargets":
        writeHttpAttachmentHeader(
            "offtargets_%s.%s" % (queryDesc, fileFormat), doDownload
        )
        skipRepetitive = fileFormat == "xls"
        otRows = list(
            iterOfftargetRows(guideData, addHeaders=True, skipRepetitive=skipRepetitive)
        )
        doReverse = not pamIsFirst
        otRows.sort(key=operator.itemgetter(4), reverse=doReverse)
        xlsWrite(
            otRows,
            "offtargets",
            sys.stdout,
            [9, 28, 28, 5],
            fileFormat,
            seq,
            org,
            pam,
            position,
            batchId,
        )

    elif fileType == "targetSeqs":
        writeHttpAttachmentHeader("targetSeqs_%s.txt" % (queryDesc), doDownload)
        minSpec = cgiGetNum(params, "minSpec", 0)
        minFusi = cgiGetNum(params, "minFusi", 0)
        writeTargetSeqs(guideData, sys.stdout, minSpec=minSpec, minFusi=minFusi)

    elif fileType == "targetLocs":
        writeHttpAttachmentHeader("targetLocs_%s.csv" % (queryDesc), doDownload)
        minSpec = cgiGetNum(params, "minSpec", 0)
        minFusi = cgiGetNum(params, "minFusi", 0)
        writeTargetLocs(
            position,
            guideData,
            sys.stdout,
            fileFormat,
            minSpec=minSpec,
            minFusi=minFusi,
        )

    elif fileType == "amplicons":
        # write amplicons of all off-targets for a single guide
        fname = makeCrispressoFname(batchName, batchId)
        writeHttpAttachmentHeader(fname, doDownload)
        pamId = cgiGetStr(params, "pamId")
        writeAmpliconFile(params, batchId, pamId, sys.stdout)

    elif fileType == "ontargetAmplicons":
        # design primers around all targets in input sequence
        writeHttpAttachmentHeader("ontargetAmplicons_%s.tsv" % (queryDesc), doDownload)
        ampLen = cgiGetNum(params, "ampLen", 140)
        tm = cgiGetNum(params, "tm", 60)
        minSpec = cgiGetNum(params, "minSpec", 0)
        minFusi = cgiGetNum(params, "minFusi", 0)
        writeOntargetAmpliconFile(
            "amplicons", batchId, ampLen, tm, sys.stdout, fileFormat, minSpec, minFusi
        )

    elif fileType == "ontargetPrimers":
        writeHttpAttachmentHeader("ontargetPrimers_%s.tsv" % (queryDesc), doDownload)
        ampLen = cgiGetNum(params, "ampLen", 140)
        tm = cgiGetNum(params, "tm", 60)
        minSpec = cgiGetNum(params, "minSpec", 0)
        minFusi = cgiGetNum(params, "minFusi", 0)
        writeOntargetAmpliconFile(
            "primers", batchId, ampLen, tm, sys.stdout, fileFormat, minSpec, minFusi
        )

    elif fileType == "satMut":
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
        writeSatMutFile(
            barcodeId, ampLen, tm, batchId, minSpec, minFusi, fileFormat, sys.stdout
        )

    elif fileType in [
        "serialcloner",
        "ape",
        "genomecompiler",
        "fasta",
        "benchling",
        "snapgene",
        "genbank",
        "vnti",
        "lasergene",
        "geneious",
    ]:
        fileFormat = params["download"]
        ext = "gb"
        if fileFormat == "serialcloner":
            ext = "xdna"
        elif fileFormat == "ape":
            ext = "ape"
        elif fileFormat == "snapgene":
            ext = "dna"
        elif fileFormat == "fasta":
            ext = "fa"
        elif fileFormat == "vnti":
            ext = "gb"
        fileName = "crispor_%s-%s.%s" % (queryDesc, fileFormat, ext)

        if fileFormat != "genomecompiler":
            writeHttpAttachmentHeader(fileName)
        else:
            print("Content-type: text/plain\n")

        if fileFormat != "fasta":
            genbankWrite(
                batchId,
                fileFormat,
                queryDesc,
                seq,
                org,
                position,
                pam,
                guideData,
                sys.stdout,
            )
        else:
            fastaWrite("crispor-" + queryDesc, seq, sys.stdout)

    else:
        errAbort("invalid value for download parameter, fileType=%s" % fileType)


def makeCrispressoFname(batchName, batchId):
    fnameDesc = ["crisporAmplicons"]
    if batchName != "":
        fnameDesc.append(batchName)
    fnameDesc.append(batchId)
    fname = "_".join(fnameDesc) + ".txt"
    return fname


def designOfftargetPrimers(
    inSeq,
    db,
    pam,
    position,
    extSeq,
    pamId,
    ampLen,
    primerLen,
    tm,
    otMatches,
    exonId=None,
    pamFullName=None,
):
    "return a list of off-target primers sorted by CFD score"
    targetChrom, targetStart, targetEnd, strand = parsePos(position)
    chromSizes = parseChromSizes(db)
    pam = setupPamInfo(pam)

    (
        guideSeq,
        pamSeq,
        pamPlusSeq,
        guideSeqWPam,
        guideStrand,
        guideSeqHtml,
        guideStart,
        guideEnd,
    ) = findGuideSeq(inSeq, pam, pamId, exonId=exonId, pamFullName=pamFullName)

    # get the coords
    coords = []
    names = []
    nameToOtScoreSeq = {}
    flankSize = 1000
    for mismCount, otMatchRows in otMatches.items():
        for otMatch in otMatchRows:
            (
                chrom,
                start,
                end,
                otSeq,
                strand,
                segType,
                segName,
                totalAlnCount,
                fromXaTag,
            ) = otMatch
            if chrom == targetChrom and start >= targetStart and end <= targetEnd:
                prefix = "ontarget_"
            else:
                prefix = ""
            segDesc = segTypeConv.get(
                segType, ""
            )  # some genomes do not have descriptions
            name = (
                "%(prefix)smm%(mismCount)d_%(segDesc)s_%(segName)s_%(chrom)s_%(start)d"
                % locals()
            )
            if start - flankSize < 0 or end + flankSize > chromSizes[chrom]:
                print(
                    (
                        "Cannot design primer for %s, too close to chromosome boundary"
                        % name
                    )
                )
            else:
                coords.append((chrom, start - flankSize, end + flankSize, name))
                otScore = calcCfdScore(guideSeq, otSeq)
                nameToOtScoreSeq[name] = (otScore, otSeq)

    # coords -> sequences
    flankSeqs = getGenomeSeqsBin(db, coords, doRepeatMask=True)
    targetSeqs = [(x[3], x[6]) for x in flankSeqs]  # strip coords, keep name+seq
    nameToSeq = dict(targetSeqs)

    # sequences -> primers
    # check the input parameters: ampLen, tm
    ampMin = ampLen - 110
    ampRange = "%d-%d" % (ampMin, ampLen)

    primers = runPrimer3(
        targetSeqs, flankSize, GUIDELEN + len(pamSeq), ampRange, tm, {}, primerLen
    )

    # sort primers by CFD off-target score
    scoredPrimers = []
    for name, primerInfo in primers.items():
        score, otSeq = nameToOtScoreSeq[name]
        scoredPrimers.append((score, name, primerInfo))
    scoredPrimers.sort(reverse=True)

    return scoredPrimers, nameToSeq, nameToOtScoreSeq, guideSeqHtml


def makeCrispressoOfftargetRows(scoredPrimers, nameToSeq, nameToOtScoreSeq):
    "yield the crispresso off-target rows"
    for score, seqName, primerInfo in scoredPrimers:
        flankSeq = nameToSeq[seqName]
        score, otSeq = nameToOtScoreSeq[seqName]
        lSeq, lTm, lPos, rSeq, rTm, rPos = primerInfo
        if lSeq is None:
            continue
        ampSeq = flankSeq[lPos : rPos + 1]
        row = [seqName, ampSeq, otSeq, "NA", "NA"]
        yield row


def writeAmpliconFile(params, batchId, pamId, outFh):
    "create the table of off-target amplicons for crispressoPooled and write to outFh"
    ampLen = cgiGetNum(params, "ampLen", 140)
    tm = cgiGetNum(params, "tm", 60)
    primerLen = cgiGetNum(params, "primerLen", 20)

    (
        inSeq,
        db,
        pam,
        position,
        extSeq,
        multiseq,
        koMethod,
        geneModel,
        koGeneId,
        multipam,
    ) = readBatchParams(batchId)

    if multipam:
        pamFullName = pamId.split(".")[0]
        pam = setupPamInfo(pamFullName)

    pamOtMatches = parseOfftargets(db, batchId)
    otMatches = pamOtMatches[pamId]

    scoredPrimers, nameToSeq, nameToOtScoreSeq, guideSeqHtml = designOfftargetPrimers(
        inSeq,
        db,
        pam,
        position,
        extSeq,
        pamId,
        ampLen,
        primerLen,
        tm,
        otMatches,
        pamFullName=pamFullName,
    )

    for row in makeCrispressoOfftargetRows(scoredPrimers, nameToSeq, nameToOtScoreSeq):
        outFh.write("\t".join(row))
        outFh.write("\n")


def otPrimerPage(params):
    "print a page with all primers for all off-targets"
    # initialize everything
    batchId, pamId = params["batchId"], params["pamId"]

    ampLen = cgiGetNum(params, "ampLen", 150)
    tm = cgiGetNum(params, "tm", 60)
    primerLen = cgiGetNum(params, "primerLen", 20)

    (
        inSeq,
        db,
        pam,
        position,
        extSeq,
        multiseq,
        koMethod,
        geneModel,
        koGeneId,
        multipam,
    ) = readBatchParams(batchId)

    if multipam:
        pamFullName = pamId.split(".")[0]
        exonId = None
        pam = setupPamInfo(pamFullName)
    elif multiseq:
        pamFullName = None
        exonId = int(pamId.split(".")[0])
        batchInfo = readBatchAsDict(batchId)
        exonSeqs = batchInfo["exonSeqs"]
        for (exonIdx, exonSeq), (posStrIdx, exonPosStr) in zip(exonSeqs, multiseq):
            if exonId == exonIdx and exonId == posStrIdx:
                inSeq = exonSeq
                position = exonPosStr
                break
    else:
        exonId = None
        pamFullName = None

    pamOtMatches = parseOfftargets(db, batchId)
    otMatches = pamOtMatches[pamId]

    onlyExons = params.get("onlyExons")
    onlyChrom = params.get("onlyChrom")
    if onlyExons or onlyChrom:
        targetChrom, _, _, _ = parsePos(position)
        newOtMatches = {}
        for editDist, matches in otMatches.items():
            filtered = []
            for m in matches:
                otChrom = m[0]
                segType = m[5]
                if onlyExons and segType != "ex":
                    continue
                if onlyChrom and otChrom != targetChrom:
                    continue
                filtered.append(m)
            if filtered:
                newOtMatches[editDist] = filtered
        otMatches = newOtMatches

    scoredPrimers, nameToSeq, nameToOtScoreSeq, guideSeqHtml = designOfftargetPrimers(
        inSeq,
        db,
        pam,
        position,
        extSeq,
        pamId,
        ampLen,
        primerLen,
        tm,
        otMatches,
        exonId=exonId,
        pamFullName=pamFullName,
    )

    # primers -> table rows
    primerTable = []
    for otScore, seqName, primerInfo in sorted(scoredPrimers, reverse=True):
        if otScore is None:
            otScore = "noOfftargetScore"
        else:
            otScore = "%.3f" % otScore
            otScore = otScore[:4]
        lSeq, lTm, lPos, rSeq, rTm, rPos = primerInfo
        if batchName != "":
            baseName = batchName + "_" + seqName
        else:
            baseName = seqName

        # Nextera handle sequences, from Matt Canver
        prefixForw = "<b>TCGTCGGCAGCGTC</b>"
        prefixRev = "<b>GTCTCGTGGGCTCGG</b>"

        if lSeq is None:
            primerTable.append(
                (baseName + "_F", "Primer3: not found at this Tm", "N.d.", otScore)
            )
        else:
            primerTable.append((baseName + "_F", prefixForw + lSeq, lTm[:4], otScore))

        if rSeq is None:
            primerTable.append(
                (baseName + "_R", "Primer3: not found at this Tm", "N.d.", otScore)
            )
        else:
            primerTable.append((baseName + "_R", prefixRev + rSeq, rTm[:4], otScore))

    prefix = ""
    if batchName != "":
        prefix = "<strong>%s :</strong>" % batchName

    printBackLink()
    print("Contents:<br>")
    print("<ul>")
    print("<li><a href='#primers'>PCR Primers</a>")
    print("<li><a href='#ampliconSeqs'>Amplicon sequences</a></li>")
    print("<li><a href='#crispresso'>Crispresso support</a></li>")
    print("</ul>")
    print("<hr>")

    print(
        (
            "<h2 id='primers'>%sPCR primers for off-targets of %s</h2>"
            % (prefix, guideSeqHtml)
        )
    )

    if onlyExons and onlyChrom:
        print(
            "<p>Only primers to amplify off-targets in exons within the chromosome of the target sequence are displayed</p>"
        )
    elif onlyExons:
        print("<p>Only primers to amplify off-targets in exons are displayed</p>")
    elif onlyChrom:
        print(
            "<p>Only primers to amplify off-targets in the chromosome of the target sequence are displayed</p>"
        )

    print(
        "<p>In the table below, Illumina Nextera Handle sequences have been added and highlighted in bold. Primers for the on-target have been added for convenience. The table below is sorted by the CFD off-target score. Sites with very low CFD scores &lt; 0.02 are unlikely to be cleaved, see our study <a href='http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2'>Haeussler et al. 2016</a>, Figure 2.</p>"
    )
    print(
        "<p>In the protocol by Matthew Canver, Harvard, two PCRs are run: one PCR to amplify the potential off-target, then a second PCR to extend the handles with Illumina barcodes. Please <a href='downloads/prot/canverProtocol.pdf'>click here</a> to download the protocol. Alternatively, you can have a look at <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3988262/#S8'>Fu et al, 2014</a>.</p>"
    )
    print(
        "<p>If a primer was not found, the reason is usually that the region around the off-target is too repetitive. "
        "To avoid unspecific primers, all repeats are masked for the primer design (not for off-target search). "
        "If you think that we should change the parameters here or should use different primer3 settings, please let "
        "us know.</p>"
    )

    if pamId not in pamOtMatches:
        errAbort("pamId %s not valid" % pamId)

    print(("""<form id="paramForm" action="%s" method="GET">""" % basename(__file__)))

    printAmpLenAndTm(ampLen, primerLen, tm)
    printHiddenFields(
        params,
        {
            "batchId": batchId,
            "pamId": pamId,
            "otPrimers": "1",
            "ampLen": None,
            "tm": None,
            "primerLen": None,
        },
    )
    print("""<input type="submit" name="submit" value="Update">""")
    print("</form>")
    print("<p>")

    printPrimerTable(primerTable, withTm=True, withScore=True)

    print("<h2 id='ampliconSeqs'>Off-target amplicon sequences with primers</h2>")
    print(
        "<p>These only list off-targets that have primers in the table above. Primers underlined, off-targets in bold.</p>"
    )

    print("<table>")
    for score, seqName, primerInfo in scoredPrimers:
        print("<tr>")
        flankSeq = nameToSeq[seqName]
        lSeq, lTm, lPos, rSeq, rTm, rPos = primerInfo
        if lSeq is None:
            continue
        ampSeq = flankSeq[lPos : rPos + 1]
        ulCoords = [(0, len(lSeq)), (rPos - lPos - len(rSeq), rPos - lPos + 1)]
        boldCoords = [(1000 - lPos, 1000 - lPos + GUIDELEN + len(pam))]
        print(("<td>%s</td> " % seqName))
        print("<td><tt>")
        print(markupSeq(ampSeq, ulCoords, boldCoords))
        print("</tt></td>")
        print("</tr>")
    print("<table>")

    # print("<h2>Input file for <a href='http://crispresso.rocks/'>Crispresso</a></h2>")
    print("<h2 id='crispresso'>Input file for Crispresso</h2>")
    print(
        "<p><a href='http://crispresso.rocks/'>Crispresso</a>, written by Luca Pinello, is a software package to quantify the Cas9-induced mutations on off- or on-targets.</p>"
    )

    downloadUrl = cgiGetSelfUrl({"download": "amplicons"})
    print(
        (
            "<p><a href='%s'>Click here</a> to download an amplicon input file for Crispresso. For each off-target, it includes the off-target name, its PCR amplicon and the guide sequence. Keep a copy of this file.</p>"
            % downloadUrl
        )
    )
    print(
        "After sequencing, run CRISPRessoPooled. The tool will map the reads to the amplicons and analyse the mutations:<br>"
    )
    fname = makeCrispressoFname(batchName, batchId)
    print(
        (
            "<tt>CRISPRessoPooled -r1 Reads1.fastq.gz -r2 Reads2.fastq.gz -f %s --name MY_EXPERIMENT</tt></p>"
            % fname
        )
    )


def printBackLink(toDonorPage=False, returnUrl=False):
    """print a link back to the main batch page
    optionally, print a link back to the donor design page
    or returns the url
    """

    newParams = {}
    newParams["batchId"] = cgiParams.get("batchId", "")
    if "batchName" in cgiParams:
        newParams["batchName"] = cgiParams["batchName"]
    if toDonorPage:
        newParams["pamId"] = urllib.parse.quote(str(cgiParams["pamId"]))
        newParams["pam"] = cgiParams["pam"]
        newParams["doRecoding"] = cgiParams["doRecoding"]
        newParams["cutUpstream"] = cgiParams["cutUpstream"]
        newParams["insertDistance"] = cgiParams["insertDistance"]
        linkText = "return to donor DNA design"
    else:
        linkText = "return to the list of all guides"

    paramStrs = ["%s=%s" % (key, val) for key, val in newParams.items()]
    paramStr = "&".join(paramStrs)
    url = basename(__file__) + "?" + paramStr

    if returnUrl:
        return url
    else:
        print("<p><a href='%s'>&larr; %s </a></p>" % (url, linkText))


def microHomPage(params):
    "show the Bae et al microhomology sequences"
    printBackLink()
    batchId, pamId = params["batchId"], params["pamId"]
    (
        inSeq,
        db,
        pam,
        position,
        extSeq,
        multiseq,
        koMethod,
        geneModel,
        koGeneId,
        multipam,
    ) = readBatchParams(batchId)
    pam = setupPamInfo(pam)

    if multiseq:
        batchInfo = readBatchAsDict(batchId)
        # get the sequence corresponding to exonId
        exonSeqs = batchInfo["exonSeqs"]
        pamExonId = int(pamId.split(".")[0])
        for exonId, seq in exonSeqs:
            if exonId == pamExonId:
                inSeq = seq
                (
                    guideSeq,
                    pamSeq,
                    pamPlusSeq,
                    guideSeqWPam,
                    guideStrand,
                    guideSeqHtml,
                    guideStart,
                    guideEnd,
                ) = findGuideSeq(inSeq, pam, pamId, exonId=exonId)
                break
    else:
        (
            guideSeq,
            pamSeq,
            pamPlusSeq,
            guideSeqWPam,
            guideStrand,
            guideSeqHtml,
            guideStart,
            guideEnd,
        ) = findGuideSeq(inSeq, pam, pamId)

    targetChrom, targetStart, targetEnd, strand = parsePos(position)
    # gStart = int(pamId.replace("s", "").replace("+","").replace("-",""))
    # gEnd = gStart+GUIDELEN
    strand = pamId[-1]

    scoreCode = params["showMh"]

    print("<h2>DNA break repair outcome predictions</h2>")
    print("<strong>Guide:</strong> %s<p>" % guideSeqHtml)

    # mutScores = crisporEffScores.calcMutSeqs(pamIds, longSeqs, enz, scoreNames=[scoreCode])
    outcome = readOutcomeData(batchId, scoreCode)[pamId]
    # extend guide to get +- 50 bp around it
    longSeq = getExtSeq(
        inSeq, guideStart, guideEnd, strand, 50 - GUIDELEN, 50, extSeq=extSeq
    )
    if scoreCode == "oof":
        oof, mhSeqs = outcome
        mhSeqs.sort(reverse=True)
        print("""The following table lists possible deletions. """)

    elif scoreCode == "lindel":
        fsProb, mhSeqs = outcome
        print(
            """The following table lists possible deletions and insertions, scored by the <i>Lindel</i> repair model.<br> """
        )

    else:
        errAbort("Invalid score code")

    print(
        """Each sequence below represents the context around the guide's target, with deleted nucleotides shown as "-".<br>"""
    )

    if scoreCode == "oof":
        print(
            """Sequences  are ranked by micro-homology score. This score is correlated to the likelihood of finding a particular deletion, as predicted by the method of <a target=_blank href='http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html'>Bae et al. 2014</a>.<p>"""
        )
        scoreName = "mh- Score"

    elif scoreCode == "lindel":
        print(
            """Sequences are ranked by their probability, as predicted by the method of <a target=_blank href='https://academic.oup.com/nar/article/47/15/7989/5511473'>Chen et al. 2019</a>.<p>"""
        )
        scoreName = "Probability"

    else:
        errAbort("Unknown score code")

    print("<table>")
    print("<tr>")
    print("<th>%s</th>" % scoreName)
    print("<th>Sequence</th>")
    print("<th>Effect</th>")
    print("</th>")

    if scoreCode == "oof":
        print("<tr>")
        print("<td></td>")
        print("<td><tt>%s</tt></td>" % longSeq[10:-10].upper())
        print("<td>(Wild-type sequence)</td>")
        print("</tr>")

    for row in mhSeqs:
        score, seq = row[:2]
        print("<tr>")
        if score == 0:
            print("<td></td>")
        else:
            if scoreCode == "oof":
                print("<td>%d</td>" % score)
            else:
                print("<td>%.2f%%</td>" % score)

        print("<td><tt>%s</tt></td>" % seq)
        if scoreCode == "oof":
            delCount = seq.count("-")
            if delCount % 3 == 0:
                resStr = "no frameshift"
            else:
                resStr = "frameshift"
            print("<td><tt>%d bp deleted &rarr; %s</tt></td>" % (delCount, resStr))
        else:
            if score == 0.0:
                print("<td>wild-type</td>")  # for lindel
            else:
                print("<td>%s</td>" % (row[2]))
        print("</tr>")
    print("</table>")


def printSatMutPage(params):
    "special form for sat mutagenesis downloads, for our Nat Prot paper"
    defBarcode = "1"
    barcode = params.get("barcode", defBarcode)
    batchId = params["batchId"]

    batchInfo = readBatchAsDict(batchId)
    pam = batchInfo["pam"]

    printBackLink()

    print("<h3>Oligonucleotides for a Lentiviral Saturating Mutagenesis Screen</h3>")
    print(
        """
    <p>This page allows you to download all files for the ordering and the analysis of your oligonucleotide pool.<br>"""
    )
    # <br> You can use Excel filters or command line programs like AWK to filter the list of oligos below
    # , e.g. to include only guides with a specificy score > 50.</p>

    print(
        (
            """<p><strong>1 - Subpool barcodes:</strong> The minimum order is
usually several thousand guides and it is cheaper to order more oligos. To reduce the cost per oligo, subsets can be
tagged with a "barcode" (unrelated to Illumina sequencing index barcodes) so they
can be selectively amplified from the pool.<br>
    <form id="paramForm" action="%(action)s" method="GET">
    Subpool barcode:
    """
            % {"action": basename(__file__)}
        )
    )

    printDropDown("barcode", satMutBarcodes, barcode)
    print("</p>")

    print(
        "<p><strong>2 - PCR primers:</strong> Output file C includes two primers per target, for cleavage analysis with high-throughput sequencing. Primers are prefixed with Illumina Adapters. The forward primer prefix is TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG and the reverse prefix is GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG. These prefixes are already added to the primers in the Excel/Tab-sep tables.<br>"
    )
    ampLen = "150"
    tm = "60"
    primerLen = "20"
    printAmpLenAndTm(ampLen, primerLen, tm)
    print("</p>")

    print("<strong>3 - Filters</strong><br>Minimum specificity score: ")
    print(
        """<input id="minSpec" type="text" size="5" name="minSpec" value="30" /><br>"""
    )
    print(
        """<small>A minimal value of 30 will remove only the guides in repeated regions. For screens, many researchers do not care a lot about off-targets. Increase this threshold if you want to more aggressively remove guides with many predicted off-targets.</small>"""
    )
    print("<br>")

    effScoreName = "fusi"
    if pamIsSaCas9(pam):
        effScoreName = "najm"
    scoreLabel = scoreDescs.get(
        effScoreName, ["- Error: NO KNOWN SCORE FOR THIS ENZYME - "]
    )[0]

    print(("Minimum %s efficiency score: " % scoreLabel))
    print(
        """<input id="minFusi" type="text" size="5" name="minFusi" value="10" /><br>"""
    )
    print(
        """<small>A minimal value of 10 will only remove the least efficient guides. Increase this if want to enrich more for predicted high efficiency guides.</small>"""
    )
    print("<br>")

    outTypes = [
        ("satMut", "A - saturating mutagenesis screen oligonucleotides"),
        (
            "targetSeqs",
            "B1 - Selection Analysis: guide sequences (CrispressoCount, .txt format)",
        ),
        (
            "targetLocs",
            "B2 - Selection Analysis: guide locations (CrisprSurf, .csv format)",
        ),
        (
            "ontargetPrimers",
            "C - On-target sequencing: PCR primers, one pair for each guide target",
        ),
        (
            "ontargetAmplicons",
            "D - Cleavage Validation/Analysis: PCR Amplicons and guide sequence (CrispressoPooled)",
        ),
    ]

    formats = [
        ("html", "Do not download, display as webpage to get an idea"),
        ("xls", "Default: Excel for A and text for B/C/D"),
        ("tsv", "For programmers: always text"),
    ]

    print(
        "<p><strong>4 - File format:</strong> Excel tables include a header with information how the oligonucleotides were constructed.<br>Text files have no header and are easier to process with other software.<br>Cleavage analysis files for Crispresso need to be in text format and are therefore not available as Excel files.<br>"
    )

    print("File format:")
    printDropDown("format", formats, "html")
    print("<br>")

    printHiddenFields(params, {"satMut": None, "download": "useGet"})

    # print("The output table consists of the guide target sequences with their scores, the oligonucleotides to order and the possible sequencing primers for targeted PCR amplification of the target.")

    print("<p><strong>Output files:</strong><ul>")
    print(
        "<li>A - Saturating Mutagenesis Oligonucleotides: the oligonucleotides to order from your Custom Oligonucleotide Array Supplier"
    )
    print(
        "<li>B - Selection Analysis: for every oligo from A, its sequence and genome location. For CrisprSurf quantification</li>"
    )
    print(
        "<li>C - On-target sequencing primers: one forward and one reverse primer for every target from B in the input sequence</li>"
    )
    print(
        "<li>D - Cleavage Analysis: for each pair of primers from C, a table with the PCR amplicon and the guide. For CrispressoPooled, to analyze DNA cleavage induced by this guide</ul>"
    )

    print("You can click the four buttons below and save all four files.")

    print("</p>")
    print("""<p>""")
    print(
        """<input id="submitForm" style="" type="submit" name="get-satMut" value="Get A - oligo list">"""
    )
    print(
        """<input id="submitForm" style="" type="submit" name="get-targetLocs" value="Get B - oligo locations">"""
    )
    print(
        """<input id="submitForm" style="" type="submit" name="get-ontargetPrimers" value="Get C - target primers">"""
    )
    print(
        """<input id="submitForm" style="" type="submit" name="get-ontargetAmplicons" value="Get D - target amplicons">"""
    )

    print("</form>")


def printLibForm(params, returnLink=True):
    """ """
    sampleGenes = "PITX2\nMTOR\nTP53\nABO\n3661\nNM_134261"
    print(("""<form id="libForm" action="%s" method="POST">""" % basename(__file__)))

    # print("Organism: ")
    # printDropDown("org", [("human", "Human"),("mouse", "Mouse")], "human")
    # print("<br>")

    if returnLink:
        url = "crispor.py"
        print(("<p><a href='%s'>&larr; return to the CRISPOR main page</a></p>" % url))
    print('<div style="margin: 32px;">')
    print(
        "<h2>Batch Gene Targeting Assistant: Paste a list of genes to download a list of guides</h2>"
    )

    print(
        "Note: if you are planning a saturating mutagenesis screen, e.g. of a non-coding sequence, this is not the right tool. Submit your sequence on the normal <a href='crispor.py'>CRISPOR</a> page, then use the link 'Saturating mutagenesis' at the top of the guide table to get oligonucleotides of all guides in the input sequence."
    )

    print(
        "Both tools are best used in conjunction with our wet-lab <a href='https://www.nature.com/articles/nprot.2018.005'>protocol</a> (<a href='http://biorxiv.org/content/early/2017/04/07/125245'>preprint</a>)<p>"
    )

    print("<strong>Lentiviral screen library:</strong>")
    printDropDown("libName", libLabels, "geckov2")
    print("<br>")

    print("<strong>Number of guides per gene:</strong> ")
    printDropDown("guideCount", [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6)], "3")
    print("<br>")

    print("<strong>Number of non-targeting control guides (max: 1000): </strong>")
    print(
        """<input id="ctrlCount" type="text" size="5" name="ctrlCount" value="10" />"""
    )
    print("<br>")

    print("<strong>Subpool barcode: </strong>")
    printDropDown("barcode", satMutBarcodes, "1")
    print("<br>")

    print("<strong>Custom oligo prefix and suffix : </strong>")
    print(
        """<input id="custPrefix" type="text" size="20" name="custPrefix" value="" />"""
    )

    # print("<strong>Custom oligo suffix (empty for defaults): </strong>")
    print(
        """<input id="custSuffix" type="text" size="20" name="custSuffix" value="" />"""
    )
    print("both empty = defaults for cloning into pLentiGuide-Puro<p>")

    print(
        (
            """
    Enter a list of gene symbols, Entrez Gene IDs or Refseq IDs, one per line (case-insensitive):<br>
    <small>Type 'all' below to get all guides in the library and no gene filtering.</small>
    <small>You may need to <a href='https://discover.nci.nih.gov/matchminer/MatchMinerLookup.jsp'>convert old symbols</a>.</small>
    <small>To convert from UniProt IDs, try <a href="http://www.uniprot.org/uploadlists/">the UniProt converter.</a></small>
    <textarea tabindex="1" style="width:100%%" name="geneIds" rows="25" \
    placeholder="Paste a list of gene symbols here">%s</textarea><p>"""
            % sampleGenes
        )
    )

    print(
        """<div style="text-align: center;">
                <input style="height:40px; width:100px;" id="submitGenes" type="submit" name="SUBMIT" value="Submit">
          </div>
          """
    )

    print('<input type="hidden" name="libDesign" value="1">')
    print("""</form>""")
    print("</div>")


def isInPar(db, chrom, start, end):
    """return None if not in PAR or "1" or "2" if genome is hg19 or hg38 and chrom:start-end is in a PAR1/2 region"""
    if db not in ("hg19", "hg38"):
        return None
    if not chrom in ("chrX", "chrY"):
        return None

    # all coordinates are from https://en.wikipedia.org/wiki/Pseudoautosomal_region
    # and look like they're 1-based
    if db == "hg38":
        if chrom == "chrX":
            if start >= 10001 and end < 2781479:
                return "1"
            if start >= 155701383 and end < 156030895:
                return "2"
        elif chrom == "chrY":
            if start >= 10001 and end < 2781479:
                return "1"
            if start >= 56887903 and end < 57217415:
                return "2"
    elif db == "hg19":
        if chrom == "chrX":
            if start >= 60001 and end < 2699520:
                return "1"
            if start >= 154931044 and end < 155260560:
                return "2"
        elif chrom == "chrY":
            if start >= 10001 and end < 2649520:
                return "1"
            if start >= 59034050 and end < 59363566:
                return "2"

    return None


def getControls(org):
    "return controls as a list"
    dbFname = "screenData/%s_controls.sqlite" % (org)
    conn = sqlite3.Connection(dbFname, timeout=10)
    cur = conn.cursor()
    sql = "SELECT guideSeq from guides"
    try:
        cur.execute(sql)
    except SQLITEERROR:
        errAbort("Cannot open the file %s" % dbFname)

    rows = cur.fetchall()
    guideList = []
    for guideSeq in rows:
        guideList.append(guideSeq[0])
    return list(set(guideList))


def getLibGuides(org, libName, geneIdStr):
    "return a dict with geneId -> list of guide sequences"
    dbFname = "screenData/%s.sqlite" % (libName)
    conn = sqlite3.Connection(dbFname, timeout=10)
    cur = conn.cursor()

    guideData = defaultdict(list)  # geneId -> guideRows
    orderedGenes = []  # ordered list of guideIds in guideData
    inSearchIds = geneIdStr.split("\n")
    notFoundGenes = []
    for inId in inSearchIds:
        inId = inId.strip().split()[0]
        if inId.lower() == "all":
            sql = "SELECT geneSym, entrezId, refseqId, guideSeq, pam from guides"
            try:
                cur.execute(sql)
            except SQLITEERROR:
                errAbort("Cannot open the file %s" % dbFname)
        else:
            if inId.isdigit():
                searchField = "entrezId"
            elif inId.startswith("NM_"):
                searchField = "refseqId"
            else:
                searchField = "geneSym"

            sql = (
                "SELECT geneSym, entrezId, refseqId, guideSeq, pam from guides WHERE %s = ? COLLATE NOCASE"
                % searchField
            )
            try:
                cur.execute(sql, (inId,))
            except SQLITEERROR:
                errAbort("Cannot open the file %s" % dbFname)

        rows = cur.fetchall()

        if len(rows) == 0:
            notFoundGenes.append(inId)
            continue

        # reformat to dict gene -> seq
        doneGenes = set()
        for geneSym, entrezId, refseqId, guideSeq, pam in rows:
            guideData[geneSym].append((entrezId, refseqId, guideSeq, pam))
            if geneSym not in doneGenes:
                orderedGenes.append(geneSym)
                doneGenes.add(geneSym)

    return guideData, orderedGenes, notFoundGenes


def createGuideTable(
    lentiTemp,
    geneIdStr,
    guideCount,
    org,
    libName,
    barcodeId,
    controlCount,
    custPrefix,
    custSuffix,
):
    "write a table with lenti lib guides and oligos to temp/lenti/ and return its filename"

    keyStr = (
        geneIdStr
        + str(guideCount)
        + org
        + libName
        + str(barcodeId)
        + str(controlCount)
        + custPrefix
        + custSuffix
    )
    hasher = hashlib.sha1(keyStr.encode("latin1"))
    digest = hasher.digest()[0:20]
    lentiJobId = base64.urlsafe_b64encode(digest)
    lentiJobId = lentiJobId.decode("latin1").translate(transTab)

    outFname = join(lentiTemp, lentiJobId + ".tsv")

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
            row = [
                "%s_%d" % (geneSym, guideId + 1),
                str(entrezId),
                refseqId,
                guideSeq,
                "%s%s%s" % (oligoPrefix, guideSeq, oligoSuffix),
                pam,
            ]
            ofh.write("\t".join(row))
            ofh.write("\n")

            if guideId + 1 >= guideCount:
                break

    ctrls = getControls(org)
    for ctrlId, guideSeq in enumerate(ctrls[:controlCount]):
        row = [
            "control_%d" % (ctrlId),
            "",
            "",
            guideSeq,
            "%s%s%s" % (oligoPrefix, guideSeq, oligoSuffix),
            "",
        ]
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
    # org = cgiGetStr(params, "org", "human")
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

    tabFname = createGuideTable(
        lentiTemp,
        geneIdStr,
        guideCount,
        org,
        libName,
        barcodeId,
        controlCount,
        custPrefix,
        custSuffix,
    )

    for line in open(tabFname):
        if "Not found genes" in line:
            notFoundGenes = line.strip("\n").split("\t")[1].split(",")

    url = "crispor.py?libDesign=1"
    print(
        ("<p><a href='%s'>&larr; return to the CRISPOR Batch input page</a></p>" % url)
    )

    # print("Organism: %s<br>" % org)
    libLabel = dict(libLabels)[libName].split("(")[0]  # strip the 'recommended' note
    print(("<strong>Library:</strong> %s<br>" % libLabel))
    print(("<strong>Number of guides per gene:</strong> %d<br>" % guideCount))
    print(
        (
            "<strong>Number of non-targeting controls:</strong> %d (all controls are from the GeCKOV2 library)<br> "
            % controlCount
        )
    )

    barcodeDict = dict(satMutBarcodes)
    satMutOpt, optFields = buildPoolOptions(barcodeId, custPrefix, custSuffix)
    oligoPrefix, oligoSuffix = satMutOpt[:2]

    for label, val in optFields.items():
        print(("<strong>%s:</strong> %s<br>\n" % (label, val)))
    print("<p>")

    print(
        "For details on these sequences, see our <a href='https://www.nature.com/articles/nprot.2018.005'>protocol</a> (<a href='https://www.biorxiv.org/content/early/2017/09/11/125245'>preprint</a>).<p>"
    )

    if len(notFoundGenes) != 0 and notFoundGenes != [""]:
        print(
            (
                "Input gene identifiers that were not found: %s<p>"
                % ",".join(notFoundGenes)
            )
        )

    if guideCount > 4 and "gecko" not in libName:
        print(
            (
                "Note: you asked for %d guides per gene but this library includes only four guides per gene, so the maximum number of guides per genes below is four.<p>"
                % guideCount
            )
        )

    print(
        (
            "<a href='%s'>Download table</a><p>"
            % relpath(tabFname, dirname(abspath(__file__)))
        )
    )

    print('<table class="libTable">')
    print("<tr><th style='width:10em'>ID of guide</th>")
    print("<th style='width:10em'>Target Entrez ID</th>")
    print("<th style='width:10em'>Target Refseq ID</th>")
    print("<th style='width:14em'>Guide RNA<br>(click to show in CRISPOR)</th>")
    print(
        (
            "<th style='width:40em'>Full oligonucleotide including guide RNA<br><small>%s</small></th>"
            % optFields["Oligonucl. structure"]
        )
    )
    print("</tr>")

    genomeDbs = {"human": "hg19", "mouse": "mm10"}
    genomeDb = genomeDbs.get(org)

    for row in lineFileNext(open(tabFname)):
        print("<tr>")
        print(("<td>%s</td>" % (row.guideId)))
        print(("<td>%s</td>" % (row.entrezId)))
        print(("<td>%s</td>" % (row.refseqId)))
        if row.pam != "":
            print(
                (
                    '<td><tt><a target=_blank href="crispor.py?org=%s&seq=%s&pam=NGG">%s</a></tt></td>'
                    % (genomeDb, row.guideSeq + row.pam, row.guideSeq)
                )
            )
        else:
            print(("<td><tt>%s</tt></td>" % row.guideSeq))
        print(("<td><tt>%s</tt></td>" % (row.oligoSeq)))
        print("</tr>")
    print("</table>")


def printKoForm(params):
    """form for knock-out mode"""

    genomes = readGenomes()
    annGenomes = readAnnGenomes()
    scriptName = basename(__file__)

    haveHuman = False
    for g in genomes:
        if g[0] == "hg19":
            haveHuman = True

    cookies = http.cookies.SimpleCookie(os.environ.get("HTTP_COOKIE"))
    if "lastKOorg" in cookies and "lastKOpam" in cookies and "lastKOmethod" in cookies:
        lastorg = cookies["lastKOorg"].value
        lastpam = cookies["lastKOpam"].value
    else:
        if not haveHuman:
            global DEFAULTORG
            DEFAULTORG = ALTORG
        lastorg = DEFAULTORG
        lastpam = DEFAULTPAM

    koBeDesc = []
    for bePam in pamVariantModels.keys():
        if bePam == "NGG":
            desc = "20bp-NGG - SpCas9 based editors"
        else:
            desc = "20bp-%s - SpCas9 variant %s based editors" % (bePam, bePam)
        koBeDesc.append((bePam, desc))

    print(
        """
    <script>

    function populatePamDropdown(desc) {

        // populates the PAM dropdown with input PAM options list

        $("#pamDropDown").empty();

        // from https://stackoverflow.com/questions/740195

        for (let i = 0; i < desc.length; i++) {

            // console.log(desc[i][1]);

            let pamItem = desc[i];
            let pamOpt = new Option(pamItem[1], pamItem[0]);

            $(pamOpt).html(pamItem[1]);
            $("#pamDropDown").append(pamOpt);

        };
    }

    function toggleMethod() {
    const koMethods = document.getElementsByName('koMethod');
    const flankLen = document.getElementById('flankLen');
    const promoterLen = document.getElementById('promoterLen');
    const exonSelect = document.getElementById('exonSelect');
    const beDesc = %s;
    const pamDesc = %s;
    const lastpam = %s;

    let selectedValue;
    for (const method of koMethods){
        if (method.checked) {
            selectedValue = method.value;
            break; }
        }
        if (selectedValue === 'excision') {
            flankLen.style.display = 'flex';
        } else {
            flankLen.style.display = 'none'
        }
        if (selectedValue === 'promoter') {
            promoterLen.style.display = 'flex'
        } else {
            promoterLen.style.display = 'none'
        }
        if (selectedValue === 'splicing') {
            $('#exonSelectContainer').show();
        } else {
            $('#exonSelectContainer').hide();
        }
        if (selectedValue === 'stop') {

            // this method uses base editors so the PAM dropdown needs to be updated
            populatePamDropdown(beDesc)

        } else {

            // rebuild the original PAM dropdown
            populatePamDropdown(pamDesc)
        }
    }



// get the exon number for the selected geneID
$(document).ready(function() {
    $('.js-select-gene').on('select2:select', function (e) {
        var data = e.params.data;
        var exonSelect = $('#exonSelect');
        exonSelect.empty();

        // get exon frames
        var exFrames = data.exFrames.split(',').map(s => s.trim());

        if (data.exonCount != undefined) {

            const allExonsOpt = new Option('target all exons ', 'allExons', false, false);
            exonSelect.append(allExonsOpt);

            for (var i = 0; i < data.exonCount; i++) {
                j = i+1;

                // if the current frame and the next frame are the same, removing the exon won't destroy the reading frame
                // test version, not sure if labeled exons are really out-of-frame

                var frame = exFrames[i];
                var nextFrame = exFrames[j];

                // don't take UTRs into account
                if (frame === nextFrame || frame === -1) {
                    oofText = ""
                } else {oofText = " (out of frame exon)"};

                var exonText = 'target exon ' + j + oofText;
                var option = new Option(exonText, i, false, false);
                exonSelect.append(option);
            }
            exonSelect.trigger('change');
        }
    });
});
</script>

    """
        % (json.dumps(koBeDesc), json.dumps(pamDesc), json.dumps(lastpam))
    )

    print(
        """
    <form id="KoForm", method="get">

        <input type=hidden name="assist" value="1">
        <input type=hidden name="expType" value="ko">
        <div style="display:grid; clear:both; width: 100%%; min-width: 1650px; grid-template-columns: 42% 58%; grid-template-rows: auto auto; place-self:center; justify-self:center; space:20px; padding:12px;">
        <div class="windowstep subpanel" style="width:90%; grid-column:1; grid-row:1;">

            <details id="ko1" open>
            <summary><small>Show / Hide step 1</small></summary>
            <div class="title" style="cursor:pointer" onclick="$('#helpstep3').toggle('fast')">
                Step 1
            </div>
            <div class="substep" style="margin-bottom:20px;">
                Select a genome
            <br>
            """
    )
    printOrgDropDown(lastorg, genomes)

    print(
        """
                <div id="trackHubNote" style="margin-bottom:12px; margin-top:12px;">
                    <small>Note: pre-calculated exonic guides for this species are on the <a id='hgTracksLink' target=_blank href="">UCSC Genome Browser</a>.</small>
                </div>
            <small style="float:left">We have %d genomes, but not yours? Search <a href="https://www.ncbi.nlm.nih.gov/assembly">NCBI assembly</a> and send a GCF_/GCA_ ID to <a href="mailto:%s">CRISPOR support</a>.</small><br>
            </div>
        </div>
        </details>
    """
        % (len(genomes), contactEmail)
    )
    print(
        """
        <div class="windowstep subpanel" style="display:flex; width:90%%; flex-direction:column; grid-column:1; grid-row:2;">

        <details id="ko2" open>
        <summary><small>Show / Hide step 2</small></summary>
        <div>
            <div class="title" style="cursor:pointer;">
                Step 2
            </div>
            <div style="margin-bottom:35px; margin-top:12px;">
                Select a Protospacer Adjacent Motif (PAM)
                <img src="%simage/info-small.png" title="The most common system uses the NGG PAM recognized by Cas9 from S. <i>pyogenes</i>. The VRER and VQR mutants were described by <a href='http://www.nature.com/nature/journal/vaop/ncurrent/abs/nature14592.html' target='_blank'>Kleinstiver et al</a>, Cas9-HF1 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4851738/'>Kleinstiver 2016</a>, eSpCas1.1 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4714946/'>Slaymaker 2016</a>, Cpf1 by <a href='http://www.cell.com/abstract/S0092-8674(15)01200-3'>Zetsche 2015</a>, SaCas9 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/25830891/'>Ran 2015</a> and KKH-SaCas9 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/26524662/'>Kleinstiver 2015</a>, modified As-Cpf1s by <a href='http://biorxiv.org/content/early/2016/12/04/091611'>Gao et al. 2017</a>." class="tooltipsterInteract">
                <br>
            </div>
            <div style="margin-bottom:15px;">

    """
        % HTMLPREFIX
    )

    printPamDropDown(lastpam)

    print(
        """
            </div>
            <div style="margin-bottom:15px;">
                <br>See <a target=_blank href="manual/manual.html#enzymes">notes on enzymes</a> in the manual.<br>
            </div>
        </div>
       </div>
       </details>
        """
    )

    print(
        """
        <div class="windowstep subpanel" style="width:100%; grid-column:2; grid-row:1/-1;">
            <div class="title" style="cursor:pointer; margin-bottom:20px;" onclick="$('#helpstep3').toggle('fast')">
                Step 3
            </div>
            <div style="margin-bottom:12px;"> Select a gene and choose to target a specific transcript or exons common to all transcripts</div>
          """
    )

    printGeneSelection("koGeneId", commonExons=True)

    print(
        """
                <div style="margin-top:12px; margin-bottom: 24px;">
                    <small>Currently, %d out of %d genomes are annotated with genes. If yours isn't included, use CRISPOR classic.</small><br>
                </div>

            <small>Planning a lentiviral gene knockout screen? Use <a href="crispor.py?libDesign=1">CRISPOR Batch</a><br></small>

            <p style="margin-top:24px">Choose one of the following approaches to inactivate your gene</p>
        """
        % (len(annGenomes), len(genomes))
    )

    # Knock out methods [(id, text, additional inputs)]
    methods = [
        (
            "frameshift",
            " Frameshift mutation in the first third of the coding sequence",
        ),
        (
            "stop",
            " Introduce a premature STOP codon in the first half of the coding sequence with base editing",
        ),
        ("excision", " Excision of the gene locus"),
        ("promoter", " Removal of the promoter"),
        ("splicing", " Disruption of splicing by targeting a splice site"),
    ]

    for methodId, methodDesc in methods:
        if methodId == "frameshift":
            methodChecked = "checked"
        else:
            methodChecked = ""
        print(
            """
        <input type="radio" %(methodChecked)s name="koMethod" id="%(methodId)s" value="%(methodId)s" onchange="toggleMethod()"/>%(methodDesc)s<br>
        """
            % locals()
        )

        if methodId == "excision":
            print(
                """
                <small id="flankLen" style="text-align:left; display:none; align-items:center; margin-left:25px; margin-top:12px; margin-bottom:12px;">
                <input type="range" name="flankLen" value="500" min="100" max="1000" oninput="this.nextElementSibling.value = this.value"/>
                &nbsp&nbsp target&nbsp<output style="">500</output> bp upstream of the transcription start site (TSS) and downstream of the transcription end site (TES)
                </small>
            """
            )
        elif methodId == "promoter":
            print(
                """
                <small id="promoterLen" style="text-align:left; display:none; align-items:center; margin-left:25px; margin-top:12px; margin-bottom:12px;">
                <input type="range" name="promoterLen" value="500" min="50" max="2000" oninput="this.nextElementSibling.value = this.value"/>
                &nbsp&nbsp Remove a&nbsp<output style="">500</output> bp region upstream of the transcription start site (TSS)
                </small>
            """
            )

    print(
        """
            <div id="exonSelectContainer" style="display: none; margin-top: 12px; margin-bottom: 12px; margin-left: 25px;">
            <select class="js-select-hidden" name="exonSelect" id="exonSelect" style="width:20%;">
            </select>
            </div>
            """
    )

    print(
        """
    <div style="align-items:center; text-align:center;margin-top:20px;">
        <br><input id="submitKoGeneID" type="submit" name="submit" value="SUBMIT" style="height:35px; width:100px;">
    </div>
    </div>
 """
    )
    print(
        """
    </div>
    </div>
    </div>
    </form>
    """
    )

    # print this script at the end so that it is executed last
    print(
        """
    <script>
        /* hide the track hub note if genome is not hg19 */
        ucscTrackDbs=['hg19', 'hg38', 'rn5', 'mm10', 'mm9', 'ci2', 'danRer7', 'sacCer3', 'dm6'];
        function showHideHubNote() {
            const valSel = $("#genomeDropDown").val();
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

<script>
// save the states of detail elements on page reload
(function() {
    var $details = $('details[id]');
    $details.each(function() {
        var savedState = localStorage.getItem('details-' + this.id);
        if (savedState !== null) {
            this.open = savedState === 'true';
        }
    });

    $details.on('toggle', function() {
        localStorage.setItem('details-' + this.id, this.open);
    });
})();

// set PAM dropdown on page reload based on the selected method
$(document).ready(function() {
    const koMethods = document.getElementsByName('koMethod');
    const beDesc = %s;
    const pamDesc = %s;

    for (method of koMethods) {
        if (!(method.checked)) {
            continue;
        }
        if (method.id === "stop") {
            populatePamDropdown(beDesc);
        } else {
            populatePamDropdown(pamDesc);
        }}
})

</script>

          """ % (json.dumps(koBeDesc), json.dumps(pamDesc))
    )


def printTagsAndLinkers(tag=True, qTAG=True, tagNames=None):
    """
    prints the dropdown menus for tags and linkers
    Optionally, the tags in tagNames (list) can be preselected
    """

    print(
        """
    <script>
        $(document).ready(function() {
        $('.js-select-tag').select2({
            placeholder: 'select a tag',
            width: '200px'
            });
        });
        $(document).ready(function() {
        $('.js-select-linker').select2({
            placeholder: 'select a linker',
            width: '200px'
            });
        });
        $(document).ready(function() {
        $('.js-select-marker').select2({
            placeholder: 'a selectable marker',
            width: '200px'
            });
        });
        $(document).ready(function() {
        $('.js-select-expression').select2({
            placeholder: 'an expression method',
            width: '250px'
            });
        });
       $(document).ready(function() {
       $('.js-select-qtag').select2({
           placeholder: 'Choose a tag',
           width: '200px'
           });
       });

    </script>
    """
    )

    # tags and linkers options
    tags = {
        "Fluorescent proteins": [
            "eGFP",
            "mStrayGold",
            "mNeon",
            "moxGFP",
            "mScarlet",
            "mCherry",
            "sTagRFP",
            "miRFP670nano3",
        ],
        "Proximity Biotinylation": ["miniTurbo", "ultraID"],
        "Targeted degradation": ["dTAG"],
        "Epitopes": [
            "3XFLAG",
            "3XHA",
            "V5",
            "SBP",
            "SBP3Flag",
            "3FlagSBP",
            "Streptavidin",
        ],
    }

    linkers = {
        "flexible linker": ["GGGGS", "GSGGG", "(GGGGS)x2", "XTEN"],
        "rigid linkers": [
            "",
        ],
        "multimer linkers": [""],
    }

    # qTAG cassette elements
    markers = {
        "Mammalian selection": ["Blast", "Puro", "Zeo"],
        "Fluorescent selection": ["moxGFP", "mScarlet"],
    }

    expressionSeqs = ["2A ribosomal skipping peptide", "EF1α promoter"]

    qTags = {
        "Fluorescent proteins": [
            "mStrayGold",
            "mNeon",
            "moxGFP",
            "mScarlet",
            "sTagRFP",
            "miRFP670nano3",
        ],
        "Proximity Biotinylation": ["miniTurbo", "ultraID"],
        "Targeted degradation": ["dTAG"],
        "Epitopes": ["3XFLAG", "3XHA", "V5"],
    }

    print(
        """<div class="windowstep subpanel" id="tagPanel" style="width:100%; height:75px; display:flex; flex-direction:row; margin-bottom: 18px;"> """
    )

    if tag:

        print(
            """<div id="tagLinkerDisplay" style="display:flex; flex-direction:row; align-items:center; gap:10px; padding: 12px;">"""
        )
        print(
            """<div>
             <select name="linkerseq" id="linkerseq" class="js-select-linker" style="width:100%;" autocomplete="off">
             <option selected="selected"></option>
           """
        )
        for linkerType in linkers:
            print("""<optgroup label="%s">""" % linkerType)
            for linker in linkers[linkerType]:
                if tagNames and linker in tagNames:
                    selected = "selected"
                else:
                    selected = ""
                print(
                    """<option %s value="%s">%s</option>""" % (selected, linker, linker)
                )
            print("</optgroup>")

        print("</select>")

        print("</div>")

        print("<div>and</div>")

        print(
            """<div>
              <select name="tagseq" id="tagseq" class="js-select-tag" style="width:100%;" autocomplete="off">
              <option selected="selected"></option>
            """
        )
        for tagType in tags:
            print("""<optgroup label="%s">""" % tagType)
            for tag in tags[tagType]:
                if tagNames and tag in tagNames:
                    selected = "selected"
                else:
                    selected = ""

                print("""<option %s value="%s">%s</option>""" % (selected, tag, tag))
            print("</optgroup>")

        print(
            """</select>
              </div>"""
        )

        print("</div>")

    # qTAG options
    if qTAG:
        if qTAG and not tag:
            display = "display: flex"
        else:
            display = "display: none"

        print(
            """<div id="qTagDisplay" style="%s; flex-direction:row; align-items:center; gap:18px; padding: 12px;">"""
            % display
        )

        # TAG sequence
        print(
            """
        <div>
        <select name="qTag" id="qTag" class="js-select-qtag" style="width:100%; margin-right: 24px;" autocomplete="off">"""
        )
        print("<option></option>")
        for tagType in qTags:
            print("""<optgroup label="%s">""" % tagType)
            for tag in qTags[tagType]:
                if tagNames and tag in tagNames:
                    selected = "selected"
                else:
                    selected = ""
                print("""<option %s value="%s">%s</option>""" % (selected, tag, tag))
            print("</optgroup>")
        print(
            """
        </select>
        </div>
        """
        )

        print("<div>and</div>")

        # Markers
        print(
            """
        <div>
        <select name="markerseq" id="markerseq" class="js-select-marker" style="width:100%;" autocomplete="off">"""
        )
        print("<option></option>")
        print("""<option value="none">None</option>""")
        for markerType in markers:
            print("""<optgroup label="%s">""" % markerType)
            for marker in markers[markerType]:
                if tagNames and marker in tagNames:
                    selected = "selected"
                else:
                    selected = ""

                print(
                    """<option %s value="%s">%s</option>""" % (selected, marker, marker)
                )
            print("</optgroup>")
        print(
            """
        </select>
        </div>
        """
        )

        print("<div>and</div>")

        # Expression method
        print(
            """
        <div>
        <select name="expressionSeq" id="expressionSeq" class="js-select-expression" style="width:100%;" autocomplete="off">"""
        )
        print("<option></option>")
        print("""<option value="none">in-frame fusion to target gene</option>""")
        for expressionSeq in expressionSeqs:
            if tagNames and expressionSeq in tagNames:
                selected = "selected"
            else:
                selected = ""

            if expressionSeq == "EF1":
                expressionSeqHtml = expressionSeq + "&alpha;"
            else:
                expressionSeqHtml = expressionSeq
            print(
                """<option %s value="%s">%s</option>"""
                % (selected, expressionSeq, expressionSeqHtml)
            )
        print(
            """
        </select>
        </div>
        """
        )

        print("</div>")
    print("</div>")


def printKiForm(params):
    # form for knock-in mode in the assistant
    # inspired by the options from protoSpaceJam (https://protospacejam.sf.czbiohub.org/)

    genomes = readGenomes()
    annGenomes = readAnnGenomes()
    scriptName = basename(__file__)

    haveHuman = False
    for g in genomes:
        if g[0] == "hg19":
            haveHuman = True

    cookies = http.cookies.SimpleCookie(os.environ.get("HTTP_COOKIE"))
    if "lastKIorg" in cookies:
        lastorg = cookies["lastKIorg"].value
        lastmultipam = cookies["lastKIpam"].value
    else:
        if haveHuman is False:
            global DEFAULTORG
            DEFAULTORG = ALTORG
        lastorg = DEFAULTORG
        lastmultipam = "20bp-NGG"

    print(
        """
    <script>

        function toggleTargetRegion() {
            const targetRegions = document.getElementsByName('targetRegions')
            const geneTarget = document.getElementById('geneTarget')

            const seqTarget = document.getElementById('seqTarget')
            const seqTargetText = document.getElementById('seqTargetText')

            const endSeqDisplay = document.getElementById('endSeqDisplay')
            const seqEditText = document.getElementById('seqEditText')

            const tagInsertDisplay = document.getElementById('tagInsertDisplay')
            const geneSelection = document.getElementById('geneSelection')

            let selectedValue;
            for (const target of targetRegions) {
                if (target.checked) {
                    selectedValue = target.value;
                    break; }
                }
            if (selectedValue === 'seq') {
                seqTarget.style.display = 'block';
                endSeqDisplay.style.display = 'block';
                seqTargetText.style.display = 'block';
                seqEditText.style.display = 'block';
            } else {
                seqTargetText.style.display = 'none';
                seqEditText.style.display = 'none';
            }

            if (selectedValue === 'gene') {
                geneTarget.style.display = 'block';
                tagInsertDisplay.style.display = 'block';
                seqTarget.style.display = 'none';
                endSeqDisplay.style.display = 'none';
            } else {
                seqTarget.style.display = 'block';
                endSeqDisplay.style.display = 'block';
                geneTarget.style.display = 'none';
                tagInsertDisplay.style.display = 'none';
                $('#geneSelection').val(null).trigger('change');

            }

        }

        function toggleInsertpos() {
            const insertpos = document.getElementsByName('insertpos');
            const ntext = document.getElementById('Ntext')
            const ctext = document.getElementById('Ctext')
            let selectedValue;
            for (const pos of insertpos) {
                if (pos.checked) {
                    selectedValue = pos.value;
                    break; }
            }
            if (selectedValue === 'Nter') {
                ntext.style.display = 'block';
            } else {
                ntext.style.display = 'none';
            }
            if (selectedValue === 'Cter') {
                ctext.style.display = 'block';
            } else {
                ctext.style.display = 'none';
            }
        }

        function toggleInsertseq() {
          const insertype = document.getElementsByName('insertype')
          const taglist = document.getElementById('taglist')
          const insertSeq = document.getElementById('insertSeq')
          const qTagDisplay = document.getElementById('qTagDisplay')
          const tagLinkerDisplay = document.getElementById('tagLinkerDisplay')
          const tagPanel = document.getElementById('tagPanel')

          const tagseq = document.getElementById('tagseq')
          const linkerseq = document.getElementById('linkerseq')
          const markerseq = document.getElementById('markerseq')
          const expressionSeq = document.getElementById('expressionSeq')
          const qTag = document.getElementById('qTag')

          let selectedValue;
          for (const type of insertype){
            if (type.checked) {
                selectedValue = type.value;
                break; }
            }
            if (selectedValue === 'tagLinker') {
                taglist.style.display = 'block';
                tagLinkerDisplay.style.display = 'flex';
                tagPanel.style.display = 'flex';
           } else {
                tagLinkerDisplay.style.display = 'none';
                $('#tagseq').val(null).trigger('change');
                $('#linkerseq').val(null).trigger('change');
               };
            if (selectedValue === 'qTag') {
                taglist.style.display = 'block';
                qTagDisplay.style.display = 'flex';
                tagPanel.style.display = 'flex';
           } else {
                qTagDisplay.style.display = 'none';
                $('#markerseq').val(null).trigger('change');
                $('#expressionSeq').val(null).trigger('change');
                $('#qTag').val(null).trigger('change');
               };
            if (selectedValue === 'custom') {
                insertSeq.style.display = 'block';
                tagPanel.style.display = 'none';
            } else {
                insertSeq.style.display = 'none';
            }
        }

// prevents the form from submitting if the enter key is pressed
function handleEnter(event) {
    if (event.keyCode === 13) {
        event.preventDefault();
        event.target.blur();
    }
}

function changeSeqCase(value) {
    const textarea = document.getElementById('endSeq');
    const start = textarea.selectionStart;
    const end = textarea.selectionEnd;
    const selectedText = textarea.value.substring(start, end);

    let modText;
    if (value === 'uppercase') {
        modText = selectedText.toUpperCase();
        }
    if (value === 'lowercase') {
        modText = selectedText.toLowerCase();
        }
    textarea.setRangeText(modText, start, end, 'select');
}

</script>
    """
    )

    print(
        """
    <script>
/* set the dropbox to hg19 and paste the example sequence into the input box. */
function resetToExample() {
    $("textarea[name='startSeq']").val("%s");
    $("#genomeDropDown").val("%s").trigger("chosen:updated");
    }

function insertExample() {
    $("textarea[name='endSeq']").val("%s");
    }
function delExample() {
    $("textarea[name='endSeq']").val("%s");
    }
function substExample() {
    $("textarea[name='endSeq']").val("%s");
    }
function replExample() {
    $("textarea[name='endSeq']").val("%s");
    }

/* clear the startSeq and endSeq input boxes
should pass element name as a parameter, but it
doesn't work ? */

function clearStartSeq() {
    $('textarea[name="startSeq"]').val("");
    }
function clearEndSeq() {
    $('textarea[name="endSeq"]').val("");
    }

    </script>
    """
        % (
            DEFAULTKISEQ,
            DEFAULTORG,
            DEFAULTINSERT,
            DEFAULTDEL,
            DEFAULTSUBST,
            DEFAULTREPL,
        )
    )

    print(
        """
    <script>
        $(document).ready(function() {
        $('.js-example-basic-single').select2();
        });
    </script>
    """
    )

    print(
        """
    <form id="KiForm" method="GET">
        <input type=hidden name="assist" value="1">
        <input type=hidden name="expType" value="ki">

       <div style="display:flex; clear:both; padding:12px; width: 100%; min-width: 1550px;">
       <div style="width: 50% ;flex:0 0 41%; display:flex; flex-direction:column; gap:14%;">

        <div class="windowstep subpanel" style="width:90%; grid-column:1; grid-row:1; height:30%;">

            <details id="ki1" open>
            <summary><small>Show / Hide step 1</small></summary>
            <div class="title" style="cursor:pointer;" onclick="$('#helpstep3').toggle('fast')">
                Step 1
            </div>
            <div class="substep" style="margin-bottom:20px;">
                Select a genome
            <br>
            """
    )

    printOrgDropDown(lastorg, genomes)

    print(
        """
                <div id="trackHubNote" style="margin-bottom:12px; margin-top:12px;">
                    <small>Note: pre-calculated exonic guides for this species are on the <a id='hgTracksLink' target=_blank href="">UCSC Genome Browser</a>.</small>
                </div>
            <small style="float:left">We have %d genomes, but not yours? Search <a href="https://www.ncbi.nlm.nih.gov/assembly">NCBI assembly</a> and send a GCF_/GCA_ ID to <a href="mailto:%s">CRISPOR support</a>.</small><br>
            </div>
        </div>
        </details>
    """
        % (len(genomes), contactEmail)
    )
    print(
        """
        <div class="windowstep subpanel" style="width:90%%; grid-column:1; grid-row:2; height: 30%%;">

        <details id="ki2" open>
        <summary><small>Show / Hide step 2</small></summary>
            <div class="title" style="cursor:pointer; margin-bottom:12px;" onclick="$('#helpstep3').toggle('fast')">
                Step 2
            </div>
        Select a list of Protospacer Adjacent Motifs (PAMs)
    <img src="%simage/info-small.png" title="The most common system uses the NGG PAM recognized by Cas9 from S. <i>pyogenes</i>. The VRER and VQR mutants were described by <a href='http://www.nature.com/nature/journal/vaop/ncurrent/abs/nature14592.html' target='_blank'>Kleinstiver et al</a>, Cas9-HF1 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4851738/'>Kleinstiver 2016</a>, eSpCas1.1 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4714946/'>Slaymaker 2016</a>, Cpf1 by <a href='http://www.cell.com/abstract/S0092-8674(15)01200-3'>Zetsche 2015</a>, SaCas9 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/25830891/'>Ran 2015</a> and KKH-SaCas9 by <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/26524662/'>Kleinstiver 2015</a>, modified As-Cpf1s by <a href='http://biorxiv.org/content/early/2016/12/04/091611'>Gao et al. 2017</a>." class="tooltipsterInteract">
        <br>
        <div style="margin-top:12px;"> """
        % HTMLPREFIX
    )

    print(
        """ <select class=js-example-basic-single style="width:85%;" name="multipam"> """
    )
    for pamKey in multiPamDict:
        if pamKey == lastmultipam:
            pamSelected = "selected"
        else:
            pamSelected = ""
        pamText = multiPamDict[pamKey][1]
        print(
            """ <option %s value="%s">%s</option> """ % (pamSelected, pamKey, pamText)
        )
    print(""" </select> """)

    print(
        """
        </div>
        <div style="margin-bottom: 6px;">
            <br>See <a target=_blank href="manual/manual.html#enzymes">notes on enzymes</a> in the manual.<br>
        </div>
        </div>
        </div>
        </details>
        """
    )

    print(
        """
        <div style="flex:0 0 59%%; display:flex; flex-direction:column; width: 100%%;">
        <div class="windowstep subpanel" style="width:100%%; grid-column:2; grid-row:1;">
            <div class="title" style="cursor:pointer; margin-bottom:4px;" onclick="$('#helpstep3').toggle('fast')">
                Step 3
            </div>
            <div style="display:flex; align-items:center; justify-contents: space-around; gap: 20%%;">
                <div style="align-self: left;">
                Choose one of the following methods :
                </div>
                <div style="border:dashed 0.5px grey; border-radius:6px; padding:8px;">
                    <input type="radio" checked name="targetRegions" value="seq" onchange="toggleTargetRegion()" autocomplete="off"/>Enter a sequence with the desired modifications<br>
                    <input type="radio" name="targetRegions" value="gene" onchange="toggleTargetRegion()" autocomplete="off"/>Select a transcript to tag a protein in Nter or Cter<br>
                </div>
            </div>
            <div id="seqTarget" style="margin-top:20px; display: flex; flex-direction: column;">
                <div id="seqTargetText">
                    Enter the target sequence manually here<br>
                </div>
                <div style="display:flex; flex-direction: row; margin-top: 6px; gap: 4px; align-items: center;">
                    <small><a href="javascript:clearStartSeq()">Clear Box</a> - </small>
                    <small><a href="javascript:resetToExample()">Set a default example</a></small>
                    <small>
                        <input style="margin-left: 24px;" type="checkbox" name="noPerfectMatch" value=1 />Sequence not in reference genome
                        <img src=" %s image/info-small.png" class="tooltipsterInteract" title="By default, CRISPOR only uses the coordinates of input sequences that are identical to the reference genome, to differentiate off-target and on-target sites. However, the coordinates are needed to design donor DNA. You can select this option to use the coordinates of the best match in the genome as the on-target site instead, allowing to search in sequences that differ from the genome. For example, this option can be used to correct pathogenic mutations in disease models.">
                    </small>
                </div>
                <textarea name="startSeq" style="display: block;" rows="8" cols="108" placeholder="Paste the target sequence here (max. 2300bp). The sequence should be identical to the selected genome." autocorrect="off" autocapitalize="off" spellcheck="false"></textarea>
            </div>
        <div id="geneTarget" style="display: none;">
            <div style="margin-bottom:15px; margin-top:20px;">Select a transcript</div>
          """ % HTMLPREFIX
    )

    printGeneSelection("koGeneId")

    print(
        """
                <div style="margin-top:20px;">
                    <small>Currently, %d out of %d genomes are annotated with genes. If yours insn't included, select "Enter a sequence" above.</small><br>
                </div>
                <div style="margin-top: 37px; margin-bottom: 16px; display:flex; flex-direction:row; align-items:center;">
                Insert :
                    <input type="radio" checked name="insertpos" value="Nter" onchange="toggleInsertpos()" autocomplete="off"/>After the START codon
                    <input type="radio" name="insertpos" value="Cter" onchange="toggleInsertpos()" autocomplete="off"/>Before the STOP codon
            <small id='Ntext' style="display: block;">&nbsp&nbsp(tagging the N-terminal end of the protein)</small>
            <small id='Ctext' style="display: none;">&nbsp&nbsp(tagging the C-terminal end of the protein)</small>
                </div>

            </div>
            </div>"""
        % (len(annGenomes), len(genomes))
    )

    print(
        """
            <div class="windowstep subpanel" style="display:flex; flex-direction:column; width:100%%;">
                <div class="title" style="cursor:pointer;" onclick="$('#helpstep3').toggle('fast')">
                    Step 4
                </div>
                <div id="endSeqDisplay" style="display: block; margin-bottom:12px; margin-top:12px;">
                    <div style="display: flex; flex-direction: row;">
                        <div id="seqEditText" style="margin-right:20px;">
                            Enter the edited sequence (target sequence with edits)<br>
                            <small>Modified bases in UPPERCASE, with the rest in lowercase</small>
                        </div>
                        <div style="display: flex; flex-direction: row; justify-content: space-around; width:50%%;">
                            <button type="button" onclick="changeSeqCase('uppercase')" style="width: 30%%; justify-self: center; background: #ffffff; color: #0480be; box-shadow: 0 2px 10px 2px #9bdcfd; webkit-box-shadow: 0 2px 10px 2px #9bdcfd; moz-box-shadow: 0 2px 10px 2px #9bdcfd;"><small>Change selection to uppercase</small></button>
                            <button type="button" onclick="changeSeqCase('lowercase')" style="width: 30%%; justify-self: center; background: #ffffff; color: #0480be; box-shadow: 0 2px 10px 2px #9bdcfd; webkit-box-shadow: 0 2px 10px 2px #9bdcfd; moz-box-shadow: 0 2px 10px 2px #9bdcfd;"><small>Change selection to lowercase</small></button>
                        </div>
                    </div><br>
                    <div style="display: flex; flex-direction: column;">
                        <div style="display: flex; flex-direction: row; gap: 4px; align-items: center;">
                            <small><a href="javascript:clearEndSeq()">Clear Box</a> -</small>
                            <small>Set an example for :</small>
                            <small class="tooltipsterInteract" title="In this example, the CDS of the mCherry fluorescent protein (%d bp) is inserted after the START codon in the human homeobox D9 (HOXD9) gene."><a href="javascript:insertExample()"> Insertion</a></small>
                            <small>/</small>
                            <small class="tooltipsterInteract" title="In this example, a 3bp deletion is introduced to delete the START codon in the human homeobox D9 (HOXD9) gene."><a href="javascript:delExample()"> Deletion</a></small>
                            <small>/</small>
                            <small class="tooltipsterInteract" title="In this example, a A to G substitution is introduced in the first base of the START codon in the human homeobox D9 (HOXD9) gene."><a href="javascript:substExample()"> Substitution</a></small>
                            <small>/</small>
                            <small class="tooltipsterInteract" title="In this example, a 3bp replacement is introduced to change the START codon in the human homeobox D9 (HOXD9) gene into a STOP codon."><a href="javascript:replExample()"> Replacement</a></small>
                        </div>
                        <textarea name="endSeq" id="endSeq" rows="8" cols="108" autocorrect="off" autocapitalize="off" spellcheck="false" placeholder="Paste the edited sequence here. Edits should be in uppercase (except for deletions), with the rest of the sequence in lowercase. Types of modification supported are insertion, deletion, single substitution and replacements (up to 10 bp, including e.g. with two substitions 10 bp apart)."></textarea>
                    </div>
                </div>

                <div id="tagInsertDisplay" style="display: none; margin-bottom:12px; margin-top:24px;">
                    Enter the sequence to insert<br>
                    <input type="radio" checked style="margin-top:24px;" name="insertype" value="tagLinker" onchange="toggleInsertseq()" autocomplete="off"/>Choose from a list of linkers and tags
<input type="radio" style="margin-top:13px;" name="insertype" value="qTag" onchange="toggleInsertseq()" autocomplete="off"/>qTAG <img src=" %s image/info-small.png" class="tooltipsterInteract" title="The qTAG system combines the tagging sequence with a marker (fluorescent protein or antibiotic resistance gene). The marker is flanked by loxP sites to allow its removal when successfully edited cells have been selected. For more information, see <a href='https://doi.org/10.1038/s44318-024-00337-5' target='blank'>Philip et al. 2025</a>">

                    <input type="radio" name="insertype" value="custom" onchange="toggleInsertseq()" autocomplete="off"/>Paste a custom sequence
                    <textarea spellcheck="false" autocorrect="false" id="insertSeq" name="insertSeq" style="display: none;" rows="6" cols="100" placeholder="Paste the sequence you want to insert here (case insensitive). Please keep the sequence in frame."></textarea>

                <div style="width:95%%; display: block; margin-top: 24px; margin-bottom: 24px;" id="taglist">
          """
        % (len(DEFAULTINSERTSEQ), HTMLPREFIX)
    )

    printTagsAndLinkers()

    print(
        """
                </div>
        </div>
    <input id="submitKoGeneID" type="submit" name="submit" value="SUBMIT" style="align-self:center; height:40px; width:100px;">
          """
    )
    print(
        """
        </div>
        </div>
    </div>
    </form>
    """
    )

    print(
        """
    <script>
        /* hide the track hub note if genome is not hg19 */
        ucscTrackDbs=['hg19', 'hg38', 'rn5', 'mm10', 'mm9', 'ci2', 'danRer7', 'sacCer3', 'dm6'];
        function showHideHubNote() {
            const valSel = $("#genomeDropDown").val();
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

<script>
// save the states of detail elements on page reload
(function() {
    var $details = $('details[id]');
    $details.each(function() {
        var savedState = localStorage.getItem('details-' + this.id);
        if (savedState !== null) {
            this.open = savedState === 'true';
        }
    });

    $details.on('toggle', function() {
        localStorage.setItem('details-' + this.id, this.open);
    });
})();
</script>

    """
    )


def wrongInputRedirect(msg):
    """Show a warning message and exits"""

    printCrisporBodyStart()

    print(
        """<div style="justify-content: center; display: flex; gap: 12px; margin-top: 24px; margin-bottom: 24px;">"""
    )
    print("<p>Sorry, the input you entered can't be used by CRISPOR : %s</p><br>" % msg)
    print("</div>")

    printTeforBodyEnd()
    sys.exit(0)


def printBody(params):
    "main dispatcher function"

    # TODO: first, if batchId is the only parameters,
    # we will check for a flag file to see if the job is running and output the status file if it is.

    org = params.get("org")
    pam = params.get("pam")
    koGeneId = params.get("koGeneId")
    expType = params.get("expType")
    submit = params.get("submit")

    # need a different way to handle errors returned by crisprSearch()

    errMsg = (
        "<p>Something unexpected occured. This is probably a bug, please contact us at %s and send us the url of this page.</p>"
        % (contactEmail)
    )

    printTeforBodyStart()
    if (
        submit is None
        and "batchId" not in params
        and "geneIds" not in params
        and "warnMsg" not in params
    ):
        printAssistant(params)
    global doCfdFix
    if "fixCfd" in params:
        doCfdFix = True

    if submit and "newSearch" in params:
        printCrisporBodyStart()
        crisprSearch(params)
    elif submit and (koGeneId or ("startSeq" in params and "endSeq" in params)):
        if expType == "ko":
            # Knock-out mode
            if koGeneId is not None:
                pam = params.get("pam")
                koMethod = params.get("koMethod")
                if koMethod == "excision":
                    targetLen = int(params.get("flankLen"))
                elif koMethod == "promoter":
                    targetLen = int(params.get("promoterLen"))
                else:
                    targetLen = None

                # select a base editor anyway to prevent bugs
                if koMethod == "stop" and pam not in ["NGN", "NRN"]:
                    global MAXSEQLEN
                    MAXSEQLEN = 1e4

                exonSelect = params.get("exonSelect")
                multiseq, geneModel = getExonsFromID(
                    koGeneId, org, pam, koMethod, targetLen, exonSelect
                )
                params["multiseq"] = multiseq
                params["geneModel"] = geneModel
                params["exonSelect"] = exonSelect
                geneModel = params.get("geneModel")  # ???
                if multiseq:
                    printCrisporBodyStart()
                    if geneModel and len(geneModel) == 0:
                        print(
                            "this is a non-coding transcript. Please choose another method to perform a knock-out on it, such as removing the promoter."
                        )
                        return
                    else:
                        crisprSearch(params)

        elif expType == "ki":
            targetRegion = params["targetRegions"]

            # Knock-in : "manual editing" mode
            if targetRegion == "seq":
                startSeq = params.get("startSeq")
                startSeq = re.sub("[\t\n\s]", "", startSeq)
                endSeq = params.get("endSeq")
                endSeq = re.sub("[\t\n\s]", "", endSeq)
                if startSeq and endSeq:
                    kiType, insertIdx, startSeq, insertSeq, clippedSeq = (
                        processCustomInsertSeq(startSeq, endSeq, targetRegion)
                    )
                    if kiType is None:
                        msg = "Insertion type currently not supported"
                        wrongInputRedirect(msg)

                    elif kiType == "noEdits":
                        msg = "No edits were found in the sequence"
                        wrongInputRedirect(msg)

                    elif kiType == "longReplacement":
                        msg = "replacement of sequences longer than 10bp is currently not supported. If you want to replace a large sequence (eg. a CDS), please refer to LINK"
                        wrongInputRedirect(msg)

                    elif kiType == "multiInsert":
                        msg = "Multiple insertions are currently not supported"
                        wrongInputRedirect(msg)

                    elif kiType == "multiDel":
                        msg = "Multiple deletions are currently not supported"
                        wrongInputRedirect(msg)

                    elif kiType == "longDel":
                        msg = "Deletions larger than 500bp are currently not supported"
                        wrongInputRedirect(msg)

                    else:
                        # don't forget to cleanSeq()
                        params["kiType"] = kiType
                        params["insertIdx"] = insertIdx
                        params["insertSeq"] = insertSeq
                        params["seq"] = startSeq
                        params["clippedSeq"] = clippedSeq

            # Knock-in : protein tagging mode
            elif (
                targetRegion == "gene"
                and koGeneId
                and (
                    ("tagseq" in params and "linkerseq" in params)
                    or "insertSeq" in params
                    or (
                        "markerseq" in params
                        and "expressionSeq" in params
                        and "qTag" in params
                    )
                )
            ):
                insertPos = params["insertpos"]

                linkerseq = params.get("linkerseq")
                tagseq = params.get("tagseq")

                markerseq = params.get("markerseq")
                expressionSeq = params.get("expressionSeq")
                qTag = params.get("qTag")
                if qTag is not None:
                    params["kiType"] = "qTag"
                else:
                    params["kiType"] = "tagging"

                try:
                    targetSeq, targetPos, insertIdx, geneModel = getTargetSeq(params)
                    if (linkerseq and tagseq) or (markerseq and expressionSeq and qTag):
                        tagNames, insertSeq = getInsertSeq(
                            linkerseq, tagseq, markerseq, expressionSeq, qTag, insertPos
                        )
                        params["insertSeq"] = insertSeq
                        params["tagNames"] = tagNames

                    params["seq"] = targetSeq
                    params["pos"] = targetPos
                    params["insertIdx"] = insertIdx
                    params["geneModel"] = geneModel
                except ValueError as err:
                    error = str(err)
                    if error == "frameErr":
                        print(
                            """<p>The insert sequence you entered is not in frame. Note that any charater other than "ATGCN" (case insensitive) is automatically removed.</p>"""
                        )
                    elif error == "insertErr":
                        print(
                            """<p>The target sequence you entered contains multiple insertion sites "//". Please use only one.</p>"""
                        )
            if "seq" in params:
                printCrisporBodyStart()
                crisprSearch(params)
            else:
                printAssistant(params)
                printKiForm(params)

        # selection of a transcriptId / exon in classic mode
        elif koGeneId and "exonSelect" in params:
            exonSelect = int(params["exonSelect"])
            if "seq" in params:
                params.pop("seq")
            try:
                multiseq, _ = getExonsFromID(koGeneId, org, pam, method="allExons")
                exonIds = [seqInfo[0] for seqInfo in multiseq]
                if exonSelect in exonIds:
                    for seqInfo in multiseq:  # seqInfo is (exonId, posStr)
                        if seqInfo[0] == exonSelect:
                            params["pos"] = seqInfo[1]
                            break
                else:
                    raise ValueError
            except ValueError:
                print(
                    "<p>The exon you selected is too small to be processed (< 23bp). Please select another exon</p>"
                )
            if "pos" in params:
                printCrisporBodyStart()
                try:
                    crisprSearch(params)
                except ValueError:
                    print(errMsg)

    elif "batchId" in params and "satMut" not in params:
        printCrisporBodyStart()

        if "donorType" in params and submit:
            (
                HA5,
                HA3,
                insertSeq,
                recodedArmSeq,
                mutEvents,
                noModel,
                recodeArm,
                HA5repeats,
                HA3repeats,
            ) = writeDonorSeq(params)
            showDonor(
                HA5,
                HA3,
                insertSeq,
                recodedArmSeq,
                mutEvents,
                noModel,
                recodeArm,
                HA5repeats,
                HA3repeats,
                params,
            )
        elif "doRecoding" in params:
            donorDesignPage(params)
        elif "pamId" in params and "doRecoding" not in params:
            if "guideSeq" in params:
                showSecondaryStructure(params)
            elif "pam" in params:
                primerDetailsPage(params)
            elif "otPrimers" in params:
                otPrimerPage(params)
            elif "showMh" in params:
                microHomPage(params)
            else:
                errAbort("Unrecognized CGI parameters.")
        elif "donorType" not in params:
            try:
                crisprSearch(params)
            except ValueError:
                print(errMsg)

    elif "satMut" in params:
        printCrisporBodyStart()
        printSatMutPage(params)

    elif (
        ("seq" in params or "pos" in params)
        and "org" in params
        and ("pam" in params or "customPAM" in params)
    ):
        printCrisporBodyStart()
        try:
            crisprSearch(params)
        except ValueError:
            print(errMsg)

    # need to rewrite this in a cleaner way

    elif "batchId" not in params and submit is None:
        if expType == "ko":
            printKoForm(params)
        elif expType == "ki":
            printKiForm(params)
        elif expType == "classic":
            printForm(params)
        elif "geneIds" in params:
            printCrisporBodyStart()
            printLibGuides(params)
        elif "libDesign" in params:
            # printCrisporBodyStart()
            printLibForm(params, returnLink=False)

        else:
            printForm(params)

        printReleaseNote()

    elif "batchId" not in params:
        printAssistant(params)
        if "libDesign" in params:
            printLibForm(params)
        elif params.get("expType") == "ki":
            printKiForm(params)
        elif params.get("expType") == "ko":
            printKoForm(params)
        else:
            printForm(params)

        printReleaseNote()


def processCustomInsertSeq(startSeq, endSeq, targetRegion):
    """from the starting sequence and the edited sequence, returns the knock-in
    type (insertion, replacement or substitution), the insert site position index
    and the insert sequence"""

    if startSeq.lower() == endSeq.lower():
        kiType = "noEdits"
        return kiType, None, None, None, None

    noEditEndSeq = "".join([base for base in endSeq if base.islower()])

    # check if all the edits are in uppercase
    doProcessInsert = False
    replace = False
    if startSeq.lower() == noEditEndSeq:
        doProcessInsert = True
    elif len(startSeq) - len(noEditEndSeq) == 1 and len(startSeq) == len(endSeq):
        doProcessInsert = True
    elif noEditEndSeq == endSeq and len(endSeq) < len(startSeq):
        kiType = "deletion"

        # new code : using difflib

        matcher = difflib.SequenceMatcher(
            a=startSeq.lower(), b=endSeq.lower(), autojunk=False
        )
        deletions = []
        for tag, i1, i2, j1, j2 in matcher.get_opcodes():
            if tag == "delete":
                deletions.append({"start": i1, "seq": startSeq[i1:i2]})

        if len(deletions) == 1:
            deletion_info = deletions[0]
            insertIdx = deletion_info["start"]
            insertSeq = deletion_info["seq"]
            # for deletions, don't clip the sequence to a maximum length
            clippedSeq = False

            # insertSeq = ''.join(deletions)

            # will never happen now

            if len(insertSeq) > 500:
                return "longDel", None, None, None, None
            else:
                return kiType, insertIdx, startSeq, insertSeq, clippedSeq
        elif len(deletions) > 1:
            return "multiDel", None, None, None, None
        else:
            return None, None, None, None, None

    elif len(startSeq) - len(noEditEndSeq) > 1:
        replace = True
        doProcessInsert = True
    else:
        return None, None, None, None, None

    if doProcessInsert:
        editSeqs = []
        targetPos = 0  # position on the non edited sequence
        for basePos, base in enumerate(endSeq):

            if base.islower():
                targetPos += 1
            stretchStart = (basePos == 0 and base.isupper()) or (
                basePos > 0 and base.isupper() and endSeq[basePos - 1].islower()
            )
            # find stretches of edited sequences
            # new stretch
            if stretchStart:
                editStretch = []
                for editBase in endSeq[basePos:]:
                    if editBase.islower():
                        editSeqs.append((basePos, "".join(editStretch)))
                        break
                    else:
                        editStretch.append(editBase)
        if (
            len(editSeqs) == 0
        ):  # multiple insertions need to be close enough to use a single donor DNA
            return None, None, None, None, None
        if len(editSeqs) > 1:
            editDists = [e[0] for e in editSeqs]
            minDist, maxDist = min(editDists), max(editDists)

            # get the last edit
            for pos, edit in editSeqs:
                if pos != maxDist:
                    continue
                lastEdit = edit

            # if multiple edits are < 10bp apart, merge them and treat it as a replacement
            if maxDist - minDist < 10:
                kiType = "replacement"
                insertIdx = editDists[0]
                insertSeq = endSeq[minDist: maxDist + len(lastEdit)].upper()

                if len(insertSeq) > 10:
                    return "multiInsert", None, None, None, None
            else:
                return "multiInsert", None, None, None, None

        elif len(editSeqs) == 1 and replace:
            if len(editSeqs[0][1]) > 10:
                kiType = "longReplacement"
            else:
                kiType = "replacement"
        elif (
            len(editSeqs) == 1
            and len(editSeqs[0][1]) == 1
            and len(startSeq) > len(noEditEndSeq)
        ):
            kiType = "substitution"
        elif len(editSeqs) == 1 and len(startSeq) == len(noEditEndSeq):
            kiType = "insertion"
        else:
            return None, None, None, None, None
        insertIdx, insertSeq = editSeqs[0]

    # clip the sequence 60bp in 5' and 3' of the editing site
    clippedSeq = False
    if insertIdx > 101:
        clippedSeq = True
        newInsertIdx = 100
        seqStart = insertIdx - newInsertIdx
    else:
        seqStart = 0
        newInsertIdx = insertIdx

    if len(startSeq) - insertIdx > 101:
        clippedSeq = True
        seqEnd = insertIdx + 100
    else:
        seqEnd = len(startSeq)

    newStartSeq = startSeq[seqStart:seqEnd]

    return kiType, newInsertIdx, newStartSeq, insertSeq, clippedSeq


def getInsertSeq(linkerSeq, tagSeq, markerSeq, expressionSeq, qTag, insertPos):
    "from a tag and a linker sequence, return the insert sequence to be used for the HDR donor"

    if (linkerSeq and tagSeq) or (
        linkerSeq and tagSeq and markerSeq and expressionSeq and qTag
    ):
        if insertPos == "Nter":
            insertSeq = taggingSeqs[tagSeq] + taggingSeqs[linkerSeq]
            tagNames = [tagSeq, linkerSeq]
        else:
            insertSeq = taggingSeqs[linkerSeq] + taggingSeqs[tagSeq]
            tagNames = [linkerSeq, tagSeq]

    elif markerSeq and expressionSeq and qTag:
        if insertPos == "Nter":
            insertSeq = (
                taggingSeqs["lox71"]
                + taggingSeqs[markerSeq]
                + taggingSeqs[expressionSeq]
                + taggingSeqs["lox66"]
                + taggingSeqs[qTag]
            )
            tagNames = ["lox71", markerSeq, expressionSeq, "lox66", qTag]
        else:
            insertSeq = (
                taggingSeqs[qTag]
                + taggingSeqs["lox71"]
                + taggingSeqs[expressionSeq]
                + taggingSeqs[markerSeq]
                + taggingSeqs["lox66"]
            )
            tagNames = [qTag, "lox71", expressionSeq, markerSeq, "lox66"]

    return tagNames, insertSeq


def getExonsFromID(geneId, org, pam, method, targetLen=None, exonSelect=None):
    """from a geneId, returns a list of tuples [(exonId, exonPosStr)]
    and the gene model (only within the cds)"""

    if targetLen is None and (method == "excision" or method == "promoter"):
        raise ValueError("targetLen needs to be set for this method")

    # set pam-dependent variables
    pam = setupPamInfo(pam)

    # few guides are found using this method : allow up to 10kb
    if method == "stop":
        maxLen = MAXSEQLEN
    elif pam in verySlowPams:
        maxLen = MAXSEQLEN3
    elif isSlowPam(pam):
        maxLen = MAXSEQLEN2
    else:
        maxLen = MAXSEQLEN

    chrom, strand, allExons = getGenePos(geneId, org, method, targetLen)
    if method in ["frameshift", "stop", "splicing"]:
        if len(allExons) == 0:
            return None, None
        # exons is the list to be returned
        geneModel = getGeneModel(allExons, strand)
    else:
        geneModel = None
    if method == "splicing":
        spliceDist = 10
        exons = []
        for exonNumber, (start, end) in enumerate(allExons):
            if exonSelect == "allExons" or exonNumber == int(exonSelect):
                exons.append((start - spliceDist, start + spliceDist))
                exons.append((end - spliceDist, end + spliceDist))
    elif method in ["frameshift", "stop"]:
        exons = getFirstThird(allExons, strand, GUIDELEN, maxLen, method)
    else:
        exons = allExons

    return formatExonPos(exons, chrom, strand, len(pam)), geneModel


def getGeneModel(allExons, strand):
    "returns a list of tuples (type, exonId, length) to display the gene model"

    geneModel = []

    for i, (start, end) in enumerate(allExons):
        exonLen = end - start
        if i + 1 < len(allExons):
            nextStart, nextEnd = allExons[i + 1]
            if strand == "+":
                intronLen = nextStart - end
            else:
                intronLen = start - nextEnd
            geneModel.append(("exon", i, exonLen))
            geneModel.append(("intron", i, intronLen))
        else:
            geneModel.append(("exon", i, exonLen))
    return geneModel


def getGenePos(geneID, org, method, targetLen):

    genomeDir = genomesDir
    twoBitFname = getTwoBitFname(org)
    genomePath = "%(genomeDir)s/%(org)s/" % locals()
    genomeFiles = os.listdir(genomePath)
    gpFiles = [f for f in genomeFiles if f.endswith(".gp")]

    if "SYM" in geneID:
        commonExons = True
        geneID = geneID.split("~")[0]
        allExons = []
    else:
        commonExons = False

    for gpFile in gpFiles:
        gpFilePath = os.path.join(genomePath, gpFile)
        with open(gpFilePath, "r") as genePred:

            for geneLine in genePred:
                geneLine = geneLine.split("\t")
                if (commonExons is False and geneLine[0] == geneID) or (
                    commonExons is True
                    and len(geneLine) >= 12
                    and geneLine[11] == geneID
                ):
                    geneInfo = {
                        "name": geneLine[0],
                        "chrom": geneLine[1],
                        "strand": geneLine[2],
                        "txStart": int(geneLine[3]),
                        "txEnd": int(geneLine[4]),
                        "exonStarts": [
                            int(start) for start in geneLine[8].rstrip(",").split(",")
                        ],
                        "exonEnds": [
                            int(end) for end in geneLine[9].rstrip(",").split(",")
                        ],
                        "cdsStart": int(geneLine[5]),
                        "cdsEnd": int(geneLine[6]),
                    }
                    if commonExons:
                        allExons.append(geneInfo)
                    else:
                        # for now, stop after finding one matching ID
                        # print a warning message if several transcripts are found ? (will be slower)
                        break
    # get the common exons
    if commonExons:
        if not allExons:
            raise ValueError("")

        def intersect(aStarts, aEnds, bStarts, bEnds):
            resStarts, resEnds = [], []
            for aStart, aEnd in zip(aStarts, aEnds):
                for bStart, bEnd in zip(bStarts, bEnds):
                    commonStart = max(aStart, bStart)
                    commonEnd = min(aEnd, bEnd)
                    if commonStart < commonEnd:
                        resStarts.append(commonStart)
                        resEnds.append(commonEnd)
            return resStarts, resEnds

        # Initialize with the first transcript
        resStarts = allExons[0]["exonStarts"]
        resEnds = allExons[0]["exonEnds"]
        txStart, txEnd = allExons[0]["txStart"], allExons[0]["txEnd"]
        cdsStart, cdsEnd = allExons[0]["cdsStart"], allExons[0]["cdsEnd"]

        for i in range(1, len(allExons)):
            currentInfo = allExons[i]
            nextStarts, nextEnds = intersect(
                resStarts, resEnds, currentInfo["exonStarts"], currentInfo["exonEnds"]
            )

            # Only update intersection and bounds if there are still common exons
            if nextStarts:
                res_starts, res_ends = nextStarts, nextEnds
                tx_start = max(txStart, currentInfo["txStart"])
                tx_end = min(txEnd, currentInfo["txEnd"])
                cds_start = max(cdsStart, currentInfo["cdsStart"])
                cds_end = min(cdsEnd, currentInfo["cdsEnd"])
            else:
                # No common exons with this specific transcript, skipping its contribution to bounds
                continue

        geneInfo = allExons[0].copy()
        geneInfo["exonStarts"] = resStarts
        geneInfo["exonEnds"] = resEnds
        geneInfo["txStart"] = txStart
        geneInfo["txEnd"] = txEnd
        geneInfo["cdsStart"] = cdsStart
        geneInfo["cdsEnd"] = cdsEnd

    if geneInfo is None:
        raise ValueError("")
    featureList = []
    chrom = geneInfo["chrom"]
    strand = geneInfo["strand"]

    # Get 100 bases around a targetLen region upstream of the TSS
    if method == "promoter":
        if strand == "+":
            tss = geneInfo["txStart"]
            promoterDownPos = (tss - 100, tss)
            upEnd = promoterDownPos[1] - targetLen
            promoterUpPos = (upEnd - 100, upEnd)

        else:
            tss = geneInfo["txEnd"]
            promoterUpPos = (tss + targetLen, tss + targetLen + 100)
            promoterDownPos = (tss, tss + 100)

        featureList.append(promoterUpPos)
        featureList.append(promoterDownPos)

    # Get targetLen bases upstream of the TSS and downstream of the TES
    # returns (upstreampos, downstreampos)
    if method == "excision":
        tss = geneInfo["txStart"]
        tes = geneInfo["txEnd"]
        upstreamPos = (tss - targetLen, tss)
        downstreamPos = (tes, tes + targetLen)
        if strand == "+":
            featureList = [upstreamPos, downstreamPos]
        else:
            featureList = [downstreamPos, upstreamPos]

    # get the positions of the coding exons (5' and 3' UTRs removed)
    elif method in ["frameshift", "stop", "allExons", "splicing"]:
        cdsStart = geneInfo["cdsStart"]
        cdsEnd = geneInfo["cdsEnd"]
        exonStarts = geneInfo["exonStarts"]
        exonEnds = geneInfo["exonEnds"]

        if cdsEnd - cdsStart <= 1:
            featureList = [(start, end) for start, end in zip(exonStarts, exonEnds)]

        for start, end in zip(exonStarts, exonEnds):
            if end <= cdsStart or start >= cdsEnd:
                continue
            codingExonStart = max(start, cdsStart)
            codingExonEnd = min(end, cdsEnd)
            featureList.append((codingExonStart, codingExonEnd))

        if strand == "-":
            featureList = featureList[::-1]

    return chrom, strand, featureList


def getFirstThird(exons, strand, pamlen, maxLen, method, minLen=250):
    """from a list of exon positions get the first third of the sequence (within the boundaries of minLen / maxLen)
    input is a list of tuples [(start, end)]. Returns the same format as the input list.
    """

    totalLen = sum([end - start for start, end in exons])

    if totalLen <= minLen:
        return exons

    # few guides can introduce a STOP codon : search in the first two thrids of the coding sequence
    if method == "stop":
        # use the maximum possible length in this mode
        # thirdLen = 2 * math.ceil(totalLen / 3)
        thirdLen = totalLen
    else:
        thirdLen = math.ceil(totalLen / 3)
    if thirdLen > maxLen:
        thirdLen = maxLen

    firstThirdExons = []
    currentLen = 0
    for start, end in exons:
        exonLen = end - start
        if currentLen + exonLen <= thirdLen:
            firstThirdExons.append((start, end))
            currentLen += exonLen
        # the exon can be cropped to fit the first third of the coding sequence
        else:
            remaining = thirdLen - currentLen
            if remaining > pamlen:
                if strand == "+":
                    firstThirdExons.append((start, start + remaining))
                else:
                    firstThirdExons.append((end - remaining, end))
            break

    return firstThirdExons


def formatExonPos(exons, chrom, strand, pamlen):
    """
    formats exon positions to be processed by CrisprSearch()"
    removes exons < PAMLEN and format posisions as posStr (chrom:start-end:strand)
    returns a list of tuples [(exonId, exonPosStr)] where exonId is the original exon number.
    """

    multiseq = []
    for exonId, (start, end) in enumerate(exons):
        exonLen = end - start
        if exonLen < pamlen:
            continue
        exonPosStr = "%(chrom)s:%(start)s-%(end)s:%(strand)s" % locals()
        multiseq.append((exonId, exonPosStr))

    return multiseq


def getPosAndSeq(org, seq, posStr, batchId):
    """from an insert position, writes the sequence of the donor DNA in the batch params.
    additionnaly, returns the formatted sequence and coordinates of the target region
    """

    batchInfo = readBatchAsDict(batchId)
    insertPos = batchInfo.get("insertpos")
    insertIdx = batchInfo["insertIdx"]
    codonTable = buildCodonTable(key="aa")
    kiType = batchInfo.get("kiType")

    noPerfectMatch = batchInfo.get("noPerfectMatch")

    # input is a sequence
    if posStr is None and seq:

        posStr = findPerfectMatch(batchId, seq, org, noPerfectMatch=noPerfectMatch)
        batchInfo["posStr"] = posStr

    # input is a transcriptID
    elif seq is None and posStr:
        seq = getSeq(org, posStr)
        logging.info(f"seq : {seq}")
        # Annotation of START and STOP codons (uppercase)
        if kiType in ["tagging", "qTag"]:
            if (
                insertPos == "Nter"
                and seq[insertIdx - 3: insertIdx].upper() in codonTable["M"]
            ):
                targetSeq = (
                    seq[0: insertIdx - 3].lower()
                    + seq[insertIdx - 3: insertIdx].upper()
                    + seq[insertIdx:].lower()
                )
                seq = targetSeq
            elif (
                insertPos == "Cter"
                and seq[insertIdx: insertIdx + 3].upper() in codonTable["*"]
            ):
                targetSeq = (
                    seq[0:insertIdx].lower()
                    + seq[insertIdx: insertIdx + 3].upper()
                    + seq[insertIdx + 3:].lower()
                )
                seq = targetSeq
            else:
                batchInfo["nonCoding"] = True

        batchInfo["seq"] = seq

    chrom, start, end, strand = parsePos(posStr)

    writeBatchAsDict(batchInfo, batchId)

    return seq, posStr


def writeDonorSeq(params):
    """From the parameters in the donor design from,
    returns the sequence of the donor DNA
    """

    batchId = params["batchId"]
    batchInfo = readBatchAsDict(batchId)

    guideInfo = params["guideInfo"]
    guideSeq = params["guideSeq"]
    arm5Len = int(params["arm5Len"])
    arm3Len = int(params["arm3Len"])
    donorType = params["donorType"]
    # doBarcode = params.get("doBarcode")
    donorType = params["donorType"]
    recodePam = params.get("recodePam")
    pamId = params["pamId"]
    recodeGap = params.get("recodeGap")
    recodeSeed = params.get("recodeSeed")
    selGeneModel = params.get("geneModelSelection")
    manualExStart = params.get("manualExStart")
    manualExEnd = params.get("manualExEnd")
    manualExFrame = params.get("manualExFrame")
    useManualAnnotation = params.get("useManualAnnotation")
    # trimGC = params.get("trimGC")
    # trimHomopolymers = params.get("trimHomopolymers")
    # trimRepeat = params.get("trimRepeats")
    org = batchInfo["org"]
    geneId = batchInfo.get("koGeneId")
    seq = batchInfo["seq"]
    posStr = batchInfo["posStr"]
    insertIdx = batchInfo["insertIdx"]
    kiType = batchInfo.get("kiType")
    insertSeq = batchInfo.get("insertSeq")
    chrom, start, end, strand = parsePos(posStr)
    # convert to string representation of the original tuple to a tuple again (not optimal at all)
    guideInfo = deserialize(guideInfo, inType="tuple")
    pamPat = pamId.split(".")[0]

    pamPat = setupPamInfo(pamPat)

    if geneId is None:
        transId = params.get("selTransId")
    else:
        transId = geneId

    if (
        manualExStart is not None
        and manualExEnd is not None
        and manualExFrame is not None
        and useManualAnnotation == "True"
    ):
        selExon = [
            (
                1,
                int(manualExStart),
                int(manualExEnd),
                int(manualExFrame),
                int(manualExFrame),
                0,
                strand,
            )
        ]

    elif selGeneModel:
        exonInfo, maxTransIdLen = getExonInfo(org, selGeneModel, posStr)
        # no transcript at this position : load an empty non-coding region
        if len(exonInfo) == 0:
            selExon = [(0, 0, 0, 0, 0, 0, strand)]
        # retreive the exon corresponding to the selected transcript
        for transcriptInfo in exonInfo:
            transcriptId, symbol = transcriptInfo
            if transId == transcriptId:
                selExon = exonInfo[transcriptInfo]
                break
            # if selGeneModel isn't transmitted, select a transcript anyway
            else:
                selExon = exonInfo[transcriptInfo]
    else:
        selExon = None

    if strand == "+":
        insertCoord = int(start + insertIdx)
    else:
        insertCoord = int(end - insertIdx)

    if kiType in ["deletion", "replacment"]:
        if strand == "+":
            editEnd = insertCoord + len(insertSeq)
            arm5start, arm5end = checkCoords(
                insertCoord - arm5Len, insertCoord, org, chrom
            )
            arm3start, arm3end = checkCoords(editEnd, editEnd + arm3Len, org, chrom)
        else:
            editEnd = insertCoord - len(insertSeq)
            arm5start, arm5end = checkCoords(editEnd - arm5Len, editEnd, org, chrom)
            arm3start, arm3end = checkCoords(
                insertCoord, insertCoord + arm3Len, org, chrom
            )
        newInsertSeq = ""
    else:
        arm5start, arm5end = checkCoords(insertCoord - arm5Len, insertCoord, org, chrom)
        arm3start, arm3end = checkCoords(insertCoord, insertCoord + arm3Len, org, chrom)
        newInsertSeq = insertSeq.upper()

    # for substitutions / replacments, remove the edited region in the homology arm and extend the arm by Nbp
    if kiType in ["substitution", "replacement"]:
        if strand == "+":  # à vérifier !! -> OK
            arm3start += len(insertSeq)
            arm3end += len(insertSeq)
        else:
            arm5end -= len(insertSeq)
            arm5start -= len(insertSeq)

    HA5 = getSeq(
        org,
        "%s:%s-%s:%s" % (chrom, arm5start, arm5end, strand),
        maxlen=False,
        minlen=False,
    )
    HA3 = getSeq(
        org,
        "%s:%s-%s:%s" % (chrom, arm3start, arm3end, strand),
        maxlen=False,
        minlen=False,
    )

    # annotate reapeats (= lowercase basees) now
    HA5repeats = findRepeats(HA5)
    HA3repeats = findRepeats(HA3)

    HA5 = HA5.lower()
    HA3 = HA3.lower()

    if strand == "-":
        HA5, HA3 = (HA3, HA5)

    # get the position of the PAM + guide
    # recoding
    if (recodePam or recodeGap) and (
        (selGeneModel is not None)
        or (
            manualExStart is not None
            and manualExEnd is not None
            and manualExFrame is not None
        )
    ):

        noModel = False
        # load codon frequency file
        codonFreqFname = "%s_codonFrequency.json" % org
        codonFreqFile = join(genomesDir, org, codonFreqFname)
        if isfile(codonFreqFile):
            with open(codonFreqFile) as jsonData:
                codonFreq = json.load(jsonData)
        else:
            codonFreq = None

        if strand == "+":
            annotationCoords, recodeCoords, recodeArm = getArmCoords(
                HA5,
                HA3,
                strand,
                seq,
                insertIdx,
                guideSeq,
                guideInfo,
                kiType,
                donorType,
                insertSeq,
                selExon=selExon,
                manual=useManualAnnotation,
            )
        else:
            annotationCoords, recodeCoords, recodeArm = getArmCoords(
                HA3,
                HA5,
                strand,
                seq,
                insertIdx,
                guideSeq,
                guideInfo,
                kiType,
                donorType,
                insertSeq,
                selExon=selExon,
                manual=useManualAnnotation,
            )
        if recodeArm == "HA5":
            # need to attempt recoding the PAM, then recode the seed region if no recoding could happen
            # call recodeDonor two times
            # will need to move the logic to a new function
            recodedArmSeq, mutEvents = getRecodeCodons(
                HA5,
                annotationCoords,
                recodeCoords,
                recodePam,
                recodeSeed,
                recodeGap,
                guideInfo,
                recodeArm,
                pamPat,
                codonFreq,
            )

        else:
            recodedArmSeq, mutEvents = getRecodeCodons(
                HA3,
                annotationCoords,
                recodeCoords,
                recodePam,
                recodeSeed,
                recodeGap,
                guideInfo,
                recodeArm,
                pamPat,
                codonFreq,
            )
    else:
        # if no recoding could happen because no annotation file is available, signal it
        if (recodePam or recodeSeed or recodeGap) and (
            (selGeneModel is None)
            or manualExStart is None
            or manualExEnd is None
            or manualExFrame is None
        ):
            noModel = True
        else:
            noModel = False

        mutEvents = None
        recodedArmSeq = None
        recodeArm = None

    return (
        HA5,
        HA3,
        newInsertSeq,
        recodedArmSeq,
        mutEvents,
        noModel,
        recodeArm,
        HA5repeats,
        HA3repeats,
    )


def findRepeats(seq):
    """
    in a sequence, returns the coords of repeats (flagged by repeatMasker as lowercase)
    as a list of tuples (start, end)
    """
    repeats = []
    seqLen = len(seq)
    i = 0

    while i < seqLen:
        if seq[i].islower():
            start = i
            while i < seqLen and seq[i].islower():
                i += 1
            repeats.append((start, i))
        else:
            i += 1

    return repeats


def getArmCoords(
    HA5,
    HA3,
    strand,
    seq,
    insertIdx,
    guideSeq,
    guideInfo,
    kiType,
    donorType,
    insertSeq,
    selExon=None,
    manual=None,
):
    """
    on both homology arms, get the coordinates of the codons
    and the regions to recode

    returns :
    HA5codonPos = [codon starts]
    HA3codonPos = [codon starts]
    splittedCodonPos = [codonStart(HA5), codonEnd(HA3)]
    pamCoords = [pamStart, pamEnd]
    seedCoords = [seedStart, seedEnd]
    cutPos = int
    recodeArm = "HA5" or "HA3" : the homology arm to recode
    """

    # convert the coordinates of the guide and the exon relative to the homology arms
    pamSeq, guideStart, guideStrand = guideInfo
    guideStart = int(guideStart)

    # position of the PAM
    if (guideStrand == "+" and not pamIsFirst) or (guideStrand == "-" and pamIsFirst):
        pamStart = guideStart + len(guideSeq)
    else:
        pamStart = guideStart - 3
    if pamStart < insertIdx:
        # get the position in HA5
        # 5'-----<<NGG--------------3'
        # 5'---------------\/-------3'

        pamStartPos = len(HA5) - (insertIdx - pamStart)
        pamEndPos = pamStartPos + len(pamSeq)
        recodeArm = "HA5"
        pamCoords = (pamStartPos, pamEndPos)
    else:
        # get the position in HA3
        # 5'-----\/--------------3'
        # 5'-------------<<NGG---3'

        if kiType in ["deletion", "substitution", "replacement"]:
            # for deletions, the 3' homology arm starts after insertIdx
            # for substituions / replacements, the original sequence is removed from HA3
            pamStartPos = pamStart - (insertIdx + len(insertSeq))

        else:
            pamStartPos = pamStart - insertIdx
        pamEndPos = pamStartPos + len(pamSeq)
        recodeArm = "HA3"
        pamCoords = (pamStartPos, pamEndPos)

    # position of the seed region (15 pam-proximal bases of the spacer)
    if (guideStrand == "+" and not pamIsFirst) or (guideStrand == "-" and pamIsFirst):
        seedStartPos = pamStartPos - 15
        if seedStartPos < 0:
            seedStartPos = 0
        seedEndPos = pamStartPos
    else:
        seedStartPos = pamEndPos
        seedEndPos = pamEndPos + 15
        # do not extent coord beyond the homology arm
        if recodeArm == "HA5" and seedEndPos > len(HA5):
            seedEndPos = len(HA5)
        elif recodeArm == "HA3" and seedEndPos > len(HA3):
            seedEndPos = len(HA3)

    seedCoords = (seedStartPos, seedEndPos)

    # position of the gap between cut site and insertion site
    if guideStrand == "+" or (guideStrand == "-" and pamIsFirst):
        cutPos = pamStartPos - 3
    else:
        cutPos = pamEndPos + 3

    if recodeArm == "HA5" and cutPos < len(HA5):
        gapStartPos = cutPos
        gapEndPos = len(HA5)
        gapCoords = (gapStartPos, gapEndPos)
    elif recodeArm == "HA3" and cutPos > 0:
        gapStartPos = 0
        gapEndPos = cutPos
        gapCoords = (gapStartPos, gapEndPos)
    else:
        gapCoords = None

    # convert exon coordinates
    if selExon:
        UTR5coords = []
        kozakCoords = []
        UTR3coords = []
        spliceCoords = (
            []
        )  # position of splice sites : 5bp upstream / downstream of coding exons (except upstream of the first and downstream of the last exon)
        oofCoords = []
        codonPos = []
        # splittedCodonPos = []  # for codons that overlap with the editing site
        for (
            exonNumber,
            exonStart,
            exonEnd,
            exonFrame,
            oldExonFrame,
            nextExonFrame,
            exonStrand,
        ) in selExon:

            isUTR5 = exonNumber == -1 and exonFrame == -1
            isUTR3 = exonNumber != -1 and exonFrame == -1

            # number of coding exons
            codingExonLen = len(
                [
                    exonStart
                    for exonNumber, exonStart, exonEnd, exonFrame, oldExonFrame, nextExonFrame, exonStrand in selExon
                    if exonNumber != -1
                ]
            )

            # whole exon in the 5' homology arm
            if exonStart < insertIdx and exonEnd < insertIdx:
                if recodeArm == "HA3":
                    continue
                exonStartPos = len(HA5) - (insertIdx - exonStart)
                exonEndPos = len(HA5) - (insertIdx - exonEnd)

                if isUTR5:
                    UTR5coords.append((exonStartPos, exonEndPos - 6))
                    kozakCoords.append((exonEndPos - 6, exonEndPos))
                    continue

                elif isUTR3:
                    UTR3coords.append((exonStartPos + 5, exonEndPos))
                    continue

                # position of the splicing donor site
                if exonNumber > 1:
                    spliceDon = (
                        exonStartPos - 5 if exonStartPos - 5 > 0 else 0,
                        exonStartPos,
                    )
                    spliceCoords.append(spliceDon)

                # need a way to figure out the last exon : len(selExon) when exonFrame != -1 probably
                # no, because selExon only contains the exons overlapping the input sequence
                # just remove UTR coords from splicing coords

                # position of the splicing acceptor site
                if exonNumber < codingExonLen:
                    spliceAcc = (
                        exonEndPos,
                        exonEndPos + 5 if exonEndPos + 5 < len(HA5) else len(HA5),
                    )
                    spliceCoords.append(spliceAcc)

                exonOffset = (3 - exonFrame) % 3
                for i in range(exonStartPos + exonOffset, exonEndPos, 3):
                    if i + 3 > exonEndPos:
                        break
                    codonPos.append(i)

            # whole exon in the 3' homology arm
            elif exonStart > insertIdx and exonEnd > insertIdx:
                if recodeArm == "HA5":
                    continue

                if kiType in [
                    "deletion",
                    "substitution",
                    "replacement",
                ] and exonEnd < insertIdx + len(insertSeq):
                    continue
                exonStartPos = exonStart - insertIdx
                exonEndPos = exonEnd - insertIdx
                print(exonStartPos)
                # correct codon positions for replacements and deletions
                if kiType in ["deletion", "substitution", "replacement"]:
                    exonStartPos = exonStartPos - len(insertSeq)
                    print(exonStartPos)

                if isUTR5:
                    UTR5coords.append((exonStartPos, exonEndPos - 6))
                    kozakCoords.append((exonEndPos - 6, exonEndPos))
                    continue

                elif isUTR3:
                    UTR3coords.append((exonStartPos, exonEndPos))
                    continue

                # position of the splicing donor site
                if exonNumber > 1:
                    spliceDon = (
                        exonStartPos - 5 if exonStartPos - 5 > 0 else 0,
                        exonStartPos,
                    )
                    spliceCoords.append(spliceDon)

                # position of the splicing acceptor site
                if exonNumber < codingExonLen:
                    spliceAcc = (
                        exonEndPos,
                        exonEndPos + 5 if exonEndPos + 5 < len(HA3) else len(HA3),
                    )
                    spliceCoords.append(spliceAcc)

                exonOffset = (3 - exonFrame) % 3
                for i in range(exonStartPos + exonOffset, exonEndPos, 3):
                    if i + 3 > exonEndPos:
                        break
                    else:
                        codonPos.append(i)

            # edit site overlaps the exon
            else:
                if recodeArm == "HA5":
                    exonStartPos = len(HA5) - (insertIdx - exonStart)
                    exonEndPos = len(HA5)

                    if exonNumber > 1:
                        spliceDon = (
                            exonStartPos - 5 if exonStartPos - 5 > 0 else 0,
                            exonStartPos,
                        )
                        spliceCoords.append(spliceDon)

                else:
                    # don't process exons that end within the deletion / replacement
                    if kiType in [
                        "deletion",
                        "substitution",
                        "replacement",
                    ] and exonEnd < insertIdx + len(insertSeq):
                        continue

                    exonStartPos = 0
                    # here, exonFrame changes according to the "missing" part of the exon!
                    exonMissingLen = len(HA5) - (len(HA5) - (insertIdx - exonStart))
                    exonOffset = exonMissingLen % 3
                    exonFrame = (exonFrame + exonMissingLen) % 3

                    # for deletions, replacement or substitutions the homology arm has been shifted, so the exonFrame needs to be adjusted
                    if kiType in ["deletion", "substitution", "replacement"]:
                        exonFrame = (exonFrame + len(insertSeq) % 3) % 3
                        exonEndPos = (
                            exonStartPos + exonEnd - (insertIdx + len(insertSeq))
                        )
                    else:
                        exonEndPos = exonEnd - insertIdx
                    if exonNumber < codingExonLen:
                        spliceAcc = (
                            exonEndPos,
                            exonEndPos + 5 if exonEndPos + 5 < len(HA3) else len(HA3),
                        )
                        spliceCoords.append(spliceAcc)

                if isUTR5:
                    UTR5coords.append((exonStartPos, exonEndPos - 6))
                    kozakCoords.append((exonEndPos - 6, exonEndPos))
                    continue

                elif isUTR3:
                    UTR3coords.append((exonStartPos, exonEndPos))
                    continue
                # get the start position of codons in the homology arm to recode
                exonOffset = (3 - exonFrame) % 3
                for i in range(exonStartPos + exonOffset, exonEndPos, 3):
                    if i + 3 > exonEndPos:
                        # end of exon is not in frame : don't recode
                        oofCoords.append((i, exonEndPos))
                        break
                    codonPos.append(i)

    annotationCoords = (
        codonPos,
        UTR3coords,
        UTR5coords,
        kozakCoords,
        spliceCoords,
        oofCoords,
    )
    recodeCoords = (pamCoords, seedCoords, gapCoords)

    return annotationCoords, recodeCoords, recodeArm


def getRecodeCodons(
    HA,
    annotationCoords,
    recodeCoords,
    recodePam,
    recodeSeed,
    recodeGap,
    guideInfo,
    recodeArm,
    pamPat,
    codonFrequency,
):
    """
    /!\ Needs refactoring, alongside recodeDonor()

    from the coordinates in the homology arm, returns :
    - the positions of codons to recode
    - the positions of the regions where no recoding is allowed
    """

    codonPos, UTR3coords, UTR5coords, kozakCoords, spliceCoords, oofCoords = (
        annotationCoords
    )
    pamCoords, seedCoords, gapCoords = recodeCoords

    _, _, guideStrand = guideInfo
    recodeRegions = []
    # get the regions to mutate and combine overlaps

    if recodeGap and gapCoords:
        recodeRegions.append(gapCoords)

    UTRpos = set()
    for start, end in UTR5coords:
        for i in range(start, end):
            UTRpos.add(i)
    for start, end in UTR3coords:
        for i in range(start, end):
            UTRpos.add(i)

    if len(kozakCoords) > 0:
        kozakEnd = kozakCoords[0][1]
        startPos = {kozakEnd, kozakEnd + 1, kozakEnd + 2}
    else:
        startPos = {}

    # don't mutate at these positions
    kozakPos = set()
    for start, end in kozakCoords:
        for i in range(start, end):
            kozakPos.add(i)
    splicePos = set()
    for start, end in spliceCoords:
        for i in range(start, end):
            # avoid flagging splicing sites in UTRs
            if i in kozakPos or i in UTRpos:
                continue
            else:
                splicePos.add(i)

    if recodePam and pamCoords is not None:
        pamPos = set()
        # the positions of the N bases of the PAM motif : to be removed from gapCoords
        pamNpos = set()
        seedPos = set()
        pamStart, pamEnd = pamCoords
        seedStart, seedEnd = seedCoords
        # make sure to mutate only the non N bases of the PAM
        for i in range(pamStart, pamEnd):
            # NGG / CCN
            patPos = i - pamStart if guideStrand == "+" else pamEnd - i - 1
            if pamPat[patPos] == "N":
                pamNpos.add(i)
            elif pamEnd < len(HA):
                pamPos.add(i)
        for i in range(seedStart, seedEnd):
            seedPos.add(i)

    # recode the region between the cut site and edit site
    if recodeGap and gapCoords is not None:
        gapPos = set()
        gapStart, gapEnd = gapCoords
        for i in range(gapStart, gapEnd):
            if i not in pamPos and i not in seedPos and i not in pamNpos:
                gapPos.add(i)
    if len(pamPos) == 0 and len(seedPos) == 0 and len(gapPos) == 0:
        return HA, None
    noRecodePos = (kozakPos, splicePos, UTRpos, startPos)

    # first, attempt to recode the PAM motif
    if recodePam and pamCoords is not None:
        recodedCodons = []
        mutHA, mutEvents, recodedCodons = recodeDonor(
            HA, pamPos, noRecodePos, codonPos, codonFrequency, maxMut=1, motif="pam"
        )

        # couldn't recode the PAM motif : try to recode the seed region
        if mutHA.lower() == HA.lower() or (recodeGap and gapCoords is not None):
            # when recoding between the cut site and insertion site, try to recode all codons
            if recodeGap:
                maxMut = len(seedPos)
            else:
                maxMut = 2
            mutHA, seedMutEvents, recodedSeedCodons = recodeDonor(
                mutHA,
                seedPos,
                noRecodePos,
                codonPos,
                codonFrequency,
                excludeCodons=recodedCodons,
                maxMut=maxMut,
                motif="seed",
            )
            recodedCodons.extend(recodedSeedCodons)
            mutEvents.update(seedMutEvents.items())

        if recodeGap and gapCoords is not None:
            mutHA, gapMutEvents, _ = recodeDonor(
                mutHA,
                gapPos,
                noRecodePos,
                codonPos,
                codonFrequency,
                excludeCodons=recodedCodons,
                motif="gap",
            )
            mutEvents.update(gapMutEvents.items())
    else:
        mutHA = HA
        mutEvents = {}

    return mutHA, mutEvents


def recodeDonor(
    HA,
    recodePos,
    noRecodePos,
    codonPos,
    codonFrequency,
    excludeCodons=[],
    maxMut=None,
    motif="PAM",
):
    """
    introduce silent mutations in the homology arm based on the codons to recode
    """

    codonTable = buildCodonTable()
    revCodonTable = buildCodonTable(key="aa")
    mutHA = list(HA)

    kozakPos, splicePos, UTRpos, startPos = noRecodePos

    if maxMut is None:
        maxMut = len(recodePos)

    # dict to store the mutation events
    mutEvents = {}
    # print(f"RECODEPOS : {recodePos} <br> POSINCODON: {posInCodon}")

    # list to store successfully recoded codons (to avoid recoding them later)
    recodedCodons = []

    # don't recode previously recoded codons
    # e.g, if this codon was recoded to modify the PAM motif but it overlaps the region between the cut site and edit site,
    # don't recode it to avoid breaking the PAM-blocking mutation
    if excludeCodons:
        for excludeStart in excludeCodons:
            for i in range(excludeStart, excludeStart + 3):
                if i in recodePos:
                    recodePos.remove(i)

    # list of codons to mutate
    mutCodons = set()
    posInCodon = set()
    for codonStart in codonPos:
        for pos in range(codonStart, codonStart + 3):
            if pos in recodePos and codonStart + 3 <= len(HA):
                for i in range(pos, pos + 3):
                    posInCodon.add(i)
                mutCodons.add(codonStart)
                break

    mutCount = 0
    # positions to mutate in non-coding regions
    # where to mutate ??
    if len(posInCodon) < len(recodePos):
        mutNonCoding = sorted(set([pos for pos in recodePos if pos not in posInCodon]))
        transitions = {"A": "G", "G": "A", "T": "C", "C": "T"}
        # transversions = {"A": "C", "G": "T", "T": "G", "C": "A"}
        # bases = ["A", "T", "G", "C"]

        # for now, mutate every 3 bases
        for posIdx in (
            range(0, len(mutNonCoding), 3)
            if motif != "Seed"
            else range(0, len(mutNonCoding), -3)
        ):
            pos = mutNonCoding[posIdx]

            if mutCount > maxMut:
                break

            if pos >= len(HA) or pos < 0:
                mutEvents[("out of range", None, pos)] = ("", None, motif)
                break
            base = HA[pos].upper()
            if pos in kozakPos:
                mutEvents[(base, None, pos)] = (
                    "No (in kozak consensus sequence)",
                    None,
                    motif,
                )
            elif pos in splicePos:
                mutEvents[(base, None, pos)] = ("No (in a splicing site)", None, motif)
            else:
                # choice = [b for b in bases if b != base]
                # newBase = choice[pos % 3]  # select a "random" base for now
                mutCount += 1
                newBase = transitions[base]
                mutHA[pos] = newBase.upper()
                mutEvents[(base, None, pos)] = (newBase, None, motif)

    # in the seed region, iterate over the list of codons to prioritize recoding near the PAM
    for codonNb, codonStart in (
        enumerate(list(sorted(mutCodons)))
        if motif != "seed"
        else enumerate(reversed(list(sorted(mutCodons))))
    ):
        codon = HA[codonStart : codonStart + 3].upper()
        if codonFrequency:
            codonFreq = codonFrequency[codon][2]
        else:
            codonFreq = None
        aa = codonTable[codon]
        if not aa:
            continue
        synCodons = revCodonTable.get(aa)
        # check codons that introduce a mutation at recodePos

        if codonNb >= maxMut and len(mutEvents) >= maxMut:
            break

        # show a special message for the START codon
        if {codonStart, codonStart + 1, codonStart + 2} == startPos and HA[
            codonStart : codonStart + 3
        ].upper() == "ATG":
            mutEvents[(codon, "1", codonStart)] = ("No (START codon)", None, motif)
            continue

        lastFreq = 0
        for synCodon in synCodons:
            # in case we want to codon with the closest frequency as the original one
            # codonFreq = codonFrequency[codon][2]
            keep = False

            if codonFrequency:
                synCodonFreq = codonFrequency[synCodon][2]
                # keep the codon with the highest frequency
                # print(codon, synCodon, synCodonFreq, lastFreq)
                if lastFreq != 0 and synCodonFreq <= lastFreq:
                    continue
            else:
                synCodonFreq = None

            if synCodon == codon:
                continue
            # keep the codon that introduces a mutation at recodePos
            for i in range(3):
                mutPos = codonStart + i
                if mutPos in recodePos and codon[i] != synCodon[i]:
                    recodedCodons.append(codonStart)
                    keep = True
                    # print(f"Attempt {codonNb} : {codon} -> {synCodon}")
                    if codonFrequency is None:
                        mutEvents[(codon, "No codon frequency table", codonStart)] = (
                            synCodon,
                            "No codon frequency table",
                            motif,
                        )
                        break
                    else:
                        mutEvents[(codon, round(codonFreq, 2), codonStart)] = (
                            synCodon,
                            round(synCodonFreq, 2),
                            motif,
                        )
                        lastFreq = synCodonFreq
            # mutate the homolgy arm with synCodon and flag the mutation in uppercase (repeated until no codon with a higher frequency is found)
            if keep is True:
                for i in range(3):
                    mutHA[codonStart + i] = (
                        synCodon[i].upper()
                        if codon[i].upper() != synCodon[i].upper()
                        else synCodon[i].lower()
                    )
                # if no codon frequency table is availble, keep the first synonymous codon
                if codonFrequency is None:
                    break

    mutDonor = "".join(mutHA)

    # print(mutEvents)
    # store codons that were not muted
    mutList = [mutEvent[0] for mutEvent in mutEvents]
    for codonNb, codonStart in (
        enumerate(list(sorted(mutCodons)))
        if motif != "seed"
        else enumerate(reversed(list(sorted(mutCodons))))
    ):
        codon = HA[codonStart : codonStart + 3].upper()
        # TO VERIFY
        if codonNb >= maxMut:
            break
        if codon not in mutList:
            if codonFrequency:
                codonFreq = round(codonFrequency[codon][2], 2)
            else:
                codonFreq = None
            mutEvents[(codon, codonFreq, codonStart)] = (
                "could not recode",
                "NA",
                motif,
            )

    return mutDonor, mutEvents, recodedCodons


def getTargetSeq(params):
    """Used in knock-in mode : from a geneID or a sequence,
    returns the target sequence or coordinates (str),
    the insert position index on the target sequence (int, 0 based),
    and the gene model
    """

    insertpos = params.get("insertpos")
    org = params.get("org")
    geneId = params.get("koGeneId")
    targetHalfLen = (
        60  # extend N bp in 5' and 3' from the insert site to get the target sequence
    )

    if insertpos == "custom":
        # Custom sequence template (not used anymore)

        customseq, warnMsgCustomPos = cleanSeq(params.get("customseq"), org, "split")
        splitSeq = customseq.split("//")
        if len(splitSeq) == 2:
            preseq = splitSeq[0].lower()
            postseq = splitSeq[1].lower()
        else:
            raise ValueError("insertErr")

        targetSeq = preseq + postseq
        targetPos = None
        insertIdx = len(preseq)

    elif geneId:
        chrom, strand, exons = getGenePos(geneId, org, method="allExons", targetLen=0)
        geneModel = getGeneModel(exons, strand)
        start, end = exons[0] if insertpos == "Nter" else exons[-1]
        isStart = (insertpos == "Nter" and strand == "+") or (
            insertpos == "Cter" and strand == "-"
        )

        if isStart:
            insertPoint = start + 3  # after START (strand +) or before STOP (strand -)
        else:
            insertPoint = end - 3  # before STOP or after START
        targetStart, targetEnd = checkCoords(
            insertPoint - targetHalfLen, insertPoint + targetHalfLen, org, chrom
        )

        targetPos = "%s:%s-%s:%s" % (chrom, targetStart, targetEnd, strand)
        targetSeq = None
        insertIdx = targetHalfLen

    return targetSeq, targetPos, insertIdx, geneModel


def checkCoords(start, end, org, chrom):
    """clip coordinates to chromosome boundaries"""

    chromSize = parseChromSizes(org)[chrom]
    if start < 0:
        start = 0
    if end > chromSize:
        end = chromSize

    return start, end


def processDonor(DonorSeq):
    """
    flags and process the donor DNA to be "sythesis-ready"
    returns a dict with a list of tuples containing bin coordinates, GC content and free energy of 50bp bins
    inspired by protoSpaceJam (https://czbiohub-sf.github.io/protoSpaceJAM/algorithmandparameters.html)
    """

    binwidth = 20
    DonorLen = len(DonorSeq)
    bins = [
        slice(binstart, min(binstart + binwidth, DonorLen))
        for binstart in range(0, DonorLen)
    ]
    donorBins = []
    for bin in bins:
        if bin.stop - bin.start < 10:
            continue
        binseq = DonorSeq[bin]
        gc = gcContent(binseq.upper())
        donorBins.append((bin.start, bin.stop, gc))

    basesCount = {"A": 10, "T": 10, "G": 6, "C": 6}

    homopolymers = findHomopolymers(DonorSeq, basesCount)

    return donorBins, homopolymers


def iterParseBoulder(tmpOutFname):
    "parse a boulder IO style file, as output by Primer3"
    data = {}
    for line in open(tmpOutFname):
        if line == "=\n":
            if data != {}:
                yield data
                data = {}
                continue
        key, val = line.rstrip("\n").split("=")
        data[key] = val
    if data != {}:
        yield data


def runPrimer3(seqs, targetStart, targetLen, prodSizeRange, tm, addTags, primerLen):
    """
    input is a list of (seqId, seq)
    return dict seqId -> (primerseq1, tm1, pos1, primerseq2, tm2, pos2)
    addTags is a dict tag -> value
    """
    p3OutFh = makeTempFile("primer3Out", ".txt")

    minTm = tm - 3
    maxTm = tm + 3
    optTm = tm

    primerMinSize = primerLen - 5
    primerOptSize = int(primerLen)
    primerMaxSize = primerLen + 5

    p3InFh = makeTempFile("primer3In", ".txt")
    # values from https://www.ncbi.nlm.nih.gov/tools/primer-blast/
    # MAX_END_STABILITY is strange but it seems to be set that way by NCBI
    primer3ConfigDir = abspath(
        join(baseDir, "bin", "src", "primer3-2.3.6", "src", "primer3_config")
    )
    # PRIMER_MAX_POLY suggested by Yueh-Chiang.Hu@cchmc.org

    # PRIMER_MAX_POLY_X=4 suggested by Chiang.Hu@cchmc.org
    conf = (
        """PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_MIN_TM=%(minTm)d
PRIMER_OPT_TM=%(optTm)d
PRIMER_MAX_TM=%(maxTm)d
PRIMER_MAX_POLY_X=4
PRIMER_MIN_SIZE=%(primerMinSize)d
PRIMER_OPT_SIZE=%(primerOptSize)d
PRIMER_MAX_SIZE=%(primerMaxSize)d
PRIMER_MAX_END_STABILITY=9
PRIMER_MAX_POLY=4
PRIMER_PRODUCT_SIZE_RANGE=%(prodSizeRange)s
SEQUENCE_TARGET=%(targetStart)s,%(targetLen)s
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%(primer3ConfigDir)s/
"""
        % locals()
    )

    for seqId, seq in seqs:
        seqId = seqId.replace("=", " ")
        p3InFh.write("SEQUENCE_ID=%s\n" % seqId)
        p3InFh.write("SEQUENCE_TEMPLATE=%s\n" % seq)
        p3InFh.write(conf)
        for key, val in addTags.items():
            p3InFh.write("%s=%s\n" % (key, val))
        p3InFh.write("=\n")

    p3InFh.flush()

    # shutil.copy(p3InFh.name, "/tmp/primer3.tmp")

    binName = join(binDir, "primer3_core")
    cmdLine = binName + " %s > %s" % (p3InFh.name, p3OutFh.name)
    runCmd(cmdLine, ignoreExitCode=True)

    ret = {}
    # print "<br>".join(open(p3OutFh.name).read().splitlines())
    for p3 in iterParseBoulder(p3OutFh.name):
        if (
            "PRIMER_LEFT_NUM_RETURNED" not in p3
            or "PRIMER_RIGHT_NUM_RETURNED" not in p3
            or p3["PRIMER_LEFT_NUM_RETURNED"] == "0"
            or p3["PRIMER_RIGHT_NUM_RETURNED"] == "0"
        ):
            # if "PRIMER_LEFT_0_SEQUENCE" not in p3 or "PRIMER_RIGHT_0_SEQUENCE" not in p3:
            # errAbort("No primer found for this Tm. Please contact us if you think there should be a primer found.")
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
    "parse a fasta file, return list (id, seq)"
    seqs = []
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n").rstrip("\r")
        if line.startswith(">"):
            if seqId != None:
                seqs.append((seqId, "".join(parts)))
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts) != 0:
        seqs.append((seqId, "".join(parts)))
    return seqs


def parseFasta(fileObj):
    "parse a fasta file, where each seq is on a single line, return dict id -> seq"
    seqs = {}
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n").rstrip("\r")
        if line.startswith(">"):
            if seqId != None:
                seqs[seqId] = "".join(parts)
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts) != 0:
        seqs[seqId] = "".join(parts)
    return seqs


def coordsToPosStr(chrom, start, end, strand):
    "convert coords to a string"
    if chrom == None:
        return "?"
    locStr = "%s:%d-%d:%s" % (chrom, start, end, strand)
    return locStr


def findPerfectMatch(batchId, seq=None, genome=None, noPerfectMatch=None):
    """find best match for input sequence from batchId in genome and return as
    a string chrom:start-end:strand or "?" if not found "
    """
    if skipAlign:
        return "?"

    if seq is None and genome is None:
        batchInfo = readBatchAsDict(batchId)
        seq = batchInfo["seq"]
        genome = batchInfo["org"]

    # write seq to tmp file
    # tmpFaFh = tempfile.NamedTemporaryFile(dir=TEMPDIR, prefix="crisporBestMatch", suffix=".fa")
    tmpFaFh = makeTempFile(prefix="bwaswInput", suffix=".fa")
    tmpFaFh.write(">tmp\n%s" % seq)
    tmpFaFh.flush()
    logging.debug("seq: %s" % open(tmpFaFh.name).read())
    faFname = tmpFaFh.name

    # create temp SAM file
    # tmpSamFh = tempfile.NamedTemporaryFile(dir=TEMPDIR, prefix="crisporBestMatch", suffix=".sam")
    tmpSamFh = makeTempFile(prefix="bwaswOutput", suffix=".sam")
    samFname = tmpSamFh.name

    bwaIndexPath = abspath(join(genomesDir, genome, genome + ".fa"))
    remoteAddr = pipes.quote(
        os.environ.get("REMOTE_ADDR", "noIp")
    )  # so I can figure out is someone is hammering the server
    cmd = (
        "true %(batchId)s %(remoteAddr)s && $BIN/bwa bwasw -b 100 -q 100 -T 20 %(bwaIndexPath)s %(faFname)s > %(samFname)s"
        % locals()
    )
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
        logging.debug(
            "qName=%s, flag=%s, rName=%s, pos=%s, mapq=%s, cigar=%s"
            % (qName, flag, rName, pos, mapq, cigar)
        )
        if (int(flag) and 2) == 2:
            strand = "-"
        else:
            strand = "+"
        if not re.compile("[0-9]*").match(cigar):
            logging.debug("CIGAR is not number")
            continue
        if cigar == "*":
            logging.debug("CIGAR is *")
            continue
            # errAbort("Sequence not found in genome. Are you sure you have pasted the correct sequence and also selected the right genome?")
        # Todo: why do we get soft clipped sequences from BWA? repeats?
        if "S" in cigar and noPerfectMatch is None:
            logging.debug("match found, but soft-clipped, cigar: %s" % cigar)
            continue
        # allow imperfect matches
        if noPerfectMatch:
            cleanCigar = re.sub('[D/H/I/M/N/P/S/X]', "", cigar)
        else:
            cleanCigar = cigar.replace("M", "")
        if not cleanCigar.isdigit():
            logging.debug("match found, but cigar string was %s" % cigar)
            continue
        matchLen = int(cleanCigar)
        chrom, start, end = (
            rName,
            int(pos) - 1,
            int(pos) - 1 + matchLen,
        )  # SAM is 1-based
        assert (
            "|" not in chrom
        )  # We do not allow '|' in chrom name. I use this char to sep. info fields in BED.
        matchByChrom[chrom].append((chrom, start, end, strand))

    # second pass, to handle the PAR matches properly
    matches = []
    for chrom, matchList in matchByChrom.items():
        if isInPar(genome, chrom, start, end) is not None:
            # only keep matches on chrY
            if chrom == "chrX" and "chrY" in matchByChrom:
                logging.debug("In PAR region, so skipping chrX")
                continue
        for m in matchList:
            matches.append(m)

    # delete the temp files
    tmpSamFh.close()
    tmpFaFh.close()

    if len(matches) == 0:
        return "?"

    nonAltMatches = [x for x in matches if not isAltChrom(x[0])]
    if len(nonAltMatches) != 0:
        bestMatch = nonAltMatches[0]
    else:
        bestMatch = matches[0]

    logging.debug(
        "Found %d best matches, %d on non-alts. matches: %s, best match %s"
        % (len(matches), len(nonAltMatches), matches, bestMatch)
    )
    return "%s:%d-%d:%s" % (bestMatch)


def maskLowercase(seq):
    "replace all lowercase letters with 'N'"
    newSeq = []
    for c in seq:
        if c.islower():
            newSeq.append("N")
        else:
            newSeq.append(c)
    return "".join(newSeq)


def getGenomeSeqsBin(genome, coordList, doRepeatMask=False):
    """use the binary twoBitToFa to get seqs:
    coordList has format (chrom, start, end, name, score, strand)
    returns list (chrom, start, end, name, seq)
    """
    twoBitPath = join(genomesDir, "%(genome)s/%(genome)s.2bit" % locals())
    twoBitPath = abspath(twoBitPath)

    bedFh = makeTempFile("getGenomeSeqs", ".bed")
    # bedFh = open("/tmp/test.bed", "w")
    faFh = makeTempFile("getGenomeSeqs", ".fa")

    for row in coordList:
        row = list(row)
        row[3] = row[3].replace(
            " ", "_"
        )  # some weird NCBI assemblies have spaces in the chrom names, just hack around it
        row = [str(x) for x in row]
        line = "\t".join([str(x) for x in row])
        bedFh.write(line)
        bedFh.write("\n")

    bedFh.flush()

    cmd = ["$BIN/twoBitToFa", "-bed=" + bedFh.name, twoBitPath, faFh.name]
    runCmd(cmd, useShell=False)
    bedFh.close()  # delete temp file

    faSeqs = parseFastaAsList(open(faFh.name))

    assert len(faSeqs) == len(coordList)

    seqs = []
    for coordTuple, seqTuple in zip(coordList, faSeqs):
        seqId, seq = seqTuple
        if len(coordTuple) == 4:
            chrom, start, end, name = coordTuple
            strand = "+"
            score = "0"
        else:
            chrom, start, end, name, score, strand = coordTuple

        # seq = tbf[chrom][start:end]
        if strand == "-":
            seq = revComp(seq)
        if doRepeatMask:
            seq = maskLowercase(seq)
        seqs.append((chrom, start, end, name, score, strand, seq))
    return seqs


def getGenomeSeqs(genome, coordList, doRepeatMask=False):
    """DOES NOT WORK ON PYTHON3 ANYMORE - NOT USED ANYMORE, as long twobitreader has not been updated

    return dict of genome sequences,
    coordList has format (chrom, start, end, name, score, strand)
    returns list (chrom, start, end, name, seq)
    """
    twoBitPath = join(genomesDir, "%(genome)s/%(genome)s.2bit" % locals())
    twoBitPath = abspath(twoBitPath)

    tbf = twobitreader.TwoBitFile(twoBitPath)
    seqs = []
    for coordTuple in coordList:
        if len(coordTuple) == 4:
            chrom, start, end, name = coordTuple
            strand = "+"
            score = "0"
        else:
            chrom, start, end, name, score, strand = coordTuple

        seq = tbf[chrom][start:end]
        if strand == "-":
            seq = revComp(seq)
        if doRepeatMask:
            seq = maskLowercase(seq)
        seqs.append((chrom, start, end, name, score, strand, seq))
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
    hdrDist=None,
    primerLen=20,
):
    "create primer for region around chrom:start-end, write output to batch"
    " returns (leftPrimerSeq, lTm, lPos, rightPrimerSeq, rTm, rPos, amplified sequence)"

    # get 1kbp of flanking sequence
    flankSeq = getFlankSeq(genome, chrom, start, end)

    targetStart, targetLen, ampRange, addTags = getTargetForPrimerDesign(
        guideStart, ampLen, hdrDist, strand
    )

    primers = runPrimer3(
        [("seq1", flankSeq)], targetStart, targetLen, ampRange, tm, addTags, primerLen
    )
    lSeq, lTm, lPos, rSeq, rTm, rPos = list(primers.values())[0]

    if lSeq == None or rSeq == None:
        return None, None, None, None, None, None, flankSeq, ampRange, flankSeq, addTags

    targetSeq = flankSeq[lPos : rPos + 1]
    if "N" in targetSeq:
        # Get the flankSeq again but without the Ns... see maskLowercase
        flankSeq = getFlankSeq(genome, chrom, start, end, doRepeatMask=False)
        targetSeq = flankSeq[lPos : rPos + 1]

    return lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq, ampRange, flankSeq, addTags


def getFlankSeq(genome, chrom, start, end, doRepeatMask=True):
    flankStart = start - 1000
    flankEnd = end + 1000

    chromSizes = parseChromSizes(genome)

    if flankStart < 0 or flankEnd > chromSizes[chrom]:
        errAbort(
            "Not enough space on genome sequence to design primer. Need at least 1kbp on each side of the input sequence to design primers. Please design primers manually, choose a more recent genome assembly with longer contig sequences or paste a shorter input sequence (e.g. just the guide sequence alone with the PAM). Still questions? Email %s"
            % contactEmail
        )

    # get 1kbp of flanking sequence
    flankSeq = getGenomeSeqsBin(
        genome, [(chrom, flankStart, flankEnd, "seq")], doRepeatMask=doRepeatMask
    )[0][-1]
    return flankSeq


def getTargetForPrimerDesign(guideStart, ampLen, hdrDist, strand):
    if hdrDist is not None:
        # Cut for Cas9 is always between 3 and 4 nt from the PAM
        if strand == "+":  # '+' here is guide strand relative to the genome
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
        ampRange = "250-310"
    else:
        # try to get some good heuristics for the primer placement
        # primers must not overlap the target but also not be
        # too far away, especially when using next-gen sequencing
        targetStart = 1000 + guideStart
        targetLen = 23
        if ampLen <= 250:
            ampDist = 50
        else:
            ampDist = 150
        ampRange = str(ampLen - ampDist) + "-" + str((ampLen))

    # for long reads: make sure that the primers are not too close to the target
    addTags = {}
    if ampLen >= 600:
        addTags = {
            "SEQUENCE_EXCLUDED_REGION": "%d,%d" % (targetStart - 150, targetLen + 300),
            # "SEQUENCE_EXCLUDED_REGION" : "0,1500"
        }

    if hdrDist is not None:
        # Primer products should be as short as possible while staying at or above 250bp.
        # See also prodSizeRange and PRIMER_PRODUCT_SIZE_RANGE.
        addTags.update(
            {
                # The optimum size for the PCR product.
                "PRIMER_PRODUCT_OPT_SIZE": ampRange.split("-")[0],  # 250
                # Penalty weight for products longer than PRIMER_PRODUCT_OPT_SIZE.
                "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": 1,
            }
        )

    return targetStart, targetLen, ampRange, addTags


def markupSeq(seq, ulPosList, boldPosList, annots={}):
    """return seq as html with some parts underlined or in bold.
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
    openAnnots = defaultdict(int)  # current number of open spans, per cssString
    openTags = set()
    for i, nucl in enumerate(seq):
        if i in annotEnds:
            for tagStr in annotEnds[i]:
                if tagStr in openAnnots:
                    openAnnots[tagStr] -= 1
                    if openAnnots[tagStr] == 0:
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
            openAnnots[tagStr] += 1
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
        if (i + 1) % 80 == 0:
            ret.append("<br>")
    for tag in openTags:
        ret.append("</%s>" % tag)
    return "".join(ret)
    # return seq[:start]+"<u>"+seq[start:end]+"</u>"+seq[end:]


def makeHelperPrimers(guideName, guideSeq, plasmid, pam):
    "return dict with various names -> primer for primer page"
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
    # specPrimer = "TAATACGACTCACTATA%s<b>%s</b>GTTTTAGAGCTAGAAATAGCAAG" % (prefix, guideSeq)

    if pamIsSaCas9(pam):
        # See http://dx.doi.org/10.1016/j.cell.2015.08.007
        specPrimer = "GAAATTAATACGACTCACTATA%s<b>%s</b>GTTTTAGTACTCTGGAAACAGAA" % (
            prefix,
            guideSeq,
        )
        tracrPrimer = "GTTTTAGTACTCTGGAAACAGAATCTACTAAAACAAGGCAAAATGCCGTGTTTATCTCGTCAACTTGTTGGCGAGATTTTTTT"
    else:
        specPrimer = "GAAATTAATACGACTCACTATA%s<b>%s</b>GTTTTAGAGCTAGAAATAGCAAG" % (
            prefix,
            guideSeq,
        )
        tracrPrimer = "AAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC"

    primers["T7iv"].append(("guideRNA%sT7crTarget" % guideName, specPrimer))
    primers["T7iv"].append(
        ("guideRNAallT7common (constant primer used for all guide RNAs)", tracrPrimer)
    )

    # geneart primers
    primers["geneArt"].append(
        ("guideRNA%sGeneArtFw" % guideName, "TACGACTCACTATAG<b>" + guideSeq + "</b>")
    )
    primers["geneArt"].append(
        (
            "guideRNA%sGeneArtRev" % guideName,
            "TTCTAGCTCTAAAAC<b>" + revComp(guideSeq) + "</b>",
        )
    )

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
            if pi[0] == plasmid:
                plasmidLabel = pi[1][1]
        primers["mammCells"].append(
            (
                fwName + plasmidLabel,
                "%s%s<b>%s</b>%s" % (u6FwPrefix, addGPrefix, guideSeq, u6Suffix),
            )
        )
        primers["mammCells"].append(
            (
                revName + plasmidLabel,
                "%s<b>%s</b>%s%s"
                % (u6RwPrefix, revComp(guideSeq), addCSuffix, u6Suffix),
            )
        )

        if not guideSeq.lower().startswith("g"):
            primers["mammCells19"].append(
                (
                    "gN19-" + fwName + plasmidLabel,
                    "%s%s<b>%s</b>%s"
                    % (u6FwPrefix, addGPrefix, guideSeqTrunc, u6Suffix),
                )
            )
            primers["mammCells19"].append(
                (
                    "gN19-" + revName + plasmidLabel,
                    "%s<b>%s</b>%s%s"
                    % (u6RwPrefix, revComp(guideSeqTrunc), addCSuffix, u6Suffix),
                )
            )

        # if guideSeq.lower().startswith("g"):
        # primers["mammCells"].append((fwName, "ACACC<b>%s</b>G" % guideSeq))
        # primers["mammCells"].append((revName, "AAAAC<b>%s</b>G" % revComp(guideSeq)))
        # else:
        # primers["mammCells"].append((fwName, "ACACC<u>G</u><b>%s</b>G" % guideSeq))
        # primers["mammCells"].append((revName, "AAAAC<b>%s</b><u>C</u>G" % revComp(guideSeq)))
        # primers["mammCellsNote"] = True

    # add the prefix to all names
    newPrimers = {}
    for key, primList in primers.items():
        if key.endswith("Note"):
            newPrimers[key] = primList
            continue

        newList = []
        for name, seq in primList:
            if batchName != "":
                newName = batchName + "_" + name
            else:
                newName = name
            newList.append((newName, seq))
        newPrimers[key] = newList
    return newPrimers


def printPrimerTableAll(primers):
    print('<table class="primerTable">')
    for key, primerList in primers.items():
        if key.endswith("Note"):
            continue
        for name, seq in primerList:
            name = name.split()[0]
            print("<tr>")
            print("<td>%s</td>" % name)
            print("<td><tt>%s</tt></td>" % seq)
            print("</tr>")
    print("</table>")


def printPrimerTable(
    primerList, onRows=False, withTm=False, withScore=False, seqType="Primer"
):
    "given a list of (name, seq) tuples, print a table"
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
            print("<tr>")
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


# def printDropboxLink():
# print """
# <a id="dropbox-link"></a>
# <script type="text/javascript" src="https://www.dropbox.com/static/api/2/dropins.js" id="dropboxjs" data-app-key="lp7njij87kjrcxv"></script>");
# <script type="text/javascript">
# options = {
# success: function(files) {
# $('textarea[name=hgct_customText]').val(files[0].link);
# alert('Here is the file link: ' + files[0].link);
# },
# cancel: function() {},
# linkType: 'direct',
# multiselect: true,
# };
# var button = Dropbox.createChooseButton(options);
# document.getElementById('dropbox-link').appendChild(button);
# </script>
# """


def mergeParamDicts(params, changeParams):
    """changeParams is a dict that can override elements in params.
    if value==None in changeParams, the whole element will get removed.
    if onlyParams is set, only copy over the keys in onlyParams (a list)
    """
    newParams = {}
    newParams.update(params)
    newParams.update(changeParams)
    for key, val in changeParams.items():
        if val == None:
            del newParams[key]
    return newParams


def printHiddenFields(params, changeParams, form=None, excludeParams=None):
    " "

    formStr = ""
    if form:
        formStr = """ form="%s" """ % form

    if excludeParams and isinstance(excludeParams, list):
        for param in excludeParams:
            params.pop(param, None)

    newParams = mergeParamDicts(params, changeParams)
    for key, val in newParams.items():
        print(('<input type="hidden" name="%s" value="%s" %s>' % (key, val, formStr)))


def cgiGetSelfUrl(changeParams, anchor=None, onlyParams=None):
    """create a URL to myself with different parameters than what we have.
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
    # urllib crashes if val is not a string
    paramStrs = [
        "%s=%s" % (key, urllib.parse.quote(val))
        for key, val in newParams.items()
        if isinstance(val, str)
    ]
    paramStr = "&".join(paramStrs)
    url = basename(__file__) + "?" + paramStr
    if anchor is not None:
        url += "#" + anchor
    return url


def printDropDown(name, nameValList, default, onChange=None, style=None, form=None):
    """print a dropdown box and set a default"""
    addStr = ""
    if onChange is not None:
        addStr = """ onchange="%s" """ % onChange
    addStr2 = ""
    if style is not None:
        addStr2 = """ style='%s' """ % style
    addStr3 = ""
    if form is not None:
        addStr3 = """ form='%s' """ % form

    print(('<select id="dropdown" name="%s"%s%s%s>' % (name, addStr, addStr2, addStr3)))
    for name, desc in nameValList:
        name = str(name)
        addString = ""
        if default is not None and str(name) == str(default):
            addString = ' selected="selected"'
        print(('<option value="%s"%s>%s</option>' % (name, addString, desc)))
    print("</select>")


def findGuideSeq(inSeq, pam, pamId, exonId=None, pamFullName=None):
    """given the input sequence and the pamId, return the guide sequence,
    the sequence with the pam and its strand.
    """
    startDict, endSet = findAllPams(inSeq, pam)
    pamInfo = list(
        flankSeqIter(
            inSeq, startDict, len(pam), False, exonId=exonId, pamFullName=pamFullName
        )
    )
    for (
        guidePamId,
        pamStart,
        guideStart,
        guideStrand,
        guideSeq,
        pamSeq,
        pamPlusSeq,
    ) in pamInfo:
        if guidePamId != pamId:
            continue
        guideSeqWPam = concatGuideAndPam(guideSeq, pamSeq)
        # prettify guideSeqWPam to highlight the PAM
        if pamIsFirst:
            guideSeqHtml = "<i>%s</i> %s" % (
                guideSeqWPam[: len(pam)].upper(),
                guideSeqWPam[len(pam) :].upper(),
            )
        else:
            guideSeqHtml = "%s <i>%s</i>" % (
                guideSeqWPam[: -len(pam)].upper(),
                guideSeqWPam[-len(pam) :].upper(),
            )
        guideEnd = guideStart + GUIDELEN
        return (
            guideSeq,
            pamSeq,
            pamPlusSeq,
            guideSeqWPam,
            guideStrand,
            guideSeqHtml,
            guideStart,
            guideEnd,
        )
    errAbort("pamId %s not found? This is a bug." % pamId)


def findOntargetPos(otMatches, pamId, position, absentOk=False):
    "find position of guide sequence in genome at MM0"
    if pamId not in otMatches or 0 not in otMatches[pamId]:
        if absentOk:
            return None, None, None, None, None, False
        else:
            errAbort(
                "No perfect match found for guide sequence in the genome. Cannot design primers for a non-matching guide sequence.<p>Are you sure you have selected the right genome? <p> If you have selected the right genome and entered a cDNA as the query sequence, please note that sequences that overlap a splice site are not part of the genome and cannot be used as guide sequences."
            )

    matchList = otMatches[pamId][0]  # = get all matches with 0 mismatches
    isUnique = True
    if len(matchList) > 1:
        # this only gets executed if there are multiple exact matches for the target
        # we iterate over all offs and try find the correct one.
        targetChrom, targetStart, targetEnd, strand = parsePos(position)

        filtMatch = None
        # search for off-target that is the on-target
        for match in matchList:
            # example match: ('scaffold_1', 578, 601, 'TATTGGATTGGTCCAATCGTTGG', '-', 'ex', 'GSADVT00000001001', 293)
            chrom, start, end = match[:3]
            if chrom == targetChrom and start >= targetStart and end <= targetEnd:
                filtMatch = match
                break

        if filtMatch is None:
            errAbort(
                "Multiple matches for this guide, but no single match is within the target sequence? This can happen when your input is not part of the genome, in which case bulk-primer design makes no sense, as far as we can tell. If you think this would make sense for you, contact us at %s."
                % contactEmail
            )

        chrom, start, end = filtMatch[:3]
        isUnique = False

        matchList = [filtMatch]

    global batchName
    batchName = batchName.replace(" ", "_")

    chrom, start, end, seq, strand, segType, segName, alnCount, hasXa = matchList[0]
    start = int(start)
    end = int(end)
    return chrom, start, end, strand, segName, isUnique


def printAmpLenAndTm(ampLen, primerLen, tm):
    "print form fields for amplicon length and TM"
    print("Maximum amplicon length:")
    dropDownSizes = [
        ("100", "100 bp - for >= 75bp paired reads"),
        ("150", "150 bp - for >= 100bp paired reads "),
        ("200", "200 bp - for >= 125bp paired reads"),
        ("300", "300 bp - for >= 200bp paired reads"),
        ("400", "400 bp - for >= 250bp paired reads"),
        ("500", "500 bp - for >= 300bp paired reads"),
        ("600", "600 bp - for Sanger reads"),
        ("800", "800 bp - for Sanger reads"),
    ]

    printDropDown(
        "ampLen", dropDownSizes, ampLen, onChange="""$('#submitPcrForm').click()"""
    )

    print("&nbsp;&nbsp;&nbsp; Primer length:")
    primerLenList = [
        ("15", "15 bases"),
        ("16", "16 bases"),
        ("17", "17 bases"),
        ("18", "18 bases"),
        ("19", "19 bases"),
        ("20", "20 bases"),
        ("21", "21 bases"),
        ("22", "22 bases"),
        ("23", "23 bases"),
        ("24", "24 bases"),
        ("25", "25 bases"),
        ("26", "26 bases"),
        ("27", "27 bases"),
        ("28", "28 bases"),
        ("29", "29 bases"),
        ("30", "30 bases"),
    ]
    printDropDown(
        "primerLen",
        primerLenList,
        primerLen,
        onChange="""$('#submitPcrForm').click()""",
    )

    print("&nbsp;&nbsp;&nbsp; Primer Tm:")
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


def printValidationPcrSection(
    batchId,
    genome,
    pamId,
    position,
    params,
    guideStart,
    guideEnd,
    primerGuideName,
    guideSeq,
):
    "print the PCR section of the primer page"

    # check the input parameters: ampLen, tm
    ampLen = params.get("ampLen", "400")
    if not ampLen.isdigit():
        errAbort("ampLen parameter must be a number")
    ampLen = int(ampLen)

    primerLen = params.get("primerLen", "22")
    if not primerLen.isdigit():
        errAbort("primerLen parameter must be a number")
    primerLen = int(primerLen)

    tm = params.get("tm", "64")
    if not tm.isdigit():
        errAbort("tm parameter must be a number")
    tm = int(tm)

    # Optional, for HDR knock-in experiments
    hdrDist = params.get("hdrDist", None)
    if hdrDist == "":
        hdrDist = None
    if hdrDist is not None:
        try:
            hdrDist = int(hdrDist)
        except (ValueError, TypeError):
            errAbort("hdrDist parameter must be a number")

    print("<h2 id='ontargetPcr'>PCR to amplify the on-target site</h2>")

    otMatches = parseOfftargets(genome, batchId)
    chrom, start, end, strand, gene, isUnique = findOntargetPos(
        otMatches, pamId, position, absentOk=True
    )
    if chrom is None:
        print("Not found in genome, cannot design primer")
        return None, None, None

    if not isUnique:
        print(
            "<div class='substep' style='border: 1px black solid; padding:5px; background-color: aliceblue'>"
        )
        print(
            (
                "<strong>Warning</strong>: Found multiple perfect matches for this guide sequence in the genome. For the PCR, we are using the on-target match in the input sequence %s:%d-%d (gene: %s), but this guide will not be specific. Is this a polyploid organism? Try selecting another guide sequence or email %s to discuss your strategy or modifications to this software.<p>"
                % (chrom, start + 1, end, gene, contactEmail)
            )
        )
        print("</div>")

    lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq, ampRange, flankSeq, addTags = (
        designPrimer(
            genome,
            chrom,
            start,
            end,
            strand,
            0,
            batchId,
            ampLen,
            tm,
            hdrDist,
            primerLen,
        )
    )

    primerPosList = []
    if lSeq is not None:
        primerPosList.append((0, len(lSeq)))
        primerPosList.append(((len(targetSeq) - len(rSeq)), len(targetSeq)))

    guideStartOnTarget = None
    guideEndOnTarget = None
    targetHtml = ""
    if lPos != None:
        guideStartOnTarget = 1000 - lPos
        guideEndOnTarget = guideStartOnTarget + GUIDELEN + PAMLEN
        annots = defaultdict(dict)
        annots[(guideStartOnTarget, guideEndOnTarget)]["css"] = {
            "font-weight": "bold",
            "background-color": "yellow",
        }
        targetHtml = markupSeq(targetSeq, primerPosList, [], annots)

    allPrimersFound = True

    if batchName == "":
        primerPrefix = ""
    else:
        primerPrefix = batchName + "_"

    print(
        "Use these primers to amplify a genomic fragment around the on-target site:<br>"
    )

    print('<table class="primerTable">')
    print("<tr>")
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
    print("</tr></table><p>")

    print(
        "<h3>Genome fragment with validation primers (underlined) and guide sequence (yellow)</h3>"
    )

    print(("""<form id="ampLenForm" action="%s" method="GET">""" % basename(__file__)))

    printAmpLenAndTm(ampLen, primerLen, tm)

    printHiddenFields(params, {"ampLen": None, "tm": None, "primerLen": None})

    print(
        """<input id="submitPcrForm" style="display:none" type="submit" name="submit" value="submit">"""
    )
    print("</form><p>")

    if strand == "-":
        print(
            "Your guide sequence is on the reverse strand relative to the genome sequence, so it is reverse complemented in the sequence below.<p>"
        )

    if not chrom.startswith("ch"):
        chromLong = "chr" + chrom
    else:
        chromLong = chrom

    print("""<div style='word-wrap: break-word; word-break: break-all;'>""")
    if allPrimersFound:
        print(
            "<strong>Genomic sequence %s:%d-%d including primers, genomic forward strand:</strong>"
            % (chromLong, start, end)
        )
        print("<br><tt>%s</tt><br>" % (targetHtml))
    else:
        print(
            "<strong>Warning: No primers were found at this Tm, please design them manually e.g. with NCBI PrimerBlast.</strong><br>"
        )
        print("<br><tt>%s</tt><br>" % (targetHtml))

    print("""</div>""")
    if rPos is not None:
        print("<strong>Sequence length:</strong> %d<p>" % (rPos - lPos))

    # primer3 may have used some special tags, add them to the description as a string
    p3ArgStr = ""
    if len(addTags) != 0:
        primer3Tags = []
        for key, val in addTags.items():
            primer3Tags.append("%s=%s" % (key, val))
        p3ArgStr = ", ".join(primer3Tags)

    print(
        "<small>Method: Primer3.2 with default settings, target length %s bp, %s</small>"
        % (ampRange, p3ArgStr)
    )

    return targetSeq, guideStartOnTarget, guideEndOnTarget


def printEnzymeSection(
    mutEnzymes, targetSeq, guideSeqWPam, guideStartOnTarget, guideEndOnTarget
):
    "print the section about restriction enzymes in the target seq"
    print("<h2 id='restrSites'>Restriction Sites for PCR product validation</h2>")

    print("Cas9 induces mutations, usually 3bp 5' of the PAM site.")
    print(
        "If a mutation is induced, then it is very likely that one of the following enzymes no longer cuts your PCR product amplified from the mutant sequence."
    )
    print(
        "For each restriction enzyme, the guide sequence with the restriction site underlined is shown below.<p>"
    )

    print("<table>")
    print("<tr>")
    print(
        "<th>Enzyme</th><th>Pattern</th><th>Guide with Restriction Site</th><th>Suppliers</th>"
    )
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
            annots[(start, end)]["css"] = {"font-weight": "bold"}
            guideHtmls.append(markupSeq(guideSeqWPam.upper(), [], [], annots))
            # print guideSeqWPam
            # if strand=="-":
            # allSitePos.add( (guideEnd-end, guideEnd-start) )
            # else:
            allSitePos.add((guideStartOnTarget, guideEndOnTarget))

        print(", ".join(guideHtmls))
        print("</tt></td>")

        supplNames = [rebaseSuppliers.get(x, x) for x in suppliers]
        print("<td>%s</td>" % ", ".join(sorted(supplNames)))
        print("</tr>")
        # print "<br>"
    print("</table>")

    print("<h3>All restriction enzyme sites on the amplicon sequence</h3>")
    print(
        "Restriction sites are shown in yellow, the guide sequence is highlighted in bold. Use this schema to check if the sites are unique enough to give separate bands on a gel:<p>"
    )

    for enzName, pat in patList:
        annots = defaultdict(dict)
        fragLens = []
        lastEnd = 0
        for pos in findPat(targetSeq, pat):
            start = pos
            end = pos + len(pat)
            annots[(start, end)]["css"] = {"background-color": "yellow"}

            fragLens.append(str(start - lastEnd) + "bp")
            lastEnd = end
        # print markupSeq(targetSeq, [], [(guideStart, guideEnd)], annots)
        # annots[(guideStart, guideEnd)]["css"] = {"font-weight":"bold"}
        fragLens.append(str(len(targetSeq) - lastEnd) + "bp")

        print(
            (
                "<div style='margin-top: 6px'>Enzyme: <strong>%s</strong>, Site: %s, Restriction fragment lengths: %s<br>"
                % (enzName, pat, ", ".join(fragLens))
            )
        )
        print("<tt>")
        print(
            markupSeq(targetSeq, [], [(guideStartOnTarget, guideEndOnTarget)], annots)
        )
        print("</tt><br></div>")


def printCloningSection(batchId, primerGuideName, guideSeq, params, pam):
    "print the cloning/expression section of the primer page"
    print("<h2 id='cloning'>Cloning and expression of guide RNA</h2>")

    plasmid = params.get("plasmid", addGenePlasmids[0][0])
    plasmidToName = dict(addGenePlasmids)
    if plasmid not in plasmidToName:
        errAbort("Invalid value for the parameter 'plasmid'")

    primers = makeHelperPrimers(primerGuideName, guideSeq, plasmid, pam)

    # T7 from plasmids
    if not pamIsCpf1(pam):
        print("<h3 id='t7plasmid'>T7 <i>in vitro</i> expression from a plasmid</h3>")
        print(
            'To produce guide RNA by in vitro transcription with T7 RNA polymerase, the guide RNA sequence can be cloned into a variety of plasmids (see <a href="http://addgene.org/crispr/empty-grna-vectors/">AddGene website</a>).<br>'
        )

        print(
            "For the guide sequence %s, the following primers should be ordered for cloning into the BsaI-digested plasmid <a href='https://www.addgene.org/42250/'>DR274</a> generated by the Joung lab.<p>"
            % guideSeq
        )

    printPrimerTable(primers["T7"])

    # T7 from primers, in vitro
    print(
        "<h3 id='t7oligo'>T7 <i>in vitro</i> expression from overlapping oligonucleotides</h3>"
    )

    primerType = "spCas9"
    if pamIsSaCas9(pam):
        primerType = "saCas9"

    if not pamIsSaCas9(pam) and not pamIsSpCas9(pam):
        printWarning(
            "The T7 primer information below only applies to %s enzymes. Make sure that you understand the importance of the scaffold sequence in the primers below (both primers depend on the scaffold sequence) and how it is specific to the tracrRNA sequence of your enzyme. Contact us if in doubt, we can quickly adapt the primer for this particular enzyme. "
            % primerType
        )

    print(
        "For %s, template for <i>in vitro</i> synthesis of guide RNA with T7 RNA polymerase can be prepared by annealing and primer extension of the following primers:<p>"
        % primerType
    )

    print(
        "<p>Note: this sequence should only be used to synthetize guides for spCas9. The tracr sequence will vary if you are using a different enzyme.</p>"
    )

    printPrimerTable(primers["T7iv"])

    print(
        """T7 RNA polymerase starts transcription most efficiently if the first two nucleotides to be transcribed are GG. A common recommendation is to add the prefix GG- if our guide does not start with G (5'-N20-(NGG)-3'), to add G- if your guide starts with a single G (5'-GN19-(NGG)-3') and to not add anything if your guide starts with GG already (5'-GGN18-(NGG)-3').<p>"""
    )

    print(
        'One protocol for template preparation from oligonucleotides and in-vitro transcription can be found in <a href="http://www.cell.com/cell-reports/abstract/S2211-1247(13)00312-4">Bassett et al. Cell Rep 2013</a>. We also provide our own <a href="downloads/prot/sgRnaSynthProtocol.pdf">optimized protocol</a> for T7 guide expression.<p>'
    )

    print(
        """<a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4038517/?report=classic">Gagnon et al. PLoS ONE 2014</a> prefixed guides with GG to ensure high efficiency in vitro transcription by T7 RNA polymerase. It has been shown by other authors that the 5' nucleotides of the guide have little or no role in target specificity and it is therefore generally accepted that prefixing guides with GG should not affect activity.<br>"""
    )

    print(
        'However, in our lab, we found that in vitro transcription with T7 RNA polymerase is efficient enough when the sequence starts with a single G rather than with GG. This took some optimization of the reaction conditions including using large amounts of template DNA and running reactions overnight. <a href="downloads/prot/sgRnaSynthProtocol.pdf">Click here</a> to download our optimized protocol for T7 guide expression.<br>'
    )

    print(
        "Do not use G-prefixing with high-fidelity Cas9 Variants like HF1 and eSpCas9 1.1 when this adds a mismatch in the genome as the efficiency will most likely be very low."
    )

    # T7 from primers, in vitro
    print("<h3 id='geneArt'>T7 <i>in vitro</i> expression with the GeneArt kit</h3>")
    print("Use these two primers for the Invitrogen GeneArt kit:<p>")

    printPrimerTable(primers["geneArt"])

    # MAMMALIAN CELLS
    print("<h3 id='u6plasmid'>U6 expression from an Addgene plasmid</h3>")
    if "tttt" in guideSeq.lower():
        print(
            "The guide sequence %s contains the motif TTTT, which terminates RNA polymerase. This guide sequence cannot be transcribed in mammalian cells."
            % guideSeq
        )
    else:
        print(
            "The guide sequence %s does not contain the motif TTTT, which terminates RNA polymerase, so it can be transcribed in mammalian cells."
            % guideSeq
        )

        print("<br>")
        if not pamIsCpf1(pam):
            print(
                (
                    """<p><form style="margin-bottom: 0px" id="plasmidForm" action="%s#u6plasmid" method="GET">"""
                    % basename(__file__)
                )
            )
            print("Select your Addgene plasmid: ")

            # we need a separate form here (not PCR form), as the target anchor
            # to jump to after a submit is different
            plasmidNames = [(x, y) for x, (y, z) in addGenePlasmids]
            printDropDown(
                "plasmid",
                plasmidNames,
                plasmid,
                onChange="""$('#submitPlasmidForm').click()""",
            )
            printHiddenFields(cgiParams, {"plasmid": None, "submit": None})
            print(
                """<input id="submitPlasmidForm" style="display:none" type="submit" name="submit" value="submit">"""
            )
            print("""</form></p>""")

            print(
                (
                    "<p>To clone the guide into <i><a href='https://www.addgene.org/%s/'>%s</a></i>, use these primers:"
                    % (plasmid, plasmidToName[plasmid][0])
                )
            )
        else:
            print(
                "To express guide RNA for Cpf1 in mammalian cells, two plasmids are available. To clone the guide RNA sequence into the plasmids <a href='https://www.addgene.org/78956/'>pU6-As-crRNA</a> or <a href='https://www.addgene.org/78957/'>pU6-Lb-crRNA</a>, where guide RNA expression is driven by a human U6 promoter, the following primers should be used :"
            )

        print("<br>")
        if "mammCellsNote" in primers:
            print(
                "<p><strong>Note:</strong> Efficient transcription from the U6 promoter requires a 5' G. This is not the case for this guide. Several options are possible, you can either add an additional G- prefix to the N20 guide sequence, called  gN20 guides here, or replace the first with a G and create a gN19 guide. For users of HF1 and eSpCas9: G- prefixing with the high-fidelity variants may reduce efficiency, as it introduces a mismatch.</p>"
            )

            print("<strong>Primers for gN20 guides:</strong>")
            printPrimerTable(primers["mammCells"])

            print("<p><strong>Primers for gN19 guides:</strong><br>")
            print(
                "<a href='https://www.nature.com/articles/s41551-019-0505-1'>Kim et al 2020</a>. showed that changing "
                "the first nucleotide to 'G' is slightly more efficient.</p>"
            )

            printPrimerTable(primers["mammCells19"])
        else:
            printPrimerTable(primers["mammCells"])

        _, _, _, enzyme, protoUrl = addGenePlasmidInfo[plasmid]
        print(("<p>The plasmid has to be digested with: <i>%s</i><br>" % enzyme))
        print(
            (
                "<a href='%s'>Click here</a> to download the cloning protocol for <i>%s</i>"
                % (protoUrl, plasmidToName[plasmid][0])
            )
        )

    if pamIsCpf1(pam):
        print("<h3 id='ciona'>Direct PCR for <i>C. intestinalis</i></h3>")
        print(
            """A method only used at the moment by <i>Ciona intestinalis</i> (alias <i>Ciona robusta</i>) labs. The DNA construct is assembled during the PCR reaction; expression cassettes are generated with One-Step Overlap PCR (OSO-PCR) <a href="http://www.sciencedirect.com/science/article/pii/S0012160616306455">Gandhi et al., Dev Bio 2016</a> (<a href="http://biorxiv.org/content/early/2017/01/01/041632">preprint</a>) following <a href="downloads/prot/cionaProtocol.pdf">this protocol</a>. The resulting unpurified PCR product can be directly electroporated or injected into Ciona eggs.<br>"""
        )
        if batchName != "":
            primerStart = batchName
        else:
            primerStart = "sg"
        ciPrimers = [
            (
                batchName + ".%s.sgF" % primerGuideName,
                "g<b>" + guideSeq[1:] + "</b>gtttaagagctatgctggaaacag",
            ),
            (
                batchName + ".%s.U6R" % primerGuideName,
                "<b>" + revComp(guideSeq[1:]) + "</b>catctataccatcggatgccttc",
            ),
        ]
        printPrimerTable(ciPrimers)

    print("<h3 id='gibson'>Lentiviral vectors: cloning with Gibson assembly</h3>")
    print(
        """Order the following oligonucleotide to clone with Gibson assembly into the vector <a href='https://www.addgene.org/52963/'>pLentiGuide-puro</a>. See the <a href="https://www.nature.com/articles/nprot.2018.005">protocol by Matt Canver</a>.<br>"""
    )
    print(
        """To clone with restriction enzymes into this vector, see the section <a href="#u6plasmid">U6 expression from an AddGene plasmid</a> and choose pLentiGuide-puro from the list of AddGene plasmids.<br>"""
    )

    satMutUrl = cgiGetSelfUrl({"satMut": "1"}, onlyParams=["batchId"])
    print(
        (
            """If you use lentiviral vectors, you may be interested in our tools for <a href="%s">saturating mutagenesis</a>"""
            """ and for <a href="crispor.py?libDesign=1">gene knockout libraries</a>."""
            % satMutUrl
        )
    )

    lentiPrimers = [
        (
            "batchOligo%s" % primerGuideName,
            "<i>GGAAAGGACGAAACACCG</i>"
            + guideSeq
            + "<i>GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC</i>",
        ),
    ]

    printPrimerTable(lentiPrimers, seqType="Oligonucleotide")

    print("<h3 id='primerSummary'>Summary of main cloning/expression primers</h3>")
    printPrimerTableAll(primers)


def primerDetailsPage(params):
    """create primers with primer3 around site identified by pamId in batch
    with batchId. Output primers as html
    """
    # retrieve batch information
    batchId, pamId, pam = params["batchId"], params["pamId"], params["pam"]

    pam = setupPamInfo(pam)

    (
        inSeq,
        genome,
        pamSeq,
        position,
        extSeq,
        multiseq,
        koMethod,
        geneModel,
        koGeneId,
        multipam,
    ) = readBatchParams(batchId)
    batchInfo = readBatchAsDict(batchId)
    exonSeqs = batchInfo.get("exonSeqs")

    # get the pam in multipam mode
    if multipam:
        pamFullName = pamId.split(".")[0]
        exonId = None
        pam = setupPamInfo(pamFullName)
    # get the exonID and its sequence if in multiseq mode
    elif multiseq:
        pamFullName = None
        exonId = int(pamId.split(".")[0])
        batchInfo = readBatchAsDict(batchId)
        exonSeqs = batchInfo["exonSeqs"]
        for (exonIdx, exonSeq), (posStrIdx, exonPosStr) in zip(exonSeqs, multiseq):
            if exonId == exonIdx and exonId == posStrIdx:
                inSeq = exonSeq
                position = exonPosStr
                break
    else:
        exonId = None
        pamFullName = None

    seqLen = len(inSeq)
    batchBase = join(batchDir, batchId)

    (
        guideSeq,
        pamSeq,
        pamPlusSeq,
        guideSeqWPam,
        guideStrand,
        guideSeqHtml,
        guideStart,
        guideEnd,
    ) = findGuideSeq(inSeq, pam, pamId, exonId=exonId, pamFullName=pamFullName)

    # search for restriction enzymes that overlap the mutation site
    allEnzymes = readRestrEnzymes()
    mutEnzymes = matchRestrEnz(
        allEnzymes, guideSeq.upper(), pamSeq.upper(), pamPlusSeq, pam
    )

    # create a more human readable name of this guide
    if exonSeqs or pamFullName:
        guidePos = int(pamId.split(".")[1].strip("s+-")) + 1
    else:
        guidePos = int(pamId.strip("s+-")) + 1
    guideStrand = pamId[-1]
    if guideStrand == "+":
        primerGuideName = str(guidePos) + "fw"
    else:
        primerGuideName = str(guidePos) + "rv"

    # primer helper
    print(
        """
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
    """
    )

    # output the page header
    print(
        """<div style='width: 80%; margin-left:10%; margin-right:10%; text-align:left;'>"""
    )
    printBackLink()
    print("<h2>")
    if batchName != "":
        print(batchName + ":")
    print("Guide sequence: %s</h2>" % (guideSeqHtml))

    print("Contents:<br>")
    print("<ul>")
    print("<li><a href='#cloning'>Cloning or expression of guide RNA</a>")
    print(
        "<ul><li><a href='#t7plasmid'>T7 <i>in vitro</i> expression from a plasmid</a></li></ul>"
    )
    print(
        "<ul><li><a href='#t7oligo'>T7 <i>in vitro</i> expression from overlapping oligonucleotides</a></li></ul>"
    )
    print("<ul><li><a href='#geneArt'>T7 expression with the GeneArt kit</a></li></ul>")
    print(
        "<ul><li><a href='#u6plasmid'>U6 expression from an Addgene plasmid</a></li></ul>"
    )
    print(
        "<ul><li><a href='#ciona'>Direct PCR for <i>C. intestinalis</i></a></li></ul>"
    )
    print(
        "<ul><li><a href='#gibson'>Lentiviral vectors: Cloning with Gibson assembly</a></li></ul>"
    )
    print(
        "<ul><li><a href='#primerSummary'>Summary of main cloning/expression primers</a></li></ul>"
    )
    print("<li><a href='#ontargetPcr'>PCR to amplify the on-target site</a></li>")
    if len(mutEnzymes) != 0:
        print("<li><a href='#restrSites'>Restriction sites for PCR validation</a></li>")
    print("<li><a href='#offtargetPcr'>PCR to amplify off-target sites</a></li>")
    print(
        "<li><a href='#donorGuide'>Guide mutations that minimize on-target activity</a></li>"
    )
    print("<li><a href='#satMut'>Saturating mutagenesis using all guides</a></li>")
    print("</ul>")
    print("<hr>")

    printCloningSection(batchId, primerGuideName, guideSeq, params, pam)
    print("<hr>")

    targetSeq, guideStartOnTarget, guideEndOnTarget = printValidationPcrSection(
        batchId,
        genome,
        pamId,
        position,
        params,
        guideStart,
        guideEnd,
        primerGuideName,
        guideSeq,
    )
    print("<hr>")

    if len(mutEnzymes) != 0 and targetSeq is not None:
        printEnzymeSection(
            mutEnzymes, targetSeq, guideSeqWPam, guideStartOnTarget, guideEndOnTarget
        )
        print("<hr>")

    print("<h2 id='offtargetPcr'>PCR to amplify off-target sites</h2>")
    offtUrl = cgiGetSelfUrl({"otPrimers": "1"}, onlyParams=["batchId", "pamId"])
    print(
        (
            "<p>Primers for all off-targets can be downloaded from the <a href='%s'>Off-target PCR</a> page.</p>"
            % offtUrl
        )
    )

    print(
        "<h2 id='donorGuide'>BETA: Guide mutations to minimize on-target activity</h2>"
    )
    doDonorGuide = cgiGetStr(params, "donorGuide", "off")
    if doDonorGuide == "on":
        # print("Guide sequence without PAM is: <tt>%s</tt><p>" % guideSeq)
        # inSeq, genome, pamSeq, position, extSeq, multiseq, koMethod, geneModel, koGeneId, multipam = readBatchParams(batchId)
        seq, org, pam, position, guideData = readBatchAndGuides(batchId)
        for guideRow in guideData:
            (
                guideScore,
                guideCfdScore,
                effScores,
                pamStart,
                guideStart,
                strand,
                rowPamId,
                guideSeq,
                pamSeq,
                otData,
                otDesc,
                last12Desc,
                mutEnzymes,
                ontargetDesc,
                repCount,
                gcFrac,
                freeEnergy,
                doRecoding,
                cutUpstream,
                mainScore,
                beScoring,
                beOutcomes
            ) = guideRow
            if rowPamId != pamId:
                continue

            cfdSums = defaultdict(float)
            for (
                otSeq,
                score,
                cfdScore,
                editDist,
                pos,
                gene,
                alnHtml,
                inLinkage,
            ) in otData:
                # print("offtarget=%s, cfd=%f, " % (otSeq, cfdScore))
                lastThree = otSeq[: -len(pamSeq)][-3:]
                # print("last 3 bp=%s<br>" % lastThree)
                cfdSums[lastThree] += cfdScore

            sumList = [(y, x) for (x, y) in list(cfdSums.items())]
            sumList.sort()
            print(
                "Sums of CFD scores for all possible mutations of the last three nucleotides of the guide:<br>"
            )
            for cfdSum, seq in sumList:
                print(("<li>%s: %f" % (seq, cfdSum)))

    else:
        url = cgiGetSelfUrl({"donorGuide": "on"})
        print(
            (
                "<a href='%s#donorGuide'>Click here to list mutated guides "
                "sorted by off-targetactivity</a>" % url
            )
        )

    print("<h2 id='satMut'>Saturating mutagenesis using all guides</h2>")
    satMutUrl = cgiGetSelfUrl({"satMut": "1"}, onlyParams=["batchId"])
    print(
        (
            "<p>Oligonucleotides of all guides for pooled cloning into a lentiviral vector can be downloaded from the <a href='%s'>Saturating mutagenesis page</a>.</p>"
            % satMutUrl
        )
    )

    print("<hr>")

    print("</div>")


def donorDesignPage(params):

    batchId = params["batchId"]
    pamId = params["pamId"]
    doRecoding = params["doRecoding"]
    cutUpstream = params["cutUpstream"]
    insertDistance = params["insertDistance"]
    manualExStart = params.get("manualExStart")

    # must fix this
    if manualExStart == "None":
        manualExStart = None
    manualExEnd = params.get("manualExEnd")
    if manualExEnd == "None":
        manualExEnd = None
    manualExFrame = params.get("manualExFrame")
    if manualExFrame == "None":
        manualExFrame = None
    selGeneModel = params.get("geneModelSelection")
    useManualAnnotation = None
    if selGeneModel == "None":
        selGeneModel = None
    elif selGeneModel == "manual":
        useManualAnnotation = True
    selTransId = params.get("selTransId")
    if selTransId in ["allTrans", "None"]:
        selTransId = None
    pamFullName = pamId.split(".")[0]
    pam = setupPamInfo(pamFullName)

    batchInfo = readBatchAsDict(batchId)
    org = batchInfo["org"]
    insertSeq = batchInfo["insertSeq"]
    inSeq = batchInfo["seq"]
    geneId = batchInfo.get("koGeneId")
    insertPos = batchInfo.get("insertpos")
    geneModel = batchInfo.get("geneModel")
    kiType = batchInfo["kiType"]
    tagNames = batchInfo.get("tagNames")
    donorSeq = batchInfo.get("donorSeq")
    posStr = batchInfo["posStr"]
    chrom, start, end, strand = parsePos(posStr)
    htmlprefix = HTMLPREFIX

    # double nicking params
    doubleNicking = params.get("doubleNicking")
    if doubleNicking:
        revPamId = params["revPamId"]
        revDoRecoding = params["revDoRecoding"]
        revCfd = params["revCfd"]

        (revGuideSeq,
         revPamSeq, _, _,
         revGuideStrand, _,
         revGuideStart, _) = findGuideSeq(inSeq, "NGG", revPamId, pamFullName="NGG")

        fwPamId = params["fwPamId"]
        fwDoRecoding = params["fwDoRecoding"]
        fwCfd = params["fwCfd"]
        (fwGuideSeq,
         fwPamSeq, _, _,
         fwGuideStrand, _,
         fwGuideStart, _) = findGuideSeq(inSeq, "NGG", fwPamId, pamFullName="NGG")

        doubleNickingParams = {
            "doubleNicking": doubleNicking,
            "revPamId": revPamId,
            "revGuideInfo": (revGuideSeq, revPamSeq, revGuideStrand, revGuideStart),
            "revDoRecoding": revDoRecoding,
            "revCfd": revCfd,
            "fwPamId": fwPamId,
            "fwGuideInfo": (fwGuideSeq, fwPamSeq, fwGuideStrand, fwGuideStart),
            "fwDoRecoding": fwDoRecoding,
            "fwCfd": fwCfd
            }

        if any((revDoRecoding, fwDoRecoding)) is False:
            doRecoding = False
        else:
            doRecoding = True

        # recode the guide with the highest CFD score
        if revCfd >= fwCfd:
            pamId = revPamId
        else:
            pamId = fwPamId

    # check if an ssODN can be used
    maxssLen = 200 - len(insertSeq)
    if maxssLen < 40 and kiType != "deletion":
        donorTypeDisplay = "none"
        donorTypeMsg = "the insert sequence is too long to use a single stranded oligonucleotide as the HDR template"
        ssChecked = ""
        dsChecked = "checked"
    else:
        donorTypeDisplay = "block"
        donorTypeMsg = ""
        ssChecked = "checked"
        dsChecked = ""

    # from CRISPRdesigner : use the genomic strand if the cut site is downstream of the insertion site,
    # use the reverse complement if the cut site is upstream of the insertion site
    # if the cut site is within 10bp of the insertion site, use the target strand as a template

    if (cutUpstream == "upstream" and not doubleNicking) or (("onPos" in cutUpstream or doubleNicking) and strand == "-"):
        antisenseChecked = "checked"
        senseChecked = ""
    else:
        antisenseChecked = ""
        senseChecked = "checked"

    (
        guideSeq,
        pamSeq,
        pamPlusSeq,
        guideSeqWPam,
        guideStrand,
        guideSeqHtml,
        guideStart,
        guideEnd,
    ) = findGuideSeq(inSeq, pam, pamId, pamFullName=pamFullName)

    # save guide Info in params for recoding
    guideInfo = (pamSeq, guideStart, guideStrand)

    printKiSteps(batchId, step=2)

    print(
        """
    <script>
    function updateValues(source, val, cursor, text) {
        let numVal = parseInt(val);
        const donorType = document.querySelector('input[name="donorType"]:checked').value;
        let minVal = 50;
        let maxVal = 2000;
        let maxTotal = 100000;

        if (donorType === 'ss') {
            minVal = 20;
            maxVal = %(maxssLen)s - minVal;
            maxTotal = %(maxssLen)s;
        }

        if (isNaN(numVal)) return;

        // Clamp
        if (numVal < minVal) numVal = minVal;
        if (numVal > maxVal) numVal = maxVal;

        if (donorType === 'ss') {
            let isArm5 = (cursor.indexOf('arm5') !== -1);
            let otherTextId = isArm5 ? 'arm3LenText' : 'arm5LenText';
            let otherRangeId = isArm5 ? 'arm3Len' : 'arm5Len';

            let otherVal = parseInt(document.getElementById(otherTextId).value);
            if (isNaN(otherVal)) otherVal = minVal;

            if (numVal + otherVal > maxTotal) {
                let newOtherVal = maxTotal - numVal;
                if (newOtherVal < minVal) {
                    newOtherVal = minVal;
                    numVal = maxTotal - newOtherVal;
                }
                document.getElementById(otherTextId).value = newOtherVal;
                document.getElementById(otherRangeId).value = newOtherVal;
            }
        }

        document.getElementById(cursor).value = numVal;
        document.getElementById(text).value = numVal;
    }

    function clampValue(el) {
        let val = parseInt(el.value);
        const donorType = document.querySelector('input[name="donorType"]:checked').value;
        let minVal = 50;
        let maxVal = 2000;

        if (donorType === 'ss') {
            minVal = 20;
            maxVal = %(maxssLen)s - minVal;
        }

        if (isNaN(val)) {
            val = (donorType === 'ss') ? Math.floor(maxVal/2) : 800;
        }
        val = Math.max(minVal, Math.min(maxVal, val));
        el.value = val;

        let rangeId = el.id.replace('Text', '');
        updateValues('number', val, rangeId, el.id);
    }

    // prevents the form from submitting if the enter key is pressed
    function handleEnter(event) {
        if (event.keyCode === 13) {
            event.preventDefault();
            event.target.blur();
        }
    }

    function toggleTemplateStrand() {
        const donorType = document.getElementsByName('donorType')
        const strandDisplay = document.getElementById('templateStrandDisplay')
        const ssODNmsg = document.getElementById('ssODNmsg')

        const arm5Range = document.getElementById('arm5Len');
        const arm5Number = document.getElementById('arm5LenText');
        const arm3Range = document.getElementById('arm3Len');
        const arm3Number = document.getElementById('arm3LenText');
        const trimOptions = document.getElementById('trimOptions');

        let selectedValue
        for (const type of donorType) {
            if (type.checked) {
                selectedValue = type.value;
                break;
                }
        }

        if (selectedValue === 'ss') {
            strandDisplay.style.display = 'block';
            ssODNmsg.style.display = 'block';

            const minVal = 20;
            const maxVal = %(maxssLen)s - minVal;

            [arm5Range, arm5Number, arm3Range, arm3Number].forEach(el => {
                el.min = minVal;
                el.max = maxVal;
            });

            let halfVal = maxVal / 2
            if (halfVal < 45) {
                resetVal = halfVal;
                } else {
                resetVal = 45;
                };

            updateValues('number', resetVal, 'arm5Len', 'arm5LenText');
            updateValues('number', resetVal, 'arm3Len', 'arm3LenText');

        } else {
            strandDisplay.style.display = 'none';
            ssODNmsg.style.display = 'none';

            const minVal = 50;
            const maxVal = 2000;

            [arm5Range, arm5Number, arm3Range, arm3Number].forEach(el => {
                el.min = minVal;
                el.max = maxVal;
            });

            updateValues('number', 800, 'arm5Len', 'arm5LenText');
            updateValues('number', 800, 'arm3Len', 'arm3LenText');
        }
    }

    function toggleRecodeGap() {
        const recodePam = document.getElementById('recodePam');
        const recodeGapDisplay = document.getElementById('recodeGapDisplay');

        if (recodePam.checked) {
            recodeGapDisplay.style.display = 'block';
        } else {
            recodeGapDisplay.style.display = 'none';
        };
    }

    $(document).ready(function() {
        if ($('input[name="donorType"][value="ss"]').is(':checked')) {
            toggleTemplateStrand();
        }
    });

    </script>
    """
        % locals()
    )
    print(
        """<div style='width: 80%; margin-top: 54px; margin-left:10%; margin-right:10%; text-align:left;'>"""
    )
    print("""<h2>guide Sequence : %s</h2>""" % guideSeqHtml)
    if insertDistance != "None":
        if insertDistance == "0":
            print("<small>DSB at the editing site</small>")
        else:
            print(
                "<small>DSB %s bp  %s of the editing site</small>"
                % (insertDistance, cutUpstream.lstrip("onPos"))
            )

    # showDonor(donorSeq, armLen, insertPos, geneId, inSeq, kiType)

    # custom track test
    """
    ctFname = join(batchDir, batchId + ".txt")  # custom track settings
    ofh = open(ctFname, "w")
    ofh.write("browser position %s\n" % posStr)
    ofh.write(
        'track type=bigBed name="CRISPOR %(batchId)s" description="CRISPOR Results %(batchId)s %(batchId)s" itemRgb=On visibility=pack\n'
        % locals()
    )
    ofh.close()
    ctUrl = batchDir + "/%s.txt" % batchId
    """
    if geneId:
        if "ENST" in geneId:
            transcriptUrl = (
                """<a href="https://www.ensembl.org/Multi/Search/Results?q=%s;site=ensembl;page=1" target="blank">%s</a>"""
                % (geneId, geneId)
            )
        else:
            transcriptUrl = (
                """<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s/" target="blank">%s</a>"""
                % (geneId, geneId)
            )
        print("<h2>Transcript : %s</h2>" % transcriptUrl)
        if insertPos == "Nter":
            insertText = "N-terminal"
        elif insertPos == "Cter":
            insertText = "C-terminal"

        print(
            "<p>Insertion of a %sbp sequence in %s</p>" % (len(insertSeq), insertText)
        )
    if doRecoding == "True":
        print(
            """<div style="display:flex; flex-direction: row; align-items: center;">"""
        )
        htmlWarn("Donor needs recoding")
        print(
            """<p style="font-style: italic;">&nbsp Warning : This guide will likely cleave the edited sequence after knock-in. Additional mutations should be introduced in the donor DNA to prevent this (see recoding options below) </p>"""
        )
        print("</div>")
        recodeChecked = "checked"
        recodeGapDisplay = "block"
    else:
        recodeChecked = ""
        recodeGapDisplay = "none"

    if strand == "+":
        positiveStrandStr = "<small>(strand of the sequence)</small>"
        negativeStrandStr = ""
    else:
        positiveStrandStr = ""
        negativeStrandStr = "<small>(strand of the sequence)</small>"
    print("<hr>")
    print("""<form id="main" action="%s" method="GET"></form>""" % basename(__file__))
    print("""<form id="updateModel"></form>""")
    if useManualAnnotation is None and (selTransId is None or selGeneModel is None):
        geneModels, selGeneModel, selTransId = getSelGeneModel(
            org, noGenes=False, noAllTrans=True
        )
        if geneId and geneModels:
            found = False
            for model, modelStr in geneModels:
                exonInfo, _ = getExonInfo(org, model, posStr)
                for tId, sym in list(exonInfo.keys()):
                    if tId == geneId:
                        selGeneModel = model
                        selTransId = tId
                        found = True
                        break
                if found:
                    break
    else:
        geneModels = None

    commonChanges = {
        "batchId": batchId,
        "guideSeq": guideSeq,
        "guideInfo": guideInfo,
        "pam": pam,
        "pamId": pamId,
        "doRecoding": doRecoding,
        "insertDistance": insertDistance,
        "cutUpstream": cutUpstream,
        "doubleNicking": doubleNicking
    }

    updateChanges = commonChanges.copy()
    updateChanges.update({"geneModelSelection": None, "selTransId": selTransId})

    if doubleNicking:
        updateChanges.update(doubleNickingParams)

    printHiddenFields(params, updateChanges, form="updateModel")
    mainChanges = commonChanges.copy()

    mainChanges.update({"geneModelSelection": selGeneModel})

    if doubleNicking:
        mainChanges.update(doubleNickingParams)

    if geneId is None:
        mainChanges["selTransId"] = None
    else:
        mainChanges["selTransId"] = selTransId
    printHiddenFields(params, mainChanges, form="main")

    print(
        """
    <h2>Options to design the donor DNA</h2>
    <small> %(donorTypeMsg)s </small>
    <div style="margin-bottom:24px; margin-top:12px;">
        Specify the fasta header for the donor sequence (optional) : &nbsp<input form="main" name="donorName"/><br>
    </div>
    <div style="display: flex; flex-direction: column; margin-bottom: 24px;">
        <div style="display: flex; margin-bottom: 24px;">
            <div style="display: %(donorTypeDisplay)s; border: 0.5px dashed; border-color: grey; padding:8px; border-radius: 8px;">
                Select donor type  <img src=" %(htmlprefix)s image/info-small.png" title="Choose between a single stranded oligodeoxyribonucleotide (ssODN), recommended for small (<50bp) edits, and a double stranded donor DNA, recommended for knock-in of large sequences. The maximum length of the ssODN can't exceed 200 bases" class="tooltipsterInteract"><br>


                <input type="radio" form="main" %(dsChecked)s name="donorType" value="ds" autocomplete="off" onchange="toggleTemplateStrand()"/>Double-stranded donor<br>
                <input type="radio" form="main" %(ssChecked)s name="donorType" value="ss" autocomplete="off" onchange="toggleTemplateStrand()"/>Single-stranded donor<br>
            </div>
            <div id="templateStrandDisplay" style="margin-left: 5%%; margin-right:5%%; border: 0.5px dashed; border-color: grey; padding:8px; border-radius: 8px; display: none;">
                Select which strand to use as template <img src=" %(htmlprefix)s image/info-small.png" title="By default, the positive strand is used as a template for guides that introduce a DSB downstream of the editing site, and the negative strand is used if the DSB occurs upstream of this position.<br>
                If the distance between the cut site and insertion site is less than ~10bp, both strands can be used as a template. In this case, the strand of the input sequence is selected by default.<br> Otherwise, selecting a template strand ensures that the 3' homology arm is complementary to the 3' end at site of the DSB. For more information, see <a href='https://doi.org/10.1073/pnas.1711979114' target='blank'>Paix et al. 2017</a>.<br>
                For designs relying on a double nicking strategy with a pair of guides, there is no evidence for strand preference (<a href='https://doi.org/10.1038/s41598-021-98965-y' target='blank'>Schubert et al. 2021</a>), so the strand of the target sequence is selected by default." class="tooltipsterInteract"><br>

                <input type="radio" form="main" %(senseChecked)s name="polarity" value="positive" autocomplete="off"/>positive strand %(positiveStrandStr)s <br>
                <input type="radio" form="main" %(antisenseChecked)s name="polarity" value="negative" autocomplete="off"/>negative strand %(negativeStrandStr)s
            </div>
            <div id="ssODNmsg" style="display: none;">
                The maximum length of a single strand donor DNA is 200bp.<br>
                With this edit, the combined length of the homology arms can be %(maxssLen)s bp at maximum.<br>
            </div>
        </div>
    <div style="display:flex; margin-bottom:12px; flex-direction: column;">
        5' homology arm length
        <div style="margin-top:12px; display:flex; gap:12px;">
            <input type="range" form="main" id="arm5Len" name="arm5Len" value="800" min="50" max="2000" oninput="updateValues('range', this.value, this.id, 'arm5LenText')" style="width:80%%;">
            <input type="number" form="main" id="arm5LenText" value="800" min="50" max="2000" oninput="updateValues('number', this.value, this.id, 'arm5Len')" onblur="clampValue(this)" onkeydown="handleEnter(event)"> bp<br>
        </div>
        <br>
        3' homology arm length
       <div style="margin-top:12px; display:flex; gap:12px;">
            <input type="range" form="main" id="arm3Len" name="arm3Len" value="800" min="50" max="2000" oninput="updateValues('range', this.value, this.id, 'arm3LenText')" style="width:80%%;">
            <input type="number" form="main" id="arm3LenText" value="800" min="50" max="2000" oninput="updateValues('number', this.value, this.id, 'arm3Len')" onblur="clampValue(this)" onkeydown="handleEnter(event)"> bp<br>
        </div>
    </div>
          """
        % locals()
    )

    if doRecoding == "True":
        recodingMsg = "If you use this guide, at least 15 bp of the target sequence is present in the donor DNA with the default design. The donor will likely be cleaved or the locus re-cleaved after insertion by the nuclease, so recoding is strongly recommended."

    else:
        recodingMsg = "If you use this guide, the target sequence between the genome and the donor DNA will differ. The donor will not be cleaved (or re-cleaved after insertion) by the nuclease, so recoding is not needed in this case."

    print(
        """
    <div id="recodingOptions">
        <hr>
        <h2>Recoding Options</h2>
        <small style="margin-bottom:24px;"> %(recodingMsg)s </small><br>
        """
        % locals()
    )

    # Transcript model selection

    transcriptDisplay = "block"
    if (
        selGeneModel == "manual"
        and manualExStart is not None
        and manualExEnd is not None
        and manualExFrame is not None
    ):
        transcriptDisplay = "None"
        print(
            """<input type="hidden" form="main" value=True id="useManualAnnotation" name="useManualAnnotation" style="margin-top: 12px;"/>
              <p>The manual annotation you specified will be used as a model for recoding.</p>"""
        )
    elif selGeneModel is not None and selTransId is not None:
        print(
            """<p>The transcript you selected (%s) will be used as a model for recoding"""
            % selTransId
        )
        transcriptDisplay = "None"

    elif geneModels and (selGeneModel is None or selTransId is None):
        exonInfo, maxTransIdLen = getExonInfo(org, selGeneModel, posStr)
        print(
            """<div id="transcriptSelection" style="display: %s;">"""
            % transcriptDisplay
        )
        if geneId is None:
            print(
                """<p>Select a transcript ID to use as a model for recoding <img src=" %s image/info-small.png" title="This step will attempt to introduce silent mutations, so a gene model needs to be selected to get the position of codons.<br> To visualize the sequence of each transcript, you can go back to the previous step by clicking on 'Select guide sequences' above. Then, select a gene model and a transcript using the dropdown menu on top of the sequence." class="tooltipsterInteract"></p>"""
                % HTMLPREFIX
            )
            print("Gene model:")
            printDropDown(
                "geneModelSelection",
                geneModels,
                selGeneModel,
                style="width:20em",
                form="updateModel",
                onChange="updateModel.submit()",
            )
            print("Transcript:")
            transIdInfo = []
            for transId, sym in list(exonInfo.keys()):
                transIdInfo.append((transId, sym + " / " + transId))
            printDropDown(
                "selTransId", transIdInfo, selTransId, style="width:20em", form="main"
            )
            print("""<br>""")
        print("</div>")

    recodeTooltip = """
    Select the regions in which to introduce blocking mutations (you can check multiple boxes).<br>
    <ul>
        <li>If recoding is needed, the default is to introduce a PAM blocking silent mutation to prevent binding of the Cas protein (first option).</li>
        <li>By checking the second option, silent mutations are introduced in the 15 nucleotides at the PAM-proximal end of the guide (seed region) to prevent stable guide hybridization and DNA cleavage.</li>
        <li>Additionally, silent mutations can be introduced in the entire region between the cut site and insertion site to make sure that the RNP complex doesn't bind to the donor (third option).</li>
        </ul>
        For more information on recoding guidelines, see <a target='blank' href='https://doi.org/10.1038/s41598-021-98965-y'>Schubert et al. 2021</a>.
    """
    print(
        """
        <p>Choose below which regions of the donor DNA to recode </p>
        <div style="display: flex; gap: 25%%; align-items: center;">
            <div>
                <strong>Recode the donor to avoid re-cleavage</strong>
                <img src=" %(htmlprefix)s image/info-small.png" title="%(recodeTooltip)s" class="tooltipsterInteract"><br>
                <div style="display: flex; flex-direction: row; gap: 4px;">
                    <div style="margin-top: 8px;"><input type="checkbox" %(recodeChecked)s form="main" id="recodePam" name="recodePam" value="True" onchange="toggleRecodeGap()" autocomplete="off"/>Recode the PAM motif or the PAM-proximal end of the guide</div>
                    <img src=" %(htmlprefix)s image/info-small.png" style="width: 16px; height: 16px;" title="This step will attempt to introduce a single PAM-blocking mutation. If the PAM motif can't be recoded, two mutations will be introduced in the 15 PAM-proximal end of the guide (seed region). In this case, mutations near the PAM are prioritized." class="tooltipsterInteract"><br>
                </div>
                <div id="recodeGapDisplay" style="display: %(recodeGapDisplay)s; margin-top: 8px;"><input type="checkbox" name="recodeGap" form="main" value="True" autocomplete="off"/>Recode between the cut site and insertion site</div>
                </div>
            </div>
        </div>
    """
        % locals()
    )
    print(
        """
       <button form="main" style="align-self:center; margin-top: 50px; width:250px; height:50px;" type="submit" id="submit" name="submit" value="SUBMIT">Design Donor DNA</button>
       </div>
    </div>
     """
    )


def runTests():
    guideSeq = "CTCTTTACGCAGAGGGATGT"
    testRes = {
        "ATTTTTATGCAGAGTGATGT": 0.4,
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
        "CTCTTTCCACAGGGGAATGT": 0,
    }

    testRes2 = {
        "GAGTCTAAGCAGAAGAAGAA": 2.2,
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
        "GGGTAAGAGTAGAAGAAGAA": 0.8,
    }

    guideSeq = "GAGTCCGAGCAGAAGAAGAA"
    for seq, expScore in testRes2.items():
        score = calcHitScore(guideSeq, seq)


def parseArgs():
    "parse command line options into args and options"
    parser = optparse.OptionParser(
        """usage: %prog [options] org inFile guideOutFile

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
    """
    )

    parser.add_option(
        "-d",
        "--debug",
        dest="debug",
        action="store_true",
        help="show debug messages, do not delete temp directory",
    )
    parser.add_option(
        "-t", "--test", dest="test", action="store_true", help="run internal tests"
    )
    pamNames = ",".join([x for x, y in pamDesc])
    parser.add_option(
        "-p",
        "--pam",
        dest="pam",
        action="store",
        help="PAM-motif to use, default %default. TTTN triggers special Cpf1 behavior: the PAM is assumed to be 5' of the guide. Common PAMs are: "
        + pamNames,
        default="NGG",
    )
    parser.add_option(
        "-o",
        "--offtargets",
        dest="offtargetFname",
        action="store",
        help="write offtarget info to this filename",
    )
    parser.add_option(
        "-m",
        "--maxOcc",
        dest="maxOcc",
        action="store",
        type="int",
        help="MAXOCC parameter, guides with more matches are not even processed, default %default",
        default=MAXOCC,
    )
    parser.add_option(
        "",
        "--mm",
        dest="mismatches",
        action="store",
        type="int",
        help="maximum number of mismatches, default %default",
        default=4,
    )
    parser.add_option(
        "",
        "--guideLen",
        dest="guideLen",
        type="int",
        action="store",
        help="Lenght of the guide. Default is: 21 for PAM=NNGRRT/NNNRRT, 23 for Cpf1, 20 otherwise. Note: 19bp guides are less efficient",
    )
    parser.add_option(
        "",
        "--bowtie",
        dest="bowtie",
        action="store_true",
        help="new: use bowtie as the aligner. Careful: misses off-targets. Do not use.",
    )
    parser.add_option(
        "",
        "--skipAlign",
        dest="skipAlign",
        action="store_true",
        help="Assume that the input is not in the genome: do not align the input sequence. The on-target will be a random match with 0 mismatches. Switches off efficiency scoring as there is no sequence context.",
    )
    parser.add_option(
        "",
        "--noEffScores",
        dest="noEffScores",
        action="store_true",
        help="do not calculate any efficiency score",
    )
    parser.add_option(
        "",
        "--effScores",
        dest="effScores",
        action="store",
        help="calculate only these efficiency scores. Comma-sep list. Possible values, per enzyme: %s"
        % crisporEffScores.possibleScores,
    )
    parser.add_option(
        "",
        "--minAltPamScore",
        dest="minAltPamScore",
        action="store",
        type="float",
        help="minimum MIT off-target score for alternative PAMs, default %default",
        default=ALTPAMMINSCORE,
    )
    parser.add_option(
        "",
        "--worker",
        dest="worker",
        action="store_true",
        help="Run as worker process for web server: watches job queue and runs jobs",
    )
    parser.add_option(
        "",
        "--noFork",
        dest="noFork",
        action="store_true",
        help="for the --worker option: do not fork off a daemon",
    )
    parser.add_option(
        "",
        "--clear",
        dest="clear",
        action="store_true",
        help="clear the worker job table and exit",
    )
    parser.add_option(
        "-g",
        "--genomeDir",
        dest="genomeDir",
        action="store",
        help="directory with genomes, default %default",
        default=genomesDir,
    )
    parser.add_option(
        "",
        "--tempDir",
        dest="tempDir",
        action="store",
        help="temp directory for command line. If not specified, remove all temp files on exit",
    )
    parser.add_option(
        "",
        "--ampliconDir",
        dest="ampDir",
        action="store",
        help="For each guide, write a file with the off-target amplicons to this directory. Filename is <seqId>_<pamId>.txt. See repeat masking note under --satMutDir.",
    )
    parser.add_option(
        "",
        "--satMutDir",
        dest="satMutDir",
        action="store",
        help="write saturating mutagenesis files to this directory, one per input sequence: ontargetAmplicons.tsv, satMutOligos.tsv, ontargetPrimers.tsv and targetSeqs.txt. Repeats are masked when designing primers, so not all guides may have primers (flagged with 'None').",
    )
    # parser.add_option("", "--ontargetPrimers", dest="ontargetPrimers", \
    # action="store", help="write on-target primers and amplicons, one per guide, to this tab-sep file")
    parser.add_option(
        "",
        "--ampLen",
        dest="ampLen",
        type="int",
        default=140,
        action="store",
        help="for --ampliconDir/--satMutDir: amplicon length, default %default",
    )
    parser.add_option(
        "",
        "--ampTm",
        dest="tm",
        type="int",
        default=60,
        action="store",
        help="for --ampliconDir/--satMutDir: Tm for PCR, default %default",
    )
    # parser.add_option("", "--notInGenome", dest="notInGenome", \
    # action="store_true", help="Input is an articial sequence: do not try to find the input sequence in the genome, assume that the first perfect match of every guide is the on-target")

    (options, args) = parser.parse_args()

    if len(args) == 0 and not options.test and not options.worker and not options.clear:
        parser.print_help()
        sys.exit(0)

    if options.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    return args, options


def delBatchDir():
    "called at program exit, for command line mode"
    delTmpDirs()  # first remove any subdirs
    if not isdir(batchDir):
        return
    logging.debug("Deleting dir %s" % batchDir)
    fnames = glob.glob(join(batchDir, "*"))
    if len(fnames) > 50:
        raise Exception("cowardly refusing to remove many temp files")
    for fname in fnames:
        os.remove(fname)
    os.removedirs(batchDir)


tmpDirsDelExit = []


def delTmpDirs():
    "signal handler at program exit, to remove registered tmp dirs"
    global tmpDirsDelExit
    logging.debug("Removing tmpDirs: %s" % ",".join(tmpDirsDelExit))
    for tmpDir in tmpDirsDelExit:
        if isdir(tmpDir):
            shutil.rmtree(tmpDir)
    tmpDirsDelExit = []


def runQueueWorker(noFork):
    "in an infinite loop, take jobs from the job queue in jobs.db and finish them"
    # if userName!=None:
    # uid =  pwd.getpwnam(userName)[2]
    # os.setuid(uid)

    if not noFork:
        try:
            # Store the Fork PID
            pid = os.fork()

            if pid > 0:
                print("PID: %d" % pid)
                os._exit(0)

        except OSError as error:
            print("Unable to fork. Error: %d (%s)" % (error.errno, error.strerror))
            os._exit(1)

    print(
        (
            "%s: Worker daemon started. Waiting for jobs."
            % datetime.ctime(datetime.now())
        )
    )
    sys.stdout.flush()

    q = JobQueue()
    q.openSqlite()
    print(("Job queue: %s" % JOBQUEUEDB))

    global doAbort
    doAbort = False

    while True:
        if q.waitCount() == 0:
            # q.dump()
            time.sleep(1 + random.random() / 10)
            logging.debug("No job")
            continue

        if isfile("/tmp/stopCrispor"):
            logging.info("Stop signal received")
            sys.exit(0)

        jobType, batchId, paramStr = q.popJob()
        if jobType == "search":
            print("found job - single sequence mode")
            jobError = False
            try:
                (seq, org, pamDesc, position, extSeq, _, _, _, koGeneId, multipam) = (
                    readBatchParams(batchId)
                )
                pam = setupPamInfo(pamDesc)  # pamDesc includes info on guidelen, etc
                assert "-" not in pam
                uppSeq = seq.upper()
                startDict, endSet = findAllPams(uppSeq, pam)
                print("searching for offtargets:  ", seq, org, pam, position)
                getOfftargets(uppSeq, org, pamDesc, batchId, startDict, q)
            except:
                exStr = traceback.format_exc()
                print(" - WORKER CRASHED WITH EXCEPTION -")
                print(exStr)
                try:
                    q.startStep(batchId, "crash", exStr.replace("\n", "///"))
                except:
                    print(" - ALSO COULD NOT UPDATE DB WITH CRASH STATUS -")
                    print(traceback.format_exc())
                jobError = True

            if not jobError:
                q.jobDone(batchId)

        elif jobType in ["multipam", "multiseq"]:
            # search for multiple pams
            print("found job - multi sequence mode")
            logging.info("executed multisearch job")
            jobError = False
            (
                seq,
                org,
                pam,
                position,
                extSeq,
                multiseq,
                koMethod,
                geneModel,
                koGeneId,
                multipam,
            ) = readBatchParams(batchId)
            batchBase = join(batchDir, batchId)
            if jobType == "multipam":
                try:
                    seq, posStr = getPosAndSeq(org, seq, position, batchId)
                    processMultiPamSubmission(
                        org, seq, posStr, multipam, batchBase, batchId, q
                    )
                    logging.info("executed processMultiPamSubmission()")
                except:
                    exStr = traceback.format_exc()
                    print(" - WORKER CRASHED WITH EXCEPTION -")
                    print(exStr)
                    try:
                        q.startStep(batchId, "crash", exStr.replace("\n", "///"))
                    except:
                        print(" - ALSO COULD NOT UPDATE DB WITH CRASH STATUS -")
                        print(traceback.format_exc())
                    jobError = True

            else:
                try:
                    processMultiSeqSubmission(
                        multiseq, org, pam, batchBase, batchId, q, koMethod
                    )
                    logging.info("executed processMultiSeqSubmission()")
                except:
                    exStr = traceback.format_exc()
                    print(" - WORKER CRASHED WITH EXCEPTION -")
                    print(exStr)
                    try:
                        q.startStep(batchId, "crash", exStr.replace("\n", "///"))
                    except:
                        print(" - ALSO COULD NOT UPDATE DB WITH CRASH STATUS -")
                        print(traceback.format_exc())
                    jobError = True

            if not jobError:
                try:
                    q.jobDone(batchId)
                except:
                    print(" - COULD NOT MARK JOB AS DONE -")
                    print(traceback.format_exc())

        elif jobType is None:
            logging.debug("No job")
            time.sleep(1 + random.random() / 10)
        else:
            logging.error(
                "Illegal jobtype: %s - %s. Marking as done." % (jobType, batchId)
            )
            try:
                q.jobDone(batchId)
            except:
                print(" - COULD NOT MARK ILLEGAL JOB AS DONE -")
                print(traceback.format_exc())
    q.close()


def clearQueue():
    "empty the job queue"
    q = JobQueue()
    q.openSqlite()
    q.clearJobs()
    q.close()
    print(("Worker queue %s, table queue, now empty" % JOBQUEUEDB))


def handleOptions(options):
    "set global vars based on options"
    if options.test:
        runTests()
        import doctest

        doctest.testmod()
        sys.exit(0)

    # if options.ajax:
    # sendStatus(options.ajax)

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
    if options.maxOcc is not None:
        global MAXOCC
        MAXOCC = options.maxOcc
        # HIGH_MAXOCC = options.maxOcc

    if options.minAltPamScore is not None:
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
        GUIDELEN = options.guideLen
        logging.info(
            "Overriding guide length with %d bp as set on command line" % GUIDELEN
        )


def parseBed(inFname):
    "parse bed  file and return rows"
    rows = []
    for line in open(inFname):
        row = line.rstrip("\n").split()
        name = "noName"
        score = "0"
        strand = "+"

        if len(row) == 3:
            chrom, start, end = row
        elif len(row) == 4:
            chrom, start, end, name = row
        else:
            chrom, start, end, name, score, strand = row[:6]

        start, end = int(start), int(end)
        assert strand in ".+-"
        rows.append((chrom, start, end, name, score, strand))

    return rows


def mainCommandLine():
    "main entry if called from command line"
    global commandLineMode
    commandLineMode = True

    args, options = parseArgs()

    if options.worker:
        runQueueWorker(options.noFork)
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
            logging.info(
                "debug-mode, temporary directory %s will not be deleted" % batchDir
            )
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
        logging.info(
            " * running on sequence '%s', guideLen=%d, seqLen=%d"
            % (seqId, GUIDELEN, len(seq))
        )
        # get the other parameters and write to a new batch
        seq = seq.upper()
        pamPat = options.pam
        pamPat = setupPamInfo(pamPat)
        batchId = newBatch(seqId, seq, org, pamPat)
        logging.debug("Temporary output directory: %s/%s" % (batchDir, batchId))

        # if position=="?":
        # logging.error("no match found for sequence %s in genome %s" % (inSeqFname, org))

        startDict, endSet = findAllPams(seq, pamPat)

        getOfftargets(seq, org, pamPat, batchId, startDict, ConsQueue())

        batchInfo = readBatchAsDict(batchId)
        position = batchInfo["posStr"]

        otMatches = parseOfftargets(org, batchId)

        # Special batch primer / Crispresso mode
        if options.ampDir:
            pamSeqs = list(flankSeqIter(seq, startDict, len(pamPat), True))
            for (
                pamId,
                pamStart,
                guideStart,
                strand,
                guideSeq,
                pamSeq,
                pamPlusSeq,
            ) in pamSeqs:
                cPath = join(options.ampDir, "crispresso_%s_%s.txt" % (seqId, pamId))
                logging.info(
                    "Writing Crispresso table for seq %s, PAM %s to %s"
                    % (seqId, pamId, cPath)
                )
                ampLen = options.ampLen
                tm = options.tm
                primerLen = 20  # Need to add to the options

                (
                    seq,
                    org,
                    pam,
                    position,
                    extSeq,
                    multiseq,
                    koMethod,
                    geneModel,
                    koGeneId,
                    multipam,
                ) = readBatchParams(batchId)
                scoredPrimers, nameToSeq, nameToOtScoreSeq, guideSeqHtml = (
                    designOfftargetPrimers(
                        seq,
                        org,
                        pamPat,
                        position,
                        extSeq,
                        pamId,
                        ampLen,
                        primerLen,
                        tm,
                        otMatches[pamId],
                    )
                )

                cFh = open(cPath, "w")
                for row in makeCrispressoOfftargetRows(
                    scoredPrimers, nameToSeq, nameToOtScoreSeq
                ):
                    cFh.write("\t".join(row))
                    cFh.write("\n")

        # special saturation mutagenesis mode
        if options.satMutDir:
            seq, org, pam, position, guideData = readBatchAndGuides(batchId)
            satMutFname = join(options.satMutDir, seqId + "_satMutOligos.tsv")
            smFh = open(satMutFname, "w")
            logging.info("Writing saturating mutagenesis oligos to %s" % satMutFname)
            writeSatMutFile(
                0, options.ampLen, options.tm, batchId, None, None, "tsv", smFh
            )

            primerFname = join(options.satMutDir, seqId + "_ontargetPrimers.tsv")
            pFh = open(primerFname, "w")
            logging.info("Writing primers to %s" % primerFname)
            writeOntargetAmpliconFile(
                "primers", batchId, options.ampLen, options.tm, pFh
            )

            ampFname = join(options.satMutDir, seqId + "_ontargetAmplicons.tsv")
            ampFh = open(ampFname, "w")
            logging.info("Writing amplicons to %s" % ampFname)
            writeOntargetAmpliconFile(
                "amplicons", batchId, options.ampLen, options.tm, ampFh
            )

            guideFname = join(options.satMutDir, seqId + "_targetSeqs.tsv")
            gFh = open(guideFname, "w")
            logging.info("Writing guide sequences to %s" % guideFname)
            writeTargetSeqs(guideData, gFh)

        if not doEffScoring:
            effScores = {}
        else:
            effScores = readEffScores(batchId)
        logging.debug("Got efficiency scores: %s" % effScores)

        guideData, guideScores, hasNotFound, pamIdToSeq = mergeGuideInfo(
            seq, startDict, pamPat, otMatches, position, effScores, org=org
        )

        # write guide headers
        if guideFh is None:
            guideFh = open(join(batchDir, "guideInfo.tab"), "w")
            guideHeaders, _ = makeGuideHeaders()
            guideHeaders.insert(0, "#seqId")
            guideFh.write("\t".join(guideHeaders) + "\n")

        # write offtarget headers
        if options.offtargetFname and offtargetFh is None:
            offtargetFh = open(join(batchDir, "offtargetInfo.tab"), "w")
            offtargetHeaders.insert(0, "seqId")
            offtargetFh.write("\t".join(offtargetHeaders) + "\n")

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
    "send batch status as json"
    q = JobQueue()
    q.openSqlite()
    status = q.getStatus(batchId)
    q.close()

    if status is None:
        d = {"status": status}
    elif "Traceback" in status:
        d = {
            "status": "An error occured. Please send an email to %s and tell us that the failing batchId was %s. We can usually fix this quickly. Thanks!"
            % (contactEmail, batchId)
        }
    else:
        d = {"status": status}
    print(json.dumps(d))


def cleanJobs():
    """look for flag file cleanJobs in current dir. If present, remove jobs.db.
    this is the only way to remove the file, as the jobs.db file is owned by apache
    """
    if isfile("cleanJobs"):
        os.remove(JOBQUEUEDB)
        os.remove("cleanJobs")


def printAssistant(params):
    "prints the dispatcher menu for the different modes"

    expType = params.get("expType")
    libDesign = params.get("libDesign")

    if expType == "ko" or libDesign:
        active = "ko"
    elif expType == "ki":
        active = "ki"
    else:
        active = "classic"

    def cls(tabId):
        # `active` adds the orange underline (see style/assistant.css)
        return (
            "assistantButton active tooltipsterInteract"
            if tabId == active
            else "assistantButton tooltipsterInteract"
        )

    print(
        """
    <form action="crispor.py" name="main" method="get">
        <div class="assistantMenu" style="margin-bottom: 24px; margin-left: 18px; margin-top: -6px;">
            <div class="tabs" style="gap: 5%%;">
                <div class="title" style="font-size: 60px; font-weight: 500; font-family: Helvetica; text-shadow: 1px 1px 2px #fae4d1; align-self: center; color: #ff7f04;" class="toolstipsterInteract" title="CRISPOR is a program that helps design, evaluate and clone guide sequences for the CRISPR/Cas9 system.">CRISPOR</div>

                <button type="submit" name="mode" value="classic"
                        class="%s"
                        style="min-width: 100px;"
                        title="Original mode : enter a sequence to find guides.">
                    <span>
                        Classic<br>
                        <small style="font-size: 0.5em; margin-left: 12px;">&nbsp</small>
                    </span>
                </button>

                <button type="submit" name="expType" value="ko"
                        class="%s"
                        style="min-width: 400px;"
                        title="Assistant for knock-out experiments. Select a transcript and find guides to inactivate its product using different methods, including the introduction of indels resulting from Non-Homologous End Joining (NHEJ), substitutions with Base Edtiting (BE), or edits with Prime Editing (PE) <i>(not implemented yet)</i>.">
                    <span style="display: flex; flex-direction: row; gap: 10px;">
                        <span style="text-align: left;">
                            Knock-out<br>
                            <small style="font-size: 0.5em; margin-left: 12px;">(NHEJ / BE / PE)</small>
                        </span>
                        <span class="tabBadge assistant" style="height: 12px; align-self: center;">New</span>
                    </span>
                </button>

                <button type="submit" name="expType" value="ki"
                        class="%s"
                        style="min-width: 275px;"
                        title="Assistant to edit a sequence in multiple ways, including insertion, deletion, substitution, replacement, or protein tagging. Depending on the intended edit, multiple editing strategies are suggested, including Homology-Directed Repair (HDR) based editing with donor DNA design, Base Editing (BE) or Prime Editing (PE) <i>(not implemented yet)</i>.">
                    <span style="display: flex; flex-direction: row; gap: 10px;">
                        <span>
                            Precision Editing<br>
                           <small style="font-size: 0.5em; margin-left: 48px;">(HDR / BE / PE)</small>
                        </span>
                        <span class="tabBadge assistant" style="height: 12px; align-self: center;">New</span>
                    </span>
                </button>
                <div style="display: flex; flex-direction: row;">
                    <a style='width:150px; align-self: center;' href='crispor.py'>
                        <img style='width:150px; align-self: center;' src='%simage/2021-Logo-Do-3.jpg' alt='UCSC Logo'>
                    </a>
                    <a class="tooltipsterInteract" style="align-self: center; margin-top: 8px;" title="CELPHEDIA (The National Infrastructure for model organisms in health and biomedical research) is a national operational research infrastructure distributed over the French territory.<br>Its mission is to support academic and industrial scientific community to accelerate discoveries in biology and improve biomedical research. To this end, CELPHEDIA operates in 3 main activities with respect of ethical principles and animal welfare.<br>
                    <ul>
                        <li>Standardized service offers, in the areas of creation, functional exploration, archiving and distribution of animal models, necessary for fundamental research and preclinical approaches: rodents with the mouse as the leader, non-human primates and non-mammals including aquatic vertebrates.</li>
                        <li>Research and development activity for new technological offers.</li>
                        <li>Training courses adapted to users needs either for the use of animals in research with respect to institutional regulations or to develop specific technological skills.</li>
                    </ul>" href='https://celphedia.eu/en/' target="_blank">
                        <img style='width:150px; margin-left:25px' src='%simage/logo_Celphedia.jpg' alt='Celphedia'>
                    </a>
                </div>

            </div>
        </div>
    </form>
    """
        % (cls("classic"), cls("ko"), cls("ki"), HTMLPREFIX, HTMLPREFIX)
    )

    # image / drawing to display near each mode
    """
                    <svg class="tabIcon" width="18" height="18" viewBox="0 0 24 24" fill="none"
                         stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
                        <circle cx="12" cy="12" r="5" class="tabIconFill" stroke="none"/>
                    </svg>
    """


def callSubServer(name, data, timeout=60):
    """sends the data to the sub-server registered under <name> in SUBSERVERS
    and returns the JSON result"""

    port = SUBSERVERS[name]["port"]
    url = f"http://localhost:{port}/{name}"
    jsonData = json.dumps(data).encode("utf-8")
    request = urllib.request.Request(url, data=jsonData, method="POST")
    request.add_header("Content-Type", "application/json")

    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            result = response.read().decode("utf-8")
            return json.loads(result)

    except urllib.error.HTTPError as error:
        try:
            payload = json.loads(error.read().decode("utf-8"))
        except Exception:
            payload = {"error": str(error)}

        sys.stderr.write(
            "subserver %s crashed: %s\n%s\n"
            % (name, payload.get("error"), payload.get("trace", ""))
        )

        return payload

    except urllib.error.URLError as error:
        sys.stderr.write("subserver %s unreachable: %s\n" % (name, error))
        return {"error": str(error)}


def mainCgi():
    "main entry if called from apache"
    # XX I need a throttling system
    ip = os.environ.get("REMOTE_ADDR", "noIp")
    if ip == "18.141.51.207" or ip == "80.11.166.200":
        print("Content-type: text/html\n")
        print(
            "Your IP address is hammering crispor and has brought down the server for dozens of other users."
        )
        print("Please contact me at maxh@ucsc.edu.")
        sys.exit(0)

    # XX IS THE SCRIPT SYMLINKED ? XX
    if os.getcwd() != "/var/www/crispor":
        # only activate stackdumps if running on a development machine
        import cgitb

        cgitb.enable()

    # make all output files world-writable. Useful so we can clean the tempfiles
    os.umask(000)
    cleanJobs()

    # parse incoming parameters and clean them
    params = cgiGetParams()
    batchId = None

    if "ajax" in params:
        if params["ajax"] == "geneSearch":
            dbsearchGene(params, onlySymbol=True)
        elif params["ajax"] == "geneSearchCommon":
            dbsearchGene(params, onlySymbol=True, commonExons=True)
            return

    # print "Content-type: text/html\n"
    if "batchId" in params and "download" in params:
        downloadFile(params)
        return

    if "downloadDonor" in params:
        downloadDonor(params)
        return

    if "ajaxStatus" in params and "batchId" in params:
        sendStatus(params["batchId"])
        return

    # save seq/org/pam into a cookie, if they were provided
    if "org" in params and (
        "pam" in params or "multipam" in params or "expType" in params
    ):
        seq, org, pam, koMethod, multipam, expType = (
            params.get("seq"),
            params["org"],
            params.get("pam"),
            params.get("koMethod"),
            params.get("multipam"),
            params.get("expType"),
        )

        if seq is not None:
            seq, warnMsg = cleanSeq(seq, org)
        saveSeqOrgPamToCookies(seq, org, pam, koMethod, multipam, expType)

    # print headers
    if "downloadCrispresso" not in params:
        print("Content-type: text/html\n")
        # errAbort must know if it has to print this line again
        global contentLineDone
        contentLineDone = True

        title = "CRISPOR"
        if "org" in params:
            title = "CRISPOR: " + params["org"]

        printHeader(batchId, title)

    printBody(params)  # main dispatcher, branches based on the params dictionary

    printTeforBodyEnd()

    # some keywords for google searches
    print(
        """<div style='display:none'>CRISPR/Cas9 Guide Designer for chordate
    vertebrate ecdysozoans lophotrochozoans protostomes spongi corals plants
    butterflies metazoans genomes fruitflies insects nematodes mammals.
    </div>"""
    )
    print("</body></html>")


def handle_pdb(sig, frame):
    import pdb

    pdb.Pdb().set_trace(frame)


def main():
    # send SIGUSR1 to the process to get a pdb console. Handy for live debugging.
    signal.signal(signal.SIGUSR1, handle_pdb)

    if "REQUEST_METHOD" in os.environ and sys.argv[-1] != "--worker":
        mainCgi()
    else:
        mainCommandLine()


if __name__ == "__main__":
    main()
