#!/usr/bin/env python
import subprocess
import Cookie, time
import cgitb

cgitb.enable()

# main script of the tefor crispr tool

# cleaning things todo:
# - maybe: make temp subdirectories instead of many files with same name

import sys, cgi, re, array, random, platform, os, hashlib, base64, string, logging, operator, urllib
from collections import defaultdict, namedtuple
from sys import stdout
from os.path import join, isfile, basename, dirname

DEBUG = False
#DEBUG = True

# the segments.bed files use abbreviated genomic region names
segTypeConv = {"ex":"exon", "in":"intron", "ig":"intergenic"}

# directory where already processed batches of offtargets are stored ("cache" of bwa results)
batchDir = "temp"

def getVars():
    form = cgi.FieldStorage()
    seq = form.getfirst("seq")
    if seq == None:
        #errAbort("CGI var seq is required")
        pamId = form.getfirst("pamId")
        batchId = form.getfirst("batchId")
        if pamId!=None and batchId!=None:
            return None, None, None, pamId, batchId
        else:
            return None, None, None, None, None

    org = form.getfirst("org")
    if org is None:
        errAbort("CGI var org is required")        
    
    pam = form.getfirst("pam")
    if pam is None:
        errAbort("CGI var pam is required")
    if len(set(pam)-set("ACTGNMK"))!=0:
        errAbort("Illegal character in PAM-sequence. Only ACTGMK and N allowed.")
    return seq, org, pam, None, None

def makeTempBase(seq, org, pam):
    "create the base of temp files using a hash function and some prettyfication "
    hasher = hashlib.sha1(seq+org+pam)
    batchId = base64.urlsafe_b64encode(hasher.digest()[0:20]).translate(transTab)[:20]
    return batchId

seq, org, pam, pamId, batchId = getVars()

cookies=Cookie.SimpleCookie()

expires = 365 * 24 * 60 * 60
defaultorg = 'ensDanRer'
defaultseq = 'CCAATCAGGTCCCTCCCTACCTCAGATCGCAGCTATAATACATAGGAGTAAAGAGGCTTCTCGCATTAAGTGGCTGTGGCTTGAAGTAACGTTGTGATTTCGAGGTCAGTCTTACCTTTCGCATCCCCGCCGCAAACCTCCGATGCGTTATCAGTCGCACGTTTCCGCACCTGTCACGGTCGGGGCTTGGCGCTGCTGAGGGACACGCGTGAACCGAGGAGACGGCAAGGACATCGCCGGAGATCCGCGCCTCGACAACGAGAAACCCTGCTAGACAGACCGCTCGAGAACACCGCAGCGAGATTCAGCGTGCGGCAAAATGCGGCTTTTGACGAGAGTGCTGCTGGTGTCTCTTCTCACTCTGTCCTTGGTGGTGTCCGGACTGGCCTGCGGTCCTGGCAGAGGCTACGGCAGAAGAAGACATCCGAAGAAGCTGACACCTCTCGCCTACAAGCAGTTCATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATAACGCGCAATTCGGAGAGATTTAAAGAACTTACTCCAAATTATAATCCCGACATTATCTTTAAGGATGAGGAGAACACGGGAGCGGACAGGCTCATGACACAG'
defaultpam = 'NGG'

try:
    cookies['lastorg']
    cookies['lastseq']
    cookies['lastpam']
    cookies['lastvisit'] = str(time.time())
    cookies['lastvisit']['expires'] = expires
except KeyError:    
    cookies['lastorg'] = defaultorg
    cookies['lastorg']['expires'] = expires
    cookies['lastpam'] = defaultpam
    cookies['lastpam']['expires'] = expires
    cookies['lastseq'] = defaultseq
    cookies['lastseq']['expires'] = expires
    cookies['lastvisit'] = str(time.time())
    cookies['lastvisit']['expires'] = expires

    
if 'HTTP_COOKIE' in os.environ:
    cookie_string=os.environ.get('HTTP_COOKIE')    
    cookies.load(cookie_string)

    try:                                
        if org is not None :
            cookies['lastorg'] = org
            cookies['lastorg']['expires'] = expires
    except KeyError:
       print "the cookie was not set or has expired<br>"

    try:                                
        if seq is not None :
            cookies['lastseq'] = seq
            cookies['lastseq']['expires'] = expires
    except KeyError:
       print "the cookie was not set or has expired<br>"
    
    try:                                
        if pam is not None :
            cookies['lastpam'] = pam
            cookies['lastpam']['expires'] = expires
    except KeyError:
       print "the cookie was not set or has expired<br>"

print cookies

print "Content-type: text/html\n"

def debug(msg):
    if DEBUG:
        print msg
        print "<br>"

def errAbort(msg):
    print(msg+"<p>")
    sys.exit(0)

def matchNuc(pat, nuc):
    " returns true if pat (single char) matches nuc (single char) "
    if pat in ["A", "C", "T", "G"] and pat==nuc:
        return True
    elif pat=="M" and nuc in ["A", "C"]:
        return True
    elif pat=="K" and nuc in ["T", "G"]:
        return True
    else:
        return False

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
    if seq.startswith("random"):
        seq = rndSeq(800)
    lines = seq.strip().splitlines()
    if len(lines)>0 and lines[0].startswith(">"):
        line1 = lines.pop(0)

    if len(seq)>2000:
        errAbort("Sorry, cannot handle sequences longer than 2kbp, Please try again with a shorter sequence.")

    newSeq = []
    nCount = 0
    for l in lines:
        if len(l)==0:
            continue
        for c in l:
            if c not in "actgACTG":
                nCount +=1 
                #newSeq.append("A")
            else:
                newSeq.append(c)
    msg = ""
    if nCount!=0:
        msg = "Sequence contained %d non-ACTG letters. They were removed." % nCount

    return "".join(newSeq), msg

def revComp(seq):
    " rev-comp a string "
    revTbl = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N' , 'M' : 'K', 'K':'M'}
    newSeq = []
    for c in reversed(seq):
        newSeq.append(revTbl[c])
    return "".join(newSeq)
        
def findPams(seq, pam, strand, startDict, endSet):
    " return two values: dict with pos -> strand of PAM and set of end positions of PAMs"

    minPos = 20
    maxPos = len(seq)-(20+len(pam))

    #print "new search", seq, pam, "<br>"
    for start in findPat(seq, pam):
        # need enough flanking seq on one side
        #print "found", start,"<br>"
        if strand=="+" and start<=minPos:
            #print "no, out of bounds +", "<br>"
            continue
        if strand=="-" and start>=maxPos:
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

def showSeq(seq, lines, maxY, pam, genomeName):
    " show the sequence and the PAM sites underneath "
    print "<div class='title'>Sequence with PAMs in <div class='speciesname'>%s</div> genome</div>" % genomeName
    print "<div class='substep'>"
    print "Click on a PAM match (<div style=\'display:inline;\' class=\'linklike\'>%s</div>) to show potential guide sequences for it" % pam
    print "</div>"
    print '''<div style="overflow-x:scroll; width:100%; background:#EEEEEE; border-style: solid; border-width: 1px">'''

    print "<pre>"+rulerString(len(seq))
    print seq

    for y in range(0, maxY+1):
        texts = []
        lastEnd = 0
        for start, end, name, strand in lines[y]:
            spacer = "".join([" "]*((start-lastEnd)))
            lastEnd = end
            texts.append(spacer)
            pamId = "s"+str(start)+strand
            texts.append('''<a id="list%s" href="#%s" onmouseover="$('.hiddenExonMatch').show('fast');$('#show-more').hide();$('#show-less').show()">''' % (pamId,pamId))
            texts.append(name)
            texts.append("</a>")
        print "".join(texts)
    print("</pre><br>")

    print '''</div>'''
    
def flankSeqIter(seq, startDict, pamLen):
    """ given a seq and dictionary of pos -> strand and the length of the pamSite
    yield 20mers flanking the sites sorted by pos
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

        yield startPos, strand, flankSeq, pamSeq

def printBrowserLink(dbInfo, pos, text, alnStr):
    " print link to genome browser (ucsc or ensembl) at pos, with given text "
    if dbInfo.server.startswith("Ensembl"):
        baseUrl = "www.ensembl.org"
        if dbInfo.server=="EnsemblPlants":
            baseUrl = "plants.ensembl.org"
        elif dbInfo.server=="EnsemblMetazoa":
            baseUrl = "metazoa.ensembl.org"
        org = dbInfo.scientificName.replace(" ", "_")
        url = "http://%s/%s/Location/View?r=%s" % (baseUrl, org, pos)
    elif dbInfo.server=="ucsc":
        url = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s" % (dbInfo.name, pos)
    else:
        print "unknown genome browser server %s, please email penigault@tefor.net" % dbInfo.server

    print '''<a title="%s" target="_blank" href="%s">%s</a>''' % (alnStr, url, text)

def makeAlnStr(seq1, seq2, pam):
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
    #lines[0].append("<i>"+seq1[i:i+3]+"</i>")
    lines[0].append(" <i>"+seq1[-len(pam):]+"</i>")
    lines[1].append(" <i>"+seq2[-len(pam):]+"</i>")
    #lines[1].append("<i>"+seq2[i:i+3]+"</i>")
    lines = ["".join(l) for l in lines]
    htmlText = "<pre>guide:      %s<br>off-target: %s<br>            %s</pre>" % (lines[0], lines[1], lines[2])
    hasLast12Mm = last12MmCount>0
    return htmlText, hasLast12Mm
        
def calcMmScore(guideSeq, otSeq):
    " return mismatch score for a given off target site "
    guideSeq = guideSeq[:20]
    otSeq = otSeq[:20]
    score = 0
    for i in range(len(guideSeq)-1, 0, -1):
        if guideSeq[i]!=otSeq[i]:
            score += i
    return score

def makePosList(countDict, guideSeq, pam):
    """ for a given guide sequence, return a list of (score, posStr, geneDesc, otSeq) sorted by score and a string to describe the offtargets in the format x/y/z/w of mismatch counts"""
    # one desc in last column per OT seq
    #countDict = otMatches[pamId]
    count = 0
    otCounts = []
    posList = []
    scores = []
    last12MmCounts = []
    # for each edit distance, get the off targets and iterate over them
    for editDist in range(0, 5):
        #print countDict,"<p>"
        matches = countDict.get(editDist, [])
        #print matches

        # create a list of number of offtargets for this edit dist
        otCounts.append( str(len(matches)))
        #print otCounts,"<p>"
        last12MmOtCount = 0

        for chrom, start, end, otSeq, strand, segType, geneNameStr in matches:
            posStr = "%s:%d-%s" % (chrom, int(start)+1,end)
            segTypeDesc = segTypeConv[segType]
            geneDesc = segTypeDesc+":"+geneNameStr
            geneDesc = geneDesc.replace("|", "-")
            score = calcMmScore(guideSeq, otSeq)
            scores.append(score)
            alnHtml, hasLast12Mm = makeAlnStr(guideSeq, otSeq, pam)
            if not hasLast12Mm:
                last12MmOtCount+=1
            posList.append( (score, editDist, posStr, geneDesc, alnHtml) )
        last12MmCounts.append(str(last12MmOtCount))

    posList.sort()
    otDescStr = "/".join(otCounts)
    last12DescStr = "/".join(last12MmCounts)

    return posList, otDescStr, sum(scores), last12DescStr


def showTable(seq, startDict, pam, otMatches, dbInfo, batchId):
    " shows table of all PAM motif matches "
    # one row per guide sequence
    guideData = []

    for startPos, strand, flankSeq, pamSeq in flankSeqIter(seq, startDict, len(pam)):
        
        # position with anchor to jump to
        pamId = "s"+str(startPos)+strand
        # flank seq
        seqStr = "<tt>"+flankSeq + " <i>" + pamSeq+"</i></tt>"
        guideSeq = flankSeq + pamSeq

        #print '''<td><tt>%s</tt></td>''' % seqStr

        # matches in genome
        # one desc in last column per OT seq
        if pamId in otMatches:
            posList, otDesc, guideScore, last12Desc = makePosList(otMatches[pamId], guideSeq, pam)
            #print otDesc, "<p>"
        else:
            posList, otDesc, guideScore = None, "No match. Incorrect genome?", [0]
            last12Desc = ""
            #print "no match"
        guideData.append(( guideScore, startPos, strand, pamId, seqStr, guideSeq, posList, otDesc, last12Desc))

    guideData.sort()

    print "<br><div class='title'>Potential guide sequences for (%s) PAMs</div>" % pam
    print "<div class='substep'>(ranked from highest to lowest specificity score determined as in <a target='_blank' href='http://dx.doi.org/10.1038/nbt.2647'>Hsu et al.</a>)</div>"
    print '<table id="otTable">'
    print "<tr>"
    #print '''<a href="#" onclick="$('#otTable .hasExonMatch').hide(); alert('testRemovedAll')">hide offtargets with matches in exons</a>'''
    print '<th>Position/ <br>Strand</th><th>Guide Sequence + <i>PAM</i></th><th style="width:250px;">target and off-target matches for 0,1,2,3,4 mismatches. (<i>In italic: with no mismatches in the 12 bp adjacent to the PAM.</i>)</th><th>Genome Browser links to target and off-targets</th>'
    #print "</tr>"

    #print '''<td><a>%d/%s</a></td>''' % (startPos, strand)
    # print '''<td><a name="%s">%d/%s</a></td>''' % (pamId, startPos, strand)
    
    print '''
    <div id="show-more" class="button" 
         onclick="$('.hiddenExonMatch').show('fast');$(this).hide();$('#show-less').show()" 
         style="margin-left:auto;margin-right:auto;width:100px">
         Show More [+]
    </div>
    '''
    print '''
    <div id="show-less" class="button" 
         onclick="$('.hiddenExonMatch').hide('fast');$(this).hide();$('#show-more').show();"
         style="margin-left:auto;margin-right:auto;width:100px;display:none;">
         Show Less [-]
    </div>
    '''
    print '''<a href="http://tefor.net/crispor/download.php?batchId=%s&amp;seq=%s&amp;org=%s&amp;pam=%s&amp;pamId=%s">
                <!--<div class="button" style="margin-left:auto;margin-right:auto;width:150px;">-->
                    <img style="width:20px;vertical-align:middle;"
                         src="http://tefor.net/crispor/image/doc.png">
                    Download results
                <!--</div>-->
            </a>
            <br><br>
    ''' % (batchId,seq,org,pam,pamId)
    count = 0    
    for guideRow in guideData:
        guideScore, startPos, strand, pamId, seqStr, guideSeq, posList, otDesc, last12Desc = guideRow 
        #if count is 0:
            #print '''
            #<tr style="cursor:pointer;font-weight:bold;color:rgb(47, 129, 203);font-size:larger;">
            #    <td id="show-more" onclick="$('.hiddenExonMatch').show('fast');$(this).hide();$('#show-less').show()" colspan="4" style="text-align:center;padding:15px;">
            #        <div class="button" style="margin-left:auto;margin-right:auto;width:200px"> Show More Results [+] </div>
            #    </td>
            #    <td id="show-less" onclick="$('.hiddenExonMatch').hide('fast');$(this).hide();$('#show-more').show();" colspan="4" style="text-align:center;display:none;padding:15px;">
            #        <div class="button" style="margin-left:auto;margin-right:auto;width:200px"> Show Less Results [-] </div>
            #    </td>                
            #</tr>
            #'''                        

        print '<tr id="%s" class="hasExonMatch' % (pamId)
        if count >=10:
            print ' hiddenExonMatch'
        print '">'        
        print "<td>"
        print '<a href="#list%s">' % (pamId)
        print str(startPos)+"/"
        if strand=="+":
            print 'forward'
        if strand=="-":
            print 'reverse'
        print '</a>'
        print "</td>"

        # sequence
        print "<td>"
        print "<small>"
        print seqStr
        print "<br>"

        scriptName = basename(__file__)
        if posList!=None:
            print '<br><a href="%s?batchId=%s&pamId=%s">Primer Design</a>' % (scriptName, batchId, urllib.quote(str(pamId)))
        print "</small>"
        print "</td>"

        # mismatch description
        print "<td>"
        print otDesc
        print "<br><small><i>", last12Desc, "</i></small>"
        print "</td>"

        # links to offtargets
        print "<td><small>"
        if posList!=None:
            for score, editDist, pos, gene, alnHtml in posList:
                print '''%d:''' % (int(editDist))
                printBrowserLink(dbInfo, pos, gene, alnHtml)

        print "</small></td>"

        print "</tr>"
        count = count+1

    print "</table>"
        
def printHeader(seq,org,pam):
    
    print "<html><head>"   

    proc = subprocess.Popen("php /var/www/crispor/header.php", shell=True, stdout=subprocess.PIPE)
    script_response = proc.stdout.read()
    print script_response
    proc = subprocess.Popen("php /var/www/main/specific/googleanalytics/script.php", shell=True, stdout=subprocess.PIPE)
    script_response = proc.stdout.read()
    print script_response

    print ("""  <script>
               $(function () {
                  $(document).tooltip({
                  relative : true,
                  content: function () {
                  return $(this).prop('title');
                  }
                 });
              });
              </script>""")
    print("""<style>
        .ui-tooltip {
            background-color: #FFFFFF;
            width: 400;
            height: 80;
            position : absolute;
            text-align: left;
            border:1px solid #cccccc;
            }
            </style>""")
    print("</head>")
    print'<body id="wrapper"'
    
    if seq != None and org != None and pam != None:
        localbatchId = makeTempBase(seq, org, pam)    
    else :
        localbatchId = None

    if localbatchId is not None:
        print '''
        onload="history.pushState('crispor/crisporDev.cgi', document.title, '?batchId=%s');"
        ''' % (localbatchId)
    print'>'
    print "<div id='fb-root'></div>"
    print('<script src="http://tefor.net/crispor/facebooklike.js"></script>')    

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
    MAXLINES = 10
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
        strand = startDict[start]
        end = start+featLen
        y = firstFreeLine(lineMasks, 0, start, end)

        if y==None:
            errAbort("not enough space to plot features")

        # fill the current mask
        mask = lineMasks[y]
        for i in range(max(start-SLOP, 0), min(end+SLOP, len(seq))):
            mask[i]=1

        maxY = max(y, maxY)
        name = seq[start:end] 
        #if strand=="+":
            #name = ">"+name
        #else:
            #name = name+"<"
        ft = (start, end, name, strand) 
        ftsByLine[y].append(ft )
    return ftsByLine, maxY

def writePamFlank(seq, startDict, pam, faFname):
    " write pam flanking sequences to fasta file "
    #print "writing pams to %s<br>" % faFname
    faFh = open(faFname, "w")
    for startPos, strand, flankSeq, pamSeq in flankSeqIter(seq, startDict, len(pam)):
        faFh.write(">s%d%s\n%s\n" % (startPos, strand, flankSeq))
    faFh.close()

def runCmd(cmd):
    " run shell command, check ret code, replaces BIN and SCRIPTS special variables "
    sysId = platform.system()
    binDir = "bin/"+sysId
    scriptDir = "bin"
    if __file__.endswith("Dev.cgi"):
        scriptDir = "binDev"

    cmd = cmd.replace("BIN", binDir)
    cmd = cmd.replace("SCRIPT", scriptDir)
    cmd = "set -o pipefail; " + cmd
    debug("Running %s" % cmd)
    #print cmd
    #ret = os.system(cmd)
    ret = subprocess.call(cmd, shell=True, executable="/bin/bash")
    if ret!=0:
        print "Server error: could not run command %s.<p>" % cmd
        print "please send us an email, we will fix this error as quickly as possible. penigault@tefor.net "
        sys.exit(0)

def parseOfftargets(bedFname):
    """ parse a bed file with annotataed off target matches from overlapSelect,
    has two name fields, one with the pam position/strand and one with the
    overlapped segment 
    
    return as dict pamId -> editDist -> (chrom, start, end, strand, segType, geneNameStr)
    segType is "ex" "int" or "ig" (=intergenic)
    if intergenic, geneNameStr is two genes, split by |
    """
    # example input:
    # chrIV 9864393 9864410 s41-|-|5    chrIV   9864303 9864408 ex:K07F5.16
    # chrIV   9864393 9864410 s41-|-|5    chrIV   9864408 9864470 in:K07F5.16
    debug("reading %s" % bedFname)

    # if a offtarget overlaps an intron/exon or ig/exon boundary it will appear twice
    # in this case, we only keep the exon offtarget
    # first sort into dict (pamId,chrom,start,end,editDist,strand) -> (segType, segName)
    #pamData = defaultdict(dict)
    pamData = {}
    for line in open(bedFname):
        fields = line.rstrip("\n").split("\t")
        chrom, start, end, name, segment = fields
        pamId, strand, editDist, seq = name.split("|")
        editDist = int(editDist)
        # some gene models include colons
        segType, segName = string.split(segment, ":", maxsplit=1)
        #pamData[(pamId].setdefault(editDist, []).append( (chrom, start, end, strand, segType, segName) )
        otKey = (pamId, chrom, start, end, editDist, seq, strand)
        if otKey in pamData and segType!="ex":
            continue
        pamData[otKey] = (segType, segName)

    indexedOts = defaultdict(dict)
    for otKey, otVal in pamData.iteritems():
        pamId, chrom, start, end, editDist, seq, strand = otKey
        segType, segName = otVal
        indexedOts[pamId].setdefault(editDist, []).append( (chrom, start, end, seq, strand, segType, segName) )

    return indexedOts

def findOfftargets(faFname, genome, pam, bedFname):
    " search fasta file against genome, filter for pam matches and write to bedFName "
    pamLen = len(pam)
    cmd = "BIN/bwa aln -n 4 -o 0 -l 20 -k 4 -N -m 1000000000 genomes/%(genome)s/%(genome)s.fa %(faFname)s | BIN/bwa samse -n 100000000000 genomes/%(genome)s/%(genome)s.fa /dev/stdin %(faFname)s  | SCRIPT/xa2multi.pl | SCRIPT/samToBed %(pamLen)s | BIN/bedClip stdin genomes/%(genome)s/%(genome)s.sizes stdout | BIN/twoBitToFa genomes/%(genome)s/%(genome)s.2bit stdout -bed=stdin | SCRIPT/filterFaToBed %(pam)s | BIN/overlapSelect genomes/%(genome)s/%(genome)s.segments.bed stdin stdout -mergeOutput -selectFmt=bed -inFmt=bed | cut -f1,2,3,4,8 2> /tmp/log > %(bedFname)s " % locals()
    #cmd = "echo mainScript > /tmp/log"
    runCmd(cmd)

transTab = string.maketrans("-=/+_", "abcde")

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

def printOrgDropDown(lastorg):
    " print the organism drop down box "
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

    print '<select style="width:100%%;max-width:250px;" name="org">'
    for db, desc in genomes:
        print '<option '
        if db == lastorg :
            print 'selected '
        print 'value="%s">%s</option>' % (db, desc)
    print "</select>"

def printPamDropDown(lastpam):        
    pams = []
    pams = [ ('NGG','NGG - Streptococcus Pyogenes'),
             ('NNAGAA','NNAGAA - Streptococcus Thermophilus'),
             ('NNNNGMTT','NNNNG(A/C)TT - Neisseiria Meningitidis'),
             ('NNNNACA','NNNNACA - Campylobacter jejuni')
           ]
    
    print '<select style="width:100%%;max-width:250px;" name="pam">'
    for key,value in pams:        
        print '<option '
        if key == lastpam :
            print 'selected '
        print 'value="%s">%s</option>' % (key, value)
    print "</select>"           

def printForm(defaultorg,defaultseq,defaultpam):
    " print html input form "
    scriptName = basename(__file__)        
    
    # The returned cookie is available in the os.environ dictionary
    cookie_string = os.environ.get('HTTP_COOKIE')
        # The first time the page is run there will be no cookies
    if not cookie_string:
       #print '<p>First visit or cookies disabled</p>'
       lastorg = defaultorg
       lastseq = defaultseq
       lastpam = defaultpam
    else:
        #print '<p>The returned cookie string was "' + cookie_string + '"</p>'
        # load() parses the cookie string
        cookies.load(cookie_string)
        # Use the value attribute of the cookie to get it
        lastvisit = float(cookies['lastvisit'].value)
        lastorg   = cookies['lastorg'].value
        lastseq   = cookies['lastseq'].value
        lastpam   = cookies['lastpam'].value

        #print '<p>Your last visit was at',
        #print time.asctime(time.localtime(lastvisit)), '</p>'
        #print '<p>Your last org was: ',
        #print lastorg, '</p>'        
         
    print """
<form id="main-form" method="post" action="%s">

<div class="introtext">
    Find more about 
    <div onclick="$('#about-us').toggle('fast');" class="title" style="cursor:pointer;display:inline;font-size:large;font-style:normal">
        CRISPOR
        <img src="http://tefor.net/crispor/image/info.png" class="infopoint" style="vertical-align:text-top;">
    </div>
    <div id="about-us"> CRISPOR - CRISPr selectOR - is a program that helps design and evaluate target sites for use with the CRISPR/Cas9 system.<br>
    It uses the BWA algorithm to identify guide RNA sequences for CRISPR mediated genome editing.<br>
    It searches for off-target sites (with and without mismatches), shows them in a table and annotates them with flanking genes.<br>
    For more information on principles of CRISPR-mediated genome editing, check the <a href="https://www.addgene.org/CRISPR/guide/">Addgene CRISPR guide</a>.</div>
</div>

<div class="windowstep subpanel" style="width:60%%;">
    <div class="substep">        
        <div class="title" style="cursor:pointer;" onclick="$('#helptext1').toggle('fast')">
            Step 1
            <img src="http://tefor.net/crispor/image/info.png" class="infopoint" >
        </div>
       Submit a single sequence for guide RNA identification and analysis
    </div>

    <textarea style="width:100%%;" name="seq" rows="10"
              placeholder="Enter the sequence of the gene you want to target - example: %s">
    %s
    </textarea>
    <div id="helptext1" class="helptext">CRISPOR conserves the lowercase and uppercase format of your sequence (allowing to highlight sequence features of interest such as ATG or STOP codons)</div>
    <input style="margin-top:20px;" type="submit" name="submit" value="SUBMIT" />
</div>
<div class="windowstep subpanel">
    <div class="substep">               
        <div class="title" style="cursor:pointer;" onclick="$('#helpstep2').toggle('fast')">
            Step 2
            <img src="http://tefor.net/crispor/image/info.png" class="infopoint">
        </div>                                          
        Choose a species genome
        
    </div>
    """% (scriptName,lastseq,lastseq)

    printOrgDropDown(lastorg)

    print """<div id="helpstep2" class="helptext">More information on these species can be found on the <a href="http://www.efor.fr">EFOR</a> website.
To add your genome of interest to the list, contact CRISPOR web site manager
<a href="mailto:penigault@tefor.net">Jean-Baptiste Penigault</a>.</div>
"""
    print """        
</div>
<div class="windowstep subpanel">
    <div class="substep">
        <div class="title" style="cursor:pointer;" onclick="$('#helpstep3').toggle('fast')">
            Step 3
            <img src="http://tefor.net/crispor/image/info.png" class="infopoint">
        </div>        
        Choose a Protospacer Adjacent Motif (PAM)        
    </div>
    """
    printPamDropDown(lastpam)
    print """
    <div id="helpstep3" class="helptext">The most common system uses the NGG PAM recognized by Cas9 from S. <i>pyogenes</i></div>
</div>
</form>
    """

def crisprSearch(seq, org, pam):
    " do crispr off target search "
    # read genome info tab file
    myDir = dirname(__file__)
    genomesDir = join(myDir, "genomes")
    infoFname = join(genomesDir, org, "genomeInfo.tab")    
    dbInfo = lineFileNext(open(infoFname)).next()

    caseSeq, userMsg = cleanSeq(seq)
    seq = caseSeq.upper()

    printHeader(seq,org,pam)

    # search pams
    startDict, endSet = findPams(seq, pam, "+", {}, set())
    startDict, endSet = findPams(seq, revComp(pam), "-", startDict, endSet)
    #print startDict, endSet

    # define temp file names
    batchId = makeTempBase(seq, org, pam)
    batchBase = join(batchDir, batchId)
    print "<!-- BATCH ID %s -->" % batchId

    # save input seq, pamSeq and genome for primer design later
    inputFaFname = batchBase+".input.fa"
    open(inputFaFname, "w").write(">%s %s\n%s\n" % (org, pam, seq))

    # write guides to fasta and run bwa
    faFname = batchBase+".fa"
    otBedFname = batchBase+".bed"

    if not isfile(otBedFname):
        # write potential PAM sites to file
        writePamFlank(seq, startDict, pam, faFname)
        findOfftargets(faFname, org, pam, otBedFname)
    otMatches = parseOfftargets(otBedFname)

    featLen = len(pam)
    lines, maxY = distrOnLines(seq, startDict, featLen)

    showSeq(caseSeq, lines, maxY, pam, dbInfo.scientificName)
    showTable(seq, startDict, pam, otMatches, dbInfo, batchId)


def printSkeleton(seq, org, pam, pamId, batchId,defaultorg,defaultseq,defaultpam):
    # add logo
    print "<a href='http://tefor.net/main/'><img class='logo' src='http://tefor.net/main/images/logo/logo_tefor.png' alt='logo tefor infrastructure'></a>"
    # add menu    
    proc = subprocess.Popen("php /var/www/crispor/menu.php", shell=True, stdout=subprocess.PIPE)
    script_response = proc.stdout.read()
    print script_response

    print '<div id="bd">'

    print '<div class="centralpanel">'
    
    proc = subprocess.Popen("php /var/www/crispor/networking.php", shell=True, stdout=subprocess.PIPE)
    script_response = proc.stdout.read()
    print script_response
    
    print '<div class="subpanel">'

    # print '<div class="title">CAS9 OFF-TARGET PARSER</div>'

    # print '<div class="title">CRISPOR</div>'

    print '<div class="contentcentral">'
    
    printContent(seq, org, pam, pamId, batchId,defaultorg,defaultseq,defaultpam)
        
    print '</div>'
            
            
    print '</div>'
    print '</div>'

    proc = subprocess.Popen("php /var/www/crispor/footer.php", shell=True, stdout=subprocess.PIPE)
    script_response = proc.stdout.read()
    print script_response

    print '</div>'

def parseBoulder(tmpOutFname):
    " parse a boulder IO style file "
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

def parseFasta(fname):
    " parse a fasta file, where each seq is on a single line "
    seqs = {}
    parts = []
    seqId = None
    for line in open(fname):
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

def makePrimers(batchId, pamId):
    " create primers with primer3 around site identified by pamId in batch with batchId. Output primers as html "
    # get offtargets
    #otBedFname = join(batchDir, batchId+".bed")
    #otMatches = parseOfftargets(otBedFname)
    #if pamId not in otMatches:
        #errAbort("Cannot find pamId %s in batch %s. Internal error? Please email penigault@tefor.net." % (pamId, batchId))
    #pamOts = otMatches[pamId]

    batchBase = join(batchDir, batchId)

    faFname = batchBase+".input.fa"
    ifh = open(faFname)
    genome, pamSeq = ifh.readline().replace(">","").strip().split()
    inSeq = ifh.readline().strip()
    seqLen = len(inSeq)

    # find position of pam in inSeq
    pamFname = batchBase+".fa"
    pams = parseFasta(pamFname)
    guideSeq = pams[pamId]
    guideStrand = pamId[-1]
    if pamId.endswith("+"):
        guideSeqPlus = guideSeq
        guideStart = inSeq.find(guideSeq)
    else:
        guideSeqPlus = revComp(guideSeq)
        guideStart = inSeq.find(guideSeqPlus)

    samFname = batchBase+".input.sam"

    #print pamOts

    cmd = "BIN/bwa bwasw genomes/%(genome)s/%(genome)s.fa %(faFname)s > %(samFname)s" % locals()
    runCmd(cmd)

    chrom, start, end = None, None, None
    for l in open(samFname):
        if l.startswith("@"):
            continue
        l = l.rstrip("\n")
        fs = l.split("\t")  
        qName, flag, rName, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = fs[:11]
        if not re.compile("[0-9]*").match(cigar):
            continue
        if cigar=="*":
            errAbort("Sequence not found in genome. Are you sure you have pasted the correct sequence and also selected the right genome?")
        matchLen = int(cigar.replace("M",""))
        chrom, start, end =  rName, int(pos), int(pos)+matchLen

    # possible problem: we do not check if a match is really a 100% match.
    if chrom==None:
        errAbort("No perfect match found in genome %s. Are you sure you have selected the right genome?" % genome)

    flankStart = start - 1000
    flankEnd = end + 1000

    if flankStart<0:
        errAbort("Not enough space on genome sequence to design primer. Please design it manually")

    flankFname = join(batchDir, batchId+".inFlank.fa")
    cmd = "twoBitToFa genomes/%(genome)s/%(genome)s.2bit %(flankFname)s -seq=%(chrom)s -start=%(flankStart)d -end=%(flankEnd)d" % locals()
    runCmd(cmd)

    flankSeq = parseFasta(flankFname).values()[0]
    tmpFname = join(batchDir, batchId+".primer3.in")
    tmpOutFname = join(batchDir, batchId+".primer3.out")
    # the guide seq has to be at least 150bp away from the left PCR primer for agarose gels
    lSeq, lTm, lPos, rSeq, rTm, rPos =  runPrimer3(flankSeq, tmpFname, tmpOutFname, 1000+guideStart-150, 330, "300-600")
    targetSeq = flankSeq[lPos:rPos+1]

    guideStart = targetSeq.find(guideSeqPlus)
    guideEnd = guideStart + len(guideSeqPlus)

    targetHtml = "<i>%s</i>" % targetSeq[:len(lSeq)] + " "+ targetSeq[len(lSeq):guideStart] +"<strong>"+ targetSeq[guideStart:guideEnd] + "</strong>"+ targetSeq[guideEnd:len(targetSeq)-len(rSeq)] + \
        " <i>%s</i>" % targetSeq[-len(rSeq):] 
    print '''<div style='width: 80%; margin-left:10%; margin-right:10%; text-align:left; word-wrap: break-word; word-break: break-all;'>'''
    print "<h3>Validation primers in genome around guide sequence %s strand %s</h3>" % (guideSeq, guideStrand)
    print "<strong>Left primer:</strong> %s (Tm: %s)<br>" % (lSeq, lTm)
    print "<strong>Right primer:</strong> %s (Tm: %s)<br>" % (rSeq, rTm)
    print "<strong>Wild-type genomic sequence including primers:</strong><br> <tt>%s</tt><br>" % targetHtml
    print "<strong>Sequence length:</strong> %d<br>" % (rPos-lPos)
    print "<hr>"
    print '<small>Primer3.2 with default settings, target length 400-600 bp</small>'
    print '</div>'

def main(seq, org, pam, pamId, batchId,defaultorg,defaultseq,defaultpam):
    printHeader(seq,org,pam)
    
    printSkeleton(seq, org, pam, pamId, batchId,defaultorg,defaultseq,defaultpam)    

    print("</body></html>")

def printContent(seq, org, pam, pamId, batchId,defaultorg,defaultseq,defaultpam):
    #seq, org, pam, pamId, batchId = getVars()    
        
    if pamId!=None:
        makePrimers(batchId, pamId)    
    elif batchId!=None or pam!=None:
        crisprSearch(seq, org, pam)        
        print '<br><a class="neutral" href="http://tefor.net/crispor/crisporDev.cgi"><div class="button" style="margin-left:auto;margin-right:auto;width:70px;">Back</div></a>'    
    elif seq==None:
        printForm(defaultorg,defaultseq,defaultpam)
    else:
        errAbort("Unrecognized CGI parameters.")

main(seq, org, pam, pamId, batchId,defaultorg,defaultseq,defaultpam)