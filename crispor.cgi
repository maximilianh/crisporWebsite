#!/usr/bin/env python
import subprocess
import Cookie, time, math
import cgitb

cgitb.enable()

# main script of the tefor crispr tool
# not broken into submodules to make installation easier

# cleaning things todo:
# - maybe: make temp subdirectories instead of many files with identical name

import sys, cgi, re, array, random, platform, os, hashlib, base64, string, logging, operator, urllib
from collections import defaultdict, namedtuple
from sys import stdout
from os.path import join, isfile, basename, dirname

DEBUG = False
#DEBUG = True

# the segments.bed files use abbreviated genomic region names
segTypeConv = {"ex":"exon", "in":"intron", "ig":"intergenic"}

# directory where processed batches of offtargets are stored ("cache" of bwa results)
batchDir = "temp"

DEFAULTORG = 'ensDanRer'
DEFAULTSEQ = 'CCAATCAGGTCCCTCCCTACCTCAGATCGCAGCTATAATACATAGGAGTAAAGAGGCTTCTCGCATTAAGTGGCTGTGGCTTGAAGTAACGTTGTGATTTCGAGGTCAGTCTTACCTTTCGCATCCCCGCCGCAAACCTCCGATGCGTTATCAGTCGCACGTTTCCGCACCTGTCACGGTCGGGGCTTGGCGCTGCTGAGGGACACGCGTGAACCGAGGAGACGGCAAGGACATCGCCGGAGATCCGCGCCTCGACAACGAGAAACCCTGCTAGACAGACCGCTCGAGAACACCGCAGCGAGATTCAGCGTGCGGCAAAATGCGGCTTTTGACGAGAGTGCTGCTGGTGTCTCTTCTCACTCTGTCCTTGGTGGTGTCCGGACTGGCCTGCGGTCCTGGCAGAGGCTACGGCAGAAGAAGACATCCGAAGAAGCTGACACCTCTCGCCTACAAGCAGTTCATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATAACGCGCAATTCGGAGAGATTTAAAGAACTTACTCCAAATTATAATCCCGACATTATCTTTAAGGATGAGGAGAACACGGGAGCGGACAGGCTCATGACACAG'
DEFAULTPAM = 'NGG'

def getParams():
    " get CGI parameters and return as dict "
    form = cgi.FieldStorage()
    params = {}

    for key in ["pamId", "batchId", "pam", "seq", "org", "showAll"]:
        val = form.getfirst(key)
	if val!=None:
            params[key] = val

    if "pam" in params:
        if len(set(params["pam"])-set("ACTGNMK"))!=0:
            errAbort("Illegal character in PAM-sequence. Only ACTGMK and N allowed.")
    return params

def makeTempBase(seq, org, pam):
    "create the base of temp files using a hash function and some prettyfication "
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
    if len(seq)>2000:
        msgs.append("<strong>Sorry, this tool cannot handle sequences longer than 2kbp</strong><br>Below you find the results for the first 2000 bp of your input sequence.<br>")
        seq = seq[:2000]

    if nCount!=0:
        msgs.append("Sequence contained %d non-ACTG letters. They were removed." % nCount)

    return seq, "<br>".join(msgs)

def revComp(seq):
    " rev-comp a dna sequence with UIPAC characters "
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

def showSeqAndPams(seq, lines, maxY, pam, guideScores):
    " show the sequence and the PAM sites underneath as a sequence viewer "
    print "<div class='substep'>"
    print '<a id="seqStart"></a>'
    print "Cas9 cuts at PAM (%s) sites. Click on a PAM match to show its guide sequence.<br>" % pam
    print '''Colors <span style="color:#32cd32; text-shadow: 1px 1px 1px #ddd">green</span>, <span style="color:#ffff00; text-shadow: 1px 1px 1px #ddd">yellow</span> and <span style="text-shadow: 1px 1px 1px #ddd; color:#aa0114">red</span> indicate high, medium and low specificity of the PAM's guide sequence in the genome.'''
    print "</div>"
    print '''<div style="text-align: left; overflow-x:scroll; width:100%; background:#DDDDDD; border-style: solid; border-width: 1px">'''

    print '<pre style="display:inline; line-height: 0.9em; text-align:left">'+rulerString(len(seq))
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
            texts.append('''<a style="text-shadow: 1px 1px 1px #bbb; color: %s" id="list%s" href="#%s" onmouseover="$('.hiddenExonMatch').show('fast');$('#show-more').hide();$('#show-less').show()" onfocus="window.location.href = '#seqStart'" >''' % (color, pamId,pamId))
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
        
#def calcMmScore(guideSeq, otSeq):
    #" return mismatch score for a given off target site "
    #guideSeq = guideSeq[:20]
    #otSeq = otSeq[:20]
    #score = 0
    #for i in range(len(guideSeq)-1, 0, -1):
        #if guideSeq[i]!=otSeq[i]:
            #score += i
    #return score

def makePosList(countDict, guideSeq, pam):
    """ for a given guide sequence, return a list of (score, posStr, geneDesc,
    otSeq) sorted by score and a string to describe the offtargets in the
    format x/y/z/w of mismatch counts
    Also return the same description for just the last 12 bp and the score 
    of the guide sequence (calculated using all offtargets).
    
    """
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
            score = calcHitScore(guideSeq[:20], otSeq[:20])
            scores.append(score)
            alnHtml, hasLast12Mm = makeAlnStr(guideSeq, otSeq, pam)
            if not hasLast12Mm:
                last12MmOtCount+=1
            posList.append( (score, editDist, posStr, geneDesc, alnHtml) )
        last12MmCounts.append(str(last12MmOtCount))

    guideScore = calcMitGuideScore(sum(scores))

    posList.sort(reverse=True)
    otDescStr = "/".join(otCounts)
    last12DescStr = "/".join(last12MmCounts)

    return posList, otDescStr, guideScore, last12DescStr

# --- START OF SCORING ROUTINES 

# DOENCH SCORING 
params = [
# pasted/typed table from PDF and converted to zero-based positions
(1,'G',-0.2753771),(2,'A',-0.3238875),(2,'C',0.17212887),(3,'C',-0.1006662),
(4,'C',-0.2018029),(4,'G',0.24595663),(5,'A',0.03644004),(5,'C',0.09837684),
(6,'C',-0.7411813),(6,'G',-0.3932644),(11,'A',-0.466099),(14,'A',0.08537695),
(14,'C',-0.013814),(15,'A',0.27262051),(15,'C',-0.1190226),(15,'T',-0.2859442),
(16,'A',0.09745459),(16,'G',-0.1755462),(17,'C',-0.3457955),(17,'G',-0.6780964),
(18,'A',0.22508903),(18,'C',-0.5077941),(19,'G',-0.4173736),(19,'T',-0.054307),
(20,'G',0.37989937),(20,'T',-0.0907126),(21,'C',0.05782332),(21,'T',-0.5305673),
(22,'T',-0.8770074),(23,'C',-0.8762358),(23,'G',0.27891626),(23,'T',-0.4031022),
(24,'A',-0.0773007),(24,'C',0.28793562),(24,'T',-0.2216372),(27,'G',-0.6890167),
(27,'T',0.11787758),(28,'C',-0.1604453),(29,'G',0.38634258),(1,'GT',-0.6257787),
(4,'GC',0.30004332),(5,'AA',-0.8348362),(5,'TA',0.76062777),(6,'GG',-0.4908167),
(11,'GG',-1.5169074),(11,'TA',0.7092612),(11,'TC',0.49629861),(11,'TT',-0.5868739),
(12,'GG',-0.3345637),(13,'GA',0.76384993),(13,'GC',-0.5370252),(16,'TG',-0.7981461),
(18,'GG',-0.6668087),(18,'TC',0.35318325),(19,'CC',0.74807209),(19,'TG',-0.3672668),
(20,'AC',0.56820913),(20,'CG',0.32907207),(20,'GA',-0.8364568),(20,'GG',-0.7822076),
(21,'TC',-1.029693),(22,'CG',0.85619782),(22,'CT',-0.4632077),(23,'AA',-0.5794924),
(23,'AG',0.64907554),(24,'AG',-0.0773007),(24,'CG',0.28793562),(24,'TG',-0.2216372),
(26,'GT',0.11787758),(28,'GG',-0.69774)]

intercept =  0.59763615
gcHigh    = -0.1665878
gcLow     = -0.2026259

def calcDoenchScore(seq):
    assert(len(seq)==30)
    score = intercept

    guideSeq = seq[4:24]
    gcCount = guideSeq.count("G") + guideSeq.count("C")
    if gcCount <= 10:
        gcWeight = gcLow
    if gcCount > 10:
        gcWeight = gcHigh
    score += abs(10-gcCount)*gcWeight

    for pos, modelSeq, weight in params:
        subSeq = seq[pos:pos+len(modelSeq)]
        if subSeq==modelSeq:
            score += weight
    return 1.0/(1.0+math.exp(-score))

# MIT offtarget scoring

def calcHitScore(string1,string2):
    " see 'Scores of single hits' on http://crispr.mit.edu/about "
    # The Patrick Hsu weighting scheme
    M = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]
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
            score1 *= 1-M[pos]
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
    " Sguide defined on http://crispr.mit.edu/about "
    score = 100 / (100+hitSum)
    score = int(round(score*100))
    return score

# --- END OF SCORING ROUTINES 

def calcDoenchScoreFromSeqPos(startPos, seq, pamLen, strand):
    " extract 30 mer from seq given beginning of pam at startPos and strand, return Doench score "
    if strand=="+":
        fromPos  = startPos-24
        toPos    = startPos+6
        func     = str
    else: # strand is minus
        fromPos = startPos+pamLen-6
        toPos   = startPos+pamLen+24
        func    = revComp

    if fromPos < 0 or toPos > len(seq):
        return None
    else:
        seq30Mer = func(seq[fromPos:toPos]) 
        return int(round(100*calcDoenchScore(seq30Mer)))

def htmlHelp(text):
    " show help text with tooltip or modal dialog "
    print '''<img src="image/info.png" class="help tooltip" title="%s" />''' % text

def scoreGuides(seq, startDict, pamPat, otMatches):
    " for each pam in startDict, retrieve the guide sequence next to it and score it "
    guideData = []
    guideScores = {}
    hasNotFound = False

    for startPos, strand, flankSeq, pamSeq in flankSeqIter(seq, startDict, len(pamPat)):
        # position with anchor to jump to
        pamId = "s"+str(startPos)+strand
        # flank seq
        seqStr = "<tt>"+flankSeq + " <i>" + pamSeq+"</i></tt>"
        guideSeq = flankSeq + pamSeq

        # matches in genome
        # one desc in last column per OT seq
        if pamId in otMatches:
            posList, otDesc, guideScore, last12Desc = makePosList(otMatches[pamId], guideSeq, pamPat)
            effScore = calcDoenchScoreFromSeqPos(startPos, seq, len(pamPat), strand)
            #print otDesc, "<p>"
        else:
            posList, otDesc, guideScore = None, "No match. Incorrect genome or cDNA?", 0
            last12Desc = ""
            effScore = 0
            hasNotFound = True
            #print "no match"
        guideData.append(( guideScore, effScore, startPos, strand, pamId, seqStr, guideSeq, posList, otDesc, last12Desc))
        guideScores[pamId] = guideScore

    guideData.sort(reverse=True)
    return guideData, guideScores, hasNotFound

def printTableHead():
    " print guide score table description and columns "
    # one row per guide sequence
    print "<div class='substep'>Ranked from highest to lowest specificity score determined as in <a target='_blank' href='http://dx.doi.org/10.1038/nbt.2647'>Hsu et al.</a> and on crispr.mit.org<br></div>"
    print '<table id="otTable">'
    print '<tr style="border-left:5px solid black">'
    
    print '<th>Position/<br>Strand'
    print '</th>'

    print '<th style="width:170px">Guide Sequence + <i>PAM</i></th>'

    print '<th>Specificity Score'
    htmlHelp("The specificity score measures the uniqueness of a guide in the genome. &lt;br&gt;The higher the specificity score, the less likely are off-target effects.")
    print "</th>"

    print '<th>Efficacy Score'
    htmlHelp("The efficacy score measures the cutting efficiency of the nuclease on a sequence. &lt;br&gt; The higher the efficacy score, the more likely is cutting.")
    print "</th>"

    print '<th style="">Matches in genome for 0,1,2,3,4 mismatches.</th>'

    print '<th>No mismatches adjacent to PAM.</i>'
    htmlHelp("Like the previous column, but only counting the last 12bp of the guide sequence, those closest to the PAM.")
    print "</th>"

    print '<th>Genome Browser links to matches</th>'

def scoreToColor(guideScore):
    if guideScore > 50:
        color = "#32cd32"
    elif guideScore > 20:
        color = "#ffff00"
    else:
        color = "#aa0114"
    return color

def showGuideTable(guideData, pam, otMatches, dbInfo, batchId, org, showAll):
    " shows table of all PAM motif matches "
    print "<br><div class='title'>Predicted guide sequences for PAMs</div>" 
    printTableHead()

    count = 0
    for guideRow in guideData:
        guideScore, effScore, startPos, strand, pamId, seqStr, guideSeq, posList, otDesc, last12Desc = guideRow

        color = scoreToColor(guideScore)
        print '<tr id="%s" style="border-left: 5px solid %s">' % (pamId, color)

        # position and strand
        #print '<td id="%s">' % pamId
        print '<td>' 
        print '<a href="#list%s">' % (pamId)
        print str(startPos)+" /"
        if strand=="+":
            print 'fw'
        else:
            print 'rev'
        print '</a>'
        print "</td>"

        # sequence
        print "<td>"
        print "<small>"
        print seqStr
        print "<br>"

        scriptName = basename(__file__)
        if posList!=None:
            print '<br><a href="%s?batchId=%s&pamId=%s&pam=%s">Details</a>' % (scriptName, batchId, urllib.quote(str(pamId)), pam)
        print "</small>"
        print "</td>"

        # off-target score, aka specificity score aka MIT score
        print "<td>"
        print "%d" % guideScore
        print "</td>"

        # efficacy score
        print "<td>"
        if effScore==None:
            print "Too close to end"
        else:
            print '''%d''' % effScore
        #print '''<a href="#" onclick="alert('%s')">%0.2f</a>''' % (effScore)
        #print "<!-- %s -->" % seq30Mer
        print "</td>"

        # mismatch description
        print "<td>"
        print otDesc
        if posList==None:
            # no genome match
            htmlHelp("Sequence was not found in genome.<br>If you have pasted a cDNA sequence, note that sequences that overlap a splice site cannot be used as guide sequences<br>This warning also occurs if you have selected the wrong genome.")
        print "</td>"

        # mismatch description, last 12 bp
        print "<td>"
        print last12Desc
        print "</td>"

        # links to offtargets
        print "<td><small>"
        if posList!=None:
            i = 0
            for score, editDist, pos, gene, alnHtml in posList:
                print '''%d:''' % (int(editDist))
                printBrowserLink(dbInfo, pos, gene, alnHtml)
                i+=1
                if i==3 and not showAll:
                    break
            if not showAll and len(posList)>3:
                 print '''... <br>&nbsp;&nbsp;&nbsp;<a href="crispor.cgi?batchId=%s&showAll=1">- show all offtargets</a>''' % batchId

        print "</small></td>"

        print "</tr>"
        count = count+1

    print "</table>"

    #print '''<a style="text-align:right" href="http://tefor.net/crispor/download.php?batchId=%s&amp;seq=%s&amp;org=%s&amp;pam=%s&amp;pamId=%s">
                    #<img style="width:20px;vertical-align:middle"
                         #src="http://tefor.net/crispor/image/doc.png">
                    #Download results
                #<!--</div>-->
            #</a>
            #<br><br>
    #''' % (batchId,seq,org,pam,pamId)


def printHeader(batchId):
    " print the html header "

    print "<html><head>"   

    runPhp("header.php")
    runPhp("/var/www/main/specific/googleanalytics/script.php")

    # activate jqueryUI tooltips
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

    print '<link rel="stylesheet" type="text/css" href="style/tooltipster.css" />'
    print '<link rel="stylesheet" type="text/css" href="style/tooltipster-shadow.css" />'
    print '<script type="text/javascript" src="js/jquery.tooltipster.min.js"></script>'
    # activate tooltipster
    #print (""" <script> $(document).ready(function() { $('.tooltip').tooltipster(); }); </script> """)
       #theme: 'tooltipster-shadow',
    print (""" <script> $(document).ready(function() { $('.tooltip').tooltipster({ 
        contentAsHTML: true,
       speed : 0
        }); }); </script> """)

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
    
    if batchId is not None:
        print '''
        onload="history.pushState('crispor/crisporDev.cgi', document.title, '?batchId=%s');"
        ''' % (batchId)
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
    
    return as dict pamId -> editDist -> (chrom, start, end, seq, strand, segType, segName)
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

def printForm(params):
    " print html input form "
    #seq, org, pam = params["seq"], params["org"], params["pam"]
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
 CRISPOR (CRISPr selectOR) is a program that helps design and evaluate target sites for use with the CRISPR/Cas9 system.
    <div onclick="$('#about-us').toggle('fast');" class="title" style="cursor:pointer;display:inline;font-size:large;font-style:normal">
        <img src="http://tefor.net/crispor/image/info.png" class="infopoint" style="vertical-align:text-top;">
    </div>
    <div id="about-us"><br>
    CRISPOR uses the BWA algorithm to identify guide RNA sequences for CRISPR mediated genome editing.<br>
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

def batchParams(batchId):
    " given a batchId, return the genome, the pam and the input sequence "
    batchBase = join(batchDir, batchId)
    faFname = batchBase+".input.fa"
    ifh = open(faFname)
    genome, pamSeq = ifh.readline().replace(">","").strip().split()
    inSeq = ifh.readline().strip()
    ifh.close()
    return inSeq, genome, pamSeq

def crisprSearch(params):
    " do crispr off target search "
    if "batchId" in params:
        # if we're getting only the batchId, extract the parameters from the batch
        # this allows a stable link to a batch that is done
        seq, org, pam = batchParams(params["batchId"])
        batchId = params["batchId"]
    else:
        seq, org, pam = params["seq"], params["org"], params["pam"]
        batchId = makeTempBase(seq, org, pam)

    showAll = (params.get("showAll", 0)=="1")

    # read genome info tab file
    myDir = dirname(__file__)
    genomesDir = join(myDir, "genomes")
    infoFname = join(genomesDir, org, "genomeInfo.tab")
    dbInfo = lineFileNext(open(infoFname)).next()

    caseSeq, userMsg = cleanSeq(seq)
    seq = caseSeq.upper()
    print userMsg

    # search pams
    startDict, endSet = findPams(seq, pam, "+", {}, set())
    startDict, endSet = findPams(seq, revComp(pam), "-", startDict, endSet)

    # define temp file names
    batchBase = join(batchDir, batchId)
    #print "<!-- BATCH ID %s -->" % batchId

    # save input seq, pamSeq and genome for primer design later
    inputFaFname = batchBase+".input.fa"
    open(inputFaFname, "w").write(">%s %s\n%s\n" % (org, pam, seq))

    # write guides to fasta and run bwa
    faFname = batchBase+".fa"
    otBedFname = batchBase+".bed"
    flagFile = batchBase+".running"

    if isfile(flagFile):
       errAbort("This sequence is still being processed. Please wait for ~20 seconds "
           "and try again. If you see this message for more than 1-2 minutes, please "
           "email penigault@tefor.net")

    if not isfile(otBedFname):
        # write potential PAM sites to file
        writePamFlank(seq, startDict, pam, faFname)
        findOfftargets(faFname, org, pam, otBedFname)

    print "<div class='title'><div class='speciesname'>%s</div> sequence with PAMs</div>" % dbInfo.scientificName

    otMatches = parseOfftargets(otBedFname)
    featLen = len(pam)
    lines, maxY = distrOnLines(seq, startDict, featLen)

    guideData, guideScores, hasNotFound = scoreGuides(seq, startDict, pam, otMatches)

    if hasNotFound:
        print("<strong>Warning:</strong> At least one of the PAM-flanking sequences was not found in the genome.<br>")
        print("Did you select the right genome? <br>")
        print("If you pasted a cDNA sequence, note that sequences with score 0 are not in the genome, only in the cDNA and are not usable as CRISPR guides.<br>")

    showSeqAndPams(caseSeq, lines, maxY, pam, guideScores)

    showGuideTable(guideData, pam, otMatches, dbInfo, batchId, org, showAll)

    # XX are back buttons necessary in 2014?
    print '<br><a class="neutral" href="http://tefor.net/crispor">'
    print '<div class="button" style="margin-left:auto;margin-right:auto;width:90;">New Query</div></a>'

def runPhp(script):
    " run a file through php, write result to stdout. accepts a full or a relative path "
    if "/" in script:
        path = script
    else:
        myDir = dirname(__file__)
        path = "%s/%s" % (myDir, script)

    proc = subprocess.Popen("php "+path, shell=True, stdout=subprocess.PIPE)
    script_response = proc.stdout.read()
    print script_response

def printTeforBodyStart():
    print "<div class='logo'><a href='http://tefor.net/main/'><img src='http://tefor.net/main/images/logo/logo_tefor.png' alt='logo tefor infrastructure'></a></div>"
    runPhp("menu.php")

    print '<div id="bd">'
    print '<div class="centralpanel">'
    runPhp("networking.php")
    print '<div class="subpanel" style="background:transparent;box-shadow:none;">'
    print '<div class="contentcentral" style="background-color:transparent;">'

def printTeforBodyEnd():
    print '</div>'
    print '</div>'
    print '</div>'
    runPhp("footer.php")
    print '</div>'

def printBody(params):
    " main dispatcher function "
    if len(params)==0:
        printForm(params)
    elif "batchId" in params and "pamId" in params and "pam" in params:
        primerDetailsPage(params)
    elif ("seq" in params and "org" in params and "pam" in params) \
                or "batchId" in params:
        crisprSearch(params)
    else:
        errAbort("Unrecognized CGI parameters.")

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

def findBestMatch(genome, batchId):
    " find best match for input sequence from batchId in genome and return as (chrom, start, end) "

    batchBase = join(batchDir, batchId)
    samFname = batchBase+".input.sam"
    faFname = batchBase+".input.fa"

    cmd = "BIN/bwa bwasw genomes/%(genome)s/%(genome)s.fa %(faFname)s > %(samFname)s" % locals()
    runCmd(cmd)

    chrom, start, end = None, None, None
    for l in open(samFname):
        if l.startswith("@"):
            continue
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
            errAbort("Sequence not found in genome. Are you sure you have pasted the correct sequence and also selected the right genome?")
        matchLen = int(cigar.replace("M","").replace("S", ""))
        # XX why do we get soft clipped sequences from BWA? repeats?
        chrom, start, end =  rName, int(pos), int(pos)+matchLen
        print chrom, start, end, strand

    # possible problem: we do not check if a match is really a 100% match.
    if chrom==None:
        errAbort("No perfect match found in genome %s." \
            "Are you sure you have selected the right genome?" % genome)
    else:
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

    flankSeq = parseFasta(flankFname).values()[0]
    tmpFname = join(batchDir, batchId+".primer3.in")
    tmpOutFname = join(batchDir, batchId+".primer3.out")
    # the guide seq has to be at least 150bp away from the left PCR primer for agarose gels
    lSeq, lTm, lPos, rSeq, rTm, rPos = runPrimer3(flankSeq, tmpFname, tmpOutFname, 1000+guideStart-150, 330, "300-600")
    targetSeq = flankSeq[lPos:rPos+1]
    return lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq

def primerDetailsPage(params):
    """ create primers with primer3 around site identified by pamId in batch
    with batchId. Output primers as html
    """
    batchId, pamId, pam = params["batchId"], params["pamId"], params["pam"]
    inSeq, genome, pamSeq = batchParams(batchId)
    seqLen = len(inSeq)
    batchBase = join(batchDir, batchId)

    # find position of guide sequence in genome at MM0
    otBedFname = batchBase+".bed"
    otMatches = parseOfftargets(otBedFname)
    if pamId not in otMatches or 0 not in otMatches[pamId]:
        errAbort("No perfect match found for guide sequence in the genome. Are you sure you have selected the right genome? If you have selected the right genome but pasted a cDNA, please note that sequences that overlap a splice site are not part of the genome and cannot be used as guide sequences.")

    matchList = otMatches[pamId][0] # = get all matches with 0 mismatches
    if len(matchList)!=1:
        errAbort("Multiple perfect matches for this guide sequence. Cannot design primer. Please select another guide sequences or email penigault@tefor.net to discuss your strategy or modifications to this software.")
        # XX we could show a dialog: which match do you want to design primers for?
        # But who would want to use a guide sequence that is not unique?

    chrom, start, end, seq, strand, segType, segName = matchList[0]
    start = int(start)
    end = int(end)

    # retrieve guideSeq + PAM sequence from input sequence
    # XX guide sequence must not appear twice in there
    pamFname = batchBase+".fa"
    pams = parseFasta(pamFname)
    guideSeq = pams[pamId]
    guideStrand = pamId[-1]
    guideSeqWPam = seq

    if strand=="+":
        guideStart = inSeq.find(guideSeq)
        highlightSeq = guideSeqWPam
    else:
        guideStart = inSeq.find(revComp(guideSeq))
        highlightSeq = revComp(guideSeqWPam)

    # matchChrom, matchStart, matchEnd, matchStrand = findBestMatch(genome, batchId)

    lSeq, lTm, lPos, rSeq, rTm, rPos, targetSeq = \
        designPrimer(genome, chrom, start, end, strand, 0, batchId)

    guideStart = targetSeq.find(highlightSeq)
    guideEnd = guideStart + len(highlightSeq)

    #if guideStrand=="+":
        #guideStrand = "forward"
    #else:
        #guideStrand = "reverse"

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
    print "<h3>Genome fragment flanking guide sequence %s</h3>" % (guideSeqHtml)

    if strand=="-":
        print("Your guide sequence is on the reverse strand relative to the genome sequence, so it is reverse complemented in the sequence below.<p>")

    print "<strong>Left primer:</strong> %s (Tm: %s)<br>" % (lSeq, lTm)
    print "<strong>Right primer:</strong> %s (Tm: %s)<br>" % (rSeq, rTm)
    print '''<div style='word-wrap: break-word; word-break: break-all;'>'''
    print "<strong>Genomic sequence %s:%d-%d including primers, forward strand:</strong><br> <tt>%s</tt><br>" % (chromLong, start, end, targetHtml)
    print '''</div>'''
    print "<strong>Sequence length:</strong> %d<br>" % (rPos-lPos)
    print "<hr>"

    print "<h4>Expression of guide RNA by in vitro transcription with T7 RNA polymerase from plasmid DNA</h4>"
    print 'To produce guide RNA by in vitro transcription with T7 RNA polymerase, the guide RNA sequence can be cloned in a variety of plasmids (see <a href="http://addgene.org/crispr/empty-grna-vectors/">AddGene website</a>).<br>'
    print "For the guide sequences %s, the following primers should be ordered for cloning into the DR274 plasmid generated by the Joung lab<p>" % guideSeqWPam
    if guideSeq.lower().startswith("gg"):
        guideRnaFw = "TA %s" % guideSeq
        guideRnaRv = "AAAC %s" % revComp(guideSeq[2:])
    else:
        guideRnaFw = "TAGG %s" % guideSeq
        guideRnaRv = "AAAC %s" % revComp(guideSeq)
    print "guideRnaFw %s<br>" % guideRnaFw
    print "guideRnaRev %s<br>" % guideRnaRv

    print "<h4>Expression of guide RNA by in vitro transcription with T7 polymerase using overlapping oligonucleotides</h4>"
    print "Template for in vitro synthesis of guide RNA with T7 RNA polymerase can be prepared by annealing and primer extension of the following primers:<p>"
    print "constant primer used for all guide RNAs: AAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC<br>"

    prefix = ""
    if not guideSeq.lower().startswith("gg"):
        prefix = "GG"
    specPrimer = "TAATACGACTCACTATA-%s<strong>%s</strong>-GTTTTAGAGCTAGAAATAGCAAG" % (prefix, guideSeq)

    print "guide RNA specific primer with T7 promoter sequence: %s<br>" % specPrimer
    print 'The protocol for template preparation from oligonucleotides and in vitro transcription can be found in Gagnon et al. PLoS One 2014 May 29. <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4038517/?report=classic">Efficient mutagenesis by Cas9 protein-mediated oligonucleotide insertion and large-scale assessment of single-guide RNAs</a>'

    print "<h4>Expression of guide RNA in mammalian cells from plasmid </h4>"
    if "tttt" in guideSeq.lower():
        print "The guide sequence %s contains the motif TTTT, which terminates RNA polymerase. This guide sequence cannot be transcribed in mammalian cells." % guideSeq
    else:
        print "The guide sequence %s does not contain the motif TTTT, which terminates RNA polymerase. This guide sequence can be transcribed in mammalian cells." % guideSeq

    print "<br>"
    print "To express guide RNA in mammalian cells, a variety of plasmids are available. For example, to clone the guide RNA sequence in the plasmid pM3636, where guide RNA expression is driven by a human U6 promoter, the following primers should be used :"
    print "<br>"
    if guideSeq.lower().startswith("g"):
        guideRnaFw = "ACACC %s G" % guideSeq
        guideRnaRv = "AAAAC %s G" % revComp(guideSeq)
    else:
        guideRnaFw = "TAGG G %s" % guideSeq
        guideRnaRv = "AAAC %s CG" % revComp(guideSeq)
    print "guideRnaFw %s<br>" % guideRnaFw
    print "guideRnaRev %s<br>" % guideRnaRv

    print "<h4>Expression of guide RNA in Drosophila cells from plasmid </h4>"
    print "TODO?"
    print 
    print ""
    print "<hr>"
    print '<small>Primer3.2 with default settings, target length 400-600 bp</small>'
    print '</div>'

def main():
    # parse incoming parameters
    params = getParams()
    batchId = None

    # print headers
    # save seq/org/pam into a cookie, if they were provided
    if "seq" in params and "org" in params and "pam" in params:
        seq, org, pam = params["seq"], params["org"], params["pam"]
        saveSeqOrgPamToCookies(seq, org, pam)
        batchId = makeTempBase(seq, org, pam)
    print "Content-type: text/html\n"
    print "" # = end of http headers

    printHeader(batchId)
    printTeforBodyStart()

    printBody(params)     # main dispatcher, branches based on the params dictionary

    printTeforBodyEnd()
    print("</body></html>")

main()
