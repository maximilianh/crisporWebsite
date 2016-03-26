import sys
import glob
#import sets
import re

def openSpec(fname, mode="r"):
    """ open and return filehandle, open stdin if fname=="stdin", do nothing if none """
    if fname=="stdin":
        return sys.stdin
    elif fname=="stdout":
        return sys.stdout
    elif fname=="none" or fname==None:
        return None
    else:
        return open(fname, mode)

def writeList(fname, list):
    of = openSpec(fname, "w")
    for row in list:
        row = [str(d) for d in row]
        of.write("\t".join(row))
        of.write("\n")
    of.close()
        
def slurpdict(fname, comments=False, valField=1, doNotCheckLen=False, asFloat=False, otherFields=False, asInt=False, headers=False, keyAsInt=False):
    """ parse file with key -> value pair on each line, key/value has 1:1 relationship"""
    """ last field: set valField==-1, return as a dictionary key -> value """
    if fname==None or fname=="":
        return {}
    dict = {}
    f = openSpec(fname)
    if not f:
        return dict

    if headers:
        headers = f.readline()
    for l in f:
        fs = l.strip().split("\t")
        if comments and l.startswith("#"):
            continue
        if not len(fs)>1:
            if not doNotCheckLen:
                sys.stderr.write("info: not enough fields, ignoring line %s\n" % l)
                continue
            else:
                key = fs[0]
                val = None
        else:
            key = fs[0]

            if keyAsInt:
                key = int(key)

            if not otherFields:
                val = fs[valField]
            else:
                val = fs[1:]
            
            if asFloat:
                val = float(val)
            elif asInt:
                val = int(val)
        if key not in dict:
            dict[key] = val
        else:
            sys.stderr.write("info: file %s, hit key %s two times: %s -> %s\n" %(fname, key, key, val))
    return dict

def slurpdictlist(fname, reverse=False, filterComments=False, keyType=None, valType=None):
    """ parse file with key -> value pair on each line and return as dict -> list (1:n relationship) """
    if fname==None:
        return {}
    dict = {}
    if fname=="stdin":
        fh = sys.stdin
    else:
        fh = open(fname, "r")

    for l in fh:
        if filterComments and l.startswith("#"):
            continue
        fs = l.strip().split("\t")
        if len(fs)>1:
            # reverse?
            if reverse:
                val = fs[0]
                key = fs[1]
            else:
                key = fs[0]
                val = fs[1]
            # convert to specified types
            if keyType:
                key = keyType(key)
            if valType:
                val = valType(val)
            dict.setdefault(key, []).append(val)
    return dict

def slurpdictset(fname, reverse=False, keyType=None, valType=None):
    """ parse file with key -> value pair on each line and return as dict -> list (1:n relationship). Keytype can be e.g. types.IntType """
    if fname==None or fname=="":
        return {}
    dict = {}
    if fname=="stdin":
        fh = sys.stdin
    else:
        fh = open(fname, "r")

    for l in fh:
        fs = l.strip().split("\t")
        if len(fs)>1:
            if not reverse:
                key = fs[0]
                val = fs[1]
            else:
                val = fs[0]
                key = fs[1]

            # convert to specified types
            if keyType:
                key = keyType(key)
            if valType:
                val = valType(val)
            # add to set
            dict.setdefault(key, set()).add(val)
    return dict

def slurplist(fname, check=True, field=None, filterComments=False, valType=None, headers=False):
    """ parse a file with one string per line and return as list"""
    if fname==None:
        return []
    if fname=="stdin":
        fh = sys.stdin
    else:
        fh = open(fname, "r")
    list = []
    if headers:
        fh.readline()
    for l in fh:
        l = l.strip()
        if filterComments and l.startswith("#"):
            continue
        if len(l)==0:
            continue

        if check and l in list:
            sys.stderr.write("tabfile.py/slurplist: file=%s, duplicate key = %s, exiting\n" % (fname, l))
            sys.exit(1)

        if field==None:
            value = l
        else:
            value = l.split()[field]

        if valType:
            value = valType(value)

        list.append(value)
    return list
        
def slurplistasdict(fname, split=False, default=True):
    """ parse a file with one string per line and return as dict for faster access"""
    dict = {}
    for l in open(fname, "r"):
        l = l.strip()
        if split:
            l = l.split("\t")[0]
        if l in dict:
            sys.stderr.write("tabfile.py: key already in dict!\n")
            return None
        dict[l.strip()] = default
    return dict

def slurpdictlistlist(fname):
    """ parse file with key -> values pair on each line, many values per line, values are returned as list"""
    if fname==None:
        return []
    dict = {}
    for l in open(fname, "r"):
        if l.startswith("#"):
                continue
        fs = l.strip().split("\t")
        if len(fs)>1:
                key = fs[0]
                vals = fs[1:]
                dict.setdefault(key, []).append(vals)
    return dict

def parseTsv(fname, columnNames = None, asListOfDicts=False):
    """ retrieve only selected fields from tsv files with headers (-> like R dataframes) as a list of lists  or list of dicts"""
    """ returns a tuple (comments, headers, data) """
    comments, headers, data = [], [], []

    f = openSpec(fname)

    # parse comments and headers
    firstLine = f.readline().strip()
    #while firstLine.startswith("#"):
        #comments.append(firstLine)
        #firstLine = f.readline().strip()
    headers = firstLine.split("\t")

    if columnNames==None or len(columnNames)==0:
        columnNames=headers
    for c in columnNames:
        if c not in headers:
            sys.stderr.write("error tabfile.py: columnName %s (out of %s) not found in headers %s\n" % (c,str(columnNames), str(headers)))
            sys.exit(1)
    noList = range(0, len(headers))
    headerToNum = dict(zip(headers,noList))

    data = []
    #lno=0
    for l in f:
        if l.startswith("#"):
            continue
        #lno+=1
        #if lno==1: # ignore headers
            #continue
        fs = l.strip().split("\t")
        if asListOfDicts:
            rec = {}
        else:
            rec = []
        for c in columnNames:
            if asListOfDicts:
                rec[c]=fs[headerToNum[c]]
            else:
                rec.append(fs[headerToNum[c]])
        data.append(rec)

    return comments, headers, data

############### PSL FILES ########################
class Psl:
    def toList(self, str):
        return [int(x) for x in str.split(",") if x.strip()!=""]

    def __init__(self, line):
        (match, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, tBaseInsert, strand, qName, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, blockCount, blockSizes, qStarts, tStarts) = line.split("\t")

        self.match = int(match)
        self.misMatches = int(misMatches)
        self.repMatches = int(repMatches)
        self.nCount = int(nCount)
        self.qNumInsert = int(qNumInsert)
        self.qBaseInsert = int(qBaseInsert)
        self.tNumInsert = int(tNumInsert)
        self.tBaseInsert = int(tBaseInsert)
        self.strand = strand
        self.qName = qName
        self.qSize = int(qSize)
        self.qStart = int(qStart)
        self.qEnd = int(qEnd)
        self.tName = tName
        self.tSize = int(tSize)
        self.tStart = int(tStart)
        self.tEnd = int(tEnd)
        self.blockCount = int(blockCount)
        self.blockSizes = self.toList(blockSizes)
        self.qStarts = self.toList(qStarts)
        self.tStarts = self.toList(tStarts)

    def getQueryBlocks(self):
        regions = []
        if self.strand[0]=="+":
            for start, len in zip(self.qStarts, self.blockSizes):
                end = start+len
                regions.append( (start,end) )
        else:
            sys.stderr("regions on negative strand not supported in parser\n")
            sys.exit(1)
        return regions
                
def parsePsl(fname):
    psls = []
    for l in open(fname, "r"):
        psl = Psl(l)
        psls.append(psl)
    return psls


######################## BLAST TAB FORMAT -m8 -m9 ####################
class hit:
    """ a hit from a blast tabular output file """
    def __repr__(self):
        return self.line

def parseBlast(fname):
    hits = []
    for l in open(fname, "r"):
        if l.startswith("#"):
            continue
        fs = l.split("\t")
        if len(fs)!=12:
            print l
            print len(fs)
            sys.stderr.write("error: blast output does not have 12 fields. Make sure that you ran blast with -m 9 or -D T (for bl2seq)\n")
            sys.exit(1)
        h = hit()
        (h.qName, h.sName, h.percId, h.alnLen, h.mismatches, h.gapOpen, h.qStart, h.qEnd, h.sStart, h.sEnd, h.eVal, h.alnScore) = fs
        h.sStart = int(h.sStart)
        h.sEnd = int(h.sEnd)
        h.qStart = int(h.qStart)
        h.qEnd = int(h.qEnd)
        if h.qEnd < h.qStart:
            (h.qStart, h.qEnd) = (h.qEnd, h.qStart)
        if h.sEnd < h.sStart:
            (h.sStart, h.sEnd) = (h.sEnd, h.sStart)
        h.sStart -= 1
        h.qStart -= 1
        h.line = l.strip()
        hits.append(h)
    return hits

# ################# PARSE INPARANOID ############################################

def openParseInparanoidTable(dir, srcOrg, targetOrg, only1To1Homologs=False, reverse=False):
    def unzip(list):
        ids = []
        scores = []
        for i in range(0, len(list), 2):
           sid = list[i]
           ids.append(sid)
           #print list, i, sid, ids
           score = list[i+1]
           scores.append(score)
        return (ids,scores)
    
    fmask = dir + "/table.%s*%s*" % (srcOrg, targetOrg)
    files = glob.glob(fmask)
    if not reverse and len(files)==0:
        return openParseInparanoidTable(dir, targetOrg, srcOrg, only1To1Homologs, reverse=True)
    if reverse and len(files)==0:
        sys.stderr.write("error: Could not find a file that matches %s\n" % (fmask))
        sys.exit(1)
    if len(files)>1:
        sys.stderr.write("error: expression %s resolved to more than one file: %s\n" % (fmask, str(files)))
        sys.exit(1)
    
    sys.stderr.write("Reading %s...\n" % files[0])
    f = open(files[0], "r")
    f.readline() # skip header
    srcToTargetList={}
    #srcToClust={}
    targetParalogs={}
    srcParalogs = {}
    for l in f:
        fs = l.strip().split("\t")
        clust, score, srcIds, targetIds = fs
        if reverse:
            srcIds, targetIds = targetIds, srcIds
        #clust = int(clust)
        clust = "InpCluster"+clust
        srcIds, srcScores = unzip(srcIds.strip().split(" "))
        targetIds, targetScores = unzip(targetIds.strip().split(" "))
        if only1To1Homologs:
            if len(srcIds)>1 or len(targetIds)>1:
                #sys.stderr.write("info: dropping inparanoid cluster %s as it is not 1:1\n" % str(srcIds)+"-"+str(targetIds)+"-"+str(clust))
                continue
        for srcId in srcIds:
            #print srcId
            srcToTargetList[srcId]=zip(targetIds, [clust]*len(targetIds))
            #if srcId in srcToClust:
                #sys.stderr.write("program error: gene on source with TWO clusters??\n")
                #sys.exit(1)
            #srcToClust[srcId]=clust
        if len(targetIds) > 1:
            for targetId in targetIds:
                targetParalogs[targetId]=targetIds
        if len(srcIds) > 1:
            for srcId in srcIds:
                srcParalogs[srcId]=srcIds
    return srcToTargetList, srcParalogs, targetParalogs
        
def parsePhylipMatrix(filename):
    regex = re.compile("[0-9-]")
    f = open(filename, "r")
    data = {}
    f.readline() # consume first line
    fs = ("".join(f.readlines())).split()
    orgs = []
    for f in fs:
        if not regex.match(f):
            data[f]=[]
            orgs.append(f)
            last=f
        else:
            data[last].append(f)

    newdata = {}
    for org, values in data.iteritems():
        dict = {}
        i = 0
        for v in values:
            dict[orgs[i]]=v
            i+=1
        newdata[org]=dict

    return newdata
            
# new routine, not compat with pairCrawler
def parseInparanoid(dir, srcOrg, targetOrg, only1To1Homologs=False, reverse=False):
    def unzip(list):
        ids = []
        scores = []
        for i in range(0, len(list), 2):
           sid = list[i]
           ids.append(sid)
           #print list, i, sid, ids
           score = list[i+1]
           scores.append(score)
        return (ids,scores)

    class Data:
        pass
    
    fmask = dir + "/table.%s*%s*" % (srcOrg, targetOrg)
    files = glob.glob(fmask)
    if not reverse and len(files)==0: # use the reverse species order if other order not found
        return parseInparanoid(dir, targetOrg, srcOrg, only1To1Homologs, reverse=True)
    if reverse and len(files)==0:
        sys.stderr.write("error: Could not find a file that matches %s\n" % (fmask))
        sys.exit(1)
    if len(files)>1:
        sys.stderr.write("error: expression %s resolved to more than one file: %s\n" % (fmask, str(files)))
        sys.exit(1)
    
    sys.stderr.write("Reading %s...\n" % files[0])
    f = open(files[0], "r")
    f.readline() # skip header
    srcToTargetList={}
    #srcToClust={}
    targetParalogs={}
    srcParalogs = {}
    trgIdSet = sets.Set()
    for l in f:
        fs = l.strip().split("\t")
        clust, score, srcIds, targetIds = fs
        if reverse:
            srcIds, targetIds = targetIds, srcIds
        #clust = int(clust)
        clust = "InpCluster"+clust
        srcIds, srcScores = unzip(srcIds.strip().split(" "))
        targetIds, targetScores = unzip(targetIds.strip().split(" "))
        if only1To1Homologs:
            if len(srcIds)>1 or len(targetIds)>1:
                #sys.stderr.write("info: dropping inparanoid cluster %s as it is not 1:1\n" % str(srcIds)+"-"+str(targetIds)+"-"+str(clust))
                continue

        for targetId in targetIds:
            trgIdSet.add(targetId)
        for srcId in srcIds:
            #print srcId
            srcToTargetList[srcId]=zip(targetIds, [clust]*len(targetIds))
            #if srcId in srcToClust:
                #sys.stderr.write("program error: gene on source with TWO clusters??\n")
                #sys.exit(1)
            #srcToClust[srcId]=clust
        if len(targetIds) > 1:
            for targetId in targetIds:
                targetParalogs[targetId]=targetIds
        if len(srcIds) > 1:
            for srcId in srcIds:
                srcParalogs[srcId]=srcIds

    data = Data()
    data.srcParalogs = srcParalogs
    data.targetParalogs = targetParalogs
    data.srcToTarget = srcToTargetList
    data.targetIds = trgIdSet
    return data
        
