import sys
import os
import tempfile, logging
import logging

# ignore errors if util.py not found, usually not needed
try:
    import util
except:
    pass

def coordOverlap(start1, end1, start2, end2):
    """ returns true if two Features overlap """
    result = (( start2 <= start1 and end2 > start1) or \
            (start2 < end1 and end2 >= end1) or \
            (start1 >= start2 and end1 <= end2) or \
            (start2 >= start1 and end2 <= end1))
    return result
    # modified: end2 > end1

def overlap(f1, f2):
    """ returns true if two Features overlap """
    if f1.chrom != f2.chrom:
        return False
    result = coordOverlap(f1.start, f1.end, f2.start, f2.end)
    return result

class Features(list):
    def __repr__(self):
        buf = []
        for i in self:
            buf.append(repr(i))
        return "\n".join(buf)
    def __sort__(self):
        self.sort(sort)

    def writeToFileHandle(self, fh):
        for b in self:
            fh.write(str(b)+"\n")

    def writeToFile(self, fname):
        if fname!="stdout":
            fh = open(fname, "w")
        else:
            fh = sys.stdout
        self.writeToFileHandle(fh)
        fh.close()

    def getChromosomes(self):
        chroms = set()
        for b in self:
            chroms.add(b.chrom)
        return chroms

    def countChromosomes(self):
        return len(self.getChromosomes())

    def bedsOutsideChrom(self, chromSizes):
        """ check if all features are within the limits of the chromSizes dictionary """
        for b in self:
            if int(b.end) > int(chromSizes[b.chrom]):
                return b
        return False

    def anyFlankingOverlap(self):
        """ is there is any overlap between any two subsequent features ? """
        lastB = None
        for b in self:
            if lastB!=None and b.overlaps(lastB):
                return b
            else:
                lastB = b
        return False

    def totalLength(self):
        """ get total lengths of all features added together """
        sum = 0
        for b in self:
            sum+=(b.end-b.start)
        return sum

    def avgLength(self):
        return self.totalLength() / len(self)

class Feature:
    def __init__(self, line=None, fields=None):
        if fields==None:
            fields = line.split()
        self._initFromFields(l)

    def _initFromFields(self, fields):
        count = len(fields)
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        if count >= 4:
            self.name = fields[3]
        if count >= 5:
            self.strand= fields[4] 
        if count >= 6:
            self.score = int(fields[5])
        return self

    def __init__(self,seqid="",start=0,end=0,name="",score=0,strand="+"):
        self.chrom = seqid
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = int(score)
        self.strand = strand

    def __repr__(self):
        fields = [self.chrom,str(self.start),str(self.end)]
        if "name" in self.__dict__:
            fields.append(self.name)
        if "score" in self.__dict__:
            fields.append(str(self.score))
        if "strand" in self.__dict__:
            fields.append(self.strand)
        if "blockStarts" in self.__dict__:
            fields.append(str(self.thickStart))
            fields.append(str(self.thickEnd))
            fields.append(self.itemRgb)
            fields.append(str(self.blockCount))
            fields.append(self.blockSizes)
            fields.append(self.blockStarts)

        return "\t".join(fields)

    def includes(self, f, tol=0):
        if self.chrom!=f.chrom:
            return False
        else:
            if self.start-tol <= f.start and self.end+tol >= f.end:
                return True

    def overlaps(self, f1):
        """ returns true if two Features overlap """
        return overlap(self,f1)

    def toString(self):
        return repr(self)

    def joinBlocks(self):
        """ join all blocks with 0-spacers between them into longer ones """
        newStarts = []
        newSizes = []

        starts = self.blockStarts.strip(",").split(",")
        sizes  = self.blockSizes.strip(",").split(",")
        starts = [int(x) for x in starts]
        sizes = [int(x) for x in sizes]

        lastStart = 0
        lastEnd   = 0
        for i in range(0, len(starts)):
            start = starts[i]
            end   = starts[i]+sizes[i]
            #print "start, end", start, end

            if lastEnd != 0:
                if start!=lastEnd:
                    #print "start==lastEnd"
                    newStarts.append(lastStart)
                    newSizes.append(lastEnd-lastStart)
                    lastStart = start
                else:
                    pass
                    # do not change lastStart
            lastEnd   = end
        newStarts.append(lastStart)
        newSizes.append(lastEnd-lastStart)

        self.blockCount  = len(newSizes)
        self.blockStarts = ",".join([str(x) for x in newStarts])
        self.blockSizes  = ",".join([str(x) for x in newSizes])

class Feature:
    def __init__(self, line=None, fields=None):
        if fields==None:
            fields = line.split()
        self._initFromFields(l)

    def _initFromFields(self, fields):
        count = len(fields)
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        if count >= 4:
            self.name = fields[3]
        if count >= 5:
            self.strand= fields[4] 
        if count >= 6:
            self.score = int(fields[5])
        return self

    def __init__(self,seqid="",start=0,end=0,name="",score=0,strand="+"):
        self.chrom = seqid
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = int(score)
        self.strand = strand

    def __repr__(self):
        fields = [self.chrom,str(self.start),str(self.end)]
        if "name" in self.__dict__:
            fields.append(self.name)
        if "score" in self.__dict__:
            fields.append(str(self.score))
        if "strand" in self.__dict__:
            fields.append(self.strand)
        if "blockStarts" in self.__dict__:
            fields.append(str(self.thickStart))
            fields.append(str(self.thickEnd))
            fields.append(self.itemRgb)
            fields.append(str(self.blockCount))
            fields.append(self.blockSizes)
            fields.append(self.blockStarts)

        return "\t".join(fields)

    def includes(self, f, tol=0):
        if self.chrom!=f.chrom:
            return False
        else:
            if self.start-tol <= f.start and self.end+tol >= f.end:
                return True

    def overlaps(self, f1):
        """ returns true if two Features overlap """
        return overlap(self,f1)

    def toString(self):
        return repr(self)

    def joinBlocks(self):
        """ join all blocks with 0-spacers between them into longer ones """
        newStarts = []
        newSizes = []

        starts = self.blockStarts.strip(",").split(",")
        sizes  = self.blockSizes.strip(",").split(",")
        starts = [int(x) for x in starts]
        sizes = [int(x) for x in sizes]

        lastStart = 0
        lastEnd   = 0
        for i in range(0, len(starts)):
            start = starts[i]
            end   = starts[i]+sizes[i]
            #print "start, end", start, end

            if lastEnd != 0:
                if start!=lastEnd:
                    #print "start==lastEnd"
                    newStarts.append(lastStart)
                    newSizes.append(lastEnd-lastStart)
                    lastStart = start
                else:
                    pass
                    # do not change lastStart
            lastEnd   = end
        newStarts.append(lastStart)
        newSizes.append(lastEnd-lastStart)

        self.blockCount  = len(newSizes)
        self.blockStarts = ",".join([str(x) for x in newStarts])
        self.blockSizes  = ",".join([str(x) for x in newSizes])

def parseBedLine(line, fieldCount=None):
    f = Feature()
    fields = line.split("\t")
    if fieldCount!=None:
        fields = fields[:fieldCount]
    if not len(fields)>=3:
        logging.error("Illegal BED format, line %s" % line)
        exit(1)
    f.chrom = fields[0]
    f.start = int(fields[1])
    f.end = int(fields[2])
    if len(fields)>3:
        f.name = fields[3]
    if len(fields)>4:
        f.score = int(fields[4])
    if len(fields)>5:
        f.strand = fields[5]
    if len(fields)>6:
        f.thickStart=int(fields[6])
        f.thickEnd=int(fields[7])
        f.itemRgb=fields[8]
        f.blockCount=int(fields[9])
        f.blockSizes=fields[10]
        f.blockStarts=fields[11]
    return f

#def parseBedLine(line):
    #f = Feature()
    #fields = line.split("\t")
    #if not len(fields)>=3:
        #logging.error("Illegal BED format, line %s" % line)
        #exit(1)
    #f.chrom = fields[0]
    #f.start = int(fields[1])
    #f.end = int(fields[2])
    #if len(fields)>3:
        #f.name = fields[3]
    #if len(fields)>4:
        #f.score = int(fields[4])
    #if len(fields)>5:
        #f.strand = fields[5]
    #if len(fields)>6:
        #f.thickStart=int(fields[6])
        #f.thickEnd=int(fields[7])
        #f.itemRgb=fields[8]
        #f.blockCount=int(fields[9])
        #f.blockSizes=fields[10]
        #f.blockStarts=fields[11]
    #return f

def sort(f1, f2):
    """ sort features, feature with smaller coord's first"""
    if f1.start < f2.start:
        return f1,f2
    else:
        return f2,f1

def cmpFeatures(f1, f2):
    """ for the sort-function in python """
    return f2.start-f1.start

def parseBedFilename(fname, reqSorted=False, quiet=False, reqNoOverlaps=False, fieldCount=None):
    if not quiet:
        sys.stderr.write("Reading %s...\n" % fname)
    if fname=="stdin":
        return parseBedFile(sys.stdin, reqSorted, reqNoOverlaps, fieldCount=fieldCount)
    else:
        return parseBedFile(open(fname,"r"), reqSorted, reqNoOverlaps, fieldCount=fieldCount)

def openFilterNumber(validNames, fname, geneNameLen=25):
    """ parse fname bed features but keep only features with a name in validNames and add attribute .no to all features """
    """ remove all genes that occur two times in file """
    """ trim down gene names """
    rawbeds = parseBedFilename(fname, reqSorted=True, quiet=True)
    newbeds = []
    # find dupl names
    oldNames = set()
    dupl = set()
    for b in rawbeds:
        b.name = b.name[:25]
        if b.name in oldNames:
            dupl.add(b.name)
        oldNames.add(b.name)

    # filter genes

    if len(dupl)>0:
        sys.stderr.write("Info: file %s, dropping duplicate genes: %s\n" % (fname, ",".join(dupl)))
    for b in rawbeds:
        if "___" in b.name:
            b.name = b.name.split("___")[0]
        if b.name in dupl:
            continue
        if b.name in validNames:
            newbeds.append(b)

    # number genes
    for i in range(0, len(newbeds)):
        b = newbeds[i]
        b.no = i
    return newbeds, rawbeds

def parseBedFile(lines, reqSort=False, reqNoOverlaps=False, fieldCount=None):
    """ will return a Features() object """
    features = Features()
    lastStart = -1
    for l in lines:
        l = l.strip()
        if l.startswith("track") or l.startswith("browser") or l.startswith("#") or l=="":
            continue
        f = parseBedLine(l, fieldCount=fieldCount)
        if reqSort and lastStart > f.start:
            sys.stderr.write("error: bed file is not sorted, violating feature: %s" % str(f))
            sys.exit(1)
        features.append(f)

    if reqNoOverlaps:
        assert(features.anyFlankingOverlap()==False)

    return features

def indexBedsName(fname):
    if fname=="stdin":
        beds= parseBedFile(sys.stdin)
    else:
        beds= parseBedFile(open(fname,"r"))
    idx = {}
    for b in beds:
        idx.setdefault(b.name, []).append(b)
    return idx

def indexBedsUniqueName(fname):
    if fname=="stdin":
        beds= parseBedFile(sys.stdin)
    else:
        beds= parseBedFile(open(fname,"r"))
    idx = {}
    for b in beds:
        if b.name in idx:
            sys.stderr.write("error: feature with name %s exists two times.\n" % b.name)
            sys.exit(1)
        idx[b.name]=b
    return idx

def indexBedsChrom(fname):
    if fname=="stdin":
        beds= parseBedFile(sys.stdin, reqSort=True)
    else:
        beds= parseBedFile(open(fname,"r"), reqSort=True)
    idx = {}
    for b in beds:
        idx.setdefault(b.chrom, []).append(b)
    return idx

def indexBedsByStart(beds):
    dict = {}
    for b in beds:
        dict.setdefault(b.start, []).append(b)
    return dict

def sortBedFile(fname):
    logging.info("Sorting %s" % fname)
    cmd = "bedSort %s %s" % (fname, fname)
    ret = os.system(cmd)
    if ret !=0:
        sys.stderr.write("error: cmd %s failed with error code %d\n" % (cmd, ret))

def namesAsDict(beds):
    dict = {}
    for b in beds:
        dict[b.name]=True
    return dict

def filterByName(beds, names):
    out = []
    for b in beds:
        if b.name in names:
            out.append(b)
    return out

def show(bedList):
    for b in beds:
        print b

def writeToFile(beds,fname):
    fh = open(fname, "w")
    for b in beds:
        fh.write(str(b)+"\n")

def countChromosomes(beds):
    chroms = set()
    for b in beds:
        chroms.add(b.chrom)
    return len(chroms)

def bedAnnotateDownstream(bedFile, geneFile):
    """ annotate bed features with the gene downstream of it """

    tempfile = tempfile.NamedTemporaryFile()
    cmd = 'bedFindNeighbors %s %s --onlyDownstream > %s' % (bedFile, geneFile, tempfile.name)
    util.execCmdLine(cmd)
    beds =  parseBedFilename(tempfile.name)
    tempfile.close()
    return beds
