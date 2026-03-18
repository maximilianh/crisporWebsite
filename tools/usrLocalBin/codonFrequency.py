#!/data/www/crispor/venv/bin/python3

from os.path import abspath, join, isdir, isfile
import os
import json
import platform
import subprocess
import tempfile
import collections

"""
for all genomes, writes a json file containing the codon frequency usage based on the longest transcript
"""

# ================== globals ==================

baseDir = "/data/www/crispor/"
binDir = abspath(join(baseDir, "bin", platform.system() + "-" + platform.machine()))
genomesDir = join(baseDir, "genomes")

revTbl = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "M": "K",
    "K": "M",
    "R": "Y",
    "Y": "R",
    "g": "c",
    "a": "t",
    "c": "g",
    "t": "a",
    "n": "n",
    "V": "B",
    "v": "b",
    "B": "V",
    "b": "v",
    "W": "W",
    "w": "w",
}

# ============= utils from crispor.py =============


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


def revComp(seq):
    "rev-comp a dna sequence with UIPAC characters"
    newSeq = []
    for c in reversed(seq):
        newSeq.append(revTbl[c])
    return "".join(newSeq)


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


def getTwoBitFname(db):
    "return the name of the twoBit file for a genome"
    # at UCSC, try to use local disk, if possible
    locPath = join("/scratch", "data", db, db + ".2bit")
    if isfile(locPath):
        return locPath
    path = join(genomesDir, db, db + ".2bit")
    return path


def parsePos(text):
    """parse a string of format chr:start-end:strand and return a 4-tuple
    Strand defaults to + and end defaults to start+23
    """
    if text is not None and len(text) != 0 and text != "?":
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


def getExonPos(org):

    genomeDir = genomesDir
    twoBitFname = getTwoBitFname(org)
    genomePath = "%(genomeDir)s/%(org)s/" % locals()
    genomeFiles = os.listdir(genomePath)
    gpFiles = [f for f in genomeFiles if f.endswith(".gp")]

    if gpFiles:
        transSymbol = {}
        transNoSymbol = []
        geneInfo = None
        for gpCount, gpFile in enumerate(gpFiles):
            if gpCount > 0:
                break
            gpFilePath = os.path.join(genomePath, gpFile)
            with open(gpFilePath, "r") as genePred:

                for geneLine in genePred:
                    geneLine = geneLine.split("\t")
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
                    if len(geneLine) > 11:
                        geneInfo["altName"] = geneLine[11]

                    chrom = geneInfo["chrom"]
                    strand = geneInfo["strand"]
                    exonStarts = geneInfo["exonStarts"]
                    exonEnds = geneInfo["exonEnds"]
                    cdsStart = geneInfo["cdsStart"]
                    cdsEnd = geneInfo["cdsEnd"]
                    altName = geneInfo.get("altName")

                    # remove non-coding transcripts
                    if cdsEnd - cdsStart < 3:
                        continue
                    else:
                        currentExons = []

                    for exonStart, exonEnd in zip(exonStarts, exonEnds):
                        if exonEnd < cdsStart or exonStart > cdsEnd:
                            # remove non-coding exons
                            continue
                        elif exonStart < cdsStart and exonEnd > cdsStart:
                            # remove UTRs
                            currentExons.append((chrom, cdsStart, exonEnd, strand))
                        elif exonStart < cdsEnd and exonEnd > cdsEnd:
                            currentExons.append((chrom, exonStart, cdsEnd, strand))
                        else:
                            currentExons.append((chrom, exonStart, exonEnd, strand))

                    # check if the coding sequence of the current gene is in frame
                    cdsLen = sum(end - start for _, start, end, _ in currentExons)
                    if cdsLen % 3 != 0:
                        continue

                    if altName:
                        if altName not in transSymbol or cdsLen > transSymbol[altName][0]:
                            transSymbol[altName] = (cdsLen, currentExons)
                    else:
                        transNoSymbol.append(currentExons)

        if geneInfo is None:
            raise ValueError("")

        # Combine results
        finalExons = transNoSymbol
        for _, exons in transSymbol.values():
            finalExons.append(exons)

        return finalExons
    else:
        return None


def getExonSeq(org, uniqueExons):
    """
    Fetches sequences for a set of exons in batch using twoBitToFa.
    unique_exons: set of (chrom, start, end, strand)
    Returns: dict { (chrom, start, end, strand): sequence }
    """
    # Map exons to IDs to safely handle special characters in chrom names
    exonList = list(uniqueExons)

    # Create temp BED file
    fd, bedPath = tempfile.mkstemp(dir=join(genomesDir, org), suffix=".bed", text=True)
    with os.fdopen(fd, 'w') as f:
        for i, (chrom, start, end, strand) in enumerate(exonList):
            # BED format: chrom, start, end, name
            f.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, i))

    twoBitPath = getTwoBitFname(org)
    binPath = join(binDir, "twoBitToFa")

    # Output FASTA file
    fdOut, faPath = tempfile.mkstemp(suffix=".fa")
    os.close(fdOut)

    cmd = [binPath, "-bed=%s" % bedPath, twoBitPath, faPath]
    subprocess.check_call(cmd)

    # Parse FASTA
    seqMap = {}
    currentId = None
    currentSeq = []

    with open(faPath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if currentId is not None:
                    seq = "".join(currentSeq)
                    idx = int(currentId)
                    exon = exonList[idx]
                    chrom, start, end, strand = exon
                    if strand == "-":
                        seq = revComp(seq)
                    seqMap[exon] = seq
                currentId = line[1:]
                currentSeq = []
            else:
                currentSeq.append(line)

        # Last sequence
        if currentId is not None:
            seq = "".join(currentSeq)
            idx = int(currentId)
            exon = exonList[idx]
            chrom, start, end, strand = exon
            if strand == "-":
                seq = revComp(seq)
            seqMap[exon] = seq

    os.remove(bedPath)
    os.remove(faPath)
    return seqMap


def calcCodonFrequency():
    """ for each genome, return a list of all exon sequences """

    aaTable = buildCodonTable(key="aa")
    allCodons = {}
    allGenomes = os.listdir(genomesDir)
    for genomeDir in allGenomes:
        if isdir(join(genomesDir, genomeDir)):
            org = genomeDir
            allCodons[org] = collections.Counter()

    for org in allCodons:
        genomePath = join(genomesDir, org)

        # need to differenciate between overlapping genes and splicing variants!!
        # for now, remove duplicate exons
        transcriptExons = getExonPos(org)
        if transcriptExons is not None:
            # Collect unique exons to fetch in batch
            uniqueExons = set()
            for transcript in transcriptExons:
                for exon in transcript:
                    uniqueExons.add(exon)

            exonSeqMap = getExonSeq(org, uniqueExons)

            for transcript in transcriptExons:
                exonSeqs = []
                for exon in transcript:
                    if exon in exonSeqMap:
                        exonSeqs.append(exonSeqMap[exon])
                transcriptSeq = ''.join([seq for seq in exonSeqs])
                if len(transcriptSeq) % 3 == 0:
                    for i in range(0, len(transcriptSeq), 3):
                        codon = transcriptSeq[i:i + 3].upper()
                        allCodons[org][codon] += 1
                else:
                    pass
            # frequency of codons relative to each aa
            for aa in aaTable:
                sumCodons = 0
                refCodons = aaTable[aa]
                for refCodon in refCodons:
                    count = allCodons[org][refCodon]
                    allCodons[org][refCodon] = []
                    allCodons[org][refCodon].append(count)
                    sumCodons += count
                for refCodon in refCodons:
                    codonInfo = allCodons[org][refCodon]
                    freq = codonInfo[0] / sumCodons
                    codonInfo.append(sumCodons)
                    codonInfo.append(freq)
            jsonFname = "%s_codonFrequency.json" % org
            jsonPath = join(genomePath, jsonFname)
            with open(jsonPath, "w", encoding="utf-8") as f:
                json.dump(allCodons[org], f, indent=4, sort_keys=True)


def main():

    calcCodonFrequency()


if __name__ == "__main__":
    main()
