#!/data/www/crispor/venv/bin/python3

from os.path import abspath, join, isdir, isfile
from io import StringIO
import os
import json
import platform
import subprocess
import urllib.error
import tempfile
import collections

"""
for all genomes, writes a json file containing the codon frequency usage
"""

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


def getExonPos(org):

    genomeDir = genomesDir
    twoBitFname = getTwoBitFname(org)
    genomePath = "%(genomeDir)s/%(org)s/" % locals()
    genomeFiles = os.listdir(genomePath)
    gpFiles = [f for f in genomeFiles if f.endswith(".gp")]

    if gpFiles:
        transcriptExons = []
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
                    cdsLen = 0
                    for currentExon in currentExons:
                        chrom, start, end, strand = currentExon
                        cdsLen += (end - start)
                    if cdsLen % 3 == 0:
                        transcriptExons.append(currentExons)
                    else:
                        # CDS not in frame
                        pass

        if geneInfo is None:
            raise ValueError("")

        return transcriptExons
    else:
        return None


def fetch_exon_sequences(org, unique_exons):
    """
    Fetches sequences for a set of exons in batch using twoBitToFa.
    unique_exons: set of (chrom, start, end, strand)
    Returns: dict { (chrom, start, end, strand): sequence }
    """
    # Map exons to IDs to safely handle special characters in chrom names
    exon_list = list(unique_exons)

    # Create temp BED file
    fd, bed_path = tempfile.mkstemp(suffix=".bed", text=True)
    with os.fdopen(fd, 'w') as f:
        for i, (chrom, start, end, strand) in enumerate(exon_list):
            # BED format: chrom, start, end, name
            f.write(f"{chrom}\t{start}\t{end}\t{i}\n")

    twoBitPath = getTwoBitFname(org)
    binPath = join(binDir, "twoBitToFa")

    # Output FASTA file
    fd_out, fa_path = tempfile.mkstemp(suffix=".fa")
    os.close(fd_out)

    cmd = [binPath, f"-bed={bed_path}", twoBitPath, fa_path]
    subprocess.check_call(cmd)

    # Parse FASTA
    seq_map = {}
    current_id = None
    current_seq = []

    with open(fa_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    seq = "".join(current_seq)
                    idx = int(current_id)
                    exon = exon_list[idx]
                    chrom, start, end, strand = exon
                    if strand == "-":
                        seq = revComp(seq)
                    seq_map[exon] = seq
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        # Last sequence
        if current_id is not None:
            seq = "".join(current_seq)
            idx = int(current_id)
            exon = exon_list[idx]
            chrom, start, end, strand = exon
            if strand == "-":
                seq = revComp(seq)
            seq_map[exon] = seq

    os.remove(bed_path)
    os.remove(fa_path)
    return seq_map


def readAllExons():
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
            unique_exons = set()
            for transcript in transcriptExons:
                for exon in transcript:
                    unique_exons.add(exon)

            exon_seq_map = fetch_exon_sequences(org, unique_exons)

            for transcript in transcriptExons:
                exonSeqs = []
                for exon in transcript:
                    if exon in exon_seq_map:
                        exonSeqs.append(exon_seq_map[exon])
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

    readAllExons()


if __name__ == "__main__":
    main()
