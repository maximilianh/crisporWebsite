from __future__ import print_function
import os
from collections import defaultdict
from itertools import combinations_with_replacement, permutations, product
import RNAstructure as r
import math
import random
import numpy as np
from string import maketrans
debug = False
mfe = True
PRECISION = 10000
if mfe:
    def read(x):
        if x == "0.0":
            return 10000000
        else:
            return (-1 * math.floor(math.log(float(x)) * PRECISION + 0.5))
else:
    read = float

with(open("datafiles/normalization.txt")) as f:
    data = f.read().strip()
    pair_bonus = read(data)

allowed_ncms = tuple([(1,2),(2,1),(1,3),(3,1),(2,2),(2,3),(3,2),(2,4),(4,2),(2,5),(5,2),(3,4),(4,3),(3,5),(5,3),(3,3),(4,4),(3,),(4,),(5,),(6,)])
def strToID(s):
    return tuple([int(x) for x in s])

seqs = defaultdict(lambda: dict())
for ncm in allowed_ncms:
    filename = "datafiles/" + str(len(ncm)) + "".join(map(str,ncm)) + "_seqs.txt"
    with open(filename) as f:
        for line in f:
            seq, p = line.split()
            seqs[ncm][seq] = read(p)

junctions = dict()
with open("datafiles/junctions.txt") as j:
    for line in j:
        inner, outer, p = line.split()
        junctions[(strToID(inner)[1:], strToID(outer)[1:])] = read(p)

hinges = dict()
with open("datafiles/hinges.txt") as h:
    for line in h:
        inner, outer, pairtype, p = line.split()
        hinges[(strToID(inner)[1:], strToID(outer)[1:], pairtype)] = read(p)

pairs = dict()
with open("datafiles/pairs.txt") as pr:
    for line in pr:
        pairtype, nucs, p = line.split()
        pairs[(pairtype, nucs)] = read(p)

connect = dict()
with open("datafiles/connect.txt") as cn:
    for line in cn:
        ncm_a, ncm_b, nucs, p = line.split()
        connect[(ncm_a, ncm_b, nucs)] = read(p)

pairtypes = set(pt for (pt, nucs) in pairs.keys())

def energy(i, j, ip, jp, k, l, sequence, pairtype='cWW'):
    """the ncm "n" being added onto goes from i to ip and jp to j
    the ncm "m" being added on goes from ip to k and l to jp"""
    assert i < ip < k < l < jp < j
    id_n = (ip-i+1, jp-j+1)
    id_m = (k-ip+1, l-jp+1)
    nucs = m[0] + m[-1]
#    seq_n = sequence[i:ip] + sequence[jp:j]
    seq_m = sequence[ip:k] + sequence[l:jp]
    seq = seqs[id_m][seq_m]
    junc = junctions[(id_m, id_n)]
    hinge = hinges[(id_m, id_n, pairtype)]
    pair = pairs[(pairtype, nucs)]
    return seq * junc * hinge * pair

def energy(sequence, ncm_adding, ncm_last, pairtype, nucs):
    seq = seqs[ncm][sequence]
    junc = junctions[(ncm_adding, ncm_last)]
    hinge = hinges[(ncm_adding, ncm_last, pairtype)]
    pair = pairs[(pairtype, nucs)]
    return seq * junc * hinge * pair


def paired(i, j, structure):
    p = structure[i] and structure[i] == j
    if p: assert structure[j] == i
    return p

class SequenceEnd(Exception): pass

intab = "AUCG"
outtab ="UAGC"
complement = maketrans(intab, outtab)
def reverse_complement(s):
    return s.translate(complement)[::-1]
assert reverse_complement("AAGG") == "CCUU"

def efncm_duplex(s):
    #TODO how do I consider the sequence contribution from both ends?
    # going one way and the opposite way might not give the same answer
    rc = s.translate(complement)
    e= 1.
    # first multiply all the ncm sequences
    for i, _ in enumerate(s[:-1]):
        ip = i + 1
        j = len(s) - i
        jp = j - 1
        seq = s[i:i+2] + rc[i:i+2][::-1]
        if debug:
            print(seq)
        e *= seqs[(2,2)][seq]
    # then do all the junctions, hinges, and pairs
    for i, _ in enumerate(s[1:-1]):
        e *= junctions[((2,2), (2,2))]
        e *= hinges[((2,2), (2,2), 'cWW')]
        if debug:
            print("nucs", s[i+1] + rc[i+1])
        e *= pairs[('cWW', s[i+1] + rc[i+1])]
    return e

def efncm_hairpin(stem, loop):
    #TODO how do I consider the sequence contribution from both ends?
    # going one way and the opposite way might not give the same answer
    rc = s.translate(complement)
    e= 1.
    # first multiply all the ncm sequences
    for i, _ in enumerate(stem[:-1]):
        ip = i + 1
        j = len(stem) - i
        jp = j - 1
        seq = stem[i:i+2] + rc[i:i+2][::-1]
        if debug:
            print(seq)
        e *= seqs[(2,2)][seq]
    # then do all the junctions, hinges, and pairs
    for i, _ in enumerate(stem[1:-1]):
        e *= junctions[((2,2), (2,2))]
        e *= hinges[((2,2), (2,2), 'cWW')]
        if debug:
            print("nucs", stem[i+1] + rc[i+1])
        e *= pairs[('cWW', stem[i+1] + rc[i+1])]
    #base of hairpin
    loop_ncm_seq = stem[-1] + loop[0] + loop[-1] + rc[-1]
    loop_1_ncm_seq = loop[1:-1]
    e *= seqs[(2,2)][loop_ncm_seq]
    e *= seqs[(1,len(loop)-2)][loop_1_ncm_seq]
    e *= junctions[((2,2), (2,2))]
    e *= hinges[((2,2), (2,2), 'cWW')]
    if debug:
        print( "loop 222", loop_ncm_seq)
        print( "loop", loop_1_ncm_seq)
        print( "nucs", loop[0] + loop[-1])
    return e

def generate_duplexes(length):
    seqs = set()
    for seq in product("AUGC",repeat=length):
        s = ''.join(seq)
        if s not in seqs and reverse_complement(s) not in seqs:
            seqs.add(s)
            yield s

def duplex_dG(s):
    rna = r.HybridRNA.fromString(s,reverse_complement(s))
    rna.FoldDuplex()
    rna.RemovePairs()
    for i in range(len(s)):
        j = len(s)+3+len(s)-i
        rna.SpecifyPair(i+1, len(s)+3+len(s)-i)
    return rna.CalculateFreeEnergy()

def hairpin_dG(stem, loop):
    sequence = stem+loop+reverse_complement(stem)
    rna = r.RNA.fromString(sequence)
    rna.FoldDuplex()
    rna.RemovePairs()
    for i in range(len(stem)):
        j = rna.GetSequenceLength()-i
        rna.SpecifyPair(i+1, j)
    return rna.CalculateFreeEnergy()

def self_complementary(s):
    return s == reverse_complement(s)

def self_complementary_penalty(s):
    return 0.4 if self_complementary(s) else 0.

initiation_penalty = 4.1

def AU_end_penalty(s):
    five = 0.45 if s[0] in "AU" else 0.
    three = 0.45 if s[-1] in "AU" else 0.
    return five + three

def duplex_nearest_neighbors(s, dG=None):
    return (duplex_dG(s) if dG is None else dG) - self_complementary_penalty(s) - initiation_penalty - AU_end_penalty(s)

RT = 0.5925

def unpaired_prob(dG, paired_prob):
    return paired_prob / (math.exp(-dG/RT))

def random_duplex(length):
    return "".join([random.choice(["A","U","C","G"]) for _ in range(length)])

def random_duplex_energy(length):
    return duplex_dG(random_duplex(length))

def main():
    print("Sequence\tduplex_dG\tnearest_neighbors_dG\tNCM_probability\tunpaired_probability")
    for seq in generate_duplexes(3):
        dG = duplex_dG(seq)
        nn_dG = duplex_nearest_neighbors(seq)
        prob = efncm_duplex(seq)
        unp = unpaired_prob(nn_dG, prob)
        print("\t".join(map(str,[seq, dG, nn_dG, prob, unp])))

    print("\nSequence\tduplex_dG\tnearest_neighbors_dG\tNCM_probability\tunpaired_probability")
    for line in open("xia_duplexes.txt"):
        l = line.split()
        seq = l[0].replace('p','').replace('/','')
        dG = -1 * float(l[1])
        nn_dG = duplex_nearest_neighbors(seq, dG)
        prob = efncm_duplex(seq)
        unp = unpaired_prob(nn_dG, prob)
        print("\t".join(map(str,[seq, dG, nn_dG, prob, unp])))

if __name__ == "__main__":
    main()
