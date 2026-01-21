import os
from collections import defaultdict

allowed_ncms = tuple([(1,2),(2,1),(1,3),(3,1),(2,2),(2,3),(3,2),(2,4),(4,2),(2,5),(5,2),(3,4),(4,3),(3,5),(5,3),(3,3),(4,4),(3,),(4,),(5,),(6,)])
def strToID(s):
    return tuple([int(x) for x in s])

seqs = defaultdict(lambda: dict())
for ncm in allowed_ncms:
    filename = str(len(ncm)) + "".join(map(str,ncm)) + "_seqs.txt"
    with open(filename) as f:
        for line in f:
            seq, p = line.split()
            seqs[ncm][seq] = float(p)

junctions = dict()
with open("junctions.txt") as j:
    for line in j:
        inner, outer, p = line.split()
        junctions[(strToID(inner), strToID(outer))] = float(p)

hinges = dict()
with open("hinges.txt") as h:
    for line in h:
        inner, outer, pairtype, p = line.split()
        hinges[(strToID(inner), strToID(outer), pairtype)] = float(p)

pairs = dict()
with open("pairs.txt") as pr:
    for line in pr:
        pairtype, nucs, p = line.split()
        pairs[(pairtype, nucs)] = float(p)
