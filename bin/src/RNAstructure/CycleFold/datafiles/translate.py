from itertools import product
from collections import defaultdict
import math
import re
raw_j  = [x.split() for x in open("junction.db")]

NCMs = ["14", "15", "16", "222", "233", "223", "232", "224", "242"]
unNCMs = ["13","244","13", "212", "221", "213", "231" , "225", "252", "243", "234", "253", "235" ]
index = {"14":2, "15":3, "16":4, "222":5, "233":6, "223":7, "232":8, "224":9, "242":10}

faces = map(lambda x:"".join(x), product(['W','S','H','B'],['W','S','H','B']))
hinges = list(product(faces, ['a','p'], ['c', 't']))
nucs = ["".join(x) for x in product("AGCU","AGCU")]

def KtoG(k):
    return -0.606*math.log(k)

def GtoK(g):
    return math.exp(-g/0.606)

def GtoP(g):
    return math.exp(-g)

def convert_line(line):
    seq, g = line.split()
    return (seq,float(g))

def convert_title(title):
    t = title.split(".")[0].split("_")
    if t[0] == "1":
        return t[0] + t[2] 
    elif t[0] == "2":
        return t[0] + t[2] + t[3]
    else:
        raise ValueError

major_files = (
["1_strand_3.db.cycle", "2_strand_2_2.db.cycle",  "2_strand_3_4.db.cycle",
"1_strand_4.db.cycle", "2_strand_2_3.db.cycle",  "2_strand_3_5.db.cycle",
"1_strand_5.db.cycle", "2_strand_2_4.db.cycle",  "2_strand_4_2.db.cycle",
"1_strand_6.db.cycle", "2_strand_2_5.db.cycle",  "2_strand_4_3.db.cycle",
"2_strand_1_2.db.cycle", "2_strand_3_1.db.cycle",  "2_strand_4_4.db.cycle",
"2_strand_1_3.db.cycle", "2_strand_3_2.db.cycle",  "2_strand_5_2.db.cycle",
"2_strand_2_1.db.cycle", "2_strand_3_3.db.cycle",  "2_strand_5_3.db.cycle"]
)


hinge_aliases = {'anti':'a','para':'p','cis':'c','trans':'t'}
def hinge_line_to_code(h):
    (edges, pa, ct, counts) = h.split()
    return (edges, hinge_aliases[pa], hinge_aliases[ct])
def hinge_line_to_count(h):
    return int(h.split()[-1])
def hinge_code_to_prime(h):
    return h.split()[0]
def hinge_code_to_nucs(h):
    return "".join(h.split()[1:])

def pair_bonus():
    raw_p  = [x.split() for x in open("hinge.db")][1:]
    nucs = ("A",'C','G','U')
    pair_index = {'X':0, 'A':1, "C":2, "G":3, "U":4}
    ret = dict()
    for n,m in product(nucs,nucs):
        dg = raw_p[pair_index[n]][pair_index[m]]
        ret["".join((n,m))] = float(dg)
    return ret

def junction_bonus():
    ret = dict()
    for n, m in product(NCMs, NCMs):
        dg = raw_j[index[n]][index[m]]
        ret[(n,m)] = float(dg)
    return ret

def makeCycleVals():
    p = pair_bonus()
    cycleVals = {convert_title(f):dict(convert_line(line) for line in open(f)) for f in major_files}
    #correct the 1-3 sequences which are length 5 for some reason
    cycleVals["13"] = {seq[1:-1]:cycleVals["13"][seq] for seq in cycleVals["13"]}
    for key in cycleVals:
        for seq in cycleVals[key]:
            cycleVals[key][seq] += p[seq[0]+seq[-1]]
    print cycleVals.keys()
    print cycleVals["13"]
    return cycleVals


def parse_hinge_db(ncm,hinge_counts,hinge_totals):
    with open("%s_strand_%s%s.db.hinge" % (ncm[0], ncm[1], ("_"+ncm[2]) if len(ncm)==3 else "")) as f:
        data = f.read().replace("'","").replace('/','')
        sp = filter(None,re.split(r'\[([A-Za-z0-9 ]+)\]', data  ))
        z = zip(sp[::2], [filter(None,x.split('\n')) for x in sp[1::2]])
        for (nucs, hinges) in z:
            prime = hinge_code_to_prime(nucs)
            n = hinge_code_to_nucs(nucs)
            for h in hinges:
                hh = hinge_line_to_code(h)
                c = hinge_line_to_count(h)
                hinge_counts[(ncm, prime, n, hh)] = c
                hinge_totals[(ncm, prime, n)] += c
    return hinge_counts,hinge_totals

def make_connect(counts,totals):
    junc_bonus = junction_bonus()
    connect_table = defaultdict(lambda:0.0)
    for (a,b,n) in product(NCMs+unNCMs,NCMs+unNCMs,nucs):
        acc = 0.0
        for h in hinges:
            try:
                if a[0]=='1': continue
                # in the db, '3' means towards the ends of the sequence
                # so single NCMs have only 5' pairs and are always outer
                inner = float(counts[(a,'3',n,h)]) / float(totals[(a,'3',n)])
                outer = float(counts[(b,'5',n,h)]) / float(totals[(b,'5',n)])
                assert inner <= 1.0
                assert outer <= 1.0
            except ZeroDivisionError:
                inner = 0.0
                outer = 0.0
            acc += inner*outer
        connect_table[(a,b,n)] = (KtoG(acc)+junc_bonus[(a,b)] if acc>0.0 
                                  else float("inf"))
    return connect_table

def symmetrize_connect(connect):
    ret = dict()
    for (a,b) in product(NCMs+unNCMs, NCMs+unNCMs):
        for n,m in nucs:
            oneWay = connect[(a,b,(n,m))]
            otherWay = connect[(b,a,(m,n))]
        ret[(a,b,(n,m))] = ret[(b,a,(m,n))] = (oneWay + otherWay) / 2.

def test():
    counts = defaultdict(lambda:0)
    totals = defaultdict(lambda:0)
    for n in NCMs:
        (counts, totals) = parse_hinge_db(n, counts, totals)
    c = dict(make_connect(counts, totals))
    cycleVals = makeCycleVals()
    print "222-222-GC", c[("222","222","GC")]
    print "222-222-CG", c[("222","222","CG")]
    print "222-15-CG", c[("222","15","CG")]
    print "222 GCGC",cycleVals["222"]["GCGC"]
    print "222 CGCG",cycleVals["222"]["CGCG"]
    print "222 GCGC",cycleVals["222"]["GCGC"]
    print "15 CAAAG",cycleVals["15"]["CAAAG"]
    with open("connect.txt",'w') as f:
        pass
        # for (a,b,n) in product(NCMs+unNCMs,NCMs+unNCMs,nucs):
            # pass
            #f.write("%s\t%s\t%s\t%f\n" % (a,b,n,c[(a,b,n)]))

def dict_to_file(filename, d):
    with open(filename,'w') as f:
        for key,val in sorted(d.items()):
            if type(key) is str: sep = ""
            else: sep = "\t"
            f.write(sep.join(key) +"\t"+ str(GtoP(val)) +"\n")

def main():
    counts = defaultdict(lambda:0)
    totals = defaultdict(lambda:0)
    for n in NCMs:
        (counts, totals) = parse_hinge_db(n, counts, totals)
    connect = symmetrize_connect(make_connect(counts, totals))
    cycleVals = makeCycleVals()
    for cycleType in cycleVals:
        dict_to_file(cycleType+"_seqs.txt" ,cycleVals[cycleType])
    dict_to_file("connect.txt",connect)

if __name__ == "__main__": main()




"""
with open("junctions.txt", "w") as f:
    for n, m in product(NCMs+unNCMs, NCMs+unNCMs):
        if n in NCMs and m in NCMs:
            dg = raw_j[index[n]][index[m]]
            k = str(math.exp(-float(dg)))
            f.write( "\t".join([n, m, k]) + "\n")
        else:
            dg = 1.0
            k = str(math.exp(-float(dg)))
            f.write( "\t".join([n, m, k]) + "\n")

with open("hinges.txt","w") as f:
    for n, m, h in product(NCMs+unNCMs, NCMs+unNCMs, hinges):
        f.write( "\t".join([n, m, h, "1.0"])+"\n")

with open("pairs.txt", 'w') as f:
    raw_p  = [x.split() for x in open("hinge.db")][1:]
    nucs = ("A",'C','G','U')
    pair_index = {'X':0, 'A':1, "C":2, "G":3, "U":4}
    for n,m in product(nucs,nucs):
        for hinge in hinges:
            dg = raw_p[pair_index[n]][pair_index[m]]
            k = str(math.exp(-float(dg)))
            f.write( "\t".join([hinge, n+m, k]) + "\n")
"""



"""
-- | Known single-pair types. (Length of nucleotide string, 0-based index).

-- 4 -> (5,2)
-- 5 -> (5,3)
-- 6 -> (5,4)
--
-- found to be wrong
-- 5 -> (0,5)

knownSingleNCM :: VU.Vector (Int,Int)
knownSingleNCM = VU.fromList
  [ (4,2)
  , (5,3)
  , (6,4)
  ]

-- | Known double-pair types. ((Length of left...,Length of right nucleotide
-- string), 0-based index

-- found to be right through hacking!
-- all are (..) to next type
-- NCM type -> junctions IV.! (a,b)
-- (a+1,b+1) for the junction.db
-- (2,2) -> (5,5)
-- (2,3) -> (8,2)
-- (3,2) -> (9,2)

knownDoubleNCM :: VU.Vector ((Int,Int),Int)
knownDoubleNCM = VU.fromList
  [ ((2,2), 5)
  , ((3,3), 6)
  , ((2,3), 7)
  , ((3,2), 8)
  , ((2,4), 9)
  , ((4,2),10)
  ]
"""
