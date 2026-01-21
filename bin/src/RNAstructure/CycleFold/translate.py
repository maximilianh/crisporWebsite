from itertools import product
import math
raw_j  = [x.split() for x in open("junction.db")]

known_NCMs = ["14", "15", "16", "222", "233", "223", "232", "224", "242"]
unknown_NCMs = ["13", "212", "221", "213", "231" , "225", "252", "243", "234", "253", "235" ]
index = {"14":2, "15":3, "16":4, "222":5, "233":6, "223":7, "232":8, "224":9, "242":10}

known_hinges = map(lambda x: "".join(x), product(['c', 't'], ['W','S','H'], ['W', 'S', 'H']))

with open("junctions.txt", "w") as f:
    for n, m in product(known_NCMs, known_NCMs):
        dg = raw_j[index[n]][index[m]]
        k = str(math.exp(-float(dg)))
        f.write( "\t".join([n, m, k]) + "\n")
    for n, m in product(unknown_NCMs, unknown_NCMs):
        dg = 1.0
        k = str(math.exp(-float(dg)))
        f.write( "\t".join([n, m, k]) + "\n")

with open("hinges.txt","w") as f:
    for n, m, h in product(known_NCMs+unknown_NCMs, known_NCMs+unkknown_NCMs, known_hinges):
        f.write( "\t".join([n, m, h, "1.0"])+"\n")

with open("pairs.txt", 'w') as f:
    raw_p  = [x.split() for x in open("hinge.db")][1:]
    nucs = ("A",'C','G','U')
    pair_index = {'X':0, 'A':1, "C":2, "G":3, "U":4}
    for n,m in product(nucs,nucs):
        for hinge in known_hinges:
            dg = raw_p[pair_index[n]][pair_index[m]]
            k = str(math.exp(-float(dg)))
            f.write( "\t".join([hinge, n+m, k]) + "\n")



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
