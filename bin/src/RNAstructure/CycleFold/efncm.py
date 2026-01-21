from __future__ import print_function
import energy as e
import re
import sys
from itertools import islice, product
#TODO I need to think carefully about how to handle NCMs like 212 and 213
debug = True
pairtypes = map(lambda x: "".join(x), product('ct','HWS','HWS'))
big_hairpins=True

def parse_ct(contents):
    """parse a ct, provided as a string, into a sequence and a dict of pairs"""
    lines = filter(None, contents.split("\n")[1:])
    assert all(re.match(r'\s*[0-9]+\s+\D(\s+[0-9]+){4,4}', x) for x in lines)
    seq = "".join(map(lambda x:x.split()[1],lines))
    def lineToPair(l):
        i, j = int(l.split()[0])-1 , int(l.split()[4])-1
        return i, j if j!=-1 else None
    pairs = map(lineToPair, lines)
    return seq, pairs

def efncm(ct):
    seq, pairs = parse_ct(ct) #get structure info
    stems = find_stems((0,len(seq)), pairs) #LIFO stack of stems to score
    score = 0
    while len(stems) > 0:
        p = follow_stem(stems.pop(), pairs) #take a stem, find all of its pairs
        if debug: print("found pairs",p)
        i, j = p[-1] #last pair of stem p defines our new fragment
        inner_stems = find_stems((i+1, j-1),pairs) #put inner stems on the stack
        stems += inner_stems
        score += score_stem(p, seq, is_hairpin=not inner_stems)
    return score

def ncm_id(pairs):
    (i, j), (ip, jp) = pairs
    assert i<=ip<jp<=j
    return (ip-i+1, j-jp+1)

def ncm_seq(pairs, seq):
    (i, j), (ip, jp) = pairs
    assert i<=ip<jp<=j
    return seq[i:ip+1]+seq[jp:j+1]

def score_stem(p, seq, is_hairpin):
    if debug: print ("scoring stem starting at %d %d, length %d" % (p[0]+(len(p),)))
    ret = e.pair_bonus * len(p)
    ncms = map(ncm_id, window(p, 2))
    seqs = map(lambda x: ncm_seq(x, seq), window(p, 2))
    if is_hairpin:
        hairpin_i, hairpin_j = p[-1]
        if not (big_hairpins and hairpin_j-hairpin_i+1 > 6):
            ncms.append((hairpin_j-hairpin_i+1,))
            seqs.append(seq[hairpin_i:hairpin_j+1])
    for ncm, s in zip(ncms, seqs):
        if debug:
            print("added an ncm", ncm, s, e.seqs[ncm][s])
        ret += e.seqs[ncm][s]
    #there are either len(p)-2 junctions (stem closing multibranch loop)
    #or len(p)-1 junctions (stem closing hairpin loop
    junctions = zip(window(ncms,2), p[1:] if is_hairpin else p[1:-1])
    for junc, pair in junctions:
        #print(junc)
        #print(pair)
        #ret += e.junctions[junc]
        a,b = junc
        ncm_a = "%d%s"%(len(a),"".join(map(str,a)))
        ncm_b = "%d%s"%(len(b),"".join(map(str,b)))
        #ret += e.connect[(a,b,
        i, j = pair
        #ret +=  min(e.hinges[junc+(p,)]+e.pairs[(p,seq[i]+seq[j])] for p in pairtypes)
        ret += e.connect[(ncm_a, ncm_b,seq[i]+seq[j])]
        if debug:
        #    _, best_pairtype  = min([(e.hinges[junc+(p,)]+e.pairs[(p,seq[i]+seq[j])],p) for p in pairtypes], key = lambda x:x[0])
            print("added a connection", (ncm_a, ncm_b, seq[i]+seq[j]), e.connect[(ncm_a, ncm_b,seq[i]+seq[j])])
        #    print("added a junction", junc, e.junctions[junc])
        #    print("added a hinge", junc, best_pairtype, e.hinges[junc+(best_pairtype,)])
        #    print("added a pair", best_pairtype, seq[i]+seq[j],e.pairs[(best_pairtype, seq[i]+seq[j])])
    return ret

def window(seq, n=2):
    """
    Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
    """
    s = tuple(seq)
    it = iter(s)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def get_first(pred, seq):
    """return the first element of seq for which pred is true"""
    for i in seq:
        if pred(i): return i
    return None

def paired(p):
    """is this tuple(int,int) a pair?"""
    return p[1] is not None

def paired_together(i,j,pairs):
    """is i paired to j?"""
    result = pairs[i][1] == j
    if result: assert pairs[j][1] == i
    return result

def same_pair(a, b):
    """do a and b represent the same base pair?"""
    i,j = a
    jp, ip = b
    return (i,j)==(ip,jp) or (i,j)==(jp,ip)

def continue_stem(this_pair, pairs):
    """given a base pair, find the next pair in the stem
    or None if this_pair closes the stem"""
    i, j = this_pair
    assert pairs[i][1] == j
    left = get_first(paired, pairs[i+1:j])
    right = get_first(paired, reversed(pairs[i+1:j]))
    if left and right and same_pair(left, right):
        return left
    return None

def find_stems(fragment, pairs):
    """given a fragment (i,j), find all of the stems (ip,jp) where
    ip>=i and jp<= j and ip-1 is not paired to jp+1
    i.e. all the "closing pairs" in this fragment"""
    i, j = fragment
    if debug: print("finding stems for fragment",i,j)
    if(paired_together(i,j,pairs)): return (i,j)
    result = []
    while i < j:
        nx = get_first(paired, pairs[i:j+1])
        if nx is None: break
        ip, jp = nx
        assert jp <= j
        result.append((ip,jp))
        i = jp+1
    return result

def follow_stem(fragment, pairs):
    """return all the pairs in this stem"""
    res = [fragment]
    while True:
        fragment = continue_stem(fragment, pairs)
        if fragment is not None:
            res.append(fragment)
        else:
            return res

def is_hairpin(fragment, pairs):
    """does this pair close a hairpin loop?"""
    i,j = fragment
    return not any(map(paired, pairs[i+1:j]))

def main():
    try:
        filename = sys.argv[1]
        print(efncm(open(filename).read()) / 1)#e.PRECISION)
    except IndexError:
        print("provide a ct file")

if __name__ == "__main__":
    main()
