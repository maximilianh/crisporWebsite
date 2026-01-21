"""
with open("hinges.txt") as h:
    data = [x.strip().split() for x in h]
    out = []
    for line in data:
        line[0] = str(len(line[0])) + line[0]
        line[1] = str(len(line[1])) + line[1]
        out.append('\t'.join(line))
with open("hinges.txt",'w') as h:
    h.write("\n".join(out)+"\n")
"""
with open("junctions.txt") as j:
    data = [x.strip().split() for x in j]
    out = []
    for line in data:
        line[0] = str(len(line[0])) + line[0]
        line[1] = str(len(line[1])) + line[1]
        out.append('\t'.join(line))
with open("junctions.txt",'w') as j:
    j.write("\n".join(out)+"\n")
"""
with open("pairs.txt") as p:
    data = [x.strip().split() for x in p]
    out = []
    for line in data:
        line[1] = line[0]+line[1]
        out.append('\t'.join([line[2], line[1], line[3]]))
with open("pairs.txt",'w') as p:
    p.write("\n".join(out)+"\n")

"""
