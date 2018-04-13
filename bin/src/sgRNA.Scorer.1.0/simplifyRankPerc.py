# keep only the ranges of the SVML scores to obtain the rank percents
# up to 2 digits precision. Makes calculation of rank percent much faster.
vals = [float(x) for x in open("Hg19.RefFlat.Genes.75bp.NoUTRs.SPSites.SVMOutput.txt").read().splitlines()]
vals.sort()

lastPerc = -1 # make sure we keep the first val
ranges = []
for i, val in enumerate(vals):
    rankPerc = int(100*(float(i)/len(vals)))
    if rankPerc > lastPerc:
        ranges.append(str(val))
        lastPerc = rankPerc

ofh = open("Hg19.RefFlat.Genes.75bp.NoUTRs.SPSites.SVMOutput.ranges.txt", "w")
ofh.write("\n".join(ranges))
        
