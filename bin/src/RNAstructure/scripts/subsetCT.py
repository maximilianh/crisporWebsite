#Program reads the command line to get the name of an input ct file, the name of an output CT file,
# start nucleotide, and stop nucleotide.
#The output CT then runs from start to stop (inclusive).  Base pairs within the region are kept, base pairs spanning the
# region and other regions are discarded.

import sys

if len(sys.argv)!=5:
    print("Usage: subsetCT inputCT outputCT start stop")
    sys.exit()

#open the input CT
inputct=open(sys.argv[1])

#open the output CT
outputct=open(sys.argv[2],"w")

#process the start and stop
start=int(sys.argv[3])
stop=int(sys.argv[4])

#read and write the first line
firstline=inputct.readline()
firstspace=firstline.find(" ")
outputct.write(str(stop-start+1))
outputct.write(" of ")
outputct.write(firstline[firstspace:])

#now read the nucleotide information
for line in inputct:
    line.rstrip()
    linelist=line.split()

    if int(linelist[0])>=start and int(linelist[0])<=stop:
        #output the index
        string = "%6i" %(int(linelist[0])-start+1)
        outputct.write(string)

        #output the nucleotide
        outputct.write(" ")
        outputct.write(linelist[1])

        #output index -1
        string = "%6i" % (int(linelist[0]) - start)
        outputct.write(string)

        #output index +1
        string = "%6i" % (int(linelist[0]) - start+2)
        outputct.write(string)

        #output the pairs
        if int(linelist[4])>0:
            if  int(linelist[4])>=start and int(linelist[4])<=stop:
                string = "%6i" % (int(linelist[4]) - start + 1)
            else:
                string = "%6i" % (0)
        else:
            string = "%6i" % (0)

        outputct.write(string)

        #output the historical numbering
        string= "%6i" %(int(linelist[5]))
        outputct.write(string)

        outputct.write("\n")


inputct.close()
outputct.close()





