#fold an RNA
import RNAstructure
import sys


def checkState(rna):
    if rna.GetErrorCode() != 0:
        print rna.GetFullErrorMessage()
        sys.exit()

def checkResult(ret):
    if ret != 0:
        print RNAstructure.RNA.GetErrorMessage(ret);
        sys.exit(1)


if len(sys.argv) != 3:
	print "Usage: python example.py <input_sequence_file> <output_ct_file>"
	sys.exit(1)

#argument 1 should be an input sequence (in seq or fasta format)
#argument 2 should be the output ct file
infile,outfile = sys.argv[1:]

#one annoying thing about the RNAstructure constructors is
#that they take a lot of boolean or integer arguments
#One thing I intend to do is add sane default arguments. 
#suggestions are welcome

#also, SWIG doesn't work nicely with overloaded functions
#in a lot of cases

#constructors must be called with all arguments
#and can't use keyword arguments

print "Reading input sequence file ", infile
strand = RNAstructure.RNA(infile,2,True)

#constructor fails silently if the input is invalid!!
#have to check error codes -- not very pythonic
checkState(strand);

#likewise for most of the methods. Something like 
#FoldSingleStrand should have the same defaults
#as the Fold program
print "Folding sequence..."
strand.FoldSingleStrand(percent=0.0,
						maximumstructures=1,
						window=0,
						mfeonly=True)
checkState(strand);

print "Writing CT..."
checkResult(strand.WriteCt(outfile))
checkState(strand);

print "Done."
