# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs


#  Usage: exe/partition-cuda [options] <sequence, or file containing one>
#  options:
#  -h:        show this message
#  -b <file>: read parameters from <file>, in native binary format
#  -d:        use DNA parameters
#  -t <file>: write probability matrix as text to <file>
#  -l <file>: write -log10 probabilities as text to <file>
#  -p <file>: write ProbKnot structure in ct format to <file>
#  -m <length>: set minimum helix length for ProbKnot
#               (default: 3 base pairs)
#  -v:        show arrays
#  If none of -t, -l, -p, or -v is chosen,
#  writes ProbKnot structure in ct format to stdout


EXT=.ct  # set extension for OK files 

# Test partition-cuda_without_options.
runFullTest 'without_options' time/ivslsu.seq ---stdout # write stdout to @OUT.ct

# Test partition-cuda ProbKnot output
runFullTest 'probknot_ct' time/ivslsu.seq -p @OUT.ct ---ref='without_options' # output should be the same as the 'without_options' test.


EXT=.txt  # set extension for OK files 

# -t option
runFullTest 'probability_matrix' time/ivslsu.seq -t @OUT.txt

# -l option
runFullTest 'pairwise_probabilities' time/ivslsu.seq -l @OUT.txt

# -l option with a sequence that has lowercase bases (constrained single-stranded) at positions 17 34 37 and 58
runFullTest 'pairwise_prob_constraints' testFiles/testFile_RA7680.fasta -l @OUT.txt


endTestBlock # End a group of tests. Also cleans up orphaned files etc.
