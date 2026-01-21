# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r $EXE ProbabilityPlot # run `make` for Rsample and ProbabilityPlot
EXT=.txt  # set extension for @OKFILE

# Test partition_without_options.
initTest 'without_options' && {  # skip this block if the test is excluded.
runTest @EXE testFiles/tRNA.seq testFiles/tRNA.shape @OUT.pfs -s 1
runSubTest 'plot' ProbabilityPlot -t @OUT.pfs @OUT.txt
runDiff; endTest; }

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
