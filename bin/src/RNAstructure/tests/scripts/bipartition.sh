# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE ProbabilityPlot  # run `make` for the main exe and any other listed programs
EXT=.ps  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)
FIX_DESC='--desc @TBASE_output.pfs'  # the output ps files contain the name of the input file, which changes (e.g. partition vs partition-smp) to correct for this, we specify a description to use.

# Test bipartition_without_options.
initTest 'without_options'  && {
runTest @EXE $DOUBLESEQ  @OUT.pfs 
runSubTest 'plot' ProbabilityPlot @OUT.pfs @OUT.ps $FIX_DESC
runDiff; endTest; }

# Test bipartition_dna_option.
initTest 'dna_option'  && {
runTest @EXE $DOUBLESEQ  @OUT.pfs  -d
runSubTest 'plot' ProbabilityPlot @OUT.pfs @OUT.ps $FIX_DESC
runDiff; endTest; }

# Test bipartition_temperature_option.
initTest 'temperature_option'  && {
runTest @EXE $DOUBLESEQ  @OUT.pfs  -t 330
runSubTest 'plot' ProbabilityPlot @OUT.pfs @OUT.ps $FIX_DESC
runDiff; endTest; }

endTestBlock # End a group of tests. Also cleans up orphaned files etc.