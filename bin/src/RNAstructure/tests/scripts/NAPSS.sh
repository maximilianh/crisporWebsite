# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

TEST1=testFiles/adhr1
TEST2=testFiles/bm

runFullTest 'adhr1_d1'  $TEST1.seq $TEST1.nmr-con @OUT.ct  -d 1
runFullTest 'bm_d10'  $TEST2.seq $TEST2.nmr-con @OUT.ct  -d 10

endTestBlock # End a group of tests. Also cleans up orphaned files etc.

