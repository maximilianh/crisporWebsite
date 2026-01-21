# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.out  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Run the unit-test program, which will output to cerr if there are problems.
runFullTest 'all' ---nodiff ---exit=0

endTestBlock