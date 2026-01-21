# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs

# Verify that required input files exist.
verifyInputFiles "$SINGLECT"

EXT=.out  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test raw output
runFullTest 'raw_option'  $SINGLECT -r ---stdout

# Test file output
runFullTest 'file_option'  $SINGLECT -f @OUT.out

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
