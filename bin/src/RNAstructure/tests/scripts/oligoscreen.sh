# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.rep  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test oligoscreen_without_options.
runFullTest 'without_options'  $OLIGOLIST @OUT.rep

# Test oligoscreen_dna_option.
runFullTest 'dna_option'  $OLIGOLIST @OUT.rep -d

# Test oligoscreen_temperature_option.
runFullTest 'temperature_option'  $OLIGOLIST @OUT.rep -t 150

endTestBlock # End a group of tests. Also cleans up orphaned files etc.

