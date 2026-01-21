# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Verify that required input files exist.
verifyInputFiles "$SINGLESAV"

# Test refold_without_options.
runFullTest 'without_options'  $SINGLESAV @OUT.ct

# Test refold_max_structures_option.
runFullTest 'max_structures_option'  $SINGLESAV @OUT.ct -m 3

# Test refold_percent_difference_option.
runFullTest 'percent_difference_option'  $SINGLESAV @OUT.ct -p 5

# Test refold_window_size_option.
runFullTest 'window_size_option'  $SINGLESAV @OUT.ct -w 20

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
