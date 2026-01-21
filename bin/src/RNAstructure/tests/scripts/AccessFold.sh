# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test AccessFold_without_options.
runFullTest 'without_options'  $DOUBLESEQ @OUT.ct

# Test AccessFold_dna_option.
runFullTest 'dna_option'  $DOUBLESEQ @OUT.ct -d

# Test AccessFold_loop_option.
runFullTest 'loop_option'  $DOUBLESEQ @OUT.ct -l 10

# Test AccessFold_max_structures_option.
runFullTest 'max_structures_option'  $DOUBLESEQ @OUT.ct -m 5

# Test AccessFold_percent_difference_option.
runFullTest 'percent_difference_option'  $DOUBLESEQ @OUT.ct -p 1

# Test AccessFold_temperature_option.
runFullTest 'temperature_option'  $DOUBLESEQ @OUT.ct -t 250

# Test AccessFold_window_size_option.
runFullTest 'window_size_option'  $DOUBLESEQ @OUT.ct -w 5

endTestBlock # End a group of tests. Also cleans up orphaned files etc.

