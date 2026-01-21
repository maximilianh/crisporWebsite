# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test DuplexFold_without_options.
runFullTest 'without_options'  $DOUBLESEQ @OUT.ct

# Test DuplexFold_dna_option.
runFullTest 'dna_option'  $DOUBLESEQ @OUT.ct -d

# Test DuplexFold_loop_option.
runFullTest 'loop_option'  $DOUBLESEQ @OUT.ct -l 10

# Test DuplexFold_max_structures_option.
runFullTest 'max_structures_option'  $DOUBLESEQ @OUT.ct -m 5

# Test DuplexFold_percent_difference_option.
runFullTest 'percent_difference_option'  $DOUBLESEQ @OUT.ct -p 1

# Test DuplexFold_temperature_option.
runFullTest 'temperature_option'  $DOUBLESEQ @OUT.ct -t 150

# Test DuplexFold_window_size_option.
runFullTest 'window_size_option'  $DOUBLESEQ @OUT.ct -w 5

endTestBlock # End a group of tests. Also cleans up orphaned files etc.

