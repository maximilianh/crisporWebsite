# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Verify that required input files exist.
verifyInputFiles "$SINGLEPFS4"

# Test MaxExpect_without_options.
runFullTest 'without_options'  $SINGLEPFS4 @OUT.ct

# Test MaxExpect_dna_option.
runFullTest 'dna_option'  $SINGLESEQ4 @OUT.ct --sequence -d

# Test MaxExpect_alphabet_option.
runFullTest 'alphabet_option'  $SINGLESEQ4 @OUT.ct --sequence -a 'dna' ---ref='dna_option'

# Test MaxExpect_gamma_option.
runFullTest 'gamma_option'  $SINGLEPFS4 @OUT.ct -g 2

# Test MaxExpect_max_structures_option.
runFullTest 'max_structures_option'  $SINGLEPFS4 @OUT.ct -w 0 -s 3

# Test MaxExpect_percent_difference_option.
runFullTest 'percent_difference_option'  $SINGLEPFS4 @OUT.ct -w 0 -p 2

# Test MaxExpect_sequence_option.
runFullTest 'sequence_option' $SINGLESEQ4 @OUT.ct  --sequence  ---ref=without_options # Use the same "OK" files as the without_options test. 

# Test MaxExpect_window_size_option.
runFullTest 'window_size_option'  $SINGLEPFS4 @OUT.ct -w 0

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
