# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Verify that required input files exist.
verifyInputFiles "$SINGLESEQ" "$SINGLEPFS"

# Test ProbablePair_without_options.
runFullTest 'without_options'  $SINGLEPFS @OUT.ct

# Test ProbablePair_dna_option.
runFullTest 'dna_option'  $SINGLESEQ @OUT.ct --sequence -d

# Test ProbablePair_sequence_option.
runFullTest 'sequence_option' $SINGLESEQ  $OUT.ct  --sequence  ---ref=without_options

# Test ProbablePair_threshold_option.
runFullTest 'threshold_option'  $SINGLEPFS @OUT.ct -t 0.92

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
