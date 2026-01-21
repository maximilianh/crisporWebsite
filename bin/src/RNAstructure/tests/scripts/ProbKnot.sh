# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Verify that required input files exist.
verifyInputFiles "$SINGLESEQ2" "$SINGLEPFS2" "$ENSEMBLECT"

# Test ProbKnot_without_options.
runFullTest 'without_options'  $SINGLEPFS2 @OUT.ct

# Test ProbKnot_dna_option.
runFullTest 'dna_option'  $SINGLESEQ2 @OUT.ct --sequence -d

# Test ProbKnot_ensemble_option.
runFullTest 'ensemble_option'  $ENSEMBLECT @OUT.ct --ensemble

# Test ProbKnot_iterations_option.
runFullTest 'iterations_option'  $SINGLEPFS2 @OUT.ct -i 2

# Test ProbKnot_min_helix_option.
runFullTest 'min_helix_option'  $SINGLEPFS2 @OUT.ct -m 2

# Test the threshold option, i.e. ThreshKnot
runFullTest 'threshold_option' $SINGLEPFS2 @OUT.ct --threshold 0.4

# Test ProbKnot_sequence_option.
runFullTest 'sequence_option' $SINGLESEQ2  $OUT.ct  --sequence  ---ref=without_options

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
