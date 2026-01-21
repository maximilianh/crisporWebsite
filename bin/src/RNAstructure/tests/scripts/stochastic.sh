# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Verify that required input files exist.
verifyInputFiles "$SINGLESEQ" "$SINGLEPFS"

# Test stochastic_without_options.
runFullTest 'without_options'  $SINGLEPFS @OUT.ct

# Test stochastic_dna_option.
runFullTest 'dna_option'  $SINGLESEQ @OUT.ct --sequence -d

# Test stochastic_ensemble_option.
runFullTest 'ensemble_option'  $SINGLEPFS @OUT.ct -e 2

# Test stochastic_seed_option.
runFullTest 'seed_option'  $SINGLEPFS @OUT.ct -s 2

# Test stochastic_sequence_option.
runFullTest 'sequence_option' $SINGLESEQ  $OUT.ct  --sequence \
	---ref='without_options'

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
