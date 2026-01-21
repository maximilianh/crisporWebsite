# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs
EXT=.out  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

SEQ_name=testFiles/testFile_RA7680.fasta
#SEQ_BONUS=testFiles/testFile_5SRNA_tail2.seq

# Test ETEcalculator without options (except ti generate a file).
runFullTest 'without_options'  $SEQ_name -f @OUT.out 

# Test ETEcalculator with an alternative ensemble size
runFullTest 'alternate_ensemble_size' $SEQ_name -f @OUT.out -n 10

# Test ETEcalculator with just the raw distance output.
runFullTest 'raw_distance_option' $SEQ_name -f @OUT.out -r

# Test ETEcalculator with a different random number seed.
runFullTest 'alternate_seed'  $SEQ_name -f @OUT.out -s 43564


