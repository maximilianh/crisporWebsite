# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs
EXT=.out  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

SEQ1=$INDIR/testFile_RA7680.seq; PFS1=${SEQ1%.seq}.pfs
SEQ2=$INDIR/testFile_ivslsu.seq
SEQ3=$INDIR/testFile_ec5s.seq

# Verify that required input files exist.
verifyInputFiles $SEQ1 $SEQ2 $SEQ3 $PFS1

# Test EnsembleEnergy_without_options.
runFullTest 'without_options' $PFS1 ---stdout   # (Test output is redirected to @STDO)

# Test EnsembleEnergy_dna_option.
runFullTest 'dna_option'  $SEQ1 --sequence -d   ---stdout # (Test output is redirected to @STDO)

# Test EnsembleEnergy_sequence_option.
runFullTest 'sequence_option'  $SEQ1 --sequence    ---stdout # (Test output is redirected to @STDO)

# Alternate sequence
# runFullTest 'sequence_alt2'  $SEQ2 --sequence    ---stdout # (Test output is redirected to @STDO)

# Another Alternate sequence
# runFullTest 'sequence_alt3'  $SEQ3 --sequence    ---stdout # (Test output is redirected to @STDO)

# Test EnsembleEnergy_silent_option.
runFullTest 'silent_option'  $PFS1 --silent    ---stdout # (Test output is redirected to @STDO)

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
