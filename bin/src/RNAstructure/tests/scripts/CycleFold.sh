# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Set path to CycleFold-specific parameter files
export CYCLEFOLD_DATAPATH=$ROOT_DIR/CycleFold/datafiles/

INPUT_FASTA=$TESTS_DIR/testFiles/testFile_RA7680.fasta
INPUT_SEQ=$TESTS_DIR/testFiles/testFile_RA7680.seq

# Test CycleFold_without_options.
runFullTest 'without_options'  $INPUT_FASTA     ---stdout

# Test CycleFold_seq_option.
runFullTest 'seq_option'       $INPUT_SEQ   -s  ---stdout  ---ref='without_options'

endTestBlock # End a group of tests. Also cleans up orphaned files etc.

