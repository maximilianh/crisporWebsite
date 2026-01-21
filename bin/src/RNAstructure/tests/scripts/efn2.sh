# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs

# Verify that required input files exist.
verifyInputFiles "$SINGLECT"

EXT=.out  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

NO_COUNT=--ne
# Note: NO_COUNT (--ne) was added to all tests to disable the experimental error 
#  output (which does not work for SMP mode) so that serial and SMP test output 
#  files would be compatible.
# Additional serial-only tests were added to track experimental error output.

# Test efn2_without_options.
runFullTest 'without_options'  $SINGLECT @OUT.out  $NO_COUNT

# Test efn2_dna_option.
runFullTest 'dna_option'  $SINGLECT @OUT.out -d  $NO_COUNT

# Test efn2_alphabet_option.
runFullTest 'alphabet_option'  $SINGLECT @OUT.out -a 'dna' $NO_COUNT ---ref='dna_option'

# Test efn2_alphabet_rna_option.
runFullTest 'alphabet_option_rna'  $SINGLECT @OUT.out -a 'rna' $NO_COUNT ---ref='without_options'

# Test efn2_print_option.
runFullTest 'print_option' $SINGLECT @OUT.out -p $NO_COUNT ---ref='without_options' \
            ---add-diff=screen:@STDO==@OKFILE # also verify STDOUT output

# Test efn2_shape_option.
runFullTest 'shape_option'  $SINGLECT @OUT.out -sh testFiles/testFile_tRNA.shape  $NO_COUNT

# Test efn2_shape_intercept_option.
runFullTest 'shape_intercept_option'  $SINGLECT @OUT.out -sh testFiles/testFile_tRNA.shape -si 0.2  $NO_COUNT

# Test efn2_shape_slope_option.
runFullTest 'shape_slope_option'  $SINGLECT @OUT.out -sh testFiles/testFile_tRNA.shape -sm 1.2  $NO_COUNT
# efn2/efn2_shape_slope_option_OK.out

# Test efn2_temperature_option.
runFullTest 'temperature_option'  $SINGLECT @OUT.out -t 150  $NO_COUNT

# Test efn2_write_thermodynamic_file_option.
runFullTest 'write_thermodynamic_file_option'  $SINGLECT @OUT.out -w  $NO_COUNT

# Test efn2_pseudoknot.
runFullTest 'knotted'  testFiles/testFile_knotted.ct @OUT.out  $NO_COUNT

# Test efn2_pseudoknot with write option.
runFullTest 'knotted_write'  testFiles/testFile_knotted.ct @OUT.out -w  $NO_COUNT

if [[ ! $SMP ]]; then
    # experrs tests include calculation of experimental errors (aka 'counts')
    runFullTest 'calc_errors'  $SINGLECT @OUT.out
    runFullTest 'dna_calc_errors'  $SINGLECT @OUT.out -d  ---ref='dna_option'
    runFullTest 'alphabet_calc_errors'  $SINGLECT @OUT.out -a 'dna'  ---ref='dna_option'
    runFullTest 'alphabet_calc_errors_rna'  $SINGLECT @OUT.out -a 'rna'  ---ref='calc_errors'
    runFullTest 'temperature_calc_errors'  $SINGLECT @OUT.out -t 150
    runFullTest 'print_option_calc_errors' $SINGLECT @OUT.out -p ---ref='calc_errors' \
            ---add-diff=screen:@STDO==@OKFILE # also verify STDOUT output
    runFullTest 'calc_errors_counts'  $SINGLECT @OUT.out --count @OUT.txt ---ext=.txt
fi

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
