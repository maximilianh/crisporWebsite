# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.txt  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test scorer without_options.
runFullTest 'without_options'  $DOUBLECT2 @OUT.txt

# Test scorer exact_option.
runFullTest 'exact_option'  $DOUBLECT2 @OUT.txt -e

# Test scorer_multiple_cts.
runFullTest 'multiple_cts'  testFiles/testFile_CircleCompare_predicted.ct testFiles/testFile_knotted.ct @OUT.txt

# Test scorer_print_option.
runFullTest 'print_option'  $DOUBLECT2 @OUT.txt -p \
	---ref='without_options'

# Test scorer_print_option_screen.
runFullTest 'print_option_screen' $DOUBLECT2 @OUT.txt -p  ---stdout # use stdout instead of output file.

# Test scorer_specific_structure_option.
runFullTest 'specific_structure_option'  $DOUBLECT @OUT.txt -n 2

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
