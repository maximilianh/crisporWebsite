# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs

# Verify that required input files exist.
verifyInputFiles "$SINGLECT" "$BIMOLCT" "$KNOTSCT" "$SINGLEPFS" "$SINGLEPFS.txt"

############### PS Tests ###############
EXT=.ps  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test draw_ps_without_options.
runFullTest 'ps_without_options'  $SINGLECT @OUT.ps

# Test draw_ps_bimolecular_option.
runFullTest 'ps_bimolecular_option'  $BIMOLCT @OUT.ps

# Test draw_ps_bimolecular_circular_option.
runFullTest 'ps_bimolecular_circular_option'  $BIMOLCT @OUT.ps -c

# Test draw_ps_circular_option.
runFullTest 'ps_circular_option'  $SINGLECT @OUT.ps -c

# Test draw_ps_circular_not_specified option.
runFullTest 'ps_circular_not_specified_option'  $KNOTSCT @OUT.ps -c

# Test draw_ps_flat_option.
runFullTest 'ps_flat_option'  $SINGLECT @OUT.ps -f

# Test draw_ps_levorotatory_option.
runFullTest 'ps_levorotatory_option'  $SINGLECT @OUT.ps -L

# Test draw_ps_probability_option.
runFullTest 'ps_probability_option'  $SINGLECT @OUT.ps -p $SINGLEPFS

# Test draw_ps_probability_text_option.
runFullTest 'ps_probability_text_option'  $SINGLECT @OUT.ps -t $SINGLEPFS.txt -n 4

# Test draw_ps_shape_option.
runFullTest 'ps_shape_option'  $SINGLECT @OUT.ps -s testFiles/testFile_tRNA.shape

# Test draw_ps_specific_structure_option.
# Note that other tests require this option; so this is just a sanity check to make sure other specific structures can be selected.
runFullTest 'ps_specific_structure_option'  $SINGLECT @OUT.ps -n 3

# Test draw_ps_uncircled_option.
runFullTest 'ps_uncircled_option'  $SINGLECT @OUT.ps -u

# Test draw_ps_custom_description.
runFullTest 'ps_custom_description'  $SINGLECT @OUT.ps  --desc "A nice drawing."

# Test draw_ps_no_description.
runFullTest 'ps_no_description'  $SINGLECT @OUT.ps  --desc ''

# Test number range
runFullTest 'ps_number_range'  $SINGLECT @OUT.ps --number 2 --end 3

# Test limit maximum structures
runFullTest 'ps_limit_structures'  $SINGLECT @OUT.ps --maxstructures 2

# Test draw_ps_filename_description.
#Skip this test because the ct file doesn't have file names.
#echo '    draw_ps_filename_description testing started...'
#../exe/draw $SINGLECT draw_ps_filename_description_test_output.ps  --desc ~file  1>/dev/null 2>draw_ps_filename_description_errors.txt
#$DIFF_CMD draw_ps_filename_description_test_output.ps draw/draw_ps_filename_description_OK.ps &> draw_ps_filename_description_diff_output.txt
#checkErrorFiles draw_ps_filename_description draw_ps_filename_description_errors.txt draw_ps_filename_description_diff_output.txt
#echo '    draw_ps_filename_description testing finished.'

# Test draw_ps_multiple_descriptions.
runFullTest 'ps_multiple_descriptions'  $SINGLECT @OUT.ps  --desc '~list,Structure 1,Another One,Symbols !@#$%^&*()-=+,Hello World' 

############### SVG Tests ###############
EXT=.svg  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test draw_svg_without_options.
runFullTest 'svg_without_options'  $SINGLECT @OUT.svg -n 4 --svg

# Test draw_svg_bimolecular_option.
runFullTest 'svg_bimolecular_option'  $BIMOLCT @OUT.svg -n 4 --svg

# Test draw_svg_bimolecular_circular_option.
runFullTest 'svg_bimolecular_circular_option'  $BIMOLCT @OUT.svg -c -n 4 --svg

# Test draw_svg_circular_option.
runFullTest 'svg_circular_option'  $SINGLECT @OUT.svg -c -n 4 --svg

# Test draw_svg_circular_not_specified option.                                                                                                                                                    
runFullTest 'svg_circular_not_specified_option'  $KNOTSCT @OUT.svg -c -n 1 --svg

# Test draw_svg_flat_option.
runFullTest 'svg_flat_option'  $SINGLECT @OUT.svg -f -n 4 --svg

# Test draw_svg_levorotatory_option.
runFullTest 'svg_levorotatory_option'  $SINGLECT @OUT.svg -L -n 4 --svg

# Test draw_svg_probability_option.
runFullTest 'svg_probability_option'  $SINGLECT @OUT.svg -p $SINGLEPFS -n 4 --svg

# Test draw_svg_probability_text_option.
runFullTest 'svg_probability_text_option'  $SINGLECT @OUT.svg -t $SINGLEPFS.txt -n 4 --svg

# Test draw_svg_shape_option.
runFullTest 'svg_shape_option'  $SINGLECT @OUT.svg -s testFiles/testFile_tRNA.shape -n 4 --svg

# Test draw_svg_specific_structure_option.
# Note that other tests require this option; so this is just a sanity check to make sure other specific structures can be selected.
runFullTest 'svg_specific_structure_option'  $SINGLECT @OUT.svg -n 3 --svg

# Test draw_svg_uncircled_option.
runFullTest 'svg_uncircled_option'  $SINGLECT @OUT.svg -u -n 4 --svg

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
