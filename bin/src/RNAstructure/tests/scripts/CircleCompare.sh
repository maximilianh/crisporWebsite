# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs

# Define some input files.
PREDICTED=testFiles/testFile_CircleCompare_predicted.ct
ACCEPTED=testFiles/testFile_CircleCompare_accepted.ct
ACCEPTEDBI=testFiles/testFile_CircleCompare_accepted_bimolecular.ct

# Verify that required input files exist.
verifyInputFiles "$SINGLEPFS3" "$SINGLEPFS3.txt" "$BIMOLCT" "$PREDICTED" "$ACCEPTED" "$ACCEPTEDBI"

# ************* PS TESTS ************************
EXT=.ps  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test CircleCompare_ps_without_options.
runFullTest 'ps_without_options'  $PREDICTED $ACCEPTED  @OUT.ps

# Test CircleCompare_ps_alternative_option.
runFullTest 'ps_alternative_option'  $PREDICTED $ACCEPTED @OUT.ps -a

# Test CircleCompare_ps_bimolecular_option.
runFullTest 'ps_bimolecular_option'  $BIMOLCT $ACCEPTEDBI @OUT.ps -n 2

# Test CircleCompare_ps_exact_option.
runFullTest 'ps_exact_option'  $PREDICTED $ACCEPTED @OUT.ps -e

# Test CircleCompare_ps_file_option.
runFullTest 'ps_file_option'  $PREDICTED $ACCEPTED @OUT.ps -f

# Test CircleCompare_ps_levorotatory_option.
runFullTest 'ps_levorotatory_option'  $PREDICTED $ACCEPTED @OUT.ps -L

# Test CircleCompare_ps_probability_option.
runFullTest 'ps_probability_option'  $PREDICTED $ACCEPTED @OUT.ps -p $SINGLEPFS3

# Test CircleCompare_ps_probability_accepted_option.
runFullTest 'ps_probability_accepted_option'  $PREDICTED $ACCEPTED @OUT.ps -p2 $SINGLEPFS3

# Test CircleCompare_ps_probability_text_option.
runFullTest 'ps_probability_text_option'  $PREDICTED $ACCEPTED @OUT.ps -n 3 -t $SINGLEPFS3.txt

# Test CircleCompare_ps_shape_option.
runFullTest 'ps_shape_option'  $PREDICTED $ACCEPTED @OUT.ps -s testFiles/testFile_ivslsu_dummy.shape

# Test CircleCompare_ps_specific_structure_option.
# Note that other tests require this option; so this is just a sanity check to make sure other specific structures can be selected.
runFullTest 'ps_specific_structure_option'  $PREDICTED $ACCEPTED @OUT.ps -n 2

# Test CircleCompare_ps_uncircled_option.
runFullTest 'ps_uncircled_option'  $PREDICTED $ACCEPTED @OUT.ps -u

# ************* SVG TESTS ************************
EXT=.svg  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test CircleCompare_svg_without_options.
runFullTest 'svg_without_options'  $PREDICTED $ACCEPTED @OUT.svg -n 3 --svg

# Test CircleCompare_svg_alternative_option.
runFullTest 'svg_alternative_option'  $PREDICTED $ACCEPTED @OUT.svg -a -n 3 --svg

# Test CircleCompare_svg_bimolecular_option.
runFullTest 'svg_bimolecular_option'  $BIMOLCT $ACCEPTEDBI @OUT.svg -n 2 --svg

# Test CircleCompare_svg_exact_option.
runFullTest 'svg_exact_option'  $PREDICTED $ACCEPTED @OUT.svg -e -n 3 --svg

# Test CircleCompare_svg_file_option.
runFullTest 'svg_file_option'  $PREDICTED $ACCEPTED @OUT.svg -f -n 3 --svg

# Test CircleCompare_svg_levorotatory_option.
runFullTest 'svg_levorotatory_option'  $PREDICTED $ACCEPTED @OUT.svg -L -n 3 --svg

# Test CircleCompare_svg_probability_option.
runFullTest 'svg_probability_option'  $PREDICTED $ACCEPTED @OUT.svg -p $SINGLEPFS3 -n 3 --svg

# Test CircleCompare_svg_probability_accepted_option.
runFullTest 'svg_probability_accepted_option'  $PREDICTED $ACCEPTED @OUT.svg -p2 $SINGLEPFS3 -n 3 --svg

# Test CircleCompare_svg_probability_text_option.
runFullTest 'svg_probability_text_option'  $PREDICTED $ACCEPTED @OUT.svg -t $SINGLEPFS3.txt -n 3 --svg

# Test CircleCompare_svg_shape_option.
runFullTest 'svg_shape_option'  $PREDICTED $ACCEPTED @OUT.svg -s testFiles/testFile_ivslsu_dummy.shape -n 3 --svg

# Test CircleCompare_svg_specific_structure_option.
# Note that other tests require this option; so this is just a sanity check to make sure other specific structures can be selected.
runFullTest 'svg_specific_structure_option'  $PREDICTED $ACCEPTED @OUT.svg -n 2 --svg

# Test CircleCompare_svg_uncircled_option.
runFullTest 'svg_uncircled_option'  $PREDICTED $ACCEPTED @OUT.svg -u -n 3 --svg

endTestBlock # End a group of tests. Also cleans up orphaned files etc.