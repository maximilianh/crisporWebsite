# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs

SAVEFILE=$SINGLEPFS # a temporary file.
MINIMUM=3
MAXIMUM=7

# REPLACE_NAME causes the value of $SAVEFILE to be replaced by the text "[FILENAME]" in the output file.
# This allows the test author to change the location of SAVEFILE without changing the contents of the reference files.
# (Otherwise the path of the input file is embedded in the PS or SVG files)
REPLACE_NAME="---replace=$SAVEFILE=>[FILENAME]"

############ PS TESTS ###############
EXT=.ps  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test ProbabilityPlot_ps_without_options.
runFullTest 'ps_without_options'  $SAVEFILE @OUT.ps   "$REPLACE_NAME"

# Test ProbabilityPlot_ps_entries_option.
runFullTest 'ps_entries_option'  $SAVEFILE @OUT.ps -e 3 "$REPLACE_NAME"

# Test ProbabilityPlot_ps_log_option.
runFullTest 'ps_log_option'  ProbabilityPlot/ProbabilityPlot_text_option_OK.txt @OUT.ps --log10

# Test ProbabilityPlot_ps_matrix_option.
runFullTest 'ps_matrix_option'  testFiles/testFile_matrix.txt @OUT.ps --matrix

# Test ProbabilityPlot_ps_maximum_option.
runFullTest 'ps_maximum_option'  $SAVEFILE @OUT.ps -max $MAXIMUM  "$REPLACE_NAME"

# Test ProbabilityPlot_ps_minimum_option.
runFullTest 'ps_minimum_option'  $SAVEFILE @OUT.ps -min $MINIMUM  "$REPLACE_NAME"

############ SVG TESTS ###############
EXT=.svg  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test ProbabilityPlot_svg_without_options.
runFullTest 'svg_without_options'  $SAVEFILE @OUT.svg --svg  "$REPLACE_NAME"

# Test ProbabilityPlot_svg_entries_option.
runFullTest 'svg_entries_option'  $SAVEFILE @OUT.svg -e 3 --svg  "$REPLACE_NAME"

# Test ProbabilityPlot_svg_log_option.
runFullTest 'svg_log_option'  ProbabilityPlot/ProbabilityPlot_text_option_OK.txt @OUT.svg --log10 --svg

# Test ProbabilityPlot_svg_matrix_option.
runFullTest 'svg_matrix_option'  testFiles/testFile_matrix.txt @OUT.svg --matrix --svg

# Test ProbabilityPlot_svg_maximum_option.
runFullTest 'svg_maximum_option'  $SAVEFILE @OUT.svg -max $MAXIMUM --svg  "$REPLACE_NAME"

# Test ProbabilityPlot_svg_minimum_option.
runFullTest 'svg_minimum_option'  $SAVEFILE @OUT.svg -min $MINIMUM --svg  "$REPLACE_NAME"

############ TXT TESTS ###############
EXT=.txt  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test ProbabilityPlot_text_option.
runFullTest 'text_option'  $SAVEFILE @OUT.txt --text  "$REPLACE_NAME"

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
