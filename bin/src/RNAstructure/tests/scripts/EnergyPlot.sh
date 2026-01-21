# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs

verifyInputFiles "$SINGLESAV"
SAVEFILE=$SINGLESAV
MINIMUM="-30"
MAXIMUM="-20"

# APPEND_TEST_ARGS causes the value of $SAVEFILE to be replaced by the text "[FILENAME]" in the output file.
# This allows the test author to change the location of SAVEFILE without changing the contents of the reference files.
# (Otherwise the path of the input file is embedded in the PS or SVG files)
APPEND_TEST_ARGS="---replace=$SAVEFILE=>[FILENAME]" # runFullTest will append these arguments in its call to "runProgAndDiff"

############# PS TESTS #################
EXT=.ps  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test EnergyPlot_ps_without_options.
runFullTest 'ps_without_options'  $SAVEFILE @OUT.ps

# Test EnergyPlot_ps_entries_option.
runFullTest 'ps_entries_option'  $SAVEFILE @OUT.ps -e 3

# Test EnergyPlot_ps_maximum_option.
runFullTest 'ps_maximum_option'  $SAVEFILE @OUT.ps -max $MAXIMUM

# Test EnergyPlot_ps_minimum_option.
runFullTest 'ps_minimum_option'  $SAVEFILE @OUT.ps -min $MINIMUM

############# SVG TESTS #################
EXT=.svg  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test EnergyPlot_svg_without_options.
runFullTest 'svg_without_options'  $SAVEFILE @OUT.svg --svg

# Test EnergyPlot_svg_entries_option.
runFullTest 'svg_entries_option'  $SAVEFILE @OUT.svg -e 3 --svg

# Test EnergyPlot_svg_maximum_option.
runFullTest 'svg_maximum_option'  $SAVEFILE @OUT.svg -max $MAXIMUM --svg

# Test EnergyPlot_svg_minimum_option.
runFullTest 'svg_minimum_option'  $SAVEFILE @OUT.svg -min $MINIMUM --svg

############# TXT TESTS #################
EXT=.txt  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test EnergyPlot_text_option.
runFullTest 'text_option'  $SAVEFILE @OUT.txt --text

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
