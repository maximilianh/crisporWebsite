# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

#   The draw program produces multi-page postscript files from a single CT 
# input file (which can contain multiple structures).
#   In contrast, the Java GUI only saves a single-page postscript file from
# the currently-selected structure. In order to make the tests compatible 
# for the draw-gui tests, the OK files are copied into the working directory
# and modified to remove all pages after the first one.  
# We need to instruct the test system to look in the working directory instead 
# of in draw/ for the OK files:
OKDIR=. 

source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester  # run `make` for the Java GUI Tester and any other listed programs.

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function draw-gui() {
  local minArgs=2
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local ctFile="$1" outFile="$2" savetype='Postscript'
  BOOL_FLAGS=(c f L svg u)         # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(desc lp n p s t)    # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)
  
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiWrite "Menu 'File', 'Draw'
    TypeText '$ctFile'; TypeENTER
    PROHIBIT { DIALOG 'RNAstructure error' } 
    WAIT { Type='DrawingWindow' }, 10000"
  
  guiWriteIf -n "Menu 'Draw', 'Go to Structure*'
    Focus { FIELD 'Structure Number' }
    ClearText; TypeText '$arg_n'
    Click { BUTTON 'OK' }"

  guiWriteIf -c "Menu 'Draw', 'Render Circular*'"
  guiWriteIf -f "Menu 'Draw', 'Render Linear*'"
  guiWriteIf -L "Menu 'Draw', 'Render *Flipped*'"
  guiWriteIf -u "Menu 'Draw', 'Render *Circles*'"

  guiWriteIf -p "Menu 'Annotations', 'Add Probability Annotation'
    TypeText '$arg_p'; TypeEnter
    Prohibit {DIALOG 'RNAstructure error'}"
  guiWriteIf -s "Menu 'Annotations', 'Add SHAPE Annotation'
    TypeText '$arg_s'; TypeEnter
    Prohibit {DIALOG 'RNAstructure error'}"

  guiArg --svg && savetype=SVG
  guiWrite "Menu 'Draw', 'Write $savetype File'
    TypeText '$outFile'; TypeENTER
    PROHIBIT { DIALOG 'RNAstructure error' } 
    WAIT { Button 'OK' }, 10000
    CLICK { Button 'OK' }"
 
  #local postProcess=removeExtraPages
  #guiArg -n && postProcess=cat # Do no post processing. just copy the OK file.
  removeExtraPages <draw/${OKFILE#./}  >$OKFILE

  guiRun -c # Run the script (-c causes cleanup of $arg_* variables)
}

function removeExtraPages() {
  sed -E -e 's/^%%Page: ([0-9]+) ([0-9]+)$/%%Page: 1 1/g' -e '/showpage/q'
}

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
# Test-TODO: runFullTest 'ps_probability_text_option'  $SINGLECT @OUT.ps -t $SINGLEPFS.txt -n 4

# Test draw_ps_shape_option.
runFullTest 'ps_shape_option'  $SINGLECT @OUT.ps -s testFiles/testFile_tRNA.shape

# Test draw_ps_specific_structure_option.
# Note that other tests require this option; so this is just a sanity check to make sure other specific structures can be selected.
runFullTest 'ps_specific_structure_option'  $SINGLECT @OUT.ps -n 3

# Test draw_ps_uncircled_option.
runFullTest 'ps_uncircled_option'  $SINGLECT @OUT.ps -u

# Test draw_ps_custom_description.
# TEST-TODO runFullTest 'ps_custom_description'  $SINGLECT @OUT.ps  --desc "A nice drawing."

# Test draw_ps_no_description.
# TEST-TODO runFullTest 'ps_no_description'  $SINGLECT @OUT.ps  --desc ''

# Test draw_ps_multiple_descriptions.
# TEST-TODO runFullTest 'ps_multiple_descriptions'  $SINGLECT @OUT.ps  --desc '~list,Structure 1,Another One,Symbols !@#$%^&*()-=+,Hello World' 

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
#Test-TODO: runFullTest 'svg_probability_text_option'  $SINGLECT @OUT.svg -t $SINGLEPFS.txt -n 4 --svg

# Test draw_svg_shape_option.
runFullTest 'svg_shape_option'  $SINGLECT @OUT.svg -s testFiles/testFile_tRNA.shape -n 4 --svg

# Test draw_svg_specific_structure_option.
# Note that other tests require this option; so this is just a sanity check to make sure other specific structures can be selected.
runFullTest 'svg_specific_structure_option'  $SINGLECT @OUT.svg -n 3 --svg

# Test draw_svg_uncircled_option.
runFullTest 'svg_uncircled_option'  $SINGLECT @OUT.svg -u -n 4 --svg

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
