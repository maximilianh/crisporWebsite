# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester # run `make` for the Java GUI Tester and any other listed programs.

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function MaxExpect-gui() {
  local minArgs=2
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local inFile="$1" ctFile="$2" ACID=RNA
  BOOL_FLAGS=(d sequence)  # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(g p s w)    # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)
  
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiArg -d && ACID=DNA
  guiWrite "Menu '$ACID', 'MaxExpect:*'
    Click { BUTTON '*Save File' }
    TypeText '$inFile'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD  'CT File' }, '$ctFile'"

  # guiWriteIf writes lines to the gui script IF the argument specified by the first parameter is present.
  guiWriteIf -g "FOCUS { FIELD 'Gamma' }; ClearText; TypeText '$arg_g'"
  guiWriteIf -p "Focus { Field 'Max % Score Difference' };  ClearText; TypeText '$arg_p'"
  guiWriteIf -s "Focus { Field 'Max Number of Structures' }; ClearText; TypeText '$arg_s'"
  guiWriteIf -w "Focus { FIELD 'Window Size' }; ClearText; TypeText '$arg_w'"

  guiWrite "CLICK { BUTTON 'START' }
    Prohibit { DIALOG 'RNAstructure error' } 
    WAIT  { BUTTON 'OK' }, ${CALC_TIMEOUT}000 
    CLICK { BUTTON 'Cancel' }"

  guiRun -c # Run the script (-c causes cleanup of $arg_* variables)
}

# Verify that required input files exist.
verifyInputFiles "$SINGLEPFS4"

EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test MaxExpect_without_options.
runFullTest 'without_options'  $SINGLEPFS4 @OUT.ct

# Test MaxExpect_dna_option.
#Test-TODO:runFullTest 'dna_option'  $SINGLESEQ4 @OUT.ct --sequence -d

# Test MaxExpect_gamma_option.
runFullTest 'gamma_option'  $SINGLEPFS4 @OUT.ct -g 2

# Test MaxExpect_max_structures_option.
runFullTest 'max_structures_option'  $SINGLEPFS4 @OUT.ct -w 0 -s 3

# Test MaxExpect_percent_difference_option.
runFullTest 'percent_difference_option'  $SINGLEPFS4 @OUT.ct -w 0 -p 2

# Test MaxExpect_sequence_option.
#Test-TODO: runFullTest 'sequence_option' $SINGLESEQ4  $OUT.ct  --sequence \
#     --ref 'without_options'

# Test MaxExpect_window_size_option.
runFullTest 'window_size_option'  $SINGLEPFS4 @OUT.ct -w 0

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
