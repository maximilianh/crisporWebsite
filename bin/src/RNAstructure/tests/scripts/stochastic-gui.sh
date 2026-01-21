# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester # run `make` for the Java GUI Tester and any other listed programs.

# Verify that required input files exist.
verifyInputFiles "$SINGLESEQ" "$SINGLEPFS"

# Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function stochastic-gui() {
  local minArgs=2
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local savFile="$1" ctFile="$2" ACID=RNA
  BOOL_FLAGS=(d sequence)                    # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(e s)    # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)
  
  guiUnsetArgs # Clean up any temporary variables remaining from previous test
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiArg -d && ACID=DNA
  guiWrite "Menu '$ACID', 'Stochastic $ACID Sampling'
    Click { BUTTON 'Partition Function Save File' }
    TypeText '$savFile'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD  'CT File' }, '$ctFile'"

  # guiWriteIf writes lines to the gui script IF the argument specified by the first parameter is present.
  guiWriteIf -s "Focus { FIELD 'Random Seed' };    ClearText; TypeText '$arg_s'"
  guiWriteIf -e "Focus { FIELD 'Ensemble Size' };  ClearText; TypeText '$arg_e'"

  guiWrite "CLICK { BUTTON 'START' }
    Prohibit { DIALOG 'RNAstructure error' } 
    WAIT  { BUTTON 'OK' }, ${CALC_TIMEOUT}000 
    CLICK { BUTTON 'Cancel' }"

  guiRun -c # Run the script & clean up temp variables.
}

EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test stochastic_without_options.
runFullTest 'without_options'  $SINGLEPFS @OUT.ct

# Test stochastic_dna_option.
#Test-TODO: runFullTest 'dna_option'  $SINGLESEQ @OUT.ct --sequence -d

# Test stochastic_ensemble_option.
runFullTest 'ensemble_option'  $SINGLEPFS @OUT.ct -e 2

# Test stochastic_seed_option.
runFullTest 'seed_option'  $SINGLEPFS @OUT.ct -s 2

# Test stochastic_sequence_option.
#Test-TODO: runFullTest 'sequence_option' $SINGLESEQ  $OUT.ct  --sequence \
#	---ref='without_options'

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
