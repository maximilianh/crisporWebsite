# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester # run `make` for the Java GUI Tester and any other listed programs.

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function ProbKnot-gui() {
  local minArgs=2
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local inFile="$1" ctFile="$2" ACID=RNA
  BOOL_FLAGS=(d enseble sequence)  # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(i m)    # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)
  
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiArg -d && ACID=DNA
  guiWrite "Menu '$ACID', 'ProbKnot:*'
    Click { BUTTON '*Save File' }
    TypeText '$inFile'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD  'CT File' }, '$ctFile'"

  # guiWriteIf writes lines to the gui script IF the argument specified by the first parameter is present.
  guiWriteIf -i "FOCUS { FIELD 'Iterations' }; ClearText; TypeText '$arg_i'"
  guiWriteIf -m "Focus { Field 'Minimum Helix Length' };  ClearText; TypeText '$arg_m'"

  guiWrite "CLICK { BUTTON 'START' }
    Prohibit { DIALOG 'RNAstructure error' } 
    WAIT  { BUTTON 'OK' }, ${CALC_TIMEOUT}000 
    CLICK { BUTTON 'Cancel' }"

  guiRun -c # Run the script (-c causes cleanup of $arg_* variables)
}

EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Verify that required input files exist.
verifyInputFiles "$SINGLESEQ2" "$SINGLEPFS2" "$ENSEMBLECT"

# Test ProbKnot_without_options.
runFullTest 'without_options'  $SINGLEPFS2 @OUT.ct

# Test ProbKnot_dna_option.
#Test-TODO runFullTest 'dna_option'  $SINGLESEQ2 @OUT.ct --sequence -d

# Test ProbKnot_ensemble_option.
#Test-TODO runFullTest 'ensemble_option'  $ENSEMBLECT @OUT.ct --ensemble

# Test ProbKnot_iterations_option.
runFullTest 'iterations_option'  $SINGLEPFS2 @OUT.ct -i 2

# Test ProbKnot_min_helix_option.
runFullTest 'min_helix_option'  $SINGLEPFS2 @OUT.ct -m 2

# Test ProbKnot_sequence_option.
#Test-TODO runFullTest 'sequence_option' $SINGLESEQ2  $OUT.ct  --sequence
#    ---ref='without_options'

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
