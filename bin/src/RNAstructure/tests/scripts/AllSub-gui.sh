# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester  # run `make` for the Java GUI Tester and any other listed programs.

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function AllSub-gui() {
  local minArgs=2
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local seqFile="$1" ctFile="$2" ACID=RNA
  BOOL_FLAGS=(d)           # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(a c p t)    # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)
  
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiArg -d && ACID=DNA
  guiWrite "Menu '$ACID', 'Generate All Suboptimal $ACID Structures'
    Click { BUTTON 'Sequence File' }
    TypeText '$seqFile'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD  'CT File' }, '$ctFile'" 

  # guiWriteIf writes lines to the gui script IF the argument specified by the first parameter is present.
  guiWriteIf -c "Menu 'Force', 'Restore Constraints';  TypeText '$arg_c';  TypeENTER"
  guiWriteIf -p "Focus { Field 'Max % Energy Difference' }; ClearText; TypeText '$arg_p'"
  guiWriteIf -a "Focus { Field 'Max Absolute Energy Difference' }; ClearText; TypeText '$arg_a'"
  guiWriteIf -t "Menu 'Temperature', 'Set Temperature'
    Focus { FIELD 'input' }; 
    ClearText; TypeText '$arg_t'
    Click { BUTTON 'OK' }"

  guiWrite "CLICK { BUTTON 'START' }
    Prohibit { DIALOG 'RNAstructure error' } 
    WAIT  { BUTTON 'OK' }, ${CALC_TIMEOUT}000 
    CLICK { BUTTON 'Cancel' }"

  guiRun -c # Run the script (-c causes cleanup of $arg_* variables)
}

EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test AllSub_without_options.
runFullTest 'without_options'  $SINGLESEQ_SHORT @OUT.ct

# Test AllSub_absolute_energy_difference_option.
runFullTest 'absolute_energy_difference_option'  $SINGLESEQ_SHORT @OUT.ct -a 1

# Test AllSub_constraint_file_option.
runFullTest 'constraint_file_option'  $SINGLESEQ_SHORT @OUT.ct -c testFiles/testFile_folding5.con

# Test AllSub_dna_option.
runFullTest 'dna_option'  $SINGLESEQ_SHORT @OUT.ct -d -a 0.1

# Test AllSub_percent_difference_option.
runFullTest 'percent_difference_option'  $SINGLESEQ_SHORT @OUT.ct -p 1

# Test AllSub_temperature_option.
runFullTest 'temperature_option'  $SINGLESEQ_SHORT @OUT.ct -t 250

endTestBlock # End a group of tests. Also cleans up orphaned files etc.