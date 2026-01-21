# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

errmsg 'DuplexFold is not implemented in the Java GUI.'
exit 1

# NOTE: There is no distinction in the Java GUI between DuplexFold and bifold
#       Except that all DuplexFold tests include the "Force > Forbid Unimolecular Pairs" 
#       constraint, whereas the bifold tests only include it if the -i flag is present.

source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester  # run `make` for the Java GUI Tester and any other listed programs.

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function DuplexFold-gui() {
  local minArgs=3
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local seqFile1="$1" seqFile2="$2" ctFile="$3" ACID=RNA
  BOOL_FLAGS=(d)             # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(l m p t w)    # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)
  
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiArg -d && ACID=DNA
  guiWrite "Menu '$ACID', 'Fold $ACID Bimolecular'
    Click { BUTTON 'Sequence File 1' }
    TypeText '$seqFile1'; TypeENTER
    Click { BUTTON 'Sequence File 2' }
    TypeText '$seqFile2'; TypeENTER    
    SetText { FIELD  'CT File' }, '$ctFile'" 

  # in bifold-gui, this menu constraint requires the -i flag, but in DuplexFold, it is implicit
  guiWrite "Menu 'Force', 'Forbid Unimolecular Pairs'" # intramolecular_option

  # guiWriteIf writes lines to the gui script IF the argument specified by the first parameter is present.
  guiWriteIf -l "Menu 'Maximum Loop', 'Set Maximum Loop Size'
    Focus { FIELD 'input' }
    ClearText; TypeText '$arg_l';
    Click { BUTTON 'OK' }"
  guiWriteIf -m "Focus { Field 'Max Number of Structures' }; ClearText; TypeText '$arg_m'"
  guiWriteIf -p "Focus { Field 'Max % Energy Difference' };  ClearText; TypeText '$arg_p'"
  guiWriteIf -t "Menu 'Temperature', 'Set Temperature'
    Focus { FIELD 'input' }
    ClearText; TypeText '$arg_t'
    Click { BUTTON 'OK' }"
  guiWriteIf -w "Focus { FIELD 'Window Size' }; ClearText; TypeText '$arg_w'"

  guiWrite "CLICK { BUTTON 'START' }
    WAIT  { BUTTON 'OK' }, ${CALC_TIMEOUT}000 
    CLICK { BUTTON 'Cancel' }"

  guiRun -c  # Run the script (-c causes cleanup of $arg_* variables)
}

EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test DuplexFold_without_options.
runFullTest 'without_options'  $DOUBLESEQ @OUT.ct

# Test DuplexFold_dna_option.
runFullTest 'dna_option'  $DOUBLESEQ @OUT.ct -d

# Test DuplexFold_loop_option.
runFullTest 'loop_option'  $DOUBLESEQ @OUT.ct -l 10

# Test DuplexFold_max_structures_option.
runFullTest 'max_structures_option'  $DOUBLESEQ @OUT.ct -m 5

# Test DuplexFold_percent_difference_option.
runFullTest 'percent_difference_option'  $DOUBLESEQ @OUT.ct -p 1

# Test DuplexFold_temperature_option.
runFullTest 'temperature_option'  $DOUBLESEQ @OUT.ct -t 150

# Test DuplexFold_window_size_option.
runFullTest 'window_size_option'  $DOUBLESEQ @OUT.ct -w 5

endTestBlock # End a group of tests. Also cleans up orphaned files etc.

