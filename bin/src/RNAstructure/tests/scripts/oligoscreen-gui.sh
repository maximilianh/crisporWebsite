# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester  # run `make` for the Java GUI Tester and any other listed programs.

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function oligoscreen-gui() {
  local minArgs=2
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local  ACID=RNA list=$1 rep=$2
  BOOL_FLAGS=(d)   # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(t)  # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)
  
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiArg -d && ACID=DNA
  guiWrite "Menu 'File', 'OligoScreen'
    Click { BUTTON 'Oligomer List' }
    TypeText '$list'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD 'Report File' }, '$rep'
    Click { RADIO '$ACID' }"
  
  guiWriteIf -t "Menu 'Temperature', 'Set Temperature'
    Focus { FIELD 'input' }
    ClearText; TypeText '$arg_t'
    Click { BUTTON 'OK' }"

  guiWrite "CLICK { BUTTON 'START' }
    Prohibit { DIALOG 'RNAstructure error' } 
    WAIT  { BUTTON 'OK' }, ${CALC_TIMEOUT}000 
    CLICK { BUTTON 'OK' }"

  guiRun -c # Run the script
}
EXT=.rep  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test oligoscreen_without_options.
runFullTest 'without_options'  $OLIGOLIST @OUT.rep

# Test oligoscreen_dna_option.
runFullTest 'dna_option'  $OLIGOLIST @OUT.rep -d

# Test oligoscreen_temperature_option.
runFullTest 'temperature_option'  $OLIGOLIST @OUT.rep -t 150

endTestBlock # End a group of tests. Also cleans up orphaned files etc.

