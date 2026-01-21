# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester ProbabilityPlot
FIX_DESC='--desc @TBASE_output.pfs'  # the output ps files contain the name of the input file, which changes (e.g. partition vs partition-smp) to correct for this, we specify a description to use.

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function bipartition-gui() {
  local minArgs=3
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local seqFile1="$1" seqFile2="$2" savFile="$3" ACID=RNA
  BOOL_FLAGS=(d)                             # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(c md sh si sm t)    # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)

  guiUnsetArgs # Clean up any temporary variables remaining from previous test
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiArg -d && ACID=DNA
  guiWrite "Menu '$ACID', 'Partition Function $ACID Bimolecular'
    Click { BUTTON 'Sequence File 1' }
    TypeText '$seqFile1'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' }
    Click { BUTTON 'Sequence File 2' }
    TypeText '$seqFile2'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD  'Save File' }, '$savFile'"

  # guiWriteIf writes lines to the gui script IF the argument specified by the first parameter is present.
  guiWriteIf -c "Menu 'Force', 'Restore Constraints'
    TypeText '$arg_c';  TypeENTER
    Prohibit { DIALOG 'RNAstructure error' }"
  guiWriteIf -md "Menu 'Force', 'Maximum Pairing Distance'
    Click { OPTION 'Yes' }
    Focus { FIELD 'Maximum Distance' }
    ClearText; TypeText '$arg_md'
    Click { BUTTON 'OK' }"

  #----SHAPE----
  guiWriteIf -sh "Menu 'Force', 'Read SHAPE Reactivity -- Pseudo-Energy Constraints'
    Focus { FIELD 'SHAPE Data File' }
    TypeText '$arg_sh';"
  guiWriteIf -si   "Focus { FIELD 'Intercept *'}; ClearText;  TypeText '$arg_si'"
  guiWriteIf -sm   "Focus { FIELD 'Slope *'};     ClearText;  TypeText '$arg_sm'"
  guiWriteIf -sh "Click { BUTTON 'OK' }" # Click OK on SHAPE dialog
  #----^ SHAPE ^ ----
  
  guiWriteIf -t "Menu 'Temperature', 'Set Temperature'
    Focus { FIELD 'input' }
    ClearText; 
    TypeText '$arg_t'; 
    Click { BUTTON 'OK' };"

  guiWrite "CLICK { BUTTON 'START' }
    Prohibit { DIALOG 'RNAstructure error' } 
    WAIT  { BUTTON 'OK' }, ${CALC_TIMEOUT}000 
    CLICK { BUTTON 'Cancel' }"

  guiRun -c # Run the script (-c causes cleanup of $arg_* variables)
}

EXT=.ps  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test bipartition_without_options.
initTest 'without_options'  && {
runTest @EXE $DOUBLESEQ  @OUT.pfs 
runSubTest 'plot' ProbabilityPlot @OUT.pfs @OUT.ps  $FIX_DESC
runDiff; endTest; }

# Test bipartition_dna_option.
initTest 'dna_option'    && {
runTest @EXE $DOUBLESEQ  @OUT.pfs  -d
runSubTest 'plot' ProbabilityPlot @OUT.pfs @OUT.ps  $FIX_DESC
runDiff; endTest; }

# Test bipartition_temperature_option.
initTest 'temperature_option'    && {
runTest @EXE $DOUBLESEQ  @OUT.pfs  -t 330
runSubTest 'plot' ProbabilityPlot @OUT.pfs @OUT.ps  $FIX_DESC
runDiff; endTest; }

endTestBlock # End a group of tests. Also cleans up orphaned files etc.