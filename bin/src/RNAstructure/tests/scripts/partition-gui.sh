# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester ProbabilityPlot
FIX_DESC='--desc @TBASE_output.pfs'  # the output ps files contain the name of the input file, which changes (e.g. partition vs partition-smp) to correct for this, we specify a description to use.

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function partition-gui() {
  local minArgs=2
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local seqFile="$1" savFile="$2" ACID=RNA
  BOOL_FLAGS=(d)                   # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(c md sh si sm t)    # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)

  guiUnsetArgs # Clean up any temporary variables remaining from previous test
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiArg -d && ACID=DNA
  guiWrite "Menu '$ACID', 'Partition Function $ACID'
    Click { BUTTON 'Sequence File' }
    TypeText '$seqFile'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD  'Save File' }, '$savFile'"

  # guiWriteIf writes lines to the gui script IF the argument specified by the first parameter is present.
  guiWriteIf -c "Menu 'Force', 'Restore Constraints';  TypeText '$arg_c';  TypeENTER"
  guiWriteIf -md "Menu 'Force', 'Maximum Pairing Distance'
    Click { OPTION 'Yes' }
    Focus { FIELD 'Maximum Distance' }
    ClearText; TypeText '$arg_md'
    Click { BUTTON 'OK' }"

  #----SHAPE----
  guiWriteIf -sh "Menu 'Force', 'Read SHAPE Reactivity -- Pseudo-Energy Constraints'
    TypeText { FIELD 'SHAPE Data File' }, '$arg_sh'"
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

# Test partition_without_options.
initTest 'without_options'  
runTest  @EXE $SINGLESEQ  $OUT.pfs 
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_without_options_alternate.
# This alternate test will be used to test experimental pair bonuses later.
initTest 'without_options_alternate'  
runTest  @EXE testFiles/testFile_5SRNA_tail2.seq  $OUT.pfs 
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_without_options_alternate_2.
# This alternate test will be used to test double stranded offsets later.
initTest 'without_options_alternate_2'  
runTest  @EXE testFiles/testFile_U1a.seq  $OUT.pfs 
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
fixExponents $OUT.ps; runDiff; endTest

# Test partition_constraint_file_option.
initTest 'constraint_file_option'  
runTest  @EXE $SINGLESEQ  $OUT.pfs -c testFiles/testFile_folding3.con
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
fixExponents $OUT.ps; runDiff; endTest

# Test partition_dna_option.
initTest 'dna_option'  
runTest  @EXE $SINGLESEQ  $OUT.pfs  -d
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
fixExponents $OUT.ps; runDiff; endTest


# Test partition_double_stranded_offset_option.
# Test_TODO: initTest 'double_stranded_offset_option'  
# runTest  @EXE testFiles/testFile_U1a.seq  $OUT.pfs  -dso testFiles/testFile_double_offset_dummy.txt
# runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
# runDiff $OUT.ps $OK.ps 
# endTest

# Test partition_experimental_pair_bonus_option.
# initTest 'experimental_pair_bonus_option'  
# Test_TODO: runTest  @EXE testFiles/testFile_5SRNA_tail2.seq  $OUT.pfs  -X testFiles/testFile_bonus_matrix.txt
# runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
# runDiff $OUT.ps $OK.ps 
# endTest

# Test partition_experimental_pair_bonus_offset_option.
# initTest 'experimental_pair_bonus_offset_option'  
# Test_TODO: runTest  @EXE testFiles/testFile_5SRNA_tail2.seq  $OUT.pfs  -X testFiles/testFile_bonus_matrix.txt -xo 10
# runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
# runDiff $OUT.ps $OK.ps 
# endTest

# Test partition_experimental_pair_bonus_scaling_option.
# initTest 'experimental_pair_bonus_scaling_option'  
# # Test_TODO: runTest  @EXE testFiles/testFile_5SRNA_tail2.seq  $OUT.pfs  -X testFiles/testFile_bonus_matrix.txt -xs 0.99
# runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
# runDiff $OUT.ps $OK.ps 
# endTest

# Test partition_max_distance_option.
initTest 'max_distance_option'  
runTest  @EXE $SINGLESEQ  $OUT.pfs  -md 15
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
fixExponents $OUT.ps; runDiff; endTest

# Test partition_shape_option.
initTest 'shape_option'  
runTest  @EXE $SINGLESEQ  $OUT.pfs  -sh testFiles/testFile_tRNA.shape
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
fixExponents $OUT.ps; runDiff; endTest

# Test partition_shape_intercept_option.
initTest 'shape_intercept_option'  
runTest  @EXE $SINGLESEQ  $OUT.pfs  -sh testFiles/testFile_tRNA.shape -si 0.9
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
fixExponents $OUT.ps; runDiff; endTest

# Test partition_shape_slope_option.
initTest 'shape_slope_option'  
runTest  @EXE $SINGLESEQ  $OUT.pfs  -sh testFiles/testFile_tRNA.shape -sm 0.2
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
fixExponents $OUT.ps; runDiff; endTest

# Test partition_temperature_option.
initTest 'temperature_option'  
runTest  @EXE $SINGLESEQ  $OUT.pfs  -t 150
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC 
fixExponents $OUT.ps; runDiff; endTest

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
