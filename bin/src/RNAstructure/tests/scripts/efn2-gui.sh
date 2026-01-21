# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.
 
source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester # run `make` for the Java GUI Tester and any other listed programs.

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function efn2-gui() {
  local minArgs=2
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local  ACID=RNA
  BOOL_FLAGS=(d p s w)                             # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(sh si sm t)    # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)
  
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiArg -d && ACID=DNA
  guiWrite "Menu '$ACID', 'Efn2 $ACID'
    Click { BUTTON 'CT File' }
    TypeText '$1'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD 'Output File' }, '$2'"

  # guiWriteIf writes lines to the gui script IF the argument specified by the first parameter is present.

  guiWriteIf -w "Click { CHECK 'Write Thermodynamic Details File' }" 

  #----SHAPE----
  guiWriteIf -sh "Menu 'Force', 'Read SHAPE Reactivity -- Pseudo-Energy Constraints'
    TypeText { FIELD 'SHAPE Data File' }, '$arg_sh'"
  guiWriteIf -si "Focus { FIELD 'Intercept*'}; ClearText; TypeText '$arg_si'"
  guiWriteIf -sm "Focus { FIELD 'Slope*'    }; ClearText; TypeText '$arg_sm'"
  guiWriteIf -sh "Click { BUTTON 'OK' }" # Click OK on SHAPE dialog
  #----^ SHAPE ^ ----
  
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

# Verify that required input files exist.
verifyInputFiles "$SINGLECT"

EXT=.out  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test efn2_without_options.
runFullTest 'without_options'  $SINGLECT @OUT.out

# Test efn2_dna_option.
runFullTest 'dna_option'  $SINGLECT @OUT.out -d

# # Test efn2_print_option.
# initTest 'print_option'
# runTest $EXENAME $SINGLECT @TBASE_output.out -p
# runDiff @TBASE_output.out @OKFILE  # run diff on @OUT.out
# runDiff @STDO @OKBASE_screen_OK.txt 'screen'  # run diff on the stdout of the previous command 
# endTest

# Test efn2_shape_option.
#Test_TODO: runFullTest 'shape_option'  $SINGLECT @OUT.out -sh testFiles/testFile_tRNA.shape

# Test efn2_shape_intercept_option.
#Test_TODO: runFullTest 'shape_intercept_option'  $SINGLECT @OUT.out -sh testFiles/testFile_tRNA.shape -si 0.2

# Test efn2_shape_slope_option.
#Test_TODO: runFullTest 'shape_slope_option'  $SINGLECT @OUT.out -sh testFiles/testFile_tRNA.shape -sm 1.2
# efn2/efn2_shape_slope_option_OK.out

# Test efn2_temperature_option.
runFullTest 'temperature_option'  $SINGLECT @OUT.out -t 150

# Test efn2_write_thermodynamic_file_option.
runFullTest 'write_thermodynamic_file_option'  $SINGLECT @OUT.out -w

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
