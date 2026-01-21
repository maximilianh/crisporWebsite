# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

#EXCLUDE_TESTS+=' double_stranded_offset_option experimental_pair_bonus_option experimental_pair_bonus_offset_option experimental_pair_bonus_scaling_option unpaired_shape_intercept_option unpaired_shape_slope_option'

source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester EnergyPlot # run `make` for the Java GUI Tester and any other listed programs.

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function Fold-gui() {
  local minArgs=2
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local seqFile="$1" ctFile="$2" ACID=RNA
  BOOL_FLAGS=(d)                             # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(s t l c m w p md sh si sm)    # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)
  
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiArg -d && ACID=DNA
  guiWrite "Menu '$ACID', 'Fold $ACID Single Strand'
    Click { BUTTON 'Sequence File' }
    TypeText '$seqFile'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD  'CT File' }, '$ctFile'"

  # guiWriteIf writes lines to the gui script IF the argument specified by the first parameter is present.
  guiWriteIf -s "CLICK { CHECK 'Generate Save File' }"
  guiWriteIf -c "Menu 'Force', 'Restore Constraints';  TypeText '$arg_c';  TypeENTER"
  guiWriteIf -l "Menu 'Maximum Loop', 'Set Maximum Loop Size'
    Focus { FIELD 'input' }
    ClearText; TypeText '$arg_l';
    Click { BUTTON 'OK' }"
  guiWriteIf -md "Menu 'Force', 'Maximum Pairing Distance'
    Click { OPTION 'Yes' }
    Focus { FIELD 'Maximum Distance' }
    ClearText; TypeText '$arg_md'
    Click { BUTTON 'OK' }"
  guiWriteIf -m "Focus { Field 'Max Number of Structures' }; ClearText; TypeText '$arg_m'"
  guiWriteIf -p "Focus { Field 'Max % Energy Difference' };  ClearText; TypeText '$arg_p'"

  #----SHAPE----
  guiWriteIf -sh "Menu 'Force', 'Read SHAPE Reactivity -- Pseudo-Energy Constraints'
    TypeText { FIELD 'SHAPE Data File' }, '$arg_sh'"
  guiWriteIf -si "Focus { FIELD 'Intercept *'}; ClearText; TypeText '$arg_si'"
  guiWriteIf -sm "Focus { FIELD 'Slope *'    }; ClearText; TypeText '$arg_sm'"
  guiWriteIf -sh "Click { BUTTON 'OK' }" # Click OK on SHAPE dialog
  #----^ SHAPE ^ ----
  
  guiWriteIf -t "Menu 'Temperature', 'Set Temperature'
    Focus { FIELD 'input' }
    ClearText; TypeText '$arg_t'
    Click { BUTTON 'OK' }"
  guiWriteIf -w "Focus { FIELD 'Window Size' }; ClearText; TypeText '$arg_w'"

  guiWrite "CLICK { BUTTON 'START' }
    Prohibit { DIALOG 'RNAstructure error' } 
    WAIT  { BUTTON 'OK' }, ${CALC_TIMEOUT}000 
    CLICK { BUTTON 'Cancel' }"

  guiRun  # Run the script
  local code=$? # Store return value

  # After the script has run, move any sav files that were specified by -s (because the java GUI does not allow setting the save-file name.)
  if guiArg -s; then
    [[ -f $OUT.sav ]] && mv "$OUT.sav" "$arg_s" || 
    { errmsg "Failed to rename save-file $OUT.sav to $arg_s"; code=1; }
  fi
  
  guiUnsetArgs # Clean up temporary variables

  return $code
}

EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test Fold_without_options.
runFullTest 'without_options'  $SINGLESEQ2 @OUT.ct

# Test Fold_without_options_alternate.
# This alternate test will be used to test single and double stranded offsets later.
runFullTest 'without_options_alternate'  testFiles/testFile_U1a.seq @OUT.ct

# Test Fold_without_options_alternate_2.
# This second alternate test will be used to test experimental pair bonuses later.
runFullTest 'without_options_alternate_2'  testFiles/testFile_5SRNA_tail2.seq @OUT.ct

# Test Fold_without_options_fasta_input.
runFullTest 'without_options_fasta_input'  $SINGLESEQ2_FASTA  @OUT.ct \
    ---ref='without_options'

# Test Fold_constraint_file_option.
runFullTest 'constraint_file_option'  $SINGLESEQ2 @OUT.ct -c testFiles/testFile_folding2.con

# Test Fold_dna_option.
runFullTest 'dna_option'  $SINGLESEQ2 @OUT.ct -d

# Test Fold_dna_option_fasta_input.
runFullTest 'dna_option_fasta_input'  $SINGLESEQ2_FASTA  @OUT.ct  -d \
    ---ref='dna_option'

# Test Fold_double_stranded_offset_option.
#Test-TODO: runFullTest 'double_stranded_offset_option'  testFiles/testFile_U1a.seq @OUT.ct -dso testFiles/testFile_double_offset_dummy.txt

# Test Fold_experimental_pair_bonus_option.
#Test-TODO: runFullTest 'experimental_pair_bonus_option'  testFiles/testFile_5SRNA_tail2.seq @OUT.ct -X testFiles/testFile_bonus_matrix.txt

# Test Fold_experimental_pair_bonus_offset_option.
#Test-TODO: runFullTest 'experimental_pair_bonus_offset_option'  testFiles/testFile_5SRNA_tail2.seq @OUT.ct -X testFiles/testFile_bonus_matrix.txt -xo 10

# Test Fold_experimental_pair_bonus_scaling_option.
#Test-TODO: runFullTest 'experimental_pair_bonus_scaling_option'  testFiles/testFile_5SRNA_tail2.seq @OUT.ct -X testFiles/testFile_bonus_matrix.txt -xs 10

# Test Fold_loop_option.
runFullTest 'loop_option'  $SINGLESEQ2 @OUT.ct -l 0

# Test Fold_max_distance_option.
runFullTest 'max_distance_option'  $SINGLESEQ2 @OUT.ct -md 25

# Test Fold_max_structures_option.
runFullTest 'max_structures_option'  $SINGLESEQ2 @OUT.ct -m 2

# Test Fold_percent_difference_option.
runFullTest 'percent_difference_option'  $SINGLESEQ2 @OUT.ct -p 5

#Test-TODO: runFullTest 'minimum_free_energy_option'  $SINGLESEQ @OUT.ct -mfe

SAVEFILE=fold_save_file.sav
# Test Fold_save_file_option.
runFullTest 'save_file_option'  $SINGLESEQ2  $OUT.ct  -s $SAVEFILE \
    ---ref='without_options'

# Test Fold_sav_plot (The save plot generated by the save file option).
if [[ -z $SKIP_TESTS ]]; then
    initTest 'sav_plot' && { # skip this block if the test is excluded.
    runTest EnergyPlot $SAVEFILE @OUT.ps
    addTestOuput $SAVEFILE # mark this file as test output (for proper cleanup)
    runDiff @OUT.ps @OK.ps
    endTest; }
fi

# Test Fold_shape_option.
runFullTest 'shape_option'  $SINGLESEQ2 @OUT.ct -sh testFiles/testFile_tRNA.shape

# Test Fold_shape_intercept_option.
runFullTest 'shape_intercept_option'  $SINGLESEQ2 @OUT.ct -sh testFiles/testFile_tRNA.shape -si 0.9

# Test Fold_shape_slope_option.
runFullTest 'shape_slope_option'  $SINGLESEQ2 @OUT.ct -sh testFiles/testFile_tRNA.shape -sm -0.2

# Test Fold_single_stranded_offset_option.
#Test-TODO: runFullTest 'single_stranded_offset_option'  testFiles/testFile_U1a.seq @OUT.ct -sso testFiles/testFile_single_offset.txt

# Test Fold_temperature_option.
runFullTest 'temperature_option'  $SINGLESEQ2 @OUT.ct -t 150

# Test Fold_unpaired_shape_intercept_option.
#Test-TODO: runFullTest 'unpaired_shape_intercept_option'  $SINGLESEQ2 @OUT.ct -sh testFiles/testFile_tRNA.shape -usi 0.9

# Test Fold_unpaired_shape_slope_option.
#Test-TODO: runFullTest 'unpaired_shape_slope_option'  $SINGLESEQ2 @OUT.ct -sh testFiles/testFile_tRNA.shape -usm -0.2

# Test Fold_window_size_option.
runFullTest 'window_size_option'  $SINGLESEQ2 @OUT.ct -w 15

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
