# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

#Test-TODO: Perhaps re-write dynalign-gui() to read and parse the conf file to generate the gui-script instead of passing parameters.

source "$SCRIPT_DIR/test-tools-gui.sh"
beginTestBlock  # Begin a group of tests.
runMake -r gui-tester  # run `make` for the Java GUI Tester and any other listed programs.

CALC_TIMEOUT=45 # Wait at most 45 seconds for calculation to complete.

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function dynalign-gui() {
  local minArgs=1
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local ACID=RNA conf=$1
  
  guiUnsetArgs # clear previous argument values.

# Set defaults:
local maxtrace=750  percent=20  bpwin=2  awin=1
      #arg_fgap=.4  
      # arg_insert=1  arg_singlefold_subopt_percent=30 arg_imaxseparation=-99  
      # arg_num_processors=1  arg_optimal_only=0  arg_local=0  arg_dsv_templated=0  arg_dsvtemplatename=<template file name> arg_ct_templated=0 
      # arg_shapeslope1=1.8  arg_shapeintercept1=-0.6 arg_shapeslope2=1.8 arg_shapeintercept2=-0.6 arg_DNA=0  arg_temperature=310.15   
  guiParseConf $conf

  guiArg DNA && ACID=DNA
  guiWrite "Menu '$ACID', '$ACID Dynalign'
    Click { BUTTON 'Sequence File 1' }; TypeText '$inseq1'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD 'CT File 1' }, '$outct'
    Click { BUTTON 'Sequence File 2' }; TypeText '$inseq2'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD 'CT File 2' }, '$outct2'
    SetText { FIELD 'Alignment File' }, '$aout'" 

  [[ $savefile ]] && guiWrite "CLICK { CHECK 'Generate Save File' }"
  guiWrite "Focus { Field 'Max Number of Structures' }; ClearText; TypeText '$maxtrace'"
  guiWrite "Focus { Field 'Max % Energy Difference' };  ClearText; TypeText '$percent'"
  guiWrite "Focus { FIELD 'Structure Window Size' }; ClearText; TypeText '$bpwin'"
  guiWrite "Focus { FIELD 'Alignment Window Size' }; ClearText; TypeText '$awin'"

  guiWrite "CLICK { BUTTON 'START' }
    Prohibit { DIALOG 'RNAstructure error' } 
    WAIT  { BUTTON 'OK' }, ${CALC_TIMEOUT}000 
    CLICK { BUTTON 'Cancel' }"

  guiRun -c # Run the script (-c causes cleanup of $arg_* variables)
}


initTest 'general_test'
CONF="${TEST}_config.conf"

# Test dynalign with the general self-contained example given.
TEXT=\
'inseq1 = testFiles/testFile_RD0260.seq
inseq2 = testFiles/testFile_RD0500.seq
outct = @TEST_seq1.ct
outct2 = @TEST_seq2.ct
aout = @TEST_alignment.ali'

echo "${TEXT//'@TEST'/$TEST}" > "$CONF"
[[ $SMP ]] &&  echo "num_processors = 2" >> "$CONF"

runTest @EXE "$CONF"
runDiff @TEST_seq1.ct @OKBASE_seq1_OK.ct  'seq1'
runDiff @TEST_seq2.ct @OKBASE_seq2_OK.ct  'seq2'
runDiff @TEST_alignment.ali @OKBASE_alignment_OK.ali  'ali'
rm -f "$CONF"; endTest

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
