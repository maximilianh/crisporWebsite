# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

source "$SCRIPT_DIR/test-tools-gui.sh"

 # Define the function that will be called implicitly by `runFullTest` and explicitly by `runTest $EXE`
function OligoWalk-gui() {
  local minArgs=2
  (($# >= minArgs)) || { warnmsg "$EXEBASE requires $minArgs arguments."; return  1; }
  local seqFile="$1" repFile="$2"  ACID=RNA

  BOOL_FLAGS=(d)   # Flags with NO parameters       (e.g. -d sets variable $arg_d to "1" in guiParseArgs)
  PARAM_FLAGS=(l c u m s fi st en)  # Flags with required parameters (e.g. -s "FILE" sets variable $arg_s to "FILE" in guiParseArgs)
  
  guiParseArgs ${@:minArgs+1} # Parse remaining arguments

  guiArg -d && ACID=DNA #Change acid to DNA if -d argument is present.

  guiWrite "Menu 'RNA', 'RNA OligoWalk'
    Click { BUTTON 'CT File' }
    TypeText '$seqFile'; TypeENTER
    Prohibit { DIALOG 'RNAstructure error' } 
    SetText { FIELD  'Report File' }, '$repFile'
    SetText { SPINNER 'Oligo Length:' }, '$arg_l'
    Click { RADIO '$ACID' }"
  
  # Choose the Target Structure Mode
  local mode='Break Local Structure'
  case $arg_m in
  	1) ;; # already assigned
  	2) mode='Refold *' ;;
	3) mode='Do Not Consider *' ;;
    *) warnmsg "Mode option $arg_m not implemented."  ;;
  esac
  guiWrite "Click { RADIO '$mode' }"

  # Choose the Concentration and Unit
  local unit='uM'
  local conc="${arg_c:-'1'}"
  case $arg_u in
  	# ToDo: Allow other units and modify conc to accomodate for e.g. -2, -5 etc.
     M|0)  unit='M' ;;
    mM|-3)  unit='mM' ;;
  	uM|-6)  unit='uM' ;;
    nM|-9)  unit='nM' ;;
    pM|-12) unit='pM' ;;
	'') ;; # use default
    *) warnmsg "Unit $arg_u not implemented." ;;
  esac
  guiWrite "SetText { COMBO }, '$unit'
            SetText { FIELD  'Oligo Concentration' }, '$conc'"

  # Choose Suboptimal Mode
  case $arg_s in
  	''|0) ;; # Do nothing. Leave box unchecked.
    3) guiWrite "Click { CHECK 'Include Target Suboptimal *' }" ;;
    *) warnmsg "Suboptimal option $arg_s not implemented in GUI."  ;;
  esac

  # Run the calculation
  guiWrite "CLICK { BUTTON 'START' }
    Prohibit { DIALOG 'RNAstructure error' } 
    WAIT  { BUTTON 'Go...' }, ${CALC_TIMEOUT}000"

  guiRun -c # Run the script
}

beginTestBlock  # Begin a group of tests.
runMake -r gui-tester  # run `make` for the Java GUI Tester and any other listed programs.
EXT=.rep  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Remove <html> tags from OligoWalk output.
#function TEST_POST_PROCESS() { perl -pi -e 's/<(br|hr)>/\n/ig; s/<\/?\w+>//g' "$1" >&2 ; }

# mode 1 - break local target structure to bind oligo 
# mode 2 - refold target RNA after oligo binding
# mode 3 - no target structure considered

# suboptimal 0 - only consider optimal structure
# suboptimal 1 - like choice 3,using suboptimal structures,but the whole set from alltrace() function prediction
# suboptimal 2 - using partition function considering every possible structure  of target
# 			     suboptimal 2 can only used with mode 2
# suboptimal 3 - using suboptimal structures (heuristic method) for both oligo-free and oligo-bound target RNA
# suboptimal 4 - using stochastic sampling method to sample 1000 structures

# useprefilter 1 - using criteria to prefill functional siRNA; (-- you may not want to type -test )

STD_ARGS='-l 18 -c 10 -u uM' # Standard concentration 10 uM

runFullTest 'mode_break'	$SINGLESEQ @OUT.rep $STD_ARGS -m 1
runFullTest 'mode_refold'	$SINGLESEQ @OUT.rep $STD_ARGS -m 2
runFullTest 'mode_ignore'	$SINGLESEQ @OUT.rep $STD_ARGS -m 3

runFullTest 'dna_option'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 2 -d

runFullTest 'suboptimal_0_none'	$SINGLESEQ @OUT.rep $STD_ARGS -m 2 -s 0
#NO_GUI: runFullTest 'suboptimal_1_all'	$SINGLESEQ @OUT.rep $STD_ARGS -m 2 -s 1
#NO_GUI: runFullTest 'suboptimal_2_part'	$SINGLESEQ @OUT.rep $STD_ARGS -m 2 -s 2
runFullTest 'suboptimal_3_heur'	$SINGLESEQ @OUT.rep $STD_ARGS -m 2 -s 3

# RMW_BROKEN:  suboptimal=4 currently segfaults
# RMW_BROKEN:  runFullTest 'suboptimal_4_stoc'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 2 -s 4

#NO_GUI: runFullTest 'prefilter'	$SINGLESEQ @OUT.rep $STD_ARGS -m 2 -fi 1

runFullTest 'length_10'	$SINGLESEQ @OUT.rep $STD_ARGS -m 2 -l 10

#NO_GUI: runFullTest 'alternate_seq' $SINGLESEQ2  @OUT.rep $STD_ARGS -m 2 -st 101 -en 200
runFullTest 'alternate_seq_break' $SINGLESEQ2 @OUT.rep $STD_ARGS -m 1

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
