# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs
EXT=.rep  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Remove <html> tags from OligoWalk output.
# function TEST_POST_PROCESS() { 
# 	# Replace all <br> with newline and remove all other html tags.
# 	perl -pi -e 's/<(br|hr)>/\n/ig; s/<\/?\w+>//g' "$1" >&2 ; 
# }

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

# Standard options (unless overridden in individual tests)
STD_ARGS='-l 18 -c 10 -u uM' # Standard concentration 10 uM

runFullTest 'mode_break'	$SINGLECT  @OUT.rep $STD_ARGS -m 1 --structure
runFullTest 'mode_refold'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 2
runFullTest 'mode_ignore'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 3

runFullTest 'dna_option'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 2 -d

runFullTest 'suboptimal_0_none'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 2 -s 0
runFullTest 'suboptimal_1_all'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 2 -s 1
runFullTest 'suboptimal_2_part'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 2 -s 2
runFullTest 'suboptimal_3_heur'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 2 -s 3

# RMW_BROKEN:  suboptimal=4 currently segfaults
# RMW_BROKEN:  runFullTest 'suboptimal_4_stoc'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 2 -s 4

runFullTest 'prefilter'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 2 -fi 1

runFullTest 'length_10'	$SINGLESEQ  @OUT.rep $STD_ARGS -m 2 -l 10

runFullTest 'alternate_seq'	$SINGLESEQ2  @OUT.rep $STD_ARGS -m 2 -st 101 -en 200
#runFullTest 'alternate_seq_break'	$ENSEMBLECT @OUT.rep $STD_ARGS -m 1 -s 3 --structure

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
