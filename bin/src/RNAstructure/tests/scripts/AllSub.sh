# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.

runMake -r @EXE  # run `make` for the main exe and any other listed programs.

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