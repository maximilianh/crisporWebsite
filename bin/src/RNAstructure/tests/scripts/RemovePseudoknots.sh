# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test RemovePseudoknots_without_options.
runFullTest 'without_options'  $KNOTSCT @OUT.ct

# Test RemovePseudoknots_dna_option.
runFullTest 'dna_option'  $KNOTSCT @OUT.ct -d

# Test RemovePseudoknots_maximize_option.
runFullTest 'maximize_option'  $KNOTSCT @OUT.ct -m

# Test RemovePseudoknots_temperature_option.
runFullTest 'temperature_option'  $KNOTSCT @OUT.ct -t 150

knot_files=(knotted_4+2  knotted_4+3+1  knotted_4+4 )
for name in "${knot_files[@]}"; do
  file=RemovePseudoknots/input/"$name".dot
  runFullTest "${name}"     "$file" @OUT.ct -m
  runFullTest "${name}_MIN" "$file" @OUT.ct            ---ref="$name"
  runFullTest "${name}_MEA" "$file" @OUT.ct -m --MEA   ---ref="$name"
done


name="knotted_4n+2+2"
file=RemovePseudoknots/input/"$name".dot
runFullTest "${name}"     "$file" @OUT.ct -m
runFullTest "${name}_MIN" "$file" @OUT.ct            ---ref="$name"
# the following does not match the output of -m (without --MEA, but that is expected)
runFullTest "${name}_MEA" "$file" @OUT.ct -m --MEA   


# Test output as a dot-bracket file.
EXT=.dot  runFullTest 'save_dotbracket'  $KNOTSCT @OUT.dot -m --bracket

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
