# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

INPUT_DOT=

#TODO: # Convert a single structure to CT format
#TODO: runFullTest 'numbered_structure' $INPUT @OUT.ct -n 1

# Convert all structures to CT format
runFullTest 'without_options' testFiles/testFile_RA7680.dot @OUT.ct

# For the rest, use the ct2dot "OK" files as input to dot2ct
CT2DOT_INPUT=ct2dot/ct2dot_@TESTNAME_OK.dot

runFullTest 'all_structures'   $CT2DOT_INPUT @OUT.ct

# Read a dot file with a single structure
runFullTest 'single_structure' $CT2DOT_INPUT @OUT.ct

# Read a "standard" dot file with a title line and structure-line for each sequence
runFullTest 'format_multi'     $CT2DOT_INPUT @OUT.ct ---ref='all_structures'

# Read a dot file with just structure-lines for subsequent structures (No title line or sequence)
runFullTest 'format_simple'    $CT2DOT_INPUT @OUT.ct

# Read a dot file with side-comments (Vienna style)
runFullTest 'format_side'      $CT2DOT_INPUT @OUT.ct ---ref='all_structures'

# Read a dot file with a pseudoknot
runFullTest 'pseudoknot'       $CT2DOT_INPUT @OUT.ct


endTestBlock # End a group of tests. Also cleans up orphaned files etc.

