# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.dot  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)
DIFF_FLAGS+=' -b' # ignore differences in spaces

INPUT_CT=testFiles/testFile_RA7680.ct

# Convert a single structure to bracket notation
runFullTest 'single_structure'  $INPUT_CT 1 @OUT.dot

# Convert all structures to bracket notation
runFullTest 'all_structures'    $INPUT_CT ALL @OUT.dot

# Single-header format
runFullTest 'format_simple'     $INPUT_CT ALL @OUT.dot --format simple

# Side-title  format
runFullTest 'format_side'       $INPUT_CT ALL @OUT.dot --format side

# multiple-header  format
runFullTest 'format_multi'      $INPUT_CT ALL @OUT.dot --format multi ---ref='all_structures'

# multiple-header-and-seq  format
runFullTest 'format_full'       $INPUT_CT ALL @OUT.dot --format full

# Test reading and writing a structure with pseudoknots
runFullTest 'pseudoknot' testFiles/testFile_knotted.ct ALL @OUT.dot

endTestBlock