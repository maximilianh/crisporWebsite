# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs
EXT=.txt  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Verify that required input files exist.
verifyInputFiles "$SINGLESEQ" "$SINGLEPFS"

# Test ProbScan_helix
runFullTest 'helices' ---stdout      $SINGLEPFS -e 1

# Test ProbScan_bulge
runFullTest 'hairpins' ---stdout     $SINGLEPFS -a

# Test ProbScan_bulge
runFullTest 'bulges' ---stdout       $SINGLEPFS -b

# Test ProbScan_iloop
runFullTest 'iloops' ---stdout       $SINGLEPFS -i

# Test ProbScan_multibranch
runFullTest 'multibranch' ---stdout  $SINGLEPFS -m ProbScan/ProbScan_multibranch_test_loops.txt

# Test ProbScan_helix
runFullTest 'pairs' ---stdout        $SINGLEPFS -p 1-72,2-71,3-70

# Test ProbScan_stemloop
runFullTest 'stemloop' ---stdout     $SINGLEPFS --stemloop 53-61,52-62

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
