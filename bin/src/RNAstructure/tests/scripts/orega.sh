# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE # run `make` for the main exe and any other listed programs

# The input file
SEQ=testFiles/testFile_RA7680.seq

# Verify that required input files exist.
verifyInputFiles "$SEQ"

RSEED='--seed 123456789' # set arbitrary random seed
EXT=.out  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)
OUTSEQ=@OUTFILE.fa

# Test with 10 iterations, default population size
runFullTest 'iter_10'           $SEQ 31 20 $OUTSEQ --iter 10    $RSEED ---stdout

# Test with 20 iterations, population size = 4
runFullTest 'iter_20'           $SEQ 31 20 $OUTSEQ --iter 20 --population 4    $RSEED ---stdout

# Test with Higher mutation rate
runFullTest 'mutate_0.5'        $SEQ 31 20 $OUTSEQ --iter 20 --population 4  --mutate 0.5  $RSEED ---stdout

# Test with filteroligoA, polyA = 7
runFullTest 'filteroligoA'      $SEQ 31 20 $OUTSEQ --iter 20 --population 4  --filteroligoA 7  $RSEED ---stdout

# Test with ComplexityConstant
runFullTest 'complexity_constant'      $SEQ 31 20 $OUTSEQ --iter 20 --population 4  --ComplexityConstant 0.1 $RSEED ---stdout


# Test with filterAUG
runFullTest 'filterAUG'      $SEQ 31 20 $OUTSEQ --iter 20 --population 4  --filterAUG  $RSEED ---stdout

# Test with filterCUG
runFullTest 'filterCUG'      $SEQ 31 20 $OUTSEQ --iter 20 --population 4  --filterCUG  $RSEED ---stdout

# Test with limitG
runFullTest 'limitG'      $SEQ 31 20 $OUTSEQ --iter 20 --population 4  --limitG  $RSEED ---stdout

# Test with MutationSwitch
runFullTest 'MutationSwitch'      $SEQ 31 20 $OUTSEQ --iter 20 --population 4  --MutationSwitch  $RSEED ---stdout

# Test with Higher Recombination Rate
runFullTest 'recomb_rate_0.5'   $SEQ 31 20 $OUTSEQ --iter 20 --population 4  -rr 0.5  $RSEED ---stdout

# Test with Higher Recombination Frequency ( smaller rf )
runFullTest 'recomb_freq_3'     $SEQ 31 20 $OUTSEQ --iter 20 --population 4  -rf 3  $RSEED ---stdout

# Test with Different objective function
runFullTest 'func_simple'       $SEQ 31 20 $OUTSEQ --iter 20 --population 4  --func 1  $RSEED ---stdout

# Test with Different ranges (SEQ size is 76)
runFullTest 'range_start'       $SEQ 1  20 $OUTSEQ --iter 3 --population 4  $RSEED ---stdout
runFullTest 'range_end'         $SEQ 66 10 $OUTSEQ --iter 3 --population 4  $RSEED ---stdout
runFullTest 'range_full'        $SEQ 1  76 $OUTSEQ --iter 3 --population 4  $RSEED ---stdout

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
