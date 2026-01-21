# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

# Note: The C-libraries on Linux, Mac, and Windows do not format exponents the same way.
# For example  `std::cout << 1e-7;`   // Outputs 1e-007 on Windows/Mac and 1e-07 on Linux
# To correct for this, the function `fixExponents` (defined in test-tools.sh) converts
# 3-digit exponents with leading zeros into 2-digit exponents. e.g.:  3.4e-006 ==> 3.4e-06

beginTestBlock  # Begin a group of tests.
runMake -r $EXE ProbabilityPlot # run `make` for the main exe and any other listed programs
EXT=.ps  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)
FIX_DESC='--desc @TBASE_output.pfs'  # the output ps files contain the name of the input file, which changes (e.g. partition vs partition-smp) to correct for this, we specify a description to use.

# Test partition_without_options.
initTest 'without_options'  
runTest $EXE $SINGLESEQ  $OUT.pfs 
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_without_options.
initTest 'without_options_stdin'  
runTest $EXE - $OUT.pfs < $SINGLESEQ
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_without_options_alternate.
# This alternate test will be used to test experimental pair bonuses later.
initTest 'without_options_alternate'  
runTest $EXE testFiles/testFile_5SRNA_tail2.seq  $OUT.pfs
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_without_options_alternate_2.
# This alternate test will be used to test double stranded offsets later.
initTest 'without_options_alternate_2'  
runTest $EXE testFiles/testFile_U1a.seq  $OUT.pfs 
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps   $FIX_DESC 
fixExponents $OUT.ps; runDiff; endTest

# Test partition_constraint_file_option.
initTest 'constraint_file_option'  
runTest $EXE $SINGLESEQ  $OUT.pfs  -c testFiles/testFile_folding3.con
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_dna_option.
initTest 'dna_option'  
runTest $EXE $SINGLESEQ  $OUT.pfs  -d
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_double_stranded_offset_option.
initTest 'double_stranded_offset_option'  
runTest $EXE testFiles/testFile_U1a.seq  $OUT.pfs  -dso testFiles/testFile_double_offset_dummy.txt
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_experimental_pair_bonus_option.
initTest 'experimental_pair_bonus_option'  
runTest $EXE testFiles/testFile_5SRNA_tail2.seq  $OUT.pfs  -X testFiles/testFile_bonus_matrix.txt
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_experimental_pair_bonus_offset_option.
initTest 'experimental_pair_bonus_offset_option'  
runTest $EXE testFiles/testFile_5SRNA_tail2.seq  $OUT.pfs  -X testFiles/testFile_bonus_matrix.txt -xo 10
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_experimental_pair_bonus_scaling_option.
initTest 'experimental_pair_bonus_scaling_option'  
runTest $EXE testFiles/testFile_5SRNA_tail2.seq  $OUT.pfs  -X testFiles/testFile_bonus_matrix.txt -xs 0.99
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_max_distance_option.
initTest 'max_distance_option'  
runTest $EXE $SINGLESEQ  $OUT.pfs  -md 15
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_shape_option.
initTest 'shape_option'  
runTest $EXE $SINGLESEQ  $OUT.pfs  -sh testFiles/testFile_tRNA.shape
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_shape_intercept_option.
initTest 'shape_intercept_option'  
runTest $EXE $SINGLESEQ  $OUT.pfs  -sh testFiles/testFile_tRNA.shape -si 0.9
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_shape_slope_option.
initTest 'shape_slope_option'  
runTest $EXE $SINGLESEQ  $OUT.pfs  -sh testFiles/testFile_tRNA.shape -sm 0.2
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

# Test partition_temperature_option.
initTest 'temperature_option'  
runTest $EXE $SINGLESEQ  $OUT.pfs  -t 150
runSubTest 'plot' ProbabilityPlot   $OUT.pfs $OUT.ps  $FIX_DESC
fixExponents $OUT.ps; runDiff; endTest

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
