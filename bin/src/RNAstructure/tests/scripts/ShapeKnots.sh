# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Test ShapeKnots_shape_without_options.
runFullTest 'without_options'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct

# Test ShapeKnots_shape_without_SHAPE_file.
runFullTest 'without_SHAPE_file'  testFiles/testFile_RD0260.seq @OUT.ct

# Test ShapeKnots_intercept_option.
runFullTest 'intercept_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -si -1.5

# Test ShapeKnots_slope_option.
runFullTest 'slope_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -sm 2.8

# Test ShapeKnots_P1_option.
runFullTest 'P1_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -p1 -3

# Test ShapeKnots_P2_option.
runFullTest 'P2_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -p2 -3

# Test ShapeKnots_single_stranded_offset_option.
runFullTest 'single_stranded_offset_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -sso testFiles/testFile_single_offset.txt

# Test ShapeKnots_double_stranded_offset_option.
runFullTest 'double_stranded_offset_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -dso testFiles/testFile_double_offset_dummy.txt

# Test the warning that should be shown when an invalid offset position is encountered.
# this uses the ShapeKnots --warn ERR flag to send warnings to STDERR. 
# it also uses the test-system's --stderr flag to indicate that the STDERR should be diff'd against the reference OK file.
EXT=.out \
runFullTest 'warn_invalid_offset_position'  testFiles/testFile_RD0260.seq @OUT.ct -sso testFiles/testFile_single_offset.txt --warn ERR  ---stderr

# Test ShapeKnots_max_structures_option.
runFullTest 'max_structures_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -m 2

# Test ShapeKnots_percent_difference_option.
runFullTest 'percent_difference_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -p 0

# Test ShapeKnots_window_size_option.
runFullTest 'window_size_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -w 0

# Test ShapeKnots_internal_max_structures_option.
runFullTest 'internal_max_structures_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -im 10

# Test ShapeKnots_internal_percent_difference_option.
runFullTest 'internal_percent_difference_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -ip 0

# Test ShapeKnots_internal_window_size_option.
runFullTest 'internal_window_size_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -iw 10

# Test ShapeKnots_pseudoknotted_helices_option.
runFullTest 'pseudoknotted_helices_option'  testFiles/testFile_RD0260.seq -sh testFiles/testFile_tRNA.shape @OUT.ct -ph 5

# Test Fold_ydao_m2seq_no_data
runFullTest 'ydao_m2seq_no_data' testFiles/testFile_ydaO.seq @OUT.ct -ph 10  -ip 5

# Test Fold_ydao_m2seq_dms
runFullTest 'ydao_m2seq_dms' testFiles/testFile_ydaO.seq @OUT.ct -ph 10  -ip 5 -dms testFiles/testFile_ydaO_M2seq_1D_DMS.txt

# Test Fold_ydao_m2seq_ex
runFullTest 'ydao_m2seq_ex' testFiles/testFile_ydaO.seq @OUT.ct -ph 10  -ip 5 -x testFiles/testFile_ydaO_M2seq_EX.txt

# Test Fold_ydao_m2seq_ex_dms
runFullTest 'ydao_m2seq_ex_dms' testFiles/testFile_ydaO.seq @OUT.ct -ph 10  -ip 5 -dms testFiles/testFile_ydaO_M2seq_1D_DMS.txt -x testFiles/testFile_ydaO_M2seq_EX.txt

# Test ShapeKnots_experimentalpair_option
runFullTest 'experimental_pair_option' testFiles/testFile_bistable.seq @OUT.ct -x testFiles/testFile_bistable.pmat -xs 0.5

# Test ShapeKnots_experimentalpair_option with columnar input
runFullTest 'experimental_pair_option' testFiles/testFile_bistable.seq @OUT.ct -x testFiles/testFile_bistable.columns.pmat -xs 0.5

# Test ShapeKnots_dms_option
runFullTest 'dms_option' testFiles/testFile_bistable.seq @OUT.ct -dms testFiles/testFile_bistable.dms

# Test dmsnt_option
runFullTest 'dmsnt_option' testFiles/testFile_5Se.seq @OUT.ct -dmsnt testFiles/testFile_5S.dmsnt

# Test dmsnt_option
runFullTest 'dmsnt_ex_option' testFiles/testFile_5Se.seq @OUT.ct -dmsnt testFiles/testFile_5S.dmsnt -x testFiles/testFile_5S_bonus.txt

# Test ShapeKnots_ydao_cons
runFullTest 'ydao_cons' testFiles/testFile_ydaO.seq @OUT.ct -ph 10 -ip 5 -c testFiles/testFile_ydaO_cons.txt



endTestBlock # End a group of tests. Also cleans up orphaned files etc.
