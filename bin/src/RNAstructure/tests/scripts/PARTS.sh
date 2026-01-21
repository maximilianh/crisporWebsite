# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.


beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs

# Test PARTS in MAP mode.
initTest 'PARTS_MAP' && { #skip block if test is excluded.
echo "seq1 testFiles/testFile_RD0260.seq" > $TEST.conf
echo "seq2 testFiles/testFile_RD0500.seq" >> $TEST.conf
echo "mode map" >> $TEST.conf
echo "seq1_map_ct_op ${OUT}_seq1.ct" >> $TEST.conf
echo "seq2_map_ct_op ${OUT}_seq2.ct" >> $TEST.conf
echo "map_aln_op ${OUT}_alignment.ali" >> $TEST.conf
runTest @EXE $TEST.conf
runDiff ${OUT}_seq1.ct PARTS/PARTS_map_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct PARTS/PARTS_map_seq2_OK.ct 'seq2'
runDiff ${OUT}_alignment.ali PARTS/PARTS_map_alignment_OK.ali 'ali'
endTest; }


# Test PARTS in probability mode.
initTest 'PARTS_probability' && { #skip block if test is excluded.
echo "seq1 testFiles/testFile_RD0260.seq" > $TEST.conf
echo "seq2 testFiles/testFile_RD0500.seq" >> $TEST.conf
echo "mode pp" >> $TEST.conf
echo "seq1_pp_op ${OUT}_seq1.txt" >> $TEST.conf
echo "seq2_pp_op ${OUT}_seq2.txt" >> $TEST.conf
runTest @EXE $TEST.conf
runDiff ${OUT}_seq1.txt PARTS/PARTS_pp_seq1_OK.txt 'seq1'
runDiff ${OUT}_seq2.txt PARTS/PARTS_pp_seq2_OK.txt 'seq2'
endTest; }

# Test PARTS in stochastic mode.
initTest 'PARTS_stochastic' && { #skip block if test is excluded.
echo "seq1 testFiles/testFile_RD0260.seq" > $TEST.conf
echo "seq2 testFiles/testFile_RD0500.seq" >> $TEST.conf
echo "mode stochsample" >> $TEST.conf
echo "seq1_sample_ct_op ${OUT}_seq1.ct" >> $TEST.conf
echo "seq2_sample_ct_op ${OUT}_seq2.ct" >> $TEST.conf
echo "sample_aln_op ${OUT}_alignment.ali" >> $TEST.conf
echo "nsamp 250" >> $TEST.conf
echo "seed 24" >> $TEST.conf
runTest @EXE $TEST.conf
runDiff ${OUT}_seq1.ct PARTS/PARTS_stochastic_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct PARTS/PARTS_stochastic_seq2_OK.ct 'seq2'
runDiff ${OUT}_alignment.ali PARTS/PARTS_stochastic_alignment_OK.ali 'ali'
endTest; }

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
