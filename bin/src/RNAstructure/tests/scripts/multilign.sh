# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r $EXE DynalignDotPlot # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

initRef 'without_options' _WO # initialize variables of the form *_WO (e.g. OUT_WO, OK_WO  etc) to point to the 'without_options' test files because these are used by multiple other tests.

if [[ $SMP ]]; then SMP_PROCS="Processors = 2"; else SMP_PROCS= ;fi

# Delete intermediate and temporary files from previous tests.
function cleanTmp() { rm -f *.aout *.dsv ${EXE}_*.conf; }

# Test multilign_without_options_singly.
cleanTmp
initTest 'without_options_singly'  && {  # surround with brackets in case test is excluded
echo "
Seq1 = testFiles/testFile_RD0260.seq
Seq2 = testFiles/testFile_RD0500.seq
Seq3 = testFiles/testFile_RA7680.seq
CT1 = ${OUT}_seq1.ct
CT2 = ${OUT}_seq2.ct
CT3 = ${OUT}_seq3.ct
SequenceNumber = 3
Alignment = ${OUT}_alignment.ali
$SMP_PROCS" > $TEST.conf
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct ${OKBASE_WO}_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct ${OKBASE_WO}_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct ${OKBASE_WO}_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali ${OKBASE_WO}_alignment_OK.ali 'ali'
endTest; }

# Test multilign_without_options_groups.
cleanTmp
initTest 'without_options_groups'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct ${OKBASE_WO}_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct ${OKBASE_WO}_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct ${OKBASE_WO}_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali ${OKBASE_WO}_alignment_OK.ali 'ali'
endTest; }

# Test multilign_constraint_first_option.
cleanTmp
initTest 'constraint_first_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "Constraint1 = testFiles/testFile_folding1.con" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_constraint_all_option.
cleanTmp
initTest 'constraint_all_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "Constraint1 = testFiles/testFile_folding4.con" >> "$TEST.conf"
echo "Constraint2 = testFiles/testFile_folding4.con" >> "$TEST.conf"
echo "Constraint3 = testFiles/testFile_folding4.con" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_dna_option.
cleanTmp
initTest 'dna_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "DNA = 1" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_gap_option.
cleanTmp
initTest 'gap_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "Gap = 1.2" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_index_sequence_option.
cleanTmp
initTest 'index_sequence_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "IndexSeq = 2" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_inserts_deselected_option.
cleanTmp
initTest 'inserts_deselected_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "Insert = 0" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_iterations_option.
cleanTmp
initTest 'iterations_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "Iterations = 1" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_keep_intermediate_option.
cleanTmp
initTest 'keep_intermediate_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "KeepIntermediate = 1" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runSubTest '1-1_seq1' DynalignDotPlot 1.1_testFile_RD0260_testFile_RD0500.dsv ${OUT}_1-1_seq1.ps
runSubTest '1-1_seq2' DynalignDotPlot 1.1_testFile_RD0260_testFile_RD0500.dsv ${OUT}_1-1_seq2.ps --sequence2
runSubTest '1-2_seq1' DynalignDotPlot 1.2_testFile_RD0260_testFile_RA7680.dsv ${OUT}_1-2_seq1.ps
runSubTest '1-2_seq2' DynalignDotPlot 1.2_testFile_RD0260_testFile_RA7680.dsv ${OUT}_1-2_seq2.ps --sequence2
runSubTest '2-1_seq1' DynalignDotPlot 2.1_testFile_RD0260_testFile_RD0500.dsv ${OUT}_2-1_seq1.ps
runSubTest '2-1_seq2' DynalignDotPlot 2.1_testFile_RD0260_testFile_RD0500.dsv ${OUT}_2-1_seq2.ps --sequence2
runSubTest '2-2_seq1' DynalignDotPlot 2.2_testFile_RD0260_testFile_RA7680.dsv ${OUT}_2-2_seq1.ps
runSubTest '2-2_seq2' DynalignDotPlot 2.2_testFile_RD0260_testFile_RA7680.dsv ${OUT}_2-2_seq2.ps --sequence2
runDiff ${OUT}_1-1_seq1.ps @OKBASE_1-1_seq1_OK.ps '1-1_seq1'
runDiff ${OUT}_1-1_seq2.ps @OKBASE_1-1_seq2_OK.ps '1-1_seq2'
runDiff ${OUT}_1-2_seq1.ps @OKBASE_1-2_seq1_OK.ps '1-2_seq1'
runDiff ${OUT}_1-2_seq2.ps @OKBASE_1-2_seq2_OK.ps '1-2_seq2'
runDiff ${OUT}_2-1_seq1.ps @OKBASE_2-1_seq1_OK.ps '2-1_seq1'
runDiff ${OUT}_2-1_seq2.ps @OKBASE_2-1_seq2_OK.ps '2-1_seq2'
runDiff ${OUT}_2-2_seq1.ps @OKBASE_2-2_seq1_OK.ps '2-2_seq1'
runDiff ${OUT}_2-2_seq2.ps @OKBASE_2-2_seq2_OK.ps '2-2_seq2'
runDiff ${OUT}_seq1.ct     ${OKBASE_WO}_seq1_OK.ct     'seq1'
runDiff ${OUT}_seq2.ct     ${OKBASE_WO}_seq2_OK.ct     'seq2'
runDiff ${OUT}_seq3.ct     ${OKBASE_WO}_seq3_OK.ct     'seq3'
runDiff ${OUT}_alignment.ali ${OKBASE_WO}_alignment_OK.ali 'ali'
runDiff 1.1_testFile_RD0260_testFile_RD0500.aout @OKBASE_1-1_OK.ali '1-1-ali'
runDiff 1.2_testFile_RD0260_testFile_RA7680.aout @OKBASE_1-2_OK.ali '1-2-ali'
runDiff 2.1_testFile_RD0260_testFile_RD0500.aout @OKBASE_2-1_OK.ali '2-1-ali'
runDiff 2.2_testFile_RD0260_testFile_RA7680.aout @OKBASE_2-2_OK.ali '2-2-ali'
endTest; }

# Test multilign_local_folding_option.
cleanTmp
initTest 'local_folding_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "Local = 1" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_max_dsv_change_option.
cleanTmp
initTest 'max_dsv_change_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "MaxDsvChange = 10" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_max_pairs_option.
cleanTmp
initTest 'max_pairs_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "MaxPairs = 500" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_max_structures_option.
cleanTmp
initTest 'max_structures_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "MaxStructures = 2" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_percent_difference_option.
cleanTmp
initTest 'percent_difference_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "MaxPercent = 1" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_percent_difference_single_fold_option.
cleanTmp
initTest 'percent_difference_single_fold_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "MaxPercentSingle = 1" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# If the executable is SMP, test multilign-smp_processors_option.
if [[ $SMP ]]; then
cleanTmp
initTest 'processors_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "Processors = 3" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct ${OKBASE_WO}_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct ${OKBASE_WO}_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct ${OKBASE_WO}_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali ${OKBASE_WO}_alignment_OK.ali 'ali'
endTest; }
fi

# Test multilign_random_option.
cleanTmp
false && # SKIP the random_option test because it currently doesn't accept a seed value and the results differ on each platform.
initTest 'random_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "Random = 1" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_separation_option.
cleanTmp
initTest 'separation_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "Separation = 6" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_shape_first_option.
cleanTmp
initTest 'shape_first_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "SHAPE1 = testFiles/testFile_tRNA_dummy.shape" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_shape_all_option.
cleanTmp
initTest 'shape_all_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "SHAPE1 = testFiles/testFile_tRNA_dummy.shape" >> "$TEST.conf"
echo "SHAPE2 = testFiles/testFile_tRNA_dummy.shape" >> "$TEST.conf"
echo "SHAPE3 = testFiles/testFile_tRNA_dummy.shape" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_shape_intercept_option.
cleanTmp
initTest 'shape_intercept_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "SHAPE1 = testFiles/testFile_tRNA_dummy.shape" >> "$TEST.conf"
echo "SHAPEintercept = -0.5" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_shape_slope_option.
cleanTmp
initTest 'shape_slope_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "SHAPE1 = testFiles/testFile_tRNA_dummy.shape" >> "$TEST.conf"
echo "SHAPEslope = 2.0" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_temperature_option.
cleanTmp
initTest 'temperature_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "Temperature = 200" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_window_size_alignment_option.
cleanTmp
initTest 'window_size_alignment_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "WindowAlign = 5" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }

# Test multilign_window_size_base_pair_option.
cleanTmp
initTest 'window_size_base_pair_option'  && {  # surround with brackets in case test is excluded
echo "InSeq = {testFiles/testFile_RD0260.seq;testFiles/testFile_RD0500.seq;testFiles/testFile_RA7680.seq;}" > "$TEST.conf"
echo "OutCT = {${OUT}_seq1.ct;${OUT}_seq2.ct;${OUT}_seq3.ct;}" >> "$TEST.conf"
echo "Alignment = ${OUT}_alignment.ali" >> "$TEST.conf"
echo "WindowBP = 4" >> "$TEST.conf"
[[ $SMP ]] && echo "Processors = 2" >> "$TEST.conf"
runTest $EXE "$TEST.conf"
runDiff ${OUT}_seq1.ct @OKBASE_seq1_OK.ct 'seq1'
runDiff ${OUT}_seq2.ct @OKBASE_seq2_OK.ct 'seq2'
runDiff ${OUT}_seq3.ct @OKBASE_seq3_OK.ct 'seq3'
runDiff ${OUT}_alignment.ali @OKBASE_alignment_OK.ali 'ali'
endTest; }
cleanTmp

endTestBlock # End a group of tests. Also cleans up orphaned files etc.

