# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

SMP_SETTING=  # SMP_SETTING is blank unless this is the SMP variant of TurboFold in which case, we set the number of processors (minimum 2)
if [[ $SMP ]]; then 
	N_PROCS=$({ getconf _NPROCESSORS_ONLN || getconf NPROCESSORS_ONLN ; } 2>/dev/null) # get number of processors (getconf is POSIX, but _NPROCESSORS_ONLN is optional)
	[[ $N_PROCS -gt 1 ]] || N_PROCS=2 # If N_PROCS is empty or 0 or 1, then set it to 2.
	SMP_SETTING="Processors = $N_PROCS"
	isQuiet || warnmsg "Using '$SMP_SETTING' in SMP mode."
fi

SHAPEFILE1=testFiles/testFile_random_dummy.shape
SHAPEFILE2=testFiles/testFile_random_dummy_114.shape # remove bases over 114 due to shorter sequence
SHAPEFILE3=$SHAPEFILE1

GEN_FILES="ec5s.ct  ca5s.ct  5s_8.ct"
ALT_FILES="met-vol.ct  met-fer.ct  met-the.ct"
RD0_FILES="RA7680.ct  RD0260.ct  RD0500.ct"

################################################
# turboTest -- Run a TurboFold test.
################################################
# This process includes:
#   1. Writing the turbofold configuration file
#   2. running TurboFold
#   3. Comparing the output (alignment file and CT files) with references ("OK" files)
# Usage: 
#   turboTest <TEST_NAME>  <CONF_TEXT>  [OPTIONS]  <FILE_LIST>
#      TEST_NAME -- The name of the test, e.g. 'general_without_options'
#      CONF_TEXT -- The literal text for the configuration file passed to TurboFold.
#                   This text CAN contain test parameters (such as @TEST, @OUT, @OK etc)
#      FILE_LIST -- The list of output files to verify.
#        These are usually listed as short CT names such as 'ec5s.ct ca5s.ct 5s_8.ct' 
#	     The corresponding OK file is "guessed" to be OK_FILE=@OK_<NAME>.ct
#        (where NAME is the ct-file name without the .ct extension)
#        Custom OK files can be specified by passing the names in the 
#        following form:  <OUTPUT_FILE>==<OK_FILE>
#        Each OK_FILE *MUST* be prefixed with '@OK_', '@OKBASE_' or '@OKROOT_'.
#        The DIFF output file in each case will be named @TEST_<NAME>_diff.txt
#        After TurboFold runs, each output file is renamed from 
#           <NAME>.ct to @TEST_<NAME>_output.ct
#
#      OPTIONS -- Can be any of the following:
#         --ref <OTHER_TESTNAME>   
#               Use <OTHER_TESTNAME> to prefix the OK file names in place of 
#			    the current <TESTNAME>.
#         --stdref
#               Use general_without_options to prefix the OK file names.
#			    (same as '--ref general_without_options')
#         --altref
#               Use general_without_options_alternate to prefix the OK file names.
#			    (same as '--ref general_without_options_alternate')
function turboTest() {
	# Initilize the test with the testing system
	initTest "$1" || return  # exit if the test is excluded etc.

	# Write the conf file. (replacing special Test variables like @TEST, @OUT, @OK etc)
	replaceVarsInText "$2" > "$TEST.conf"  # Write output to $TEST.conf e.g. "TurboFold_general_without_options.conf"
	shift 2 # remove first two args.

	# Loop through remaining arguments. These will be either --options or 
	#      sets of OUTPUT==OK files.
	local files=()
	while (($#)); do
		case "$1" in
			--ref) shift; setTestRef "$1" ;;        # use the specified TESTNAME OK files as references for this test. This affects @OKBASE, @OK, and @OKFILE
			--stdref) setTestRef 'general_without_options' ;; # use the general_without_options OK files as references for this test.
			--altref) setTestRef 'general_without_options_alternate' ;; # use the general_without_options_alternate OK files as references for this test.
			*) files+=("$1") ;; # append to list of files
		esac
		shift # next argument
	done

	# Now run TurboFold with the given conf
	runTest @EXE "$TEST.conf"

	# Run a diff on the alignment file 
	# diff output is written to ${TEST}_aln_diff.txt e.g.  "TurboFold_general_without_options_aln_diff.txt"
	runDiff @OUT.aln @OK.aln 'aln'

	# Now for each of the listed input sequences, do a diff with the expected output CT
	local ctName okName
	for ctName in "${files[@]}"; do  # loop through each OUTPUT==OK file
		if [[ "$ctName" == *'=='* ]]; then
			# Custom CT and OK-file names.
			okName="${ctName##*==}" # just the part after ==
			ctName="${ctName%==*}"  # just the part before ==
		else
			okName=
		fi
		local base="${ctName%.*}"      # base filename (no extension)
		local ext="${ctName:${#base}}" # extension
		[[ $okName ]]  ||  okName="@OK{$base}$ext"

		# Iif the output file exists, do a DIFF.
		if [[ -e "$ctName" ]]; then
			# Rename the file to a standard form.
			local saveCtName="${TEST}_${base}_output$ext"
			mv "$ctName" "$saveCtName" || errmsg "Failed to rename output file: $ctName"
			# Generate name for DIFF output
			local diffName="${TEST}_${base}_diff.txt"			
			runDiff "$saveCtName"  "$okName"  "$diffName"  # diff output written to  ${TEST}_${seqName}_diff.txt
		else
			errmsg "Missing expected output file: $ctName"
		fi
	done
	endTest # perform end-of-test duties
}

# Test TurboFold_without_options_singly.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_without_options_singly' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
OutAln = @OUT.aln
SequenceNumber = 3
Mode = MEA
$SMP_SETTING"   --stdref  $GEN_FILES

# Test fasta input files.
turboTest 'general_without_options_fasta' "
Seq1 = testFiles/testFile_ec5s.fasta
Seq2 = testFiles/testFile_ca5s.fasta
Seq3 = testFiles/testFile_5s_8.fasta
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
OutAln = @OUT.aln
SequenceNumber = 3
Mode = MEA
$SMP_SETTING"   --stdref  $GEN_FILES

# Test TurboFold_without_options_groups.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_without_options_groups' "
InSeq = {testFiles/testFile_ec5s.seq;testFiles/testFile_ca5s.seq;testFiles/testFile_5s_8.seq;}
OutCT = {ec5s.ct;ca5s.ct;5s_8.ct;}
OutAln = @OUT.aln
Mode = MEA
$SMP_SETTING"   --stdref  $GEN_FILES

# TurboFold_fasta_multiple
#    This reads a single FASTA input file containing three sequences.
#    The names of the three output CT files are also specified.
turboTest 'fasta_multiple' "
InFasta = testFiles/testFile_TurboFold.fasta
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
OutAln = @OUT.aln
SequenceNumber = 3
Mode = MEA
$SMP_SETTING"   --stdref  $GEN_FILES

# TurboFold_fasta_multiple_autoname
#    This reads a single FASTA input file containing three 
#    sequences, but it does NOT specify the names of 
#    output CT files, so those names are auto-generated from
#    the sequence names (given in the FASTA file)
turboTest 'fasta_multiple_autoname' "
InFasta = testFiles/testFile_TurboFold.fasta
OutAln = @OUT.aln
SequenceNumber = 3
Mode = MEA
$SMP_SETTING" --stdref \
	'1_ecoli_5s.ct'==@OK{ec5s}.ct \
	'2_CA5SRNA.ct'==@OK{ca5s}.ct  \
	'3_THERMOCOCCUS_CELER_[XTR]_(THERMOCOCCALES).ct'==@OK{5s_8}.ct

# Test an alternate TurboFold_general_without_options.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
# Note also that these alternate test files are used to test ProbKnot mode only.
# This alternate test is needed to ensure that the ProbKnot_without_options test is different from a default test where ProbKnot mode isn't specified.
turboTest 'general_without_options_alternate' "
Seq1 = testFiles/testFile_met-vol.seq
Seq2 = testFiles/testFile_met-fer.seq
Seq3 = testFiles/testFile_met-the.seq
CT1 = met-vol.ct
CT2 = met-fer.ct
CT3 = met-the.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA" --altref $ALT_FILES

# Test TurboFold_general_gamma_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_gamma_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
Gamma = 0.9"  $GEN_FILES

# Test TurboFold_general_iterations_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_iterations_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
Iterations = 1"  $GEN_FILES

# Test TurboFold_general_shape_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_shape_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
OutAln = @OUT.aln
SHAPE1 = $SHAPEFILE1
SHAPE2 = $SHAPEFILE2
SHAPE3 = $SHAPEFILE3"  $GEN_FILES

# Test TurboFold_general_shape_intercept_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_shape_intercept_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
OutAln = @OUT.aln
SHAPE1 = $SHAPEFILE1
SHAPE2 = $SHAPEFILE2
SHAPE3 = $SHAPEFILE3
SHAPEintercept = 0.1"  $GEN_FILES


# Test TurboFold_general_shape_slope_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_shape_slope_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
OutAln = @OUT.aln
SHAPE1 = $SHAPEFILE1
SHAPE2 = $SHAPEFILE2
SHAPE3 = $SHAPEFILE3
SHAPEslope = 0.1"   $GEN_FILES


# Test TurboFold_general_temperature_option.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_temperature_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
Temperature = 330"   $GEN_FILES


[[ $SMP ]] &&  # If the executable is SMP, test the processors option.
turboTest 'processors_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
OutAln = @OUT.aln
SequenceNumber = 3
Processors = 3;
Mode = MEA"   --stdref  $GEN_FILES


# Test TurboFold_MEA-mode_without_options.
turboTest 'MEA-mode_without_options' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA"   --stdref  $GEN_FILES

# Test an alternate TurboFold_MEA-mode_without_options.
# These test files are not used for all MEA mode tests, just max structures, percent difference, and window size.
turboTest 'MEA-mode_without_options_alternate' "
Seq1 = testFiles/testFile_RA7680.seq
Seq2 = testFiles/testFile_RD0260.seq
Seq3 = testFiles/testFile_RD0500.seq
CT1 = RA7680.ct
CT2 = RD0260.ct
CT3 = RD0500.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA"  $RD0_FILES

# Test TurboFold_MEA-mode_gamma_option.
turboTest 'MEA-mode_gamma_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
MeaGamma = 0.2"   $GEN_FILES

# Test TurboFold_MEA-mode_max_structures_option.
turboTest 'MEA-mode_max_structures_option' "
Seq1 = testFiles/testFile_RA7680.seq
Seq2 = testFiles/testFile_RD0260.seq
Seq3 = testFiles/testFile_RD0500.seq
CT1 = RA7680.ct
CT2 = RD0260.ct
CT3 = RD0500.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
MaxStructures = 5
Window = 0"   $RD0_FILES

# Test TurboFold_MEA-mode_percent_difference_option.
turboTest 'MEA-mode_percent_difference_option' "
Seq1 = testFiles/testFile_RA7680.seq
Seq2 = testFiles/testFile_RD0260.seq
Seq3 = testFiles/testFile_RD0500.seq
CT1 = RA7680.ct
CT2 = RD0260.ct
CT3 = RD0500.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
MaxPercent = 4
Window = 0"   $RD0_FILES

# Test TurboFold_MEA-mode_window_size_option.
turboTest 'MEA-mode_window_size_option' "
Seq1 = testFiles/testFile_RA7680.seq
Seq2 = testFiles/testFile_RD0260.seq
Seq3 = testFiles/testFile_RD0500.seq
CT1 = RA7680.ct
CT2 = RD0260.ct
CT3 = RD0500.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = MEA
Window = 0"   $RD0_FILES

# Test TurboFold_ProbKnot-mode_without_options.
turboTest 'ProbKnot-mode_without_options' "
Seq1 = testFiles/testFile_met-vol.seq
Seq2 = testFiles/testFile_met-fer.seq
Seq3 = testFiles/testFile_met-the.seq
CT1 = met-vol.ct
CT2 = met-fer.ct
CT3 = met-the.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = ProbKnot"   $ALT_FILES

# Test TurboFold_ProbKnot-mode_iterations_option.
turboTest 'ProbKnot-mode_iterations_option' "
Seq1 = testFiles/testFile_met-vol.seq
Seq2 = testFiles/testFile_met-fer.seq
Seq3 = testFiles/testFile_met-the.seq
CT1 = met-vol.ct
CT2 = met-fer.ct
CT3 = met-the.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = ProbKnot
PkIterations = 10"   $ALT_FILES

# Test TurboFold_ProbKnot-mode_min_helix_length_option.
turboTest 'ProbKnot-mode_min_helix_length_option' "
Seq1 = testFiles/testFile_met-vol.seq
Seq2 = testFiles/testFile_met-fer.seq
Seq3 = testFiles/testFile_met-the.seq
CT1 = met-vol.ct
CT2 = met-fer.ct
CT3 = met-the.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = ProbKnot
MinHelixLength = 1"   $ALT_FILES

# Test TurboFold_Threshold-mode_without_options.
turboTest 'Threshold-mode_without_options' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
OutAln = @OUT.aln
SequenceNumber = 3
$SMP_SETTING
Mode = Threshold"   $GEN_FILES

# Test TurboFold_Threshold-mode_threshold_option.
turboTest 'Threshold-mode_threshold_option' "
Seq1 = testFiles/testFile_ec5s.seq
Seq2 = testFiles/testFile_ca5s.seq
Seq3 = testFiles/testFile_5s_8.seq
CT1 = ec5s.ct
CT2 = ca5s.ct
CT3 = 5s_8.ct
SequenceNumber = 3
$SMP_SETTING
Mode = Threshold
OutAln = @OUT.aln
Threshold = 0.92"   $GEN_FILES

# Test TurboFold_Rsample_option.
turboTest 'Rsample_option' "
InSeq = {TurboFold/input/tRNA.seq;TurboFold/input/tRNA2.seq;TurboFold/input/RD0260.seq;TurboFold/input/RD0500.seq;TurboFold/input/RA7680.seq;}
OutCT = {tRNA.ct;tRNA2.ct;RD0260.ct;RD0500.ct;RA7680.ct;}
Mode = MEA
UseRsample = 1
SHAPEFiles = {TurboFold/input/tRNA.shape;TurboFold/input/tRNA2.shape;;;;}
Seed = 1
OutAln = @OUT.aln
$SMP_SETTING" tRNA.ct tRNA2.ct RD0260.ct RD0500.ct RA7680.ct

# Should give an error about a missing shape file.
initTest 'shape_file_missing' && {
replaceVarsInText "InSeq = {testFiles/testFile_ec5s.seq;testFiles/testFile_ca5s.seq;}
OutCT = {@TEST_1.ct;@TEST_2.ct}
Mode = MEA
SHAPEFiles = {missing_shape_file.shape;}
$SMP_SETTING" > "$TEST.conf"
# run the test, but expect an exit code of 1 and perform the diff on the output from stderr
runProgAndDiff "$TEST.conf"   ---exit=1 ---stderr 
endTest; }

# SHAPEFILE1 is too long for ca5s.seq -- should give an error about an invalid nucleotide position
initTest 'shape_file_bad' && {
replaceVarsInText "InSeq = {testFiles/testFile_ec5s.seq;testFiles/testFile_ca5s.seq;}
OutCT = {@TEST_1.ct;@TEST_2.ct}
OutAln = @OUT.aln
Mode = MEA
SHAPEFiles = {$SHAPEFILE1;$SHAPEFILE1}
$SMP_SETTING" > "$TEST.conf"
# run the test, but expect an exit code of 1 and perform the diff on the output from stderr
# set RNA_WARNINGS environment variable to ERR, which causes warnings to be written to STDERR instead of STDOUT.
RNA_WARNINGS=ERR \
runProgAndDiff "$TEST.conf"   ---stderr
endTest; }

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
