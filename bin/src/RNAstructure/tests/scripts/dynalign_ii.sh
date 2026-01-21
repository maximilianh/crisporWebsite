# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs

CONFIG_TEXT='
inseq1 = %s
inseq2 = %s
outct  = @TEST_seq1.ct
outct2 = @TEST_seq2.ct
aout   = @TEST_alignment.ali

#fgap is the per nucleotide insert penalty for alignments:
fgap = .4
#slope is the per nucleotide free energy penalty for inserted domains
slope = 0.1
#intercept is the initiation free energy penalty for inserted domains
intercept = 0.5
#maxtrace is the masximum number of predicted structures:
maxtrace = 750
#percent is the maximum % change in free energy from the lowest free energy structure
percent = 20
#bpwin is the base pair window
bpwin = 2
#awin is the alignment window
awin = 1'

###############################################################
# This function will be used to template the dynalign tests.
# Usage: DynTest <TEST_NAME> <CONF_TEXT> [<REPLACEMENTS>...]
###############################################################
# TEST_NAME is the name of the test, e.g 'without options'
# CONF_TEXT is the text that should be written to the conf file. e.g. "$CONFIG_TEXT"
# REPLACEMENTS are the replacement values of printf placeholders (e.g. %s) in CONF_TEXT (if any).
function DynTest() {
    initTest "$1" || return 0 # skip the test if it is excluded

    ######### Write the dynalign configuration file ###############
    # $2 is the CONF_TEXT. Replace any occurances of "@TEST" with the actual test name (e.g. "dynalign_general_test").
    local text="${2//'@TEST'/$TEST}"
    # Write the configuration file (e.g. "dynalign_general_test.conf") for each test.
    printf "$text" "${@:3}" > "$TEST.conf" # Pass all arguments (after the first two) to printf to replace '%s' placeholders
    [[ $SMP ]] &&  printf "\n num_processors = 2" >> "$TEST.conf"

    ######### Run dynalign #########
    runTest $EXENAME $TEST.conf

    ######### Compare outputs (CT files and alignments) #########
    runDiff @TEST_seq1.ct        @OKBASE_seq1_OK.ct  'seq1'
    runDiff @TEST_seq2.ct        @OKBASE_seq2_OK.ct  'seq2'
    runDiff @TEST_alignment.ali  @OKBASE_alignment_OK.ali  'ali'

    ######### Do end-of-test processing #########
    endTest
}

#define some input files
INPUT1=testFiles/testFile_RA1540.seq
INPUT2=testFiles/testFile_RL3280.seq

# Actual tests are below:
DynTest 'general_test'         "$CONFIG_TEXT"                  $INPUT1 $INPUT2
DynTest 'dna_option'           "$CONFIG_TEXT\n DNA = 1"        $INPUT1 $INPUT2
DynTest 'alphabet_dna_option'  "$CONFIG_TEXT\n alphabet = dna" $INPUT1 $INPUT2

# Test for a segfault bug that existed in a previous version of dynalign
##### NOTE: This test is disabled because it takes a long time.
##### However, it should be re-enabled and tested once in a while
##### By setting the environment variable ENABLE_LONG_TESTS to a non-empty value.
[[ $ENABLE_LONG_TESTS ]] && DynTest 'segfault_test'  "$CONFIG_TEXT" dynalign_ii/segfault_seq1.fasta  dynalign_ii/segfault_seq2.fasta

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
