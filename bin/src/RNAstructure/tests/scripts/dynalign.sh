# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs

CONFIG_TEXT='
inseq1 = %s
inseq2 = %s
outct = @TEST_seq1.ct
outct2 = @TEST_seq2.ct
aout = @TEST_alignment.ali'

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
    [[ $SMP ]] && printf "\n num_processors = 2" >> "$TEST.conf"

    ######### Run dynalign #########
    runTest $EXENAME $TEST.conf

    ######### Compare outputs (CT files and alignments) #########
    runDiff @TEST_seq1.ct        @OKBASE_seq1_OK.ct  'seq1'
    runDiff @TEST_seq2.ct        @OKBASE_seq2_OK.ct  'seq2'
    runDiff @TEST_alignment.ali  @OKBASE_alignment_OK.ali  'ali'

    ######### Delete the conf file and do end-of-test processing #########
    #rm -f $TEST.conf
    endTest
}
#define some input files
INPUT1=testFiles/testFile_RD0260.seq
INPUT2=testFiles/testFile_RD0500.seq

# Actual tests are below:
DynTest 'general_test'         "$CONFIG_TEXT"                  $INPUT1 $INPUT2
DynTest 'dna_option'           "$CONFIG_TEXT\n DNA = 1"        $INPUT1 $INPUT2
DynTest 'dna_alphabet_option'  "$CONFIG_TEXT\n alphabet = dna" $INPUT1 $INPUT2

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
