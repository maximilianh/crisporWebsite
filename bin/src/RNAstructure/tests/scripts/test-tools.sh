#!/bin/bash
source "${BASH_SOURCE[0]%/*}"/script-utils.sh

################################################################################
#                               test-tools.sh
# This script provides several functions and for convenience in composing 
# test scripts.
# Note: These tools rely on several variables being set by the parent script:
#  EXENAME  - the test exe file name (e.g. Fold or Fold-smp etc)
#  EXEBASE  - usually the same as EXENAME, but see notes in initTest.
#             e.g. for Fold-smp EXENAME is 'Fold-smp' but EXEBASE is 'Fold'
#  EXEDIR   - relative path to test executables (e.g. ../exe)
#  MAKEROOT - relative path to the directory from which `make EXENAME` should 
#             be called (e.g. '..')
#  OKDIR    - The folder in which the test reference output ("OK") files are stored.
#             Usually the same as the EXEBASE name. e.g. Fold-smp output files 
#             are compared with those in the 'Fold' folder.
# Set Defaults:
: ${EXEDIR:=../exe}  ${MAKEROOT:=..}  ${TEST_OPTS:='DumpErrors ShowOutput'}
################################################################################


################################################################################
## checkErrorFiles - Checks for the existance of files passed as arguments. If the 
##     files are found and are not empty, they are assumed to be 
##     error files. They trigger output to stderr as well as the 
##     copying of all test-related files (not just the error files).
## Usage:   checkErrorFiles TESTNAME [FILES...]
## Example: checkErrorFiles 'efn2_without_options' efn2_without_options_errors.txt efn2_without_options_diff_output.txt
################################################################################
function checkErrorFiles() {
    local testName="$1"
    local errorLevel=0
    local file msg errNum numLines errorFiles=() 
    DETECTED_TEST_ERRORS=
    shift # remove the testName argument
    for file in "$@"; do
        if [[ -s $file ]]; then
            case $file in 
                *_errors.txt) 
                    errNum=10 #major error
                    DETECTED_TEST_ERRORS+=,'Error'
                    ;;

                *_diff_output.txt|*_diff.txt)
                    errNum=5 #possible error
                    DETECTED_TEST_ERRORS+=,'Diff'
                    ;;
                *)
                    errNum=1 #possible error
                    DETECTED_TEST_ERRORS+=,'Unknown'
                    ;;
            esac
            ((errNum>errorLevel)) && ((errorLevel=errNum))
            errorFiles+=("$file")

            # if ((errNum>errorLevel)); then
            #     ((errorLevel=errNum))
            #     numLines=$(wc -l $file | { read -r first _; echo $first; } || echo UNKNOWN )
            #     #failureReason=" (see $fileType File $file; $numLines lines)"
            #     failures+=("$file (${numLines} lines)")
            # fi
        fi
    done
    DETECTED_TEST_ERRORS=${DETECTED_TEST_ERRORS:1} #remove first comma
    if ((errorLevel)); then
        #if ((errorLevel>=10)); then
        #    msg='Test Failure'
        #else
            msg='Possible Failure'
        #fi
        errmsg "  ! ${msg}: $testName"
        if isTestOpt 'DumpErrors' || isTestOpt 'ListErrors'; then
            for file in "${errorFiles[@]}"; do 
                numLines=$(wc -l $file | { read -r first _; echo $first; } || echo UNKNOWN )
                errmsg "  ! => $file (${numLines} lines)"
                if isTestOpt 'DumpErrors'; then
                    head -n 4 $file | dumpErrorFile
                    ((numLines>4)) && echo "... ($((numLines-4)) more lines)" | dumpErrorFile
                fi
            done
        fi
    else
        okmsg  "  * Test Passed: $testName"
    fi
    ((errorLevel==0)) # return 0 (success) if no errors were found or 0 otherwise.
}

# Display stdin indented and in a special color
function dumpErrorFile() {
    [[ -t 2 ]] && echo -n "$FMT_BLD$FMT_PRP"
    sed 's/^/  |/g' # replace the start of each line of input with ' |' to align with error message
    [[ -t 2 ]] && echo -n "$FMT_RST"
}

# List of test-specific variable names. These are prefixed with '@' and can be 
# used in command-line arguments for tests.
# They are replaced at runtime with test-specific values.
# For example @OUT may expand to 'Fold_without_options'
# See initTest() for their definitions and assignment
# See fixParams() and replaceVarsInText() for their replacement at runtime.
# IMPORTANT -- This list must be updated whenever a new variable is added.
#    If a parameter name (such as 'OKDIR') starts with another name in the list
#    (such as 'OK'), the shorter name MUST be listed AFTER the longer one
#    because replacements are made using PARTIAL words, e.g. @OKDIRFOO => ${OKDIR}FOO
#    If 'OK' is listed before 'OKDIR' then '@OKDIR/file.ct' would
#    be replaced by '${OK}DIR/file.ct' instead of '${OKDIR}/file.ct'
#    ** A simple way to ensure this is to sort the list in reverse ABC order. **
TEST_PARAMNAME_LIST=(TBASE TESTNAME TEST STDO OUTFILE OUT OKROOT OKFILE OKDIR OKBASE OK INDIR EXENAME EXEBASE EXE ERRO ERRDIR)

################################################################################
# initTest - Report that the test has started.
#             Also, set the variables TEST, OUT, and OUTFILE for convenience 
#             in writing tests and use by other functions.
# Usage:   initTest TESTNAME
# Example: initTest 'without_options'
################################################################################
function initTest() {
    verifyEndTest

    shouldIncludeTest "${EXENAME}_$1" || return  # Check whether this test is included ( depends on INCLUDE_TESTS and EXCLUDE_TESTS environment variables)

    ###### Assign test-specific variables ######
    TESTNAME="$1"  # the name of the test, e.g.   'dna_option'

    # If TESTNAME contains the name of the exe, remove it. E.g.: Fold_dna_option ==> 'dna_option'
    [[ $TESTNAME == ${EXENAME}_*  ]] && TESTNAME="${TESTNAME#${EXENAME}_}"
    [[ $TESTNAME == ${EXEBASE}_*  ]] && TESTNAME="${TESTNAME#${EXEBASE}_}"

    # IMPORTANT -- update 'TEST_PARAMNAME_LIST' when adding or changing parameter names.
    TEST="${EXENAME}_$TESTNAME"   # The full test name, e.g. Fold_without_options
    TBASE="${EXEBASE}_$TESTNAME"  # The full test name using EXEBASE instead of EXENAME, e.g. Fold_without_options
    OUT="${TEST}_output"          # The output file name without extension, e.g. Fold_without_options_output
    OUTFILE="$OUT$EXT"            # The output file name WITH extension, e.g. Fold_without_options_output.ct
    OKBASE="${OKDIR}/$TBASE"      # The "OK" directory and base file name (before the OK prefix). Use this when there are multiple OK files for a single test, eg. @OKBASE-seq1_OK.ps, @OKBASE-seq2_OK.ps etc.
    OKROOT="${OKDIR}/$EXEBASE"    # The "OK" directory and EXEBASE prefix with no TESTNAME. Use this when you need a prefix that can apply to more than one test.
    OK="${OKBASE}_OK"             # The "OK" file name without extension, e.g. Fold/Fold_without_options_OK
    OKFILE="$OK$EXT"              # The "OK" file name WITH extension, e.g. Fold/Fold_without_options_OK.ct
    STDO="${TEST}_stdout.txt"     # The file that should capture the program's stdout.  e.g. Fold_without_options_stdout.txt
    ERRO="${TEST}_errors.txt"     # The file that should capture the program's stderr.  e.g. Fold_without_options_errors.txt
    CMDLOG="${TEST}_commands.log" # A log of all test-related commands run by the testing system (e.g. runTest, runDiff), for debugging and/or verifying commands.
    status "Test $EXENAME '$TESTNAME' started..."
    
    # Export some variables for sub-processes. This could be useful in conjunction with 
    # the use of a TEST_PROG_ANALYZER (e.g. valgrind) to write separate log files for each program.
    # e.g. TEST_PROG_ANALYZER='valgrind --log-file=valgrind-%q{TEST_FULL_NAME}.log'
    export TEST_EXENAME="$EXENAME"
    export TEST_FULL_NAME="$TEST" 
}

# initRef: Initialize variables to refer to another test. (used in rare cases)
# Usage: initRef TESTNAME VAR_SUFFIX
# The exact naming of files for tests is subject to change in future
# updates to the testing system, so for best future compatibility,
# one should not refer to files directly (e.g. Fold_without_options.ct)
# but rather with variables $OUT.ct or placeholders @OUT.ct
# A problem arises when trying to reference files related to a different test.
# initRef will assign the other test's values to variables with the standard 
# names appended with the given suffix. (e.g. TEST_MAIN or OK_1 etc)
# Example:    initRef 'without_options' <SUFFIX>
#    Initializes the following variables:
#      TEST<SUFFIX>,  OUT<SUFFIX>,  OKBASE<SUFFIX>,  OK<SUFFIX>,  STDO<SUFFIX>,  etc.
#    which correspond to the test named 'without_options'
function initRef() {
    local TN="$1" SX="$2"
    [[ $TN ]] || warnmsg 'initRef called with an empty TESTNAME.'
    [[ $SX ]] || warnmsg 'initRef called with an empty SUFFIX.'
    local T="${EXENAME}_$TN" R="${EXEBASE}_$TN"
    export -n   \
       TESTNAME$SX="$TN"        \
       TEST$SX="$T"             \
       TBASE$SX="$R"            \
       OUT$SX="${T}_output"     \
       OKBASE$SX="$OKDIR/$R"    \
       OKROOT$SX="$OKDIR/${EXEBASE}"    \
       OK$SX="$OKDIR/${R}_OK"   \
       OUTFILE$SX="$T$EXT"      \
       OKFILE$SX="$OKDIR/${R}_OK$EXT"   \
       STDO$SX="${T}_stdout.txt" \
       ERRO$SX="${T}_errors.txt"
}

function resetTestVars() {
    # Store previous values
    # These can be used if a test needs to refer to the previous test's output etc.
    PREV_TESTNAME="$TESTNAME"
    PREV_TEST="$TEST"
    PREV_OUT="$OUT"
    PREV_TBASE="$TBASE"
    PREV_OUTFILE="$OUTFILE"
    PREV_OKBASE="$OKBASE"
    PREV_OKROOT="$OKROOT"
    PREV_OK="$OK"
    PREV_OKFILE="$OKFILE"
    PREV_STDO="$STDO"
    PREV_ERRO="$ERRO"
    
    # Just in case runFullTest is made with e.g. $OUT instead of @OUT, the result will be the same.
    TESTNAME=@TESTNAME
    TEST=@TEST
    OUT=@OUT
    TBASE=@TBASE
    OUTFILE=@OUTFILE
    OKBASE=@OKBASE
    OKROOT=@OKROOT
    OK=@OK
    OKFILE=@OKFILE
    STDO=@STDO
    ERRO=@ERRO
    CMDLOG=

    TEST_MISC_OUTPUT_FILES=()
}

# Shows a warning if the previous test was not ended (by calling endTest)
function verifyEndTest() {
    [[ -z $TEST || $TEST == '@TEST' ]] && return 0
    
    warnmsg "Test Authoring Error: endTest was not called after test: '$TEST'."
    endTest
    return 0
}

################################################################################
# Runs checkErrorFiles, passing the in the name of the current test as well as
# the expected error and diff output files.
# Copies test-related files to ${EXENAME}_OUTPUT if the test fails.
# DELETES any test-related files UNLESS the -k or --keep flag is passed .
# Additional arguments can be passed, which specify other files to test for errors.
# Usage: endTest  [ADDITIONAL_FILES...]
################################################################################
function endTest() {
    shouldIncludeTest || return

    local keep hadErr
    [[ -z $TEST || $TEST == @TEST ]] && { warnmsg 'Called "endTest" outside of a test.'; return 1; }
    
    [[ $1 == "-k" || $1 == "--keep" ]] && { keep=1; shift; }
    fixParams "$@"; set -- "${FIXED_PARAMS[@]}" # performs parameter replacements (e.g. @OUT etc)
    
    # Check for the presence of files like ${TEST}_*_errors.txt ${TEST}_*_diff.txt
    checkErrorFiles "$TEST" "$TEST"*_{errors,diff}.txt "$@"  ||  hadErr=1 

    # If any of the following variables are empty, it could cause e.g. `rm -f $VAR*` 
    #  to delete ALL files in the directory, which would be bad. So set them
    #  to a safe value before continuing.
    : ${TBASE:=@TBASE} ${APP_PREFIX:=@APP}

    # If there was an error or if the KeepAll test option is set, then
    # Copy all test files to error output folder (e.g. "${EXENAME}_OUTPUT")
    if ((hadErr)) || isTestOpt 'KeepAll'; then
        copyTestFiles "${ERROR_OUTPUT_DIR:-${EXENAME}_OUTPUT}"
    fi

    if ((!keep)); then
        rm -f ./{"$TEST","$TBASE"}*  # remove all test-related files   (perhaps use ./{"$TEST","$TBASE"}{.,_}* to narrow results)
        #(Warning: bash expands alternates before variables, so $TEST{_,.} is NOT the same as ${TEST}{_,.}. The former becomes $TEST_ $TEST. which is NOT what we want here.)
        # If any MISC output files were listed, delete them also.
        ((${#TEST_MISC_OUTPUT_FILES[@]})) && rm -f "$dest" "${TEST_MISC_OUTPUT_FILES[@]}"
    fi

    local dt=$(date +'%Y-%m-%d %H:%M:%S') TAB=$'\t'
    ((hadErr)) && {
        echo "$dt$TAB$TEST$TAB[$DETECTED_TEST_ERRORS]" >> "!FAILED_TESTS.log"
        echo "$TEST$TAB[$DETECTED_TEST_ERRORS]" >> "$TEST_BLOCK_ERR_LOG"
        true # in case above echo failed.
    } || echo "$dt$TAB$TEST" >> "!PASSED_TESTS.log"

    #echo "  Test $EXENAME '$TESTNAME' finished."
    resetTestVars

    if isTestOpt 'ExitOnError' && ((hadErr)); then
          #exit if an error is detected and the EXIT_ON_ERROR flag is set.
          echo 'Exiting test script due to error. (ExitOnError option is set)'
          isTestOpt 'CleanExit'; exit  # return code will be 0 if CleanExit is set or 1 if not set.
    fi
}

################################################################################
# Copies all existing test-related files to the specified 
# directory, creating it if necessary.
# This includes stdout, stderr and diff output and any other 
# misc files of the form EXENAME_TESTNAME_* or EXENAME_TESTNAME.* 
# Also included are files related to the base/reference 
# program, EXEBASE: EXEBASE_TESTNAME_* or EXEBASE_TESTNAME.* 
# So both Fold-smp_TESTNAME_errors.txt and Fold_TESTNAME.conf 
# would be included for the Fold-smp program.
# Usage:  copyTestFiles  TARGET_DIRECTORY
################################################################################
function copyTestFiles() {
    local dest="$1"
    # If any of these are empty, it could cause $VAR* to match ALL files in the directory, which would be bad.
    [[ $TEST && $TBASE && $APP_PREFIX ]] || { errmsg "In $FUNCNAME: Missing some test variables." && return 1; }
    filterFiles -0 {"$TEST","$TBASE"}*_{errors,diff}.txt && rm -f "${#_FILTERED_FILES[@]}"  # Delete any error or diff files that exist but are empty.
    copyToDir -e "$dest" "$TEST"*   # Copy existing files like {EXENAME}_{TESTNAME}* to the target directory
    [[ $TBASE != $TEST ]] && copyToDir -e "$dest" "$TBASE"*   # If BASENAME is different from EXENAME, copy {BASENAME}_{TESTNAME}* also 
    copyToDir -s "$dest" "$APP_PREFIX"*   # Copy any non-empty app-related files "[{EXENAME}]"* also, like the command-log etc.
    # If any MISC output files were listed, copy them also.
    ((${#TEST_MISC_OUTPUT_FILES[@]})) && copyToDir -e "$dest" "${TEST_MISC_OUTPUT_FILES[@]}"
}

################################################################################
# List additional expected output files from tests.
#   Usage:  addTestOuput <FILES...>
# This does not perform any diffs or other test operations on the files.
# It only registers the files as belonging to the current test, so that
# when the test ends, the files are either deleted or copied to the OUTPUT 
# folder (depending on the outcome of the test and any user-specified test 
# options).
################################################################################
function addTestOuput() {
    TEST_MISC_OUTPUT_FILES+=("$@")
}

################################################################################
# Runs the test executable with whatever arguments are passed to this function. 
# Redirects output to  ${TEST}_stdout.txt and errors to ${TEST}_errors.txt
# Usage:   runTest [--name=SUBTEST]  EXE  [ARGUMENTS...]   
# Notes:
# - If the executable does not contain a path and it is not a known command
#   (as determined using `type -p COMMAND`) but $EXEDIR/COMMAND  *is*  a valid 
#   executable file, then that is used as the executable.
################################################################################
function runTest() {
    shouldIncludeTest || return

    local args=() subname i=0
    local out="$STDO" # i.e. ${TEST}_stdout.txt  (set in initTest)
    local err="$ERRO" # i.e. ${TEST}_errors.txt  (set in initTest)
    for arg; do
        ((i++))
        case "$arg" in 
            -n=*)     subname="${arg:3}"  ;;
            --name=*) subname="${arg:7}"  ;;
             *)       args+=( "${@:i}" ); break;; # add all the rest of the args and exit the loop
        esac
    done
    if [[ $subname ]]; then
      #replace the TEST name with TEST_subname in stdout and stderr file names
      out="${STDO/#$TEST/${TEST}_$subname}"
      err="${ERRO/#$TEST/${TEST}_$subname}"
    fi
    fixParams "${args[@]}"
    [[ -s $out ]] && printf '\n--- TEST %s ---\n' "$(date)" >> "$out"
    [[ -s $err ]] && {
        warnmsg 'Warning: The error output file already exists: ' "$err"
        printf '\n--- TEST %s ---\n' "$(date)" >> "$err"
    }

    # TEST_PROG_ANALYZER -- an environment variable set by the user e.g. to enable analysis by valgrind etc. Empty by default.
    verifyExe "${FIXED_PARAMS[0]}" && FIXED_PARAMS[0]="$VERIFIED_EXE"

    local harness # used for valgrind, gdb etc
    if [[ $TEST_PROG_ANALYZER ]]; then
        analyzer=$(replaceVarsInText "$TEST_PROG_ANALYZER")
        warnmsg "Using analyzer: $analyzer (stderr redirected to $err)"
    fi

    logCmd 'runTest' "$subname" $TEST_PROG_ANALYZER  "${FIXED_PARAMS[@]}"  '1>>' "$out" '2>>' "$err"  # log this to $CMDLOG


    set -o pipefail # causes the exit code of the pipeline to be the first non-zero exit code 
                    # (default is the exit code of the LAST command in the pipeline -- i.e. "tee")
    # Trap signals sent to the child process (e.g. SIGSEGV, SIGINT, SIGSEGV)
	# Note that CHLD now seems to trap even error-free exits.

    trap FAILED_runTest CHLD INT SEGV

    if isTestOpt 'ShowOutput'; then
        # Do bash redirection to capture stderr and stdout to two separate files but still output them to the screen.
        { { $analyzer "${FIXED_PARAMS[@]}"; } 2>&3 | tee -a "$out"; } 3>&1 1>&2 | tee -a "$err"
    else
        $analyzer "${FIXED_PARAMS[@]}" 1>>"$out" 2>>"$err"
    fi

    local exitcode=$? #save the exit code from the program so the following trap command doesn't alter it.
    trap - CHLD INT # remove the trap
    return $exitcode
}

# Called when the program launched by runTest encounters a segfault or other critical error.
function FAILED_runTest() {
     local code=$?
	 #exit this if the return code was 0
	 if [[ "$code" -eq 0 ]]; then
        return 0
     fi
	 #also exit if return code is 1; RNAstructure returns 1 when input was incorrect.
	 if [[ "$code" -eq 1 ]]; then
        return 0
     fi
     case $code in 
        130) errmsg "User Canceled Command (Ctrl+C): ${FIXED_PARAMS[@]}" 
             trap - INT; kill -INT $$ ;;  # the inner program was exited with Ctrl+C, so kill this script with the same signal. (see http://mywiki.wooledge.org/SignalTrap)
        139) errmsg "Segmentation Fault: ${FIXED_PARAMS[@]}" ;;
        126) errmsg "Error 126 -- Cannot Execute Command (possible permission problem or command is not an executable): ${FIXED_PARAMS[@]}" ;;
        127) errmsg "Error 127 -- Command Not Found (missing executable or possible problem with \$PATH or a typo): ${FIXED_PARAMS[@]}" ;;
        *)   errmsg "Error $code -- Fatal Error Signal: ${FIXED_PARAMS[@]}" ;;
    esac
}

################################################################################
# Runs runTest followed by runDiff
################################################################################
function runProgAndDiff() {
    shouldIncludeTest || return
    local arg val args=()  diff1=@OUTFILE  diff2=@OKFILE 
    local requireExit=0 code useStdErr replaceInOutput skipDiff osSpecific
    local diff_files=() # possible to add multiple diff-file-sets

    for arg; do
        [[ $arg == -*=* ]] && val="${arg#*=}" || val= # gets the part of arg following "=", just in case it is in the form "-name=val"
        case "$arg" in 
            ---stdout) diff1=@STDO ;;
            ---stderr) useStdErr=1 ;;
            #---ref=*) diff2=$(initRef "$val" _REF; echo "$OKFILE_REF") ;;
            ---ref=*) setTestRef "$val" ;;
            ---ext=*) diff1="@OUT$val"; diff2="@OK$val" ;;
            ---out=*) diff1="$val" ;;
            ---ok=* ) diff2="$val" ;;
            ---exit=*) requireExit="$val"  ;;
            ---replace=*) replaceInOutput="$val" ;;
            ---nodiff)  skipDiff=1 ;;
            ---add-diff=*) diff_files+=( "$val" ) ;; # add another set of files to DIFF
            ---OS ) osSpecific=1 ;;
            *)  args+=("$arg") ;;
        esac
    done
    # The following command is usually runTest, but can be modified by the 
    # calling script (e.g. in GUI tests, it is changed from runTest to runGuiTest)
    # RUNTEST_CMD -- usually set by scripts to change the way tests are run. E.g. running runTestGui instead of runTest etc.
    ${RUNTEST_CMD:-runTest} @EXE "${args[@]}"
    code=$?

    if [[ $useStdErr ]]; then
        # replace the output file with the stderr output
        cp $ERRO $OUTFILE # replace the output file with the stderr output.
        echo -n > $ERRO; # truncate the error file
    fi

    #TODO: see if we can require the exit code to always be 0 (unless specified otherwise) Do most programs have 0 exit codes by default?
    # By default, programs must exit with 0 as their exit-code.
    # --exit=ANY  means ignore the exit code (e.g. the exit code is either uninformative or is allowed to be 0 or non-zero)
    # --exit=ERR  means that a non-zero exit code is required. (In which case, 0 is considered an error)
    # --exit=<NUM>  (where <NUM> is a number 0 to 255) means that a specific exit code is required (any other number is considered an error)
    case "$requireExit" in
        ign|IGN|any|ANY|'') ;; # ignore the exit code
        err*|ERR*)  # any non-zero exit code
            if [[ $code -eq 0 ]]; then
                errmsg "Expected an Error-Code from ${FIXED_PARAMS[0]} but the exit-code was 0."
                return 1
            fi ;;
        *)  # require a specific exit code
            if [[ "$requireExit" -ne $code ]]; then
                errmsg "Unexpected Exit-Code from ${FIXED_PARAMS[0]}: $code (expected $requireExit)."
                return 1
            fi ;;
    esac
    # The exit code is OK, so continue.

    # Add the main diff pair (unless ---nodiff is given). 
    #Also add any other diff-file-sets specifically added with --add-diff
    [[ ! $skipDiff ]] && diff_files=( ":$diff1==$diff2"   "${diff_files[@]}"  )

    local diffName outName refName
    for diff_pair in "${diff_files[@]}"; do
        # Each diff_pair is in the form <NAME>:<OUT>==<OK>
        # Both <NAME> and <OK> are optional (along with the ':' and '==' separators)
        if [[ "$diff_pair" == *:* ]]; then
            diffName="${diff_pair%%:*}"
            diff_pair="${diff_pair#*:}"
        else
            diffName=
        fi
        if [[ "$diff_pair" == *'=='* ]]; then
            outName="${diff_pair%%==*}" # part before ==
            refName="${diff_pair#*==}"  # part after ==
        else
            outName="$diff_pair"
            refName= #set later
        fi
        if [[ ! $refName ]]; then
            refName=@OKFILE
            # append diffName if non-empty (this uses the special @VAR{name} syntax -- see replaceVarsInText)
            [[ $diffName ]] && refName+="{$diffName}"
        fi

        local ppexit=0
        local outfile=$(replaceVarsInText "$outName")

        # Perform text-replacement post-processing (from --replace)
        if [[ $replaceInOutput ]]; then
            replaceInFiles "${replaceInOutput%=>*}" "${replaceInOutput#*=>}" "$outfile" || 
                errmsg "Bad Post-Processing (Text Replacement) Exit-Code: $?"
        fi
        # Perform additional custom post-processing if the TEST_POST_PROCESS function is defined.
        if isFunc TEST_POST_PROCESS; then
            TEST_POST_PROCESS "$outfile" || ppexit=$?
        elif isVar TEST_POST_PROCESS; then
            $TEST_POST_PROCESS "$outfile"  || ppexit=$?
        fi
        ((ppexit)) && errmsg "Bad Post-Processing Exit-Code: $ppexit"

        # Run the actual diff.
        # runDiff <OUTPUT>  <OKFILE> <DIFF_NAME>  (note that this calls `diff <OKFILE> <OUTPUT>`` -- order is swapped)
        runDiff  "$outName" "$refName" "$diffName"
    done
}

# Change the base name for "OK" file references
function setTestRef() {
    OKBASE="${OKDIR}/${EXEBASE}_$1"   # The "OK" directory and base file name (before the OK prefix). Use this when there are multiple OK files for a single test, eg. @OKBASE-seq1_OK.ps, @OKBASE-seq2_OK.ps etc.
    OKROOT="${OKDIR}/${EXEBASE}"      # The "OK" directory and EXEBASE (with no TESTNAME). Use this when you need a prefix that can apply to more than one test.
    OK="${OKBASE}_OK"             # The "OK" file name without extension, e.g. Fold/Fold_without_options_OK
    OKFILE="$OK$EXT"              # The "OK" file name WITH extension, e.g. Fold/Fold_without_options_OK.ct
}

# Logs a command executed by the testing system to the file "$CMDLOG"
# Usage: logCmd CALLING_FUNCTION  SUB_TEST_NAME  [COMMAND_AND_PARAMETERS...]
function logCmd() {
    local out cmd="$1" sub="$2" arg
    if [[ $CMDLOG == ${TEST}_* ]]; then 
        out="$CMDLOG" 
        [[ -s $out ]] || printf "# Commands run during the test '$TEST':\n" > "$out"
    else
        out="$APP_CMDLOG"
        [[ -s $out ]] || printf "# Commands run when testing $EXENAME (NOT part of any specific test):\n" > "$out"
    fi

    [[ -z $sub ]] && sub="(MAIN)"
    [[ $cmd == runCmd ]] && cmd="CMD"
    [[ $cmd == runTest ]] && cmd="RUN"
    [[ $cmd == runDiff ]] && cmd="DIFF"
    local dt=$(date +'%F %T.%N')
    printf "#%s\t%s::%s:\n " "$dt" "$cmd" "$sub" >> "$out"
    
    for arg in "${@:3}"; do
        case "$arg" in
            '>'|'>>'|[12]">"|[12]">>"|'>&'|'&>') printf ' %s' "$arg" >> "$out" ;;
            *) printf ' %q' "$arg" >> "$out" ;;
        esac
    done
    printf '\n' >> "$out"
}

# Replace parameters that start with a special character (@)
function fixParams() {
    # e.g. @TEST_sav.txt  ==>  ${TEST}_sav.txt ==> Fold_dna_option_sav.txt
    local fixParams_ARG_RESULT
    FIXED_PARAMS=()
    for arg; do
        if [[ $arg == *@* ]]; then
            replaceVarsInText "$arg" fixParams_ARG_RESULT
            FIXED_PARAMS+=("$fixParams_ARG_RESULT")
        else
            FIXED_PARAMS+=("$arg")
        fi
    done
}

# Replace all occurances of @OUT @TEST etc in the given text.
function replaceVarsInText() {
  local text="$1"
  if [[ $text == *@* ]]; then
    local var rep # var is a variable name from TEST_PARAMNAME_LIST.  rep is the replacement text.
    local ERRDIR="${ERROR_OUTPUT_DIR:-${EXENAME}_OUTPUT}" # this variable must be set dynamically and only remain set for this scope.
    for var in "${TEST_PARAMNAME_LIST[@]}"; do
        [[ $text == *@$var* ]] || continue # if @VAR is not present, skip to next var.

        # First look for @VAR{TEXT}  (i.e. custom text in literal braces)
        while [[ $text =~ @$var'{'([^}]*)'}' ]]; do
            local match="${BASH_REMATCH[1]}" rt # rt is replacement-text
            # For some vars, this format has special meaning.
            # e.g. $OUT{name} ==> ${TEST}_name_output
            case $var in
                OUT)     rt="${TEST}_${match}_output" ;;
                OUTFILE) rt="${TEST}_${match}_output$EXT" ;;
                OK)      rt="${OKBASE}_${match}_OK" ;;
                OKFILE)  rt="${OKBASE}_${match}_OK$EXT" ;;
                STDO)    rt="${TEST}_${match}_stdout.txt" ;;
                ERRO)    rt="${TEST}_${match}_errors.txt" ;;
                *)       rt="${!var}_$match";;  # default replacement is simply $VAR_<TEXT> (same as @VAR_<TEXT>)
            esac
            local ft="@$var{$match}"
            #echo -e "\nreplacing $ft with $rt\n"
            text="${text//"$ft"/$rt}";
        done
        # Replace the simple @VAR form (i.e. without braces)
        text="${text//@"$var"/${!var}}"
    done
    text="${text//@AT/@}" #finally replace @AT with literal '@'
  fi
  _echoOrSet "$text" "$2" || true # will either echo the result or set a variable.
}

###############################################################################
# Usage: _echoOrSet <VALUE> [<VARNAME>]
#  If VARNAME is present (and non-empty) then the variable with that name is
#  set to VALUE. Otherwise VALUE is echoed to STDOUT.
# (Note that this function is simpler than `echoOrSetVar` in script-utils.sh)
###############################################################################
function _echoOrSet() {
  if [[ "$2" ]]; then
    setVar "$2" "$1"
  else
    echo "$1"
  fi
}

################################################################################
# runSubTest - Runs an executable (usually different than the main TEST one) 
#              with whatever arguments are passed to this function. 
#              Redirects output to  ${TEST}_${SUBTEST}_stdout.txt 
#                    and errors to ${TEST}_${SUBTEST}_errors.txt
# The first argument is SUBTEST - a descriptive name for what the subtest does.
# The second argument is PATH_TO_EXE.
#   if PATH_TO_EXE contains no path information, $EXEDIR (e.g. ../exe) 
#   is prepended.
# Usage:   runSubTest 'SUBTEST' PATH_TO_EXE  [ARGUMENTS...]
#           (where PATH_TO_EXE is assumed to be $EXEPATH if not specified).
################################################################################
function runSubTest() {
    runTest --name="$1"  "$2"  "${@:3}"
}


################################################################################
#   Runs the given command, but redirects output to <COMMAND>_stdout.txt and 
# errors to <COMMAND>_errors.txt, where <COMMAND> is the name of the 
# executable or command that is being run (unless the -n=NAME flag is present).
#
# Usage:    runCmd [--name=NAME] [--out=FILE] [--required]  [ARGUMENTS...]   
# Flags:   
#           -n | --name=NAME  - Use NAME instead of COMMAND to derive the 
#                               names of the output and error files.
#           -o | --out=FILE   - Redirect stdout to FILE instead of <COMMAND>_stdout.txt
#                               names of the output and error files.
#           -r | --required   - Exit the script if the command fails.
# Notes:
# - If the executable does not contain a path and it is not a known command
#   (as determined using `type -p COMMAND`) but $EXEDIR/COMMAND  *is*  a valid 
#   executable file, then that is used as the executable.
# - The error file, <COMMAND>_errors.txt will be copied into $ERROR_OUTPUT_DIR.
# - The exit code of this function will be the same as that of the underlying 
#   command.
# - This function is generally used for commands that are run in preparation
#   for tests, not during individual tests themselves. Therefore it only relies
#   on top-level variables (EXENAME, ERROR_OUTPUT_DIR) to be set, but not any 
#   test-level variables (e.g. TEST, TESTNAME, OUT, OK, OKFILE etc)
################################################################################
function runCmd() {
    local args=() name out err i=0 required=0 ret KeepAll
    for arg; do
        ((i++))
        case "$arg" in 
            -n=*|--name=*)      name="${arg#*=}"  ;;
            -o=*|--out=*)       out="${arg#*=}"; KeepAll=1   ;;
            -o|--out)           KeepAll=1   ;;
            -r|--required)      required=1        ;;
             *)       args+=( "${@:i}" ); break;; # add all the rest of the args and exit the loop
        esac
    done
    [[ -z $name ]] && name="${args[0]##*[\\/]}"
    [[ -z $out ]]  && out="[${EXENAME}]_${name}_stdout.txt"
    err="[${EXENAME}]_${name}_errors.txt"
    
    #fixParams "${args[@]}"
    [[ ! $KeepAll ]] && printf '\n#--- OUTPUT FROM %s RUN %s ---\n' "$name" "$(date)" >> "$out"
    [[ -s $err ]] && {
        warnmsg 'Warning: The error output file already exists: ' "$err"
        printf '\n--- RUN %s ---\n' $(date) >> "$err"
    }
    verifyExe "${args[0]}" && args[0]="$VERIFIED_EXE"
    
    logCmd runCmd "$name" "${args[@]}" '1>>' "$out" '2>>' "$err"
    
    #*** RUN the command  ***
    "${args[@]}" 1>>"$out" 2>>"$err"
    ret=$? # store exit code
    # copy error file to ERROR_OUTPUT_DIR
    [[ -s $err ]] && {
        copyToDir "${ERROR_OUTPUT_DIR:-${EXENAME}_OUTPUT}"  "$err" "$out" "$APP_CMDLOG" # copyToDir (defined in script-utils) is the same as `cp -t`, but creates the directory if it does not exist.
        (( ret==0 && ret++ )) # set ret to 1 (error) if it is currently 0
    }
    # remove error and output files (they have already been copied if an error occurred.)
    rm -f "$err" #; ((KeepAll)) || rm -f "$out"

    # Output error message (if appropriate)
    ((ret)) && {
        errmsg "  ! The '$name' command failed."
        ((required)) && die 'Required command failed.'
    }
    return $ret
}

function verifyExe() {
    # isExe (defined in script-utils) returns 0 if $1 is an executable file, script, function, or command (in the PATH)
    # Prioritize RNAstructure programs (in case there is a command or function with the same name)
    if [[ $1 != *[\\/]* && -x "$EXEDIR/$1" ]]; then
        VERIFIED_EXE="$EXEDIR/$1"
    elif isExe "$1"; then
        VERIFIED_EXE="$1"
    else
        return 1
    fi
}

################################################################################
# One-liner to run a test, perform diff, and check for errors.
# Usage:   runFullTest  TESTNAME  [TEST_ARGS...]
# Example: runFullTest  'without_options' $SINGLESEQ2 @OUT.ct
# See descriptions of initTest, runTest, runDiff, and endTest
################################################################################
function runFullTest() {
    shouldIncludeTest "${EXENAME}_$1" || return  # Check whether this test is included ( depends on INCLUDE_TESTS and EXCLUDE_TESTS environment variables)
    initTest "$1"
    runProgAndDiff "${@:2}" $APPEND_TEST_ARGS
    endTest
}

################################################################################
# Performs a diff on the output of a test and redirects the diff output to file.
# Usage runDiff  [FILE1] [FILE2] [DIFF_NAME]   
#  FILE1 -- First file to compare.  Usually the actual output file from a test.
#  FILE2 -- Second file to compare. Usually the "reference" or "expected output" or "original" file.
#  DIFF_NAME -- Unique name for diff-output and for logging. See notes.
#            Default: ""
#
#  All parameters are optional. Default values are:
#    FILE1=${TEST}_output$EXT  FILE2=$OKDIR/${TEST}_OK$EXT  DIFF_NAME=""
#
#  Important: runDiff swaps the file order when calling `diff` 
#    i.e. diff <FILE2> <FILE2> 
#    This is because the standard order of arguments to diff is 
#    `diff <REFERENCE> <ACTUAL>`   or  `diff <ORIGINAL> <NEW>`
#
# DIFF_NAME is a unique name for this diff operation for logging. It is also
#   used to generate the name of the diff-output file.
#   In most cases, DIFF_NAME should be empty. However if runDiff is called 
#   multiple times in one TEST, (for example if the test produces multiple 
#   output files) then DIFF_NAME must be unique for each product, so that 
#   successive diff operations do not overwrite previous ones.

#   (*) If DIFF_NAME is empty, the auto-generated diff-output file name will be
#       ${TEST}_diff.txt.
#   (*) If DIFF_NAME ends in "_diff.txt" then it is a literal file name and will 
#       be used as such without modification. The name for logging will be the
#       text before "_diff.txt", but it is stripped of test-related prefixes.
#   (*) Otherwise DIFF_NAME is assumed to be a "simple" name and is used for
#       logging (without modification). The diff-output file name will be
#       auto-generated as ${TEST}_<DIFF_NAME>_diff.txt
#       This is the usual form when multiple test-output files need to be 
#       verified during a single test.
################################################################################
function runDiff() {
    shouldIncludeTest || return

    fixParams "$@"; set -- "${FIXED_PARAMS[@]}"
    
    local fOut=${1:-$OUT$EXT} fRef=${2:-$OK$EXT} fDiff logname
    # If an OS-specific reference file exists, then use it instead of the standard one.
    #   OS-specific reference file names are the same as standard OK file names except 
    #   that _OK is replaced with _@{OS}_OK
    #   For example: ShapeKnots_window_size_option_@Mac_OK.ct  (the @Mac portion specifies that it is specific to Mac OSX)
    #   Values for OS are obtained from the getOS function in script-utils.sh and are the following:
    #       Windows, Linux, Mac, BSD, and Solaris  (abbreviations and changes in letter case will NOT be accepted!)
    local f2_OS="${fRef/_OK/_@${OS}_OK}"; [[ -f $f2_OS ]] && fRef="$f2_OS"

    # Test for the third parameter: DIFF_NAME
    if [[ "$3" == *_diff.txt ]]; then
        fDiff="$3" # DIFF_NAME ends with the proper suffix. Use it directly without modification.
        logname=$(strip-test-prefix "${3%_diff.txt}")
    elif [[ $3 ]]; then
        # DIFF_NAME is a simple name, not a full file name. Decorate it with prefix and suffix.
        fDiff="${TEST}_$3_diff.txt"
        logname="$3"
    else # DIFF_NAME is either not specified or is empty. Use default.
        fDiff="${TEST}_diff.txt"
    fi

    #  Determine if the diff output file already exists. If so, append a number to generate a unique 
    #  non-existing file name for diff output.
    if [[ -f $fDiff ]]; then
        warnmsg 'The diff output file already exists: '"$fDiff"
        local pfx="${fDiff%_diff.txt}"  i=2
        while [[ -f $fDiff ]]; do fDiff=$pfx$((i++))_diff.txt; done
    fi

    # Log the diff command regardless of whether or not we run it. (This helps with test maintenance, e.g. accept-test-output-OK.sh)
    logCmd runDiff "$logname" "$DIFF_CMD" $DIFF_FLAGS "$fRef" "$fOut"  '>&' "$fDiff" 

    if [[ ! -f $fOut ]]; then
        errmsg "      ERROR: Expected output file missing: $fOut"
    elif [[ ! -f $fRef ]]; then
        errmsg "      ERROR: Expected OK-file missing: $fRef"
    else
        "$DIFF_CMD" $DIFF_FLAGS "$fRef" "$fOut" >& "$fDiff" 
    fi
}

# Use formatting sequences (if running in a terminal).
if [[ $TERM && -t 1 ]]; then
    FMT_BLD=$(tput bold)  # Bold (Strong)
    FMT_RED=$(tput setaf 1) # Red
    FMT_GRN=$(tput setaf 2) # Green
    FMT_YLW=$(tput setaf 3) # Yellow
    FMT_BLU=$(tput setaf 4) # Blue
    FMT_PRP=$(tput setaf 5) # Purple
    FMT_CYN=$(tput setaf 6) # Cyan
    FMT_WHT=$(tput setaf 7) # White
    FMT_RST=$(tput sgr0)    # Reset
fi

# Remove test-variables from the beginning of a word
function strip-test-prefix() {
    # We should remove LONGER values before shorter ones.
    # TEST="${EXENAME}_$TESTNAME"
    # TBASE="${EXEBASE}_$TESTNAME"
    # OUT="${TEST}_output"
    # OUTFILE="$OUT$EXT"
    # OKBASE="${OKDIR}/$TBASE"
    # OKROOT="${OKDIR}/$EXEBASE"
    # OK="${OKBASE}_OK"
    # OKFILE="$OK$EXT"
    # STDO="${TEST}_stdout.txt"
    # ERRO="${TEST}_errors.txt"
    text="$1"
    for p in ERRO STDO OKFILE OK OKROOT OKBASE OUTFILE OUT TBASE TEST TESTNAME EXENAME EXEBASE OKDIR; do
        val="${!p}"
        if [[ $val && "$text" == "$val"* ]]; then
            text="${text#"$val"}"
        fi
    done
    echo "$text"
}

# Echo a message to stderr. Also echo to file $ERRO (if is empty).
function errmsg() { 
    local fmt end
    if [[ -t 2 ]]; then # if outputing to a terminal, use color.
        fmt="$FMT_BLD$FMT_RED"
        end="$FMT_RST"
    fi
    echo >&2 "$fmt""$@""$end"
    # IF this is an error during a test, we want to make sure
    #  that the test is marked as having an error. So
    #  also output to $ERRO if it doesn't yet exist.
    if [[ $ERRO && $ERRO != @ERRO ]]; then
       [[ -s $ERRO ]] || echo >"$ERRO" "$@"
    fi
}

# Show success message
function okmsg() { 
    if [[ -t 1 ]]; then
        echo "$FMT_BLD$FMT_GRN""$@""$FMT_RST"
    else
        echo "$@"
    fi
}

# Show test-authoring warning message (NOT output to stderr)
function warnmsg() { 
    isTestOpt W && return # return if the NoWarnings test option is enabled.
    if [[ -t 1 ]]; then
        echo "$FMT_BLD$FMT_YLW  ## WARNING: ""$@""$FMT_RST"
    else
        echo "  ## WARNING: $@"
    fi
}

# output test status message
function status() {
    isQuiet && return 
    if [[ $1 == -0 ]]; then
        echo "${*:2}"
    else
        echo "  $*"
    fi
}

################################################################################
# Runs `make` on the main executable as well as each additional argument.
# Make is executed in the parent directory and is run individually for each
#   listed executable.
# Usage: runMake [-s] [-r] [OTHER_PROGRAMS...]
#   -s | --silent   - Suppress messages
#   -r | --required - Call `exit 1`  if any of the make commands fail.
#   -f | --force    - Force calling `make PROG` even if PROG already exists,
#                     to pick up any recent code changes. (The default is to 
#                     use existing programs without rebuilding).
# Example:
#  `runMake Fold EnergyPlot`  is equivalent to:
#     pushd ..
#     make Fold         >& /dev/null  || errmsg "Failed to make Fold"  
#     make EnergyPlot   >& /dev/null  || errmsg "Failed to make EnergyPlot"  
#     popd
################################################################################
function runMake() {
    local silent=0   required=0  errors=0 force=
    [[ $TEST_SKIP_REMAKE ]] || force=1  #force rebuilding unless TEST_SKIP_REMAKE is set.
    targets=()
    for arg; do
        case "$arg" in 
            '') ;; # skip empty
            -s|--silent)   silent=1      ;;
            -r|--required) required=1    ;;
            -f|--force)    force=1       ;;
             @EXE) targets+=("$EXENAME") ;;
             *) targets+=("$arg")        ;;
        esac
    done

    ((silent)) || status  "Building executables for ${EXENAME} tests..." # "  Building executables: ${targets[@]} ..."
    pushd "${MAKEROOT}" > /dev/null
    local target makelog="tests/${APP_PREFIX}_MAKE.log"
    for target in "${targets[@]}"; do
        if [[ ! $force && -x "$OLDPWD/$EXEDIR/$target" ]]; then
            ((silent)) || status "  Skip making $target (already exists)."
            continue
        fi
        ((silent)) || status "  Making $target..."
        if ! make "$target" $TEST_PROG_MAKE_ARGS  >& "$makelog"; then
            ((errors++))
            errmsg "  ! Failed to make '$target'"
            isTestOpt 'DumpErrors' && tail -n 4 "$makelog" | dumpErrorFile
            ((required)) && die 'Building executables failed!'
        fi
    done
    popd > /dev/null
    ((silent)) || status "Done Building executables ($errors errors)."
    ((errors==0)) # return value is OK (0) only if there were no errors.
}

function die() {
    errmsg "  ! Aborting tests for $EXENAME: $*"$'\n\n'
    exit 1
}


################################################################################
# Verifies that each of its arguments is an existing file and exits the script 
# with an error message if the file is missing.
# If the -t (--test) flag is passed, then any missing file causes the script to 
#   output a warning, but does not cause the script to exit.
# Returns: 0 if all files exist or 1 if one or more files are missing.
# Usage: verifyInputFiles [-t] [FILES...]
################################################################################
function verifyInputFiles() {
    local file files=() errors=0 action=die # action to take if a file is missing: errmsg or die
    for arg; do
        case "$arg" in 
            '') ;; # skip empty
            -t|--test) action=errmsg   ;;
            -r|--required) action=die  ;; # the default action
             *) files+=("$arg")        ;; 
        esac
    done

    [[ -d $INDIR ]] || INDIR=testFiles
    for file in "${files[@]}"; do
        [[ -f $file ]] && continue
        #[[ -f $INDIR/$file ]] && { cp "$INDIR/$file" $file && continue; }
        # File not found. Try to make it.
        make -f "input-files.mk" "$file" && [[ -f $file ]] && continue
        # If we're still here, the file was not found and could not be made.
        $action "A required input file is missing: $file" #$'\n'"This file was not in $INDIR/ and could not be produced by \`make\`."
        ((errors++))
    done
    ((errors==0)) # return true if there were no errors.
}

################################################################################
## beginTestBlock - Announces start of test suite.
##                  Run this before any tests are run.
## Usage:  beginTestBlock [-n=NAME] [PROGRAMS_TO_BUILD...]
#        -n=NAME will change the name displayed in the message to NAME
##  
################################################################################
function beginTestBlock() {
    # Get a list of existing files. These will be compared at the end of the suite to
    # ensure that all test files are deleted.

    EXE="$EXENAME"
    OS=$(getOS) # used for OS-specific reference files.

    # Note: EXEBASE is usually the same as EXENAME except when testing a
    #   program that is implemented differently, but should have exactly 
    #   the same output as another "reference" implementation, for example Fold-smp,
    #   where EXENAME='Fold-smp' and EXEBASE='Fold'
    : ${EXEBASE:=$EXENAME} ${OKDIR:=$EXEBASE}
    APP_PREFIX="[${EXENAME}]"
    APP_CMDLOG="${APP_PREFIX}_commands.log" # E.g. [Fold-smp]_commands.log

    TEST_BLOCK_NAME="$EXENAME"
    if [[ $1 == '-n='* ]]; then
        TEST_BLOCK_NAME="${1:3}"
        shift
    fi

    if [[ $TESTING_STARTED ]]; then
        status -0 "Continuing tests: '$TEST_BLOCK_NAME'"
    else
        TESTING_STARTED=1  # guard against this being called twice in the same test-run

        status -0 "Beginning '$TEST_BLOCK_NAME' tests..."
        
        # Store the list of files in the directory, so any orphans from tests can be removed.
        filesToArray . -maxdepth 1 -type f ! -name Makefile; #FILES_BEFORE_TESTS=( "${FILES_ARRAY[@]}" )

        if listOrphanedFiles; then
            mkdir -p TRASH
            status 'Cleaning root test dir.'
            mv -f "${_ORPHANED_FILES[@]}" TRASH
        fi

        local errdir="${ERROR_OUTPUT_DIR:-${EXENAME}_OUTPUT}"
        TEST_BLOCK_ERR_LOG=$errdir/FAILED_TESTS.log
        if [[ -d ./$errdir ]]; then
            status 'Cleaning output dir.'
            rm -f ./"$errdir"/*  # clean output directory so old files don't persist
        fi
    fi
    
    resetTestVars
}

################################################################################
## endTestBlock - Announces end of test suite.
##                Run this after all tests are run.
## Usage:  endTestBlock [-n=NAME]
#        -n=NAME will change the name displayed in the message to NAME
################################################################################
function endTestBlock() {
    verifyEndTest
    status -0 "Completed '$TEST_BLOCK_NAME' tests."$'\n'

    local errdir="${ERROR_OUTPUT_DIR:-${EXENAME}_OUTPUT}"
    local app="${APP_PREFIX:-@APP}" # we don't want this to be empty, because rm -f "$APP_PREFIX"* would delete all files.

    # If the error output directory exists, (which could be due to errors or the KeepAll test option)
    # then copy the command log ($APP_CMDLOG) and any stdout/stderr output from commands
    # into the error output directory.
    [[ -d $errdir ]] && copyToDir -s "$errdir" "$app"*  
    rm -f "$app"* # delete all files related to the program being tested.

    if listOrphanedFiles; then
        warnmsg "Test Authoring error -- Orphaned test files: "
        printf '    %s\n' "${_ORPHANED_FILES[@]}"
        local dir="$errdir/~ORPHANED" 
        mkdir -p "$dir"
        mv -f "${_ORPHANED_FILES[@]}" "$dir"
    fi

    #return value should be true (0) if no error occurred. Otherwise it will be false (1) unless the CleanExit option is set.
    [[ ! -f $TEST_BLOCK_ERR_LOG ]] || isTestOpt 'CleanExit'
}

################################################################################
# Fills the array _ORPHANED_FILES with a list of files in the root test directory
#  ( RNAstructure/tests/ ) aside from Makefile and other "protected" files.
# Returns true (0) if there ARE orphaned files and false (1) otherwise.
################################################################################
function listOrphanedFiles() {
    filesToArray . -maxdepth 1 -type f  ! -name Makefile ! -name '!*' ! -name '*.mk'
    _ORPHANED_FILES=( "${FILES_ARRAY[@]}" )
    ((${#_ORPHANED_FILES[@]}>0)) # Returns true (0) if there ARE orphaned files and false (1) otherwise.
}

PREV_TEST_OPTS=-- # set to invalid value
function isTestOpt() {
    local opt abbrev
    shopt -s nocasematch # ignore case for option name comparisons
    ALL_TEST_OPTS='K-KeepAll L-ListErrors D-DumpErrors S-ShowOutput X-ExitOnError C-CleanExit G-DebugGUI Q-Quiet W-NoWarnings'
    if [[ $TEST_OPTS != $PREV_TEST_OPTS ]]; then
        # Setup option variables if the TEST_OPTS have changed.
        #echo 'Setting opts'
        PREV_TEST_OPTS="$TEST_OPTS"
        getTestOptAbbrev() { local opt;abbrev=;for opt in $ALL_TEST_OPTS; do [[ $1 == ${opt:2} || $1 == ${opt:0:1} ]] && { abbrev=${opt:0:1}; return 0; } done; return 1; }
        unset ${!__TEST_OPT__*}
        for opt in $TEST_OPTS; do 
            if getTestOptAbbrev $opt; then 
                #echo "Setting __TEST_OPT__$abbrev"
                export -n  __TEST_OPT__$abbrev=1
            else
                #echo 'Testing unknown option ' $opt
                # see if the entire option is composed of test-opt abbreviations e.g. kdgs
                local i len=${#opt}
                # just check here, don't set any options, in case the whole option is bad.
                for ((i=0; i<len; i++)); do getTestOptAbbrev "${opt:i:1}" || break; done
                if ((i<len)); then
                    errmsg "Unknown test system option: '$opt'"
                    exit 1
                else
                    for ((i=0; i<len;i++)); do
                      getTestOptAbbrev "${opt:i:1}"
                      #echo "Setting __TEST_OPT__$abbrev"
                      export -n  __TEST_OPT__$abbrev=1
                    done
                fi
            fi
        done
    fi
    if [[ $1 == --show ]]; then
        local var allset
        for opt in $ALL_TEST_OPTS; do
            var=__TEST_OPT__${opt:0:1}
            [[ ${!var} ]] && allset+=,${opt:2}
        done
        allset=${allset:1}
        echo "${allset//,/ }"
    else
        getTestOptAbbrev $1 || errmsg "Unknown test system option: '$1'"
        opt=__TEST_OPT__${abbrev:-NULL} # e.g. __TEST_OPT__k or __TEST_OPT__NULL
        shopt -u nocasematch # restore case sensitivity
        #echo "Getting $opt"
        [[ ${!opt} ]] # Return the result of the test
    fi
}
function isQuiet() { isTestOpt Q; }


# This function searches a file for all occurances of numbers having 
# 3-digit exponents with leading zeros and converts them into 
# 2-digit exponents. e.g.:  3.4e-006 ==> 3.4e-06
#
# Usage: fixExponents FILE
#
# Remarks:
# The C-libraries on Linux, Mac, and Windows do not format 
# exponents the same way. For example:
#     std::cout << 1e-7;   // Outputs 1e-007 on Windows/Mac and 1e-07 on Linux
# Note that on linux, the exponent has a minimum of two digits and does
# not print 3 digits unless necessary. On Windows it always prints 3 digits.
function fixExponents() {
    [[ -e $1 ]] && sed -i.bak -E 's/([0-9.]+e[-+])0([0-9]{2})/\1\2/g' "$1"
}

function shouldIncludeTest() {
    Y="$FMT_YLW"
    R="$FMT_RST"
    [[ $SKIP_TESTS ]] && return 1

    if (($#==0)); then
        [[ ! $EXCLUDE_CURRENT_TEST ]]
        return # return the result of the previous test
    fi

    local name="$1" list
    if [[ $EXCLUDE_TESTS ]]; then #|| ${#EXCLUDE_TESTS[@]} > 1 
        EXCLUDE_CURRENT_TEST=1 # assume excluded if exited during the loop
        set -f; list=($EXCLUDE_TESTS); set +f #set -f prevents expanding wildcards
        for expr in "${list[@]}"; do
            if [[ $name == $expr || $name == ${EXENAME}_$expr ]]; then 
                warnmsg "Skipping Excluded Test $name"
                return 1 # the test was excluded
            fi
        done
    fi
    
    EXCLUDE_CURRENT_TEST= # Assume included

    # test was NOT excluded. make sure it was included
    if [[ $INCLUDE_TESTS ]]; then #|| ${#INCLUDE_TESTS[@]} > 1
        set -f; list=($INCLUDE_TESTS); set +f  #set -f prevents expanding wildcards
        for expr in "${list[@]}"; do
            [[ $name == $expr  || $name == ${EXENAME}_$expr ]] && return 0 # the test was explicitly included
        done
        EXCLUDE_CURRENT_TEST=1
        warnmsg "Skipping NON-Included Test $name"
        return 1 # Not explicitly or implicilty included
    else
        return 0 # Implicitly included
    fi
}

function replaceInFiles() {
    fixParams "$@"
    local file find="${FIXED_PARAMS[0]}" replaceWith="${FIXED_PARAMS[1]}" file
    for file in "${FIXED_PARAMS[@]:2}"; do
        [[ -f $file ]] && sed -i.bak -e "s|$find|$replaceWith|g" "$file"
    done
}
