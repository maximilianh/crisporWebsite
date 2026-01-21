#!/bin/bash

################################################################################
#                           test-runner.sh
################################################################################
#  This script prepares the environment for test scripts and executes them.
################################################################################

################################################################################
#                        Script initialization
################################################################################
THIS_SCRIPT="${BASH_SOURCE[0]}"
THIS_SCRIPT_DIR=$(dirname "$THIS_SCRIPT")
export SCRIPT_DIR="$THIS_SCRIPT_DIR"    # path to RNAstructure/tests/scripts

source "$SCRIPT_DIR/script-utils.sh"    # General-purpose bash-script utility functions
export TESTS_DIR=$(cleanpath "$SCRIPT_DIR/..")   # path to RNAstructure/tests
export ROOT_DIR=$(cleanpath "$TESTS_DIR/..")      # path to RNAstructure
#echo "SCRIPT_DIR=$SCRIPT_DIR  TESTS_DIR=$TESTS_DIR  ROOT_DIR=$ROOT_DIR"

DEFAULT_DATAPATH="$ROOT_DIR/data_tables"

if [[ ! $DATAPATH ]]; then
    DATAPATH="$DEFAULT_DATAPATH"
    #echo "*** DATAPATH was not set. Using default: '$DATAPATH'"
elif [[ ! -d $DATAPATH ]]; then 
    echo >&2 "*** Warning: DATAPATH does not point to a valid directory: '$DATAPATH'"
    DATAPATH="$DEFAULT_DATAPATH"
    echo "*** Using default DATAPATH: '$DATAPATH'"
fi

## If the path STILL doesn't work, abort the tests. (They will all fail, so there's no point.)
if [[ ! -d $DATAPATH ]]; then 
  echo >&2 "*** ERROR: DATAPATH is not set or does not point to a valid directory: '$DATAPATH'"
  exit 1
fi

#Make sure DATAPATH is available to child processes.
export DATAPATH


################################################################################
#    Set appropriate diff flags. (Ignore carriage-returns (CR) on windows, 
#  because Windows versions of programs output CRLF instead of just LF, and we
#  need to compare these files to the reference files produced on linux.
################################################################################
#  $OSTYPE Should already be set from the shell environment.
if [[ $OSTYPE == cygwin || $OSTYPE == msys ]]; then
    [[ $DIFF_FLAGS == *'--strip-trailing-cr'* ]] || DIFF_FLAGS+=' --strip-trailing-cr'
    OSNAME=win
elif [[ $OSTYPE == darwin* ]]; then
    OSNAME=mac
fi
: ${DIFF_CMD:=diff}
export DIFF_FLAGS DIFF_CMD
export OSNAME

################################################################################
#                        Recipe Functions 
#    This section includes recipes for standard test scripts.
################################################################################

#-------------------------------------------------------------------------------
#                          runStandardScript
#    This recipe is used for commands that correspond to a stanadard shell 
#  script. "standard" means the script must be located in 
#      scripts/   or   <command>/  
#  and the script should be named one of:
#      <command>, <command>.sh, <command>.bash, or <command>_Script
#  E.g. the script for the command "bifold" could be located at 
#       "scripts/bifold.sh"   or   "bifold/bifold_Script"  etc
#-------------------------------------------------------------------------------
runStandardScript() {
    local cmd=$1;  [[ $cmd ]] || { echo >&2  'Missing command name!'; exit 1; }
    export EXENAME=$cmd 
    export EXEBASE=$EXENAME  OKDIR=$EXENAME  ERROR_OUTPUT_DIR=${EXENAME}_OUTPUT

    # Launch the first script that matches one of the expected paths. Pass a single argument -- the name of the command.
    runFirstScript  "$cmd"   {$EXENAME,scripts}/$EXENAME{,.sh,.bash,_Script}  --  $EXENAME
}

#-------------------------------------------------------------------------------
#                          runSmpScript
#    Similar to runStandardScript, but if the command is <basecmd>-smp,
#  It is assumed that the script is located in:
#      scripts/   <basecmd>-smp/   or  <basecmd>/
#  and is named:
#      <basecmd>-smp, <basecmd>-smp.sh, <basecmd>-smp.bash, <basecmd>-smp_Script
#      <basecmd>, <basecmd>.sh, <basecmd>.bash, or  <basecmd>_Script
#
#    The *-smp variants of the folder and script name are attempted first.
#  this way, an smp-specific script can override the serial analog.
#  Please see 'Note regarding variables used by individual test scripts.' below.
#-------------------------------------------------------------------------------
runSmpScript() {
    local cmd=$1;  [[ $cmd ]] || { echo >&2  'Missing command name!'; exit 1; }
    # EXEBASE: remove -smp from command name (e.g. bifold-smp => bifold)
    export EXENAME=$cmd      EXEBASE=${cmd%-smp}
    export OKDIR=$EXEBASE    ERROR_OUTPUT_DIR=${EXENAME}_OUTPUT
    export SMP=1  
    runFirstScript "$cmd"  \
        {$EXENAME,$EXEBASE,scripts}/{$EXEBASE,$EXENAME}{,.sh,.bash,_Script}  \
        --  "$EXENAME" 
}

#-------------------------------------------------------------------------------
#                          runGuiScript
#    Similar to runStandardScript, but if the command is <basecmd>-gui,
#  It is assumed that the script is located in:
#      scripts/   <basecmd>-gui/   or  <basecmd>/
#  and is named:
#      <basecmd>-gui, <basecmd>-gui.sh, <basecmd>-gui.bash, <basecmd>-gui_Script
#      <basecmd>, <basecmd>.sh, <basecmd>.bash, or  <basecmd>_Script
#
#    The *-gui variants of the folder and script name are attempted first.
#  this way, an smp-specific script can override the serial analog.
#  Please see 'Note regarding variables used by individual test scripts.' below.
#-------------------------------------------------------------------------------
runGuiScript() {
    local cmd=$1;  [[ $cmd ]] || { echo >&2  'Missing command name!'; exit 1; }
    # EXEBASE: remove -gui from command name (e.g. bifold-gui => bifold)
    export EXENAME=$cmd      EXEBASE=${cmd%-gui}
    export OKDIR=$EXEBASE    ERROR_OUTPUT_DIR=${EXENAME}_OUTPUT
    export GUI=1
    runFirstScript "$cmd"  \
        {$EXENAME,$EXEBASE,scripts}/{$EXENAME,$EXEBASE}{,.sh,.bash,_Script}  \
        --  "$EXENAME" 
}

#------------------------------------------------------------------------------
#  Note regarding variables used by individual test scripts.
#------------------------------------------------------------------------------
#      In the current RNAstructure regression test system, an
#  application is run and the output is compared to reference files to 
#  determine if the application is working as expected. 
#  The reference files usually end in "_OK.*" so we can refer to them as OK-files.
#  The OK files usually follow the naming convention  <EXENAME>/<EXENAME>_<TESTNAME>_OK.<EXT>
#  However sometimes the folder or file name differs.
#    e.g. The OK-files for `partition` are named "pfunction/..."  (i.e. different directory than usual)
#    e.g. Usually SMP versions of programs use the same OK-files as non-SMP versions, so the folder and base-name BOTH differ from the EXENAME.
#  There may be other examples where either the folder or the base-name of the OK-files differ from the EXENAME
#  To mitigate these inconsistencies, but still provide a consistent API for writing tests,
#  we set the following variables that are used to derive the names of the OK files.
#    1)  OKDIR   - the directory for OK-files
#    2)  OKNAME  - the base name for OK-files (before the TEST name)
#  I.e. an OK-File name is composed using the following:
#    OK=<OKDIR>/<OKNAME>_<TESTNAME>_OK    (The extension is usually left off for the test to add itself)

#----------------------------------------------------------------------
#  Run the first script that corresponds to a valid file/command
#  Pass any arguments to the command that appear after "--" in the arguments list.
#  Usage: runFirstScript COMMAND-NAME  [POSSIBLE_SCRIPT_FILES] [-- SCRIPT_PARAMS...]
#----------------------------------------------------------------------
runFirstScript() {  
    local cmd="$1"; shift

    for file in "$@"; do
        [[ $file == '--' ]] && break; # the command was not found. show error at end.
        if [[ -f $file ]] || type -p "$file"; then
            # We found a valid file.
            # First, remove all arguments up to and including "--"
            local arg; for arg in "$@"; do [[ $arg == '--' ]] && { shift; break; } || shift; done  
            
            # Next, include common script environment and functions.
            source "$SCRIPT_DIR/test-env.sh"        # Common environment variables for regression tests. (Variables specific to a single test should be moved out of 'test-env.sh' and into the local script.)
            source "$SCRIPT_DIR/test-tools.sh"      # Functions used in RNAstructure regression test scripts.

            #Import configuration scripts when specified by the TEST_CONFIG environment variable.
            [[ $TEST_CONFIG ]] && { readTestConfig || exit 1; }
            
            # Source the script to run it.
            isQuiet || {
                echo "Running Test Script $file"
                echo -n "Test Options: "
                isTestOpt --show
                echo
                [[ $INCLUDE_TESTS && $INCLUDE_TESTS != '*' ]] && echo "Including Tests: $INCLUDE_TESTS"
                [[ $EXCLUDE_TESTS ]] && echo "Excluding Tests: $EXCLUDE_TESTS"
            }
            source "$file" "$@"  # run the script/command with any arguments remaining after "--"
            return
        fi
    done
    echo >&2 "The test command '$cmd' failed because the required test script was not found."
    exit 1
}

################################################################################
function readTestConfig() {
  for scr in $TEST_CONFIG; do
      [[ ! -f $scr && -f .conf/$scr.conf ]] && scr=".conf/$scr.conf"
      if [[ -f $scr ]]; then
        source $scr
        isQuiet || echo "Applied test configuration from $scr"
      else
        echo "Failed to apply test configuration from $scr: File not found."
        return 1
      fi
  done
}
################################################################################


################################################################################
#                    Command Processing 
#    This section handles the actual processing of arguments to the script.
################################################################################

case "$1" in 
    --serial|--cuda) runStandardScript "$2" ;;
    --smp) runSmpScript "$2" ;;
    --gui) runGuiScript "$2" ;;
    *) echo >&2 'ERROR: Unknown script command: '"$1" ;;
esac    
