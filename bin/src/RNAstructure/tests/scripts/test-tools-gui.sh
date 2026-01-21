: ${CALC_TIMEOUT:=10} # seconds
TEST_DIR=$(winpath "$PWD") # correct paths on Windows

# Variables the calling script can change
GUI_STARTUP="StartGUI '-prefs ~test -opt draw-structures=P  -sfc -jmnu'
   WAIT { Menu 'RNA' }, 20000" #Wait for menu creation in case startup is slow.

# # GUI Scripts MUST OVERRIDE the following function.
# # The function is passed all of the arguments to the test
# # and it should output (using echo, printf etc) the content
# # of the gui-script that should be run.
# # If the return value is not 0, the test is not run.
# function guiScript() { 
#   (($# >= 2)) || { warnmsg "$EXEBASE requires 2 arguments."; return  1; }
#   local seqFile="$1" ctFile="$2"
#   PARAM_FLAGS=(s t l c m w p md sh si sm) # Flags with required parameters
#   BOOL_FLAGS=(d) # Flags with NO parameters
#   getGuiArgs ${@:3} # Parse arguments starting with $3
# }

# # Optionally override this to perform any final operations here, after the GUI test-script has run.
# # The first argument is the exit code from guiScript
# # The remaining arguments are the same as those passed to guiScript.
# function guiPostRun() { :; } 

GUI_ARG_PREFIX=arg
function guiParseArgs() {
    local pos arg flag skip param

    # Clear previous variables unless '+' is the first argument.
    if  [[ $1 == '+' ]]; then
      shift
    else
      guiUnsetArgs # Clean up any temporary variables remaining from previous test
    fi

    while ((pos++ < $#)); do
      arg="${@:pos:1}"
      if [[ ${arg:0:2} == -- ]]; then
        flag="${arg:2}" # remove --
      elif [[ ${arg:0:1} == - ]]; then
        flag="${arg:1}" # remove '-' 
      else
        # Not a flag, and not a valid argument to any other flag.
        errmsg "Unknown argument to $EXE: $arg"
        continue
      fi

      if inArr BOOL_FLAGS "$flag"; then
        # if -X is present, set variable $gui_X=1
        setVar "arg_$flag" 1
        #warnmsg "Set var arg_$flag to 1 (TRUE)"
      elif inArr PARAM_FLAGS "$flag"; then
        param="${@:1+pos++:1}"  # Get next argument AND advance position counter
        # if -X "PARAM" is present, set variable $gui_X="PARAM"
        setVar "arg_$flag" "$param" 
        #warnmsg "Set var arg_$flag to $param"
      else
        errmsg "Unknown command-line flag to $EXE: $arg"
      fi
    done
}

function guiParseConf() {
  local name eq val prefix=$2
  while IFS=' ' read name eq val; do
      [[ -z $name || $name == '#'* ]] && continue;
      if [[ $eq == '=' ]]; then
        echo "found conf  $name=$val"
        setVar "$prefix$name" "$val"
      else
        errmsg "Unknown config file parameter: $name"
      fi
  done < "$1"
}

# Unset all variables set by guiParseArgs
function guiUnsetArgs() {
    # Cleanup before next test run
    # for var in "${BOOL_FLAGS[@]}" "${PARAM_FLAGS[@]}"; do
    #   unset "_$var"
    # done
    unset ${!arg_*}
}

function guiCleanup() {
  guiUnsetArgs
}

# Write line(s) to the GUI script
function guiWrite() {
  local echoArgs arg words=()
  for arg; do
    case "$arg" in
      -e|-n) echoArgs+=" $arg" ;;
      -s) rm -f $TEST.gui ;; # -s means Start the script
      *) words+=("$arg")
    esac
  done
  [[ -s $TEST.gui ]] || echo "$GUI_STARTUP"  >> $TEST.gui
  echo $echoArgs "${words[@]}"  >> $TEST.gui
}
# A terse way to write script lines, based on the presence of arguments.
# Usage:  guiWriteIf -ARG "SCRIPT_LINES"
# This will write SCRIPT_LINES if ARG was found in the list of 
# command-line arguments (assuming that guiParseArgs was called)
function guiWriteIf() { guiArg "$1" && guiWrite "${@:2}"; }

# Returns 0 (success) if the specified argument is set.
# Usage: guiArg -X "LINES"
#  The above will write "LINES" to the gui script if -X was present
#  as a command-line argument.  The preceeding hypen in -X is optional.
function guiArg() {  
  local var=$1
  [[ ${var:0:1} == - ]] && var="${var:1}" # remove hyphen (-) if present. e.g. -d
  [[ ${var:0:1} == - ]] && var="${var:1}" # remove second hyphen (-) if present. e.g. --sequence
  local var=arg_$var
  [[ ${!var} ]]
}

function guiRun() {
    local cleanup

    for arg; do
      case "$arg" in
        -c) cleanup=1 ;;
        *) errmsg "Invalid argument to guiRun: $arg" ;;
      esac
    done

    echo 'CloseGUI' >> $TEST.gui

    # Cleanup spacing in gui-script
    sed -E -i.bak 's/^[ \t]+//g' $TEST.gui
    rm -f $TEST.gui.bak

    local runFlags='-d -t'
    if isTestOpt 'DebugGui' ; then
      echo -n "$FMT_BLD$FMT_PRP"
      nl -b a $TEST.gui
      echo -n "$FMT_RST"
    fi
    export GUITESTER_LOG_OUT_FILE="$TEST.guilog"

    nl -b a $TEST.gui >> "$TEST.guilog"

    runTest run-gui-test $runFlags $GUI_TEST_FLAGS "$TEST.gui"
    [[ ! $cleanup ]] || guiCleanup
}
