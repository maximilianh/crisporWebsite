#!/bin/bash

# build-compare.sh  -- Compile a program from source using different build 
#     flags (e.g. #defines) and compare the resulting assembly output.

SELF_FILE="${BASH_SOURCE[0]}"
SELF_NAME="${SELF_FILE##*/}"

# Set program defaults (shown in usage message)
: ${BCMP_TASKS:='BASR'}   # Operations to perform -- B(build)  A(asm-diff)  S(src+asm-diff)  R(run)
: ${DIFF_PROG:='meld'}    # Diff program
: ${CPP_FLAGS:='-fdiagnostics-color=always -O3'}     # Common compiler (e.g. g++) build flags for all prog_nums
: ${CXX:='g++'}           # c++ compiler
: ${DASM_FLAGS:='-drwCS'} # Flags for the objdump program used to produce asm+source output from disassembly
: ${OUTPUT_DIR:='bcmp-output'}
# Ensure the following are defined (but leave empty)
: ${STDIN_DATA:=''}       # input piped into each program when run
PAGER='less -R --quit-if-one-screen --no-init'

if [[ -t 1 ]]; then
  # Some formatting
  F_BR=$'\e[1;31m' # bold-red
  F_BG=$'\e[1;32m' # bold-green
  F_BY=$'\e[1;33m' # bold-yellow
  F_BB=$'\e[1;34m' # bold-blue
  F_BM=$'\e[1;35m' # bold-magenta
  F_BC=$'\e[1;36m' # bold-cyan
  F_B=$'\e[1m'    F_D=$'\e[2m'    F_U=$'\e[4m' # Set Bold, Dim, and Underline
  F_0D=$'\e[21m'  F_0B=$'\e[22m'  F_0U=$'\e[24m' # Remove Bold, Dim, and Underline
  F_0C=$'\e[39m' # Remove FG Color (i.e. use Default)
  F_0=$'\e[0m'   # Reset all formatting
else
  F_BR= F_BG= F_BY= F_BB= F_BM= F_BC= F_BW=
  F_B= F_D= F_U= F_0C= F_0D= F_0B= F_0U= F_0= 
fi
F_CMD="$F_B" # format a command
CQ='`'

function showUsage() { 
  local QS="$CQ${F_CMD}"  QE="${F_0}$CQ"
    cat <<END-OF-TEXT
${SELF_NAME} -- Compile a program from source using different build 
flags (e.g. #defines) and compare the resulting assembly output.

Usage:  ${SELF_NAME} [INPUT_FILE] [PREFIX] [OPTIONS] [-- [PROG_BUILD_FLAGS...]]
  PREFIX     -- The prefix for output files (PREFIX#.exe, PREFIX#.asm, ... etc).
      * If not specified, it is set to the basename of the first input file.
      * Important: If PREFIX starts with a hyphen or if INPUT_FILE is missing, 
        then PREFIX *must* be specified with the -p option (see below)
      * PREFIX may include a directory path. If it does not, output files will
        be created in OUTPUT_DIR ('${OUTPUT_DIR}').
  INPUT_FILE -- The first input source (cpp) file. If there are more source 
        files, list them with the -s option.
  PROG_BUILD_FLAGS -- For each program to be built, list exactly one argument
        containing the custom build flags for it. Each argument may contain 
        multiple flags, separated by spaces.
        e.g. -- "-DFRED -DALICE"  "-w -DBOB"  "-std=c++11 -D BANANA"
      * Make sure to list these AFTER the '--' separator so they aren't
        mistaken as OPTIONS.
      * If the -n option is provided, it is not necessary to list build 
        flags for each program.
      * The flag -DTEST# is always included, where '#' is the 1-based index
        of the program being built.
  OPTIONS
    -s SOURCE_FILE -- List additional source files (cpp).
    -p PREFIX      -- Specify PREFIX explicitly. (See PREFIX above.)
    -n PROG_COUNT  -- Specify the number of programs to build. If this option 
        is absent, the number of programs is taken from the number of
        PROG_BUILD_FLAGS arguments.
    -i STDIN_DATA  -- Specify a string that is piped into each program (as 
        STDIN) when it is run.
    -t BCMP_TASKS  -- list the tasks (operations) to perform. These can include:
        B(build)  A(asm diff)  S(src+asm diff)  R(run). The default is '${BCMP_TASKS}'
    -c CPP_FLAGS   -- list the common cpp build flags used for all programs.
      * The default is '${CPP_FLAGS}'.
    -x CXX         -- Set the c++ compiler. Default: '${CXX}'.
    -d DASM_FLAGS  -- list the flags for the objdump program that produces the
        asm+source output by diassembling the binary output.
      * The default is '${DASM_FLAGS}'. See ${QS}man objdump${QE} for details.
    -D DIFF_PROG   -- Use a custom DIFF program. Default is '${DIFF_PROG}'.
    -o OUTPUT_DIR  -- Set the default output directory (only used when PREFIX 
        does not include a directory component). Default: '${OUTPUT_DIR}'.
    -v  -vv  -vvv  -- Verbose output (include echoing of settings).
      * Levels 2 and 3 include even more output.
    -h             -- Show Help and Usage information.

    Environment Variables:
    The following settings take their default values from the environment:
        BCMP_TASKS  CPP_FLAGS  DASM_FLAGS  STDIN_DATA  
        DIFF_PROG   CXX        OUTPUT_DIR
END-OF-TEXT
}

#function cleanexit() { trap - EXIT; exit ${1:-}; }
function showUsageError() {
  echo >&2 "$*"
  showUsage | grep '^Usage:' >&2
  die "Run $CQ$F_CMD$SELF_FILE -h$F_0$CQ to show additional help and usage information."
}

# Determine if the required verbosity is set.
function is_verbose() { local min=${1:-'@1'}; (($VERBOSE>=${min:1})); }

# Display debug message (if VERBOSE is set to an appropriate level)
# Usage: dbg [@LEVEL] [ECHO_ARGS]...
# LEVEL: Minimum verbosity level (e.g. 1-4)
function dbg() {
  local min=@1; if [[ ${1:-} == '@'[1234] ]]; then min=$1; shift; fi
  ! is_verbose $min || echo "$@"
}
# Same as dbg, but uses printf.
function dbgp() {
  local min=@1; if [[ ${1:-} == '@'[1234] ]]; then min=$1; shift; fi
  ! is_verbose $min || printf "$@"
}

function warn() { echo >&2 "${F_BY}WARNING:${F_0} $*"; }
function die() { echo >&2 "$*"; untrap; exit 1;  }

# initial setup
PREFIX=
INPUT_FILES=()    # List of cpp input (source) files
PROG_BUILD_FLAGS=()    # Custom build flags for each program
PROG_COUNT=0      # Number of programs to build
VERBOSE=0

# Add the first argment as an input file (as long as it doesn't look like an option.)
if [[ $1 != -* ]]; then INPUT_FILES+=("$1"); shift; fi

# Add the next argment as the PREFIX  (as long as it doesn't look like an option.)
if [[ $1 != -* ]]; then PREFIX="$1"; shift; fi

# Read all options in (until '--' is encountered)
while (($#)); do
  case "$1" in
    --) shift; break ;; # exit the option loop. the rest are PROG_BUILD_FLAGS
    -s) INPUT_FILES+=("$2");  shift 2;; # List additional source file
    -n) PROG_COUNT=$2;        shift 2;; # the number of programs to build.
    -i) STDIN_DATA="$2";      shift 2;; # a string that is piped into each program
    -t) BCMP_TASKS="$2";      shift 2;; # tasks to perform.
    -c) CPP_FLAGS="$2";       shift 2;; # the common cpp build flags used for all programs.
    -x) CXX="$2";             shift 2;; # the c++ compiler
    -d) DASM_FLAGS="$2";      shift 2;; # the flags for the objdump program
    -D) DIFF_PROG="$2";       shift 2;; # custom diff program
    -v)   VERBOSE=1; shift;;
    -vv)  VERBOSE=2; shift;;
    -vvv) VERBOSE=3; shift;;
    #-v|-vv|-vvv) VERBOSE=$((${#1}-1)); shift;; # set verbosity to length(arg)-1
    -h) showUsage; exit;;
    *)  showUsageError "Invalid option/flag: '$1'"
  esac
done

# Remaining arguments are PROG_BUILD_FLAGS
while (($#)); do
  PROG_BUILD_FLAGS+=("$1"); shift
done

# Set PREFIX if unset
if [[ -z $PREFIX && ${#INPUT_FILES[@]} > 0 ]]; then 
  PREFIX="${INPUT_FILES[0]##*/}" # Remove directory
  PREFIX="${PREFIX%.*}"          # Remove suffix (..cpp etc)
fi
# Set PROG_COUNT if it has not already been set
(($PROG_COUNT)) || PROG_COUNT=${#PROG_BUILD_FLAGS[@]}
# Alternatively, use maximum: PBF_COUNT=${#PROG_BUILD_FLAGS[@]}; PROG_COUNT=$((PROG_COUNT>PBF_COUNT?PROG_COUNT:PBF_COUNT))

# Display interpretted options/settings
if is_verbose; then
  IFS="','" # array element separator for ${arr[*]}
  printf '%s\n' \
    PREFIX="$PREFIX" \
    INPUT_FILES="('${INPUT_FILES[*]}')" \
    PROG_BUILD_FLAGS="('${PROG_BUILD_FLAGS[*]}')" \
    PROG_COUNT=$PROG_COUNT \
    BCMP_TASKS="$BCMP_TASKS" \
    CPP_FLAGS="$CPP_FLAGS" \
    DASM_FLAGS="$DASM_FLAGS" \
    "STDIN_DATA(size)"="${#STDIN_DATA}" \
    OUTPUT_DIR="$OUTPUT_DIR" \
    DIFF_PROG="$DIFF_PROG"
  unset IFS
fi

# Verify arguments and options.
[[ $PREFIX ]]                 || showUsageError "Missing PREFIX."
[[ ${#INPUT_FILES[@]} != 0 ]] || showUsageError "No input files have been listed."
[[ $BCMP_TASKS ]]             || showUsageError "No tasks were specified. Set BCMP_TASKS (-t)."

if [[ $PROG_COUNT == 0 ]]; then
  echo >&2 "Unknown number of programs--assuming 1. (Specify '-n PROG_COUNT' or PROG_BUILD_FLAGS)"
  PROG_COUNT=1
fi

# Turn on improved error handling.
set -eEuo pipefail
# set -x # turn on script command echoing output for debugging
function untrap() { trap - ERR EXIT; _exited=1; }; _exited=
function reportErr() { [[ ${_exited} ]] || echo >&2 "ERROR at line $2 (rc: $1)"; exit $1; }
trap 'reportErr $? $LINENO' ERR
#trap 'rc=$?; (($rc)) && { echo >&2 "EXITING DUE TO ERROR ($rc)"; exit $rc; } || echo Done.' EXIT


prog_nums=$(seq 1 $PROG_COUNT)
# Update PREFIX to include a directory
[[ $PREFIX == */* ]] || PREFIX="$OUTPUT_DIR/$PREFIX"

function verifyDir {
  local d=$(dirname "$1")
  mkdir -p "$d"
}

function setProgVars() {
  local i=$1
  local P="$PREFIX$i"
  exe=$P.exe  asm=$P.asm  src=$P.sasm  out=$P.stdout
  defs="-DTEST$i ${PROG_BUILD_FLAGS[i-1]:-}"
  outfiles=($exe $asm $src $out)
}

function safeRm() { local f; for f; do [[ ! -f $f ]] || rm $f; done; }
function removeOldFiles() {
  local P="$PREFIX"
  if [[ $1 == '*' ]]; then
    safeRm $P?.exe  $P?.asm  $P?.sasm  $P?.stdout
  else
    setProgVars $1; safeRm "${outfiles[@]}"
  fi
}

function cleanFile() {
  if declare -f cleanOutputFile >/dev/null; then
    cleanOutputFile "$@"
  fi
}

# Removes unnecessary or noisy text from output files (e.g. asm, src+asm)
function cleanOutputFile() {
    local op="$1" f="$2"

    if [[ $op == A ]]; then
        dbg @2 "Cleaning $op $f"
        #.LFB2837:
        perl -i -pe 's/^\.L\w+://;' \
                -pe 's/\.constprop\.\d+/.constprop.0/g;' \
                "$f"
    elif [[ $op == S ]]; then
        dbg @2 "Cleaning $op $f"
        perl -i -pe 's/.*_IO_stdin_used.*//;' \
                -pe 's/^ +[a-f0-9]+://;'  \
                -pe 's/.*file format pei-x86-64.*//;'  \
                -pe 's/\.constprop\.\d+/.constprop.0/g;' \
                "$f"
    elif [[ $op == [BR] ]]; then
        :
    else
        echo >&2 "WARNING: cleanFile: Unknown operation type: $op (command: $@})"
    fi
}

function run() { 
  if is_verbose @2; then echo "Running"; printf ' %q' "$@"; fi
  "$@"
}

function testdiff() { 
    if diff -q "$1" "$2" >& /dev/null; then
      echo "${F_BG}Identical $3 Files: $A $B${F_0}"
    else
      echo >&2 "${F_BR}Files differ: $A and $B !${F_0}"
      run $DIFF_PROG $A $B
    fi
}

# Build (Compile) programs
if [[ $BCMP_TASKS == *[BC]* ]]; then
  # remove previous builds (in case of failure)
  removeOldFiles '*'
  for N in ${prog_nums[@]}; do
      setProgVars $N
      verifyDir $exe
      dbg "Compiling $exe with defs: $defs"
      $CXX $CPP_FLAGS -g $defs -o $exe "${INPUT_FILES[@]}" |& $PAGER || die "Aborting due to build errors."
      cleanFile 'B' $exe

      $CXX $CPP_FLAGS -S $defs -o $asm "${INPUT_FILES[@]}" |& $PAGER || die "Aborting due to build errors."
      cleanFile 'A' $asm

      objdump ${DASM_FLAGS:-'-drwCS'} $exe > $src
      cleanFile 'S' $src
  done
fi

# Run programs
if [[ $BCMP_TASKS == *R* ]]; then
  for N in ${prog_nums[@]}; do
    setProgVars $N
    dbg "Running Prog $N  (defs: $defs):"
    # run with "$STDIN_DATA" as STDIN
    $exe <<< "${STDIN_DATA:-}" | tee $out  || warn "Program $N exited with error: $?"
    cleanFile 'R' $out
  done
fi

# Run diff assembly
if [[ $BCMP_TASKS == *A* ]]; then
  for N in ${prog_nums[@]:1}; do
    setProgVars $((N-1)); A=$asm
    setProgVars $N;       B=$asm
    testdiff  $A $B  'ASM'
  done
fi

# Run diff objdump assembly+source
if [[ $BCMP_TASKS == *S* ]]; then
  for N in ${prog_nums[@]:1}; do
    setProgVars $((N-1)); A=$src
    setProgVars $N;       B=$src
    testdiff  $A $B  'ASM+SRC'
  done
fi

echo "Done"