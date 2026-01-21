# script-utils.sh  -  General Purpose Bash Utility Functions
# Author: Richard M. Watson


###############################################################################
## diffArr - Compares the values in ARR2 with those in ARR1
##           Each item in ARR2 that is not also in ARR1 is stored in the global 
#            variable 'DIFF_ARR'.
## Usage: diffArr ARR1 ARR2;  result=( "${DIFF_ARR[@]}" )
##    ARR1 and ARR2 are the names of the globally-available array variables.
##  ex: a1=( a b d g ); a2=( d e f g h );   diffArr a1 a2
##      echo "${DIFF_ARR[@]}"  # outputs:  e f h
###############################################################################
function diffArr() {
    local __a1=$1'[@]' # e.g. ${list1[@]}
    local __a2=$2'[@]'
    DIFF_ARR=()
    if [[ $3 == '-d' ]]; then
        echo a1: $__a1
        printf "    %s\n" "${!__a1}"
        echo a2:  $__a2
        printf "    %s\n" "${!__a2}"
        echo diff:
    fi
    
    local i j
    for i in "${!__a2}"; do
        for j in "${!__a1}"; do
            [[ $i == $j ]] && continue 2; ## continue outer loop
        done
        DIFF_ARR+=("$i")
    done
}
###############################################################################
## Get the absolute path to the specified file or folder.
## The file or folder need not exist, but its parent directory MUST exist.
## Usage: abspath PATH [--ret:OUTPUT_VAR]  or   path=$(abspath PATH)  etc.
##    If OUTPUT_VAR is specified, the result will be stored in that variable.
###############################################################################
function abspath() {
    local __file="$1" __retval
    if [[ $__file == *[\\/]* ]]; then
        __file="${1##*[\\/]}" #remove directory
        if [[ $__file == '.' || $__file == '..' ]];  then
            pushd "$1"
            __retval="$PWD"
        else
            pushd "${1%[\\/]*}"
            __retval="$PWD/$__file"
        fi
        popd
    else
        __retval="$PWD/$__file"
    fi
    echoOrSetVar "$2" "$__retval" __file __retval
}

# Fix paths (for Windows)
function winpath() {
    if [[ $OSTYPE == cygwin ]]; then
        cygpath -w "$1"
    else
        echo "$1"
    fi
}

# Determines if the given name represents a defined variable.
# The exit code is 0 (success) if the variable exists (even if it is blank) 
# or 1 if it is undefined.
function isVar() {
    # Use bash test. 
    #   -z tests to see if the value is empty.
    #   !1 uses indirection to get the value of the variable specified by $1
    #   ${var+x} is a parameter expansion which evaluates to the null if 
    #   var is unset and substitutes the string "x" otherwise
    [[ ! -z ${!1+x} ]]
}

# Determines if the given name represents a defined function.
function isFunc() { declare -f -F "$1" &>/dev/null; }

# Determines if $2 is an element in the array named by $1
function inArr() {
    # local variables won't pollute the outer scope, but if 
    # the array name is the same as the local variable name, the
    # local variable will shadow the outer one, so it is important to have
    # a unique name for the indirection variable.
    local i __varr__="$1"
    [[ $__varr__ == *'[@]' ]]  || __varr__+="[@]"
    for i in "${!__varr__}"; do
        [[ $i == $2 ]] && return 0
    done
    return 1
}

# Returns 0 (success) if the first argument is an executable file, script or 
# command (in the PATH). Returns 1 otherwise.
function isExe() {
    type -p "$1" >/dev/null
}

# Usage: copyToDir [OPTIONS] TARGET_DIRECTORY [SOURCE_FILES...]
# Same as `cp -fp`, but the directory is specified before the files, 
#          and the directory is created if it does not exist.
#          (No error is shown, and the target directory is NOT created if 
#          there are no SOURCE_FILES.)
# OPTIONS: 
#   At most one option is allowed: -e, -s, or --
#   -e  Non-existing files will be ignored (otherwise an error message is shown for each one)
#   -s  Only non-empty files will be copied (also implies -e)
#   --  Use to avoid ambiguity if TARGET_DIRECTORY could be "-e" or "-s"
function copyToDir() {
    local dir="$2" # assume an option was specified
    local _FILTERED_FILES # prevents filterFiles results from being global 
    case "$1" in
        -e) filterFiles -e "${@:3}" ;;
        -s) filterFiles -s "${@:3}" ;;
        --) _FILTERED_FILES=( "${@:3}" ) ;;
        *) dir="$1"; _FILTERED_FILES=( "${@:2}" ) ;;
    esac
    ((${#_FILTERED_FILES[@]}==0)) && return

    # If we get here, there is at least one file to copy, so create the directory
    [[ -d  $dir ]] || mkdir -p "$dir"
    cp -fp "${_FILTERED_FILES[@]}" "$dir"
}

# This function tests each file and builds an array of files that exist and 
# (optionally) are non-empty.
# Usage:  filterFiles -TEST FILES..
#  TEST is one of -e=exists, -n=doesn't exist, -s=non-empty -0=exists, but is empty
# Returns: 0 (success) if any files matched the criteria and 1 (failure) 
#          if no files do.
# Results: The matching files are stored in the global array _FILTERED_FILES
function filterFiles() {
    local file test="$1"
    shift
    _FILTERED_FILES=()
    case "$test" in
        -e) for file; do [[ -e $file ]] && _FILTERED_FILES+=($file);  done ;;
        -s) for file; do [[ -s $file ]] && _FILTERED_FILES+=($file);  done ;;
        -n) for file; do [[ ! -e $file ]] && _FILTERED_FILES+=($file);  done ;;
        -0) for file; do [[ -e $file && ! -s $file ]] && _FILTERED_FILES+=($file);  done ;;
        *) echo >&2 "$FUNCNAME requires a flag as its first argument (-e, -s, -n, -0)"; return 1 ;;
    esac 
    ((${#_FILTERED_FILES[@]}>0)) # return false (1) if _FILTERED_FILES is empty and true(0) otherwise.
}


# Uses `find` to list files and place them in a global array named "FILES_ARRAY".
# Usage:  filesToArray STARTDIR [OPTIONS];  result=( "${FILES_ARRAY[@]}" )
#   STARTDIR and OPTIONS are passed directly to find: `find STARTDIR OPTIONS...`
# Example:  List all files in the current directory:
#    filesToArray . -maxdepth 1 -type f ;  result=( "${FILES_ARRAY[@]}" )
function filesToArray() {
    FILES_ARRAY=()
    while IFS= read -r -d '' file <&4; do
        FILES_ARRAY+=("$file")
    done 4< <(find "$@" -print0)
}

#####################################################################################
## normpath - Outputs the normalized (canonical) path to a file, similar to realpath, 
##   except that by default normpath does NOT resolve symlinks (although there is a 
##   flag that causes it to do so.) Furthermore (by default) the path need not exist.
##   This function is required mainly because 'realpath' is not present on many
##   systems and 'readlink' resolves symlinks, which can be undesirable.
##
##   * This outputs the absolute path to a file, unless the --rel flag has been
##     specified, in which case the function outputs a path that is relative to the
##     indicated directory.
##   * A 'base' folder can be specified (-b|--base), and it will be prepended to 
##     the subject path when attempting to normalize it (but only if the subject 
##     path is a relative path). 
##   
###############################################################################
## Usage: normpath path [OPTIONS] 
###############################################################################
##  OPTIONS:
##
##  -m|--sym            Resolve symlinks (default is to NOT resolve)
##
##  -s|--silent         Silences error messages about non-existant paths.
##
##  -e|--exists=EXISTS  Specify whether or not the path must exist.
##      Valid values for EXISTS are:
##      *  all  --  The entire resolved path must exist.
##      *  last --  All but the last component of the path must exist.
##      *  none --  No part of the path is required to exist. (default)
##      If the test for existance fails, the output is empty and the
##      exit status will be 1.
##
##  -r|--rel=REFDIR     The resulting output path will be made relative to 
##      REFDIR (even if REFDIR is not a parent directory.)
##
##  -b|--base=BASEDIR   Specify that the original path *might be* relative 
##      to BASEDIR. This causes BASEDIR to be prepended to 'path', but only
##      if 'path' is not absolute. If the path (or a portion of it) is 
##      required to exist (due to the --exists flag), and the test for
##      existance fails when prepending the base, then the function will
##      re-attempt the test without the base.
##
##  -v|--ret=OUTPUT_VAR  Store the result in the bash variable named OUTPUT_VAR
##                      instead of echo'ing it to stdout. This can provide a
##                      substantial speed boost over executing it in a subshell.
##
###############################################################################
function normpath() {
    local __path relative_to baseDir="$PWD" 
    local silent=0 symlinks=0 exists=none
    local opt val outvar
    local argPos=0

    for arg in "$@"; do 
        # an option is -n or --name or -n=value or --name=value  or -n:value or --name:value
        if [[ $arg == -* ]]; then
            opt=${arg%%[=:]*}  #opt keeps the hyphens!
            val=${arg:${#opt}+1}
        else
            # use the index of the positional argument as the option name e.g. 1, 2, 3 etc
            opt=$((++argPos))
            val="$arg"
        fi
        case $(tolower $opt) in #options are case NOT case sensitive (,, causes the variable to be lower case)
            1)                          __path="$val"       ;;
            2|-r|--rel|--relative-to)   relative_to="$val"  ;;
            3|-b|--base)                baseDir="$val"      ;;
            -e|--exist|--exists)        exists="$val"       ;;
            -s|--silent)                silent=1            ;;
            -m|--sym)                   symlinks=1          ;;
            -v|--ret)                   outvar="$val"       ;;
            *) echo >&2 "Invalid flag passed to normpath: $arg"; return 1 ;;
        esac
    done

    while true; do #allows us to repeat the following if -base=BASEDIR fails the first time.
        # first get the absolute path
        local __abspath="$__path"
        [[ $__path == /* ]] || __abspath="$baseDir/$__path"

        if ((symlinks)); then
            __abspath=$(readlink -s -m -n "$__abspath") # -s=silent -m=allow-missing -n=no-newline
        else
            cleanpath "$__abspath" --ret:__abspath
        fi

        local error
        case "$exist" in
            all)        [[ -e $__abspath ]]; error=$?         ;;
            no|none)    error=0                             ;;
            last)       local __dir=$(dirname "$__abspath")
                        [[ ! $__dir || -e $__dir ]]; error=$?   ;;
            *) echo >&2 "Invalid flag passed to normpath: $arg"
               return 1 ;; 
        esac

        if ((error)); then
            if [[ $baseDir == "$PWD" ]]; then
                # there's no point in retrying the loop because baseDir is already PWD.
                # just show the error and exit
                ((silent)) || echo >&2 "normpath: ‘$__path’: No such file or directory"
                return 1
            else
                # the user-specified baseDir failed, so retry with PWD
                baseDir="$PWD"
            fi
        else
            # we have successfully found the path, so exit the loop
            break
        fi
    done

    if [[ $relative_to ]]; then 
        relpath "$__path"  "$relative_to" --ret:__path
    fi
 
    echoOrSetVar "--ret:$outvar" "$__path" __path relative_to baseDir silent symlinks exists opt val outvar argPos error __abspath __dir
}

###############################################################################
# Usage: relpath subjectPath referencePath [--ret:OUTPUT_VAR]
###############################################################################
#   Given a target path (which is either absolute or relative to PWD)
#   and a base path, this outputs a relative path to the target 
#   from to the base path.
#
#   Examples:
#       relpath  /foo/bar/file  /foo/ ==> bar/file
#       relpath  /foo/bar/file  /foo/etc/ ==> ../bar/file
#       relpath  /usr/bin/exe   /home/user/ ==> ../../usr/bin/exe
###############################################################################
function relpath() {
    local base="$1"
    local target="$2"
    local base_parts=()
    local targ_parts=()

    #Split path elements & canonicalize
    local IFS='/'
    local bPos=0;
    for bp in $base; do
        case "$bp" in
            ".");;
            "..") ((bPos--)) ;;
            *) base_parts[bPos++]="$bp" ;;
        esac
    done
    local tPos=0;
    for tp in $target; do
        case "$tp" in
            ".");;
            "..") ((tPos--)) ;;
            *) targ_parts[tPos++]="$bp" ;;
        esac
    done
    #unset IFS

    #Count common prefix
    common=0
    for (( i=0 ; i<bPos; i++ )); do
        if [[ ${base_parts[i]} == "${targ_parts[common]}" ]] ; then
            ((common++))
        else
            break
        fi
    done

    #Compute number of directories up
    local updir=$((bPos-common))

    #trivial case (after canonical decomposition)
    if [[ $updir -eq 0 ]]; then
        echo .
        return 0
    fi

    local __result
    #Print updirs
    for (( i=0 ; i<$updir ; i++ )); do
        __result+=../
    done

    #Print remaining path
    for (( i=common ; i<tPos ; i++ )); do
        if [[ ${targ_parts[i]} ]]; then
            ((i==common)) || __result+='/'
            __result+="${targ_parts[i]}"
        fi
    done

    echoOrSetVar "$3" "$__result" base target base_parts target_parts __result updir tPos bPos
}

###############################################################################
##  usage: cleanpath PATH  [--ret:OUTPUT_VAR]
##  Cleans up PATH by performing the following replacements:
##  * Replaces multiple adjacent slashes with a single slash:
##      PathA//pathB ==> PathA/PathB
##  * Replaces /./ with /
##      PathA/./PathB ==> /PathA/PathB
##  * Removes dir/.. 
##      PathA/PathB/../pathC  ==> /PathA/pathC
##      PathA/PathB/../../  ==> ./
###############################################################################
function cleanpath() {
  tmpopt -s extglob
  
  local __p="$1" oldPath
  local debug=0
  [[ $2 == -d || $3 == -d ]] && debug=1

  # Remove dir/.. sequences.
  #local A='*([^\/])' # pattern meaning 0 or more characters that are NOT slashes (uses extglob)
  #local B='[^\/.]'   # pattern meaning a single character that is neither "/" nor "."
  #local NAME="$A$B$A" 
  local NAME="*([^\/])[^\/.]*([^\/])" # pattern representing any file or directory name that does not consist entirely of . (e.g. not "." or "..")  (uses extglob)
  local S=/ # slash -- useful in replacement patterns to replace a literal slash
  oldPath=
  while [[ $__p != "$oldPath" ]]; do #repeat until path stops changing
    ((debug)) && echo >&2 "path=$__p"
    oldPath="$__p"
    # replace PATH/../  with  ""
    # PATH must include NO slashes, and must not be ".."
    # require the ending slash to distinguish between ../ and ..word (which IS a valid file name)

    # Replace all /./ sequences with /
    __p="${__p//\/.\//\/}"

    # Replace all // with /  (also ///, ////, ... )
    __p="${__p//\/+(\/)/\/}" #uses extglob:  +(pattern) means one or more occurrances

    #  /hello/../  ==> /
    #  /hello/../world  ==> /world
    #  /hello/..world  ==> /hello/..world (no change)
    #  /hello/world/../../  ==> /hello/../
    #  /hello/world/../../path  ==> /hello/../path
    #  ../path  ==> ../path (no change)
    #  ../../path  ==> ../../path (no change)
    #  /hello/w../  ==> (no change -- NAME must be separated from .. by a slash)
    #  /hello/..w/  ==> (no change)
    __p="${__p//$S$NAME$S..$S/$S}" # replace "/ABA/../" with "/"
  done

  # same replacement as in loop, but do not require the ending slash
  # ONLY if the pattern matches at the very END of path (%)
  #  /hello/..  ==> /
  #  /hello/world/..  ==> /hello
  #  /hello/..world  ==> (no change -- the .. must occur at the END of the path)
  __p="${__p/%$S$NAME$S../$S.}" # replace "/ABA/.." with "/."  
  ((debug)) && echo >&2 "path.2=$__p"


  # replaces dir/../  with  ./ at the beginning of the path
  # hello/../  ==> "./"  (should this be an error?)
  __p="${__p/#$NAME$S..$S/.$S}" # replace "/ABA/../" with "/"  
  ((debug)) && echo >&2 "path.3=$__p"

  # replaces dir/.. with  . (only if the full path matches that pattern)
  # hello/..  ==> "."
  # hello/../world  ==> (no change -- should have already been processed to hello/world)
  [[ $__p == $NAME$S.. ]] && __p="." # replace "ABA/.." with "."
  ((debug)) && echo >&2 "path.4=$__p"

  # Removes "./" from the beginning if followed by anything else
  # ./..  ==> ..
  # ./../  ==> ../
  # ./../NAME  ==> ../NAME
  # ./NAME  ==> NAME
  # ./NAME/  ==> NAME/
  # ./ ==> (unchanged)
  [[ $__p == ./?* ]] && __p="${__p:2}"
  ((debug)) && echo >&2 "path.5=$__p"

  # Removes "/." from the end if preceeded by anything else
  # ../.  ==> ..
  # .././  ==> (no change -- should have been already processed to ../)
  # ../NAME/.  ==> ../NAME
  # /. ==> (unchanged)
  [[ $__p == *?/. ]] && __p="${__p:0:-2}"
  ((debug)) && echo >&2 "path.6=$__p"

  # /. ==> /
  [[ $__p == /. ]] && __p="/"
  ((debug)) && echo >&2 "path.7=$__p"

  tmpopt -r extglob # reset shell option
  echoOrSetVar "$2" "$__p" __p oldPath debug S NAME
}

###############################################################################
##  usage: assign OUT_VAR command [command_args...]
##  OUT_VAR: The name of the output variable.
###############################################################################
function assign() { 
    local varName=$1 tmp # line DONE
    tmp=/tmp/$$-$RANDOM$RANDOM
    #echo tmp=$tmp
    # if [[ ! -p $FIFO_ASSIGN ]]; then
    #     FIFO_ASSIGN="/tmp/$$-$RANDOM$RANDOM"
    #     mkfifo "$FIFO_ASSIGN"
    # fi
    shift;
    #"$@" > "$FIFO_ASSIGN" &
    #echo "args=$@"
    "$@" > "$tmp"
    #echo "result=$?"
    # { 
    IFS= read -r -d '' $varName < $tmp
    # echo "val='$val' DONE='$DONE'"
    # until [[ $DONE ]]; do
    #   IFS= read line || DONE=true
    #   echo "line=$line"
    #   IFS= val+=$'\n'"$line"
    # done
    # } < $tmp
    # #done < "$FIFO_ASSIGN"
    rm $tmp
    # setVar $varName "$val"
}

###############################################################################
## Functions to deal with indirect references
###############################################################################
## isVar - return 0 if the variable name represented by $1 has been set (even
##      if it is set to NULL). If the variable is unset, 1 will be returned.
##      isVar $variableName && echo "$variableName is set!"
function isVar() { declare -p &>/dev/null "$1"; }
function setVar() { export -n "$1"="$2"; }
function getVar() { echo "${!1}"; }

###############################################################################
##  usage: echoOrSetVar RET_FLAG VALUE
##  RET_FLAG: Must be either --ret:OUT_VAR or --var:OUT_VAR where OUT_VAR
##            is a valid bash variable name. If RET_FLAG is one of the above,
##            VALUE is assigned to the variable OUT_VAR and the return code is 0.
##            (The colon can optionally be replaced by an equals sign: --ret=OUT_VAR)
##  Otherwise VALUE is echo'd to stdout and the return code is 1.
###############################################################################
function echoOrSetVar() {
  if [[ "$1" == --ret[:=]?* || "$1" == --var[:=]?* ]]; then
    local ___varName="${1:6}" ___arg
    for ___arg in "${@:3}"; do
        # verify that name of the variable is not one of the local variables of the calling function. If it is, setting it will NOT be visible to the outer function.
        [[ "$___arg" == "$___varName" ]] && { 
            echo >&2 "SCRIPT AUTHORING ERROR: In ${BASH_SOURCE[2]}::${FUNCNAME[2]} line ${BASH_LINENO[1]}: Cannot use the variable name '$___varName' as an output variable. It is a local variable in the function ${BASH_SOURCE[1]}::${FUNCNAME[1]}."
            echo "$2"
            return 1
        }
    done
    setVar "$___varName" "$2"
    return
  else
    echo "$2"
    return 1
  fi
}

###############################################################################
##  usage: var=$(tolower $var)
##  tolower, toupper: convert the input text to lower or upper case,
##     respectively and outputs the result.
##  This is required because bash < v4 lacks the ${VAR,,} or ${VAR^^} syntax.
###############################################################################
function tolower() { 
    # Re-define the function on first use, depending on bash version.
    if ((BASH_VERSINFO[0] < 4)); then
        function tolower() { echo "$*" | tr '[:upper:]' '[:lower:]'; }
    else
        function tolower() { echo "${*,,}"; }
    fi
    tolower "$*" # call the newly-defined version of the function
}
function toupper() { 
    # Re-define the function on first use, depending on bash version.
    if ((BASH_VERSINFO[0] < 4)); then
        function toupper() { echo "$*" | tr '[:lower:]' '[:upper:]'; }
    else
        function toupper() { echo "${*^^}"; }
    fi
    toupper "$*" # call the newly-defined version of the function
}

###############################################################################
##  usage: replaceInvalidVarChars <TEXT>
##  Replaces the following characters in TEXT with an unerscore (_)
##          - + @ = . ~ \ ! , ^ $ # % * " ' 
###############################################################################
function replaceInvalidVarChars() { echo "${1//[-+@ =.~\\!,^\$#%*\"\']/_}"; }

###############################################################################
##  usage: ifzero <TEST> <VALUE1> <VALUE2>
##  If TEST evaluates mathematically to 0, this outputs $VALUE1
##    Otherwise (i.e. TEST != 0) it outputs $VALUE2
###############################################################################
function ifzero() { (($1==0)) && echo $1 || echo $2; }

###############################################################################
##  usage:  trim [<TEXT>]   or   <COMMAND> | trim
##  Removes leading and trailing space from stdin or arguments.
###############################################################################
function trim() {
    local var
    if (($#)); then var="$*"; else var=$(cat); fi
    var="${var#"${var%%[![:space:]]*}"}"   # remove leading whitespace characters
    var="${var%"${var##*[![:space:]]}"}"   # remove trailing whitespace characters
    echo "$var"
}

###############################################################################
##  usage:  showVars [<VAR_NAME>...]
##  Prints out the variable names and their corresponding values, each on
##    a separate line.
###############################################################################
function showVars() {
    for var in "$@"; do echo -e "\t* \$"$var"\t"=" ""${!var}"; done
}

###############################################################################
##  usage:  joinArgs [--ret:OUTPUT_VAR]  <DELIM> [<ARGS>...]
##  Joins all ARGS together as a single string using DELIM as a separator.
##  E.G.: joinArgs '|' apple 'banana bread' '' last  ==> "apple|banana bread||last"
###############################################################################
function joinArgs() { 
    local __val __delim __OUTPUT_VAR
    [[ $1 == --ret[=:]?* ]] && { __OUTPUT_VAR="$1"; shift; }

    __delim="$1"; shift
    if ((${#__delim}==1)); then
        local IFS="$__delim"; __val="$*";
    else
        local __arg;
        __val="$1"; shift
        for __arg; do __val+="$__delim$__arg"; done
    fi
    local __RETVAL="$__val"
    echoOrSetVar "$__OUTPUT_VAR" "$__val" __val __delim __OUTPUT_VAR
}

###############################################################################
##  Usage:   if isWindows; then ...
##  The exit code is 0 (success) if the current operating system is Microsoft
##      Windows (as determined based on the OSTYPE environment variable set 
##      by the current shell).  Otherwise the exit code is 1 (failure).
###############################################################################
function isWindows() { [[ $OSTYPE == cygwin || $OSTYPE == msys ]]; }

###############################################################################
##  Usage:   if isMac; then ...
##  The exit code is 0 (success) if the current operating system is Mac OSX
##      (as determined based on the OSTYPE environment variable set 
##      by the current shell).  Otherwise the exit code is 1 (failure).
###############################################################################
function isMac() { [[ $OSTYPE == darwin* ]]; }

###############################################################################
##  Usage:   if isNix; then ...
##  The exit code is 0 (success) if the current operating system is Linux or  
##      similar *nix (as determined based on the OSTYPE environment variable set 
##      by the current shell).  Otherwise the exit code is 1 (failure).
###############################################################################
function isNix() { [[ $OSTYPE == linux* || $OSTYPE == bsd* || $OSTYPE == freebsd* || $OSTYPE == solaris* ]]; }

###############################################################################
##  Usage:  OS=$(getOS)
##  Outputs the canonicalized OS name, based on the OSTYPE bash variable.
##  Known OS list:  Windows, Linux, Mac, BSD, Solaris
###############################################################################
function getOS() {
  case $(tolower $OSTYPE) in
    cygwin|msys)    KNOWN_OS_NAME=Windows  ;;
    linux*)         KNOWN_OS_NAME=Linux    ;;
    darwin*)        KNOWN_OS_NAME=Mac      ;; 
    bsd*|freebsd*)  KNOWN_OS_NAME=BSD      ;;
    solaris*)       KNOWN_OS_NAME=Solaris  ;;
    *)              KNOWN_OS_NAME=UNKNOWN  ;;
  esac
  function getOS() { echo "$KNOWN_OS_NAME"; } # re-define for faster performance
  getOS # call re-defined function
}

###############################################################################
##  Usage: if is64Bit; then ...   OR   is64Bit && ...
##  The exit code is 0 (success) if the current shell is a 64-bit process (as 
##    reported by the HOSTTYPE environtment variable). 
##  Otherwise the exit code is 1 (failure).
###############################################################################
function is64Bit() {  [[ $HOSTTYPE == *64 ]]; }

###############################################################################
##  Usage: BITS=$(getArchBits) 
##  Outputs '32' or '64' depending on the architecture of the current shell/ENV.
###############################################################################
function getArchBits() {
  [[ $HOSTTYPE == *64 ]] && echo 64 || echo 32 
}

###############################################################################
##  Usage: if isEmptyDir <DIR>; then ...
##  Determines whether or not the specified directory is empty. The exit code
##    is 0 if the directory is empty or 1 if it contains any files or sub-
##    directories (other than . and .. )
##  If the specified directory does not exist (or is a file etc), an error
##    message is displayed and the exit code is 2
###############################################################################
function isEmptyDir() {
    [[ -d "$1" ]] && (($(ls -fq1 "$1"|wc -l)<=2))  || { printf >&2 'isEmptyDir: "%s" is not a directory.' "$1"; exit 2; }
}

###############################################################################
##  Usage: count=$(fileCount <DIR>)
##  Outputs the number of items (files or sub-directories) contained in the 
##    the specified directory. It does not count (. or ..)
##  If the specified directory does not exist (or is a file etc), the output
##    is 0.
###############################################################################
function fileCount() {
    # Count files in the directory
    [[ -d "$1" ]] && echo $(($(ls -fq1 "$1"|wc -l)-2)) || echo 0
}

###################### Shell Options (shopt) ##################################
# tmpopt - Temporarily  set(-s) or unset(-u) a shell option, first storing its 
#          previous state so it can be reset later with `tmpopt -r OPTION`.
# Usage: tmpopt [-s|-u|-r] OPTION
# Notes: Before an option is set or unset its current value is stored in a
#        stack-line manner so it can be recurively reset. 
function tmpopt() {
    local action=s  # set by default
    local val stackName stack size newval    
    [[ $1 == -[rsu] ]] && { action=${1:1}; shift; } # now action= s | u | r
    (($#)) || return 0 # no options to set
    for name in "$@"; do
        [[ -z $name ]] && continue
        shopt -q $name && val=s || val=u  # Get the current value of the option
        stackName=PREV_OPTS_$name   # The name of the stack variable containing previous values e.g. PREV_OPTS_extglob
        stack=${!stackName}         # get the current stack 
        size=${#stack}              # Get the size of the stack
        newval=$action              # The new shopt value (if action == r, newval will be changed to u or s)

        # The array does not exist, so create it (if $action is 'set' or 'unset'. Show an error if $action is 'reset')
        if [[ $action == r ]]; then
            # Action is RESET
            if ((size==0)); then
                echo >&2 "Error resetting option '$name': no previous value was saved."
                newval=u      # unset the value
            else
                newval=${stack:$size-1}  # Get the last value in the previous-value string
                stack=${stack:0:$size-1} # Remove the last character in the string
            fi
            #dbg: echo "Resetting $name to -$newval"
            shopt -$newval "$name"
        else
            # Action is SET or UNSET
            stack+=$val   # Add existing value to the end of the previous-value list
            # newval is not relevant. all options will be set together at the end
            #dbg: echo "$name: old=$val new=$action"
        fi
        export -n $stackName="$stack"
        # echo -e "$stackName=\t$stack"
    done

    if [[ $action == [us] ]]; then
        #dbg: echo "shopt -$action $@"
        shopt -$action "$@"
    fi
}
function usopt() { tmpopt -u "$@"; }
function rsopt() { tmpopt -r "$@"; }
function qopt() { for name in "$@"; do echo -n $name ''; shopt -q $name && echo on || echo off; done; }

###############################################################################
# Returns 0 (success) if any of the arguments exists as files or directories 
#  or 1 (failure) otherwise (i.e. none exist)
# Usage anyExist  [-f] [-d] [-p] [-e] [--] [FILE|GLOB]... 
# Flags: 
#   -f   -   Only match files.
#   -d   -   Only match directories.
#   -p   -   print the first match
#   -e   -   expand the glob(s)
#   --   -   interpret remaining arguments as FILES (i.e. not flags)
# Remarks: 
#  Make sure you quote the glob so it isn't expanded prematurely by the shell
function anyExist() {
    local files=() i f d p e
    for arg; do
        ((i++))
        case "$arg" in
          -[fdpe]) declare ${arg:1}=1 ;;
          --)      files+=("${@:i+1}"); break ;;
          *)       files+=("$arg") ;;
        esac
    done

    ((${#files[@]}==0)) && { echo >&2 "in $FUNCNAME: No files were specified"; return 1; }
    ((e)) && files=( ${files[@]} )
    
    for file in "${files[@]}"; do
        if [[ ( $d && -d $file ) || ( $f && -f $file ) || ( ! ( $f || $d ) && -e $file )  ]]; then
            ((p)) && echo "$file"
            return 0
        fi
    done
    #find ${files[@]} -maxdepth 0 $flags $actions -quit 2>/dev/null
    return 1
}

# Outputs the MD5 hash of a file.
# Usage: result=$(getMD5 <FILE>)
# Returns 0 on success or 1 if no MD5 tool could be found.
function getMD5() {
    if isExe md5; then
        md5 -q "$1"
    elif isExe md5sum; then
        local out=$(md5sum "$1")
        echo ${out% *}
    else
        echo >&2 "No MD5 utility found."
        return 1;
    fi
}