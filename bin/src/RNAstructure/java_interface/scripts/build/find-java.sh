#!/bin/bash
################################################################################
# This script attempts to locate the following Java-related directories:
#
#   JDK_BIN - Usually called "bin", this directory contains the Java compiler 
#       (`javac`) and archiver (`jar`).
#   JDK_INCLUDE - Usually called "include", this directory contains c++ header 
#       files required to compile Java Native Interface (JNI) libraries.
#   JAVA_HOME - On some systems JDK_BIN and JDK_INCLUDE are both contained in a
#       single parent directory (i.e. the Java installation directory).
#       (However this is not always the case.)
#   JDK_HOME - A synonymn for JAVA_HOME for systems in which JAVA_HOME is 
#       set by a third-party or refers to a version of Java that is not 
#       suitable for use by RNAstructure.
#
# Usage:   find-java.sh <TARGET> [<GUESS>]  [OPTIONS]
#   <TARGET> The target directory to find: 'home', 'bin' or 'inc'.
#   <GUESS> (optional) A potential location for TARGET. If it exists and points
#           to a suitable version of Java, it is returned (without any 
#           additional searching.)
#
# See the "Usage Information" section below or run the script with the -h flag
# to see additional usage information.
################################################################################

SCRIPT="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT")
SCRIPT_NAME=$(basename "$SCRIPT")

JAVAC_MIN_VERSION_NAME="1.8 (Java 8)" # user-friendly description of required Java version
JAVAC_MIN_VERSION=8 # An integer representing the major JDK version required.

# Store the results of the search.
FOUND_HOME=
FOUND_BIN=
FOUND_JNI=
JAVAC_VERSION_TEXT=
JAVAC_VERSION_NUM=

#set some globals
FIND_VERBOSE=
FIND_DEBUG=
unset VERIFY_CACHE VERIFY_CACHE_ASSOC; 
# create VERIFY_CACHE as an associative array if they are supported by this version of bash. Otherwise, use a standard array.
declare -A VERIFY_CACHE >&/dev/null && VERIFY_CACHE_ASSOC=1 || VERIFY_CACHE=() 


############## Show Usage information ##########################################
function showUsage() {
  # Print intro if specified
  [[ $1 == -i ]] && printf >&2 "\n%s -- Attempts to locate various Java-related directories.\n\n" "$SCRIPT_NAME"

  cat >&2 <<EOF
Usage:   find-java.sh  <TARGET>  [<GUESS>]  [OPTIONS]
  <TARGET> The target directory to find: 'home', 'bin' or 'inc'.
  <GUESS> (optional) A potential location for TARGET. If it exists and points
          to a suitable version of Java, it is returned (without any 
          additional searching.)
EOF

  # If additional information is desired, show the extended usage information below.
  [[ $1 == -i ]] && cat >&2 <<EOF

TARGET:
  home -- The JDK parent directory that contains 'bin' and 'include' sub-dirs.
  bin  -- The JDK 'bin' directory that contains javac, jar etc.
  inc  -- The JDK 'include' directory that contains the JNI header files.
              (jni.h, <OS>/jni_md.h)
  jni  -- (alias for 'inc')

OPTIONS:  (These must be specified individually, e.g. -d -v NOT -dv)
  -v  Verbose mode -- output messages (to stderr) during the process to inform 
      the user that auto-detection is in progress and to let them know how they 
      can circumvent detection. Also informs them of the result (assuming that 
      stdout is being redirected to a variable and will NOT be seen by users).
  -d  Output information useful for debugging (but not for end-users).
  -h  Show this help message.
EOF
  [[ $1 == -i ]] || printf >&2 "Run '%s -h' to show additional help and usage information.\n\n" "$SCRIPT_NAME"
}

############## Main entry-point of program #####################################
function startup() {
    detectOS

    # Parse script arguments and perform the find operation.
    local arg target guess pos=1
    for arg; do
        case "$arg" in 
            -v) FIND_VERBOSE=1	;;
            -d) FIND_DEBUG=1	;;
            *) case $((pos++)) in
                1) target="$arg" ;;
                2) guess="$arg"  ;;
                *) errmsg "Unknown argument: $arg"; showUsage; return 1 ;;
               esac ;;
        esac
    done

    case "$target" in
        home) 			findHome "$guess" && echo "$FOUND_HOME";;
        bin)  			findBin  "$guess" && echo "$FOUND_BIN";;
        jni|inc)		findJni  "$guess" && echo "$FOUND_JNI";;
        -h|--help|help) showUsage -i	;;
        # Show error for any unknown targets
        *)    			printf >&2 "Invalid target: %q\n" "$1"; showUsage ;;
    esac
}

function detectOS() {
    # OSTYPE is almost always defined by the shell
    # But if it isn't, use uname instead.
    [[ $OSTYPE ]] || OSTYPE=$(uname -s)

    #     OS     - short name of Operating System
    #     OSNAME - descriptive name of OS
    #     OS_DIR - name of jni include subfolder
    case $(tolower $OSTYPE) in
        cygwin|msys|CYGWIN*|MSYS*|MINGW*|Windows_NT*) #the upper-case names come from uname
            OS=win; OS_DIR=win32;  OSNAME=Windows ;;
        GNU*|[Ll]inux*)
            OS=lin; OS_DIR=linux;  OSNAME=Linux   ;;
        [Dd]arwin*)
            OS=mac; OS_DIR=darwin; OSNAME=MacOSX  ;;
        *) errmsg "Unkown operating system: $OSTYPE"; return 1 ;;
    esac	
}

################################################################################
# General-use functions
################################################################################
# remove whitespace from before and after the argument(s)
function trim() {
    local var="$*"
    # remove leading whitespace characters
    var="${var#"${var%%[![:space:]]*}"}"
    # remove trailing whitespace characters
    var="${var%"${var##*[![:space:]]}"}"   
    echo -n "$var"
}

# Show a message on stderr
function errmsg() {
    echo -e >&2 "$@"
    return 0
}

# Show a message (on stderr), but only if the verbose flag (-v) was given on the command-line.
function verboseMsg() {
    [[ $FIND_VERBOSE ]] && echo -e >&2 "$@"	
    return 0
}

function debug() {
    if [[ $FIND_DEBUG ]]; then
        local fn; [[ $1 == -f ]] && { fn="$(fn 1): "; shift; }
        local i=${#FUNCNAME[@]}
        local indent
        while ((i-- > 4)); do indent+='.'; done  # level 4 comes from: *main* > startup > TRUE_FUNCTION > debug
         #echo "${FUNCNAME[*]}"
         echo -e >&2 "$indent$fn$*"
    fi
    return 0
}

# outputs the path in canonical form (e.g. UNIX form if on windows)
# returns 0 if path exists or 1 otherwise.
function resolve() {
    local verbose; [[ $1 == -v ]] && { verbose=1; shift; }
    if [[ -d $1 ]] && cd "$1" 2>/dev/null && echo "$PWD"; then
        [[ $verbose && $PWD != $1 ]] && debug -f "\"$1\" = \"$PWD\""
        return 0
    else
        [[ $verbose ]] && debug -f "Directory \"$1\" does not exist."
        return 1
    fi
}

# fn [<LEVEL:=0>] -- outputs the name of the calling function (with no arguments) or the name of the function LEVEL steps back in the call stack.
function fn() {
    local i=${1:-0}
    echo "${FUNCNAME[i+1]}"
}

# Convert arguments to lower case (This is required because bash < v4 lacks the ${VAR,,} or ${VAR^^} syntax.)
function tolower() { 
    # Re-define the function on first use, depending on bash version.
    if ((BASH_VERSINFO[0] < 4)); then
        function tolower() { echo "$*" | tr '[:upper:]' '[:lower:]'; }
    else
        function tolower() { echo "${*,,}"; }
    fi
    tolower "$*" # call the newly-defined version of the function
}

function isCmd() {
    type -p >&/dev/null #return the result of this test.
}

########## Some functions to cache results of verifications, because some 
# paths will be discovered by multiple mechanisms and will be
# verified multiple times. If the result is cached, we can skip
# re-verification, which could be time consuming.

# Determine if the specified path has already failed
# Usage: prevFail <PATH> <TYPE> && ...
function prevFail() { prevResult "$1" "$2" F ; }  # F means Failed
function prevGood() { prevResult "$1" "$2" S ; }  # S means Success
# Determine if the specified path has already been tested.
# Usage: prevResult <PATH> <TARGET> <TEST_FOR_RESULT: T|S> && ...
function prevResult() {
    local key="$1@$2?$3"
    [[ $OS == win ]] && key=$(tolower "$key")  # do case insensitive comparison on Windows

    if [[ $VERIFY_CACHE_ASSOC ]]; then
        [[ ${VERIFY_CACHE[$key]} ]] && return 0 # return true if the key exists
    else
        for item in "${VERIFY_CACHE[@]}"; do
            [[ $item == $key ]] && return 0
        done
    fi
    return 1
}

function saveFail() { saveResult "$1" "$2" F ; }
function saveGood() { saveResult "$1" "$2" S ; }
# Usage: saveResult <PATH> <TARGET> <RESULT: T|S>
function saveResult() {
    local key="$1@$2?$3"
    [[ $OS == win ]] && key=$(tolower "$key")  # do case insensitive comparison on Windows
    if [[ $VERIFY_CACHE_ASSOC ]]; then
        VERIFY_CACHE[$key]=1
    else
        VERIFY_CACHE+=("$key")
    fi
}

################################################################################
# Java finding functions
################################################################################

function findBin() {
    local guess="$1"
    [[ $guess ]]     && { verifyBin "$guess"   			&& return 0; }
    [[ $JDK_BIN ]]   && { verifyBin "$JDK_BIN" 			&& return 0 || verboseMsg "JDK_BIN is set but it is invalid."; }
    [[ $JAVA_HOME ]] && { verifyBin "$JAVA_HOME/bin" 	&& return 0 || verboseMsg "JAVA_HOME is set but it does not point to a valid JDK."; }
    [[ $JDK_HOME ]]  && { verifyBin "$JDK_HOME/bin" 	&& return 0 || verboseMsg "JDK_HOME is set but it does not point to a valid JDK."; }

    verboseMsg "Detecting JDK 'bin' directory (set JDK_BIN to skip this step)."

    if findHome '' -s; then 
        verifyBin "$FOUND_HOME/bin" # sets FOUND_BIN if successful
    else
        verboseMsg "Failed to auto-detect JDK directory."
    fi

    # TODO: try other methods here, especially on Mac
    [[ $FOUND_BIN ]] && verboseMsg "Found JDK 'bin' directory: JDK_BIN=$FOUND_BIN" && return 0;

    showVersionError ||	errmsg "Failed to auto-detect the JDK 'bin' folder."
    return 1
}

# Returns true if Java was found, but is the wrong version. 
# Otherwise (Java was NOT found), this returns false (1)
function showVersionError() {
    if [[ $JAVAC_VERSION_TEXT ]]; then
        errmsg "The JDK was found, but it is the wrong version ($JAVAC_VERSION_TEXT). Version $JAVAC_MIN_VERSION_NAME or later is required."
        if [[ $JAVAC_VERSION_NUM == 7 ]]; then
            errmsg "Official public support for Java 7 (aka 1.7) ended in April 2015."
        fi
        return 0
    fi
    return 1
}

# Verifies that the directory exists, contains javac, and verifies that javac is the correct version
# Sets FOUND_BIN and returns true (0) if it is valid.
# Otherwise returns 0 (with no stdout)
function verifyBin() {
    debug -f "testing $1"
    local dir=$(resolve -v "$1") || return 1  # get cannonical form of directory
    local cached=1
     # first test for cached results for this same dir.
    prevFail "$dir" bin && debug "  FAILED (cached)" && return 1
    prevGood "$dir" bin || 
        { cached=; [[ $dir && -d "$dir" && -f "$dir/javac" ]] && verifyJavacVersion "$dir/javac" ; } ||
        { debug "  FAILED"; saveFail "$dir" bin; return 1; }

    debug "  SUCCESS${cached:+ (cached)}"
    [[ $cached ]] || saveGood "$dir" bin  # cache result in case it is repeated (verifyHome calls verifyBin)
    FOUND_BIN="$dir"
}

# if $2 is '-s' (silent), all output is suppressed and only the exit code indicates the result.
function verifyJavacVersion() {
    local path="$1"
    prevFail "$path" jxx && debug -f "FAILED \"$path\" (cached)"  && return 1
    prevGood "$path" jxx && debug -f "SUCCESS \"$path\" (cached)" && return 0

    debug -f "testing $1"
    local out=$("$path" -version 2>&1) || { debug "  javac exit: ${PIPESTATUS[0]}"; return 1; } # return if the javac command fails
    debug "  javac output: $out"
    out=$(trim "${out#javac}")

    if [[ $out ]]; then
        JAVAC_VERSION_TEXT="$out"  # set a global that is shown in showJavacVersionMessage if no paths succeed.
        JAVAC_VERSION_NUM=$(parseJavaVersion "$out")
        debug "  parsed javac version: $JAVAC_VERSION_NUM"
        if [[ $JAVAC_VERSION_NUM -ge $JAVAC_MIN_VERSION ]]; then
            debug "  javac version confirmed: $JAVAC_VERSION_NUM >= $JAVAC_MIN_VERSION"; saveGood "$path" jxx; return 0;
        elif [[ $JAVAC_VERSION_NUM == [1-7] ]]; then 
            debug "  javac is outdated: $JAVAC_VERSION_TEXT"; saveFail "$path" jxx; return 1;
        fi
        # Otherwise, we couldn't parse the version number, so continue below to test javac compilation.
    else
        JAVAC_VERSION_TEXT=
        JAVAC_VERSION_NUM=
    fi
    
    # If execution gets here, the java version test was inconclusive.
    debug "  javac version unknown. Testing compiler.";
    # if it is not exactly the required version, it might still be valid. Test an actual comilation of a test.java file:
    "$path" -d /tmp "$SCRIPT_DIR/javac-test.java" >&/dev/null && {  debug "  javac Compilation SUCCESS"; saveGood "$path" jxx; return 0; }

    debug "  javac Compilation FAILED"
    saveFail "$path" jxx
    return 1
}

function parseJavaVersion() {
    # Parse a Java version string (e.g "1.8.0_161") and return an integer 
    # representation of the Java major version (e.g. 8)
    case "$1" in 
        # Java 1 to 8
        1.[1-9]*)        echo "${1:2:1}" ;; # output 3rd character
        # Java 9
        9[-_.' ']*)      echo 9 ;;
        # Java 10 to 19
        1[0-9][-_.' ']*) echo "${1:0:2}" ;; # output first two digits
        *) echo 0 ;; # represents invalid
    esac
}

# locates the JDK_INCLUDE directory that contains c++ header files 
#   named jni.h and <OSDIR>/jni_md.h
function findJni() {
    local guess="$1"
    [[ $guess ]]       && { verifyJNI "$guess"       		&& return 0; }
    [[ $JDK_INCLUDE ]] && { verifyJNI "$JDK_INCLUDE" 		&& return 0 || verboseMsg "JDK_INCLUDE is set but it is invalid."; }
    [[ $JAVA_HOME ]]   && { verifyJNI "$JAVA_HOME/include" 	&& return 0 || verboseMsg "JAVA_HOME is set but it does not point to a valid JDK."; }

    verboseMsg "Detecting JDK 'include' directory (set JDK_INCLUDE to skip this step)."

    if findHome '' -s; then 
        verifyJNI "$FOUND_HOME/include"  # sets FOUND_JNI if successful
    else
        verboseMsg "Failed to auto-detect JDK directory."
    fi

    [[ $FOUND_JNI ]] && verboseMsg "Found JDK 'include' directory: JDK_INCLUDE=$FOUND_JNI" && return 0;


    # TODO: try other methods here, especially on Mac
    errmsg "Failed to auto-detect the JDK 'include' folder."
    return 1
}
# Verifies that the directory exists, contains jni.h and <OSDIR>/jni_md.h
# Sets FOUND_JNI and returns true (0) if it is valid.
# Otherwise returns 0 (with no stdout)
function verifyJNI() {
    debug -f "testing $1"
    local dir=$(resolve -v "$1") || return 1  # get cannonical form of directory	
    local cached=1
     # first test for cached results for this same dir.
    prevFail "$dir" jni && debug "  FAILED (cached)" && return 1
    prevGood "$dir" jni || 
        { cached=; [[ $dir && -d "$dir" && -f "$dir/jni.h" && -f "$dir/$OS_DIR/jni_md.h" ]]; } || 
        { debug "  FAILED"; saveFail "$dir" jni; return 1; }
    debug "  SUCCESS${cached:+' (cached)'}"
    [[ $cached ]] || saveGood "$dir" jni  # cache result in case it is repeated (verifyHome calls this)
    FOUND_JNI="$dir"
}

# Find JAVA_HOME, assuming that it contains 'bin' and 'include' folders.
# If $2 == -s this function will be 'silent' -- i.e. verboseMsg will be suppressed.
function findHome() {
    local home  guess="$1" silent
    [[ $2 == -s ]] && silent=1

    # First see if JAVA_HOME is already set and valid. If so verify it and return.
    [[ $guess ]]     && verifyHome "$guess"     && return 0
    [[ $JAVA_HOME ]] && verifyHome "$JAVA_HOME" && return 0
    [[ $JDK_HOME ]]  && verifyHome "$JDK_HOME"  && return 0

    [[ $silent ]] || verboseMsg "Detecting Java JDK parent directory (set JAVA_HOME or JDK_HOME to skip this step)."

    case $OS in
        win) 
            # Try the registry first.
            winGetJavaHomeRegistry ||
            winGetJavaHomeSearch   ||
            getHomeFromJava
            ;;
        lin)
            getHomeFromJava
            ;;
        mac)
            getHomeFromJava
            ;;
    esac

    if [[ $FOUND_HOME ]]; then 
        [[ $silent ]] || verboseMsg "Found JDK parent directory: JDK_HOME=$FOUND_HOME"
        return 0
    fi
    
    [[ $silent ]] || showVersionError || errmsg "Failed to auto-detect JDK directory."
    return 1
}

function verifyHome() {
    debug -f "testing $1"
    local dir=$(resolve -v "$1") || return 1  # get cannonical form of directory
    local cached=1
     # first test for cached results for this same dir.
    prevFail "$dir" home && debug "  FAILED (cached)" && return 1
    prevGood "$dir" home || 
        { cached=; [[ $dir ]] && verifyBin "$dir/bin"  && verifyJNI "$dir/include"; } ||
        { debug "  FAILED"; saveFail "$dir" home; return 1; }
    debug "  SUCCESS${cached:+' (cached)'}"
    [[ $cached ]] || saveGood "$dir" home  # cache result in case it is repeated
    # Set a global (this isn't ideal, but it is necessary to preserve changes to other globals, including the results cache etc)
    FOUND_HOME="$dir"
}

# Locate JAVA_HOME by completely resolving the full path to the specified binary
# (e.g. java or javac), which is usually a symlink (and often requires several
# levels of link resolution)
function getHomeFromBinary() {
    local path link 
    debug -f "using: $1"
    path=$(type -p $1)
    [[ -z $path ]] && { debug "  $1 is invalid or inaccessible."; return 1; }

    debug "  path: $path"
    link=$(readlink -f "$path" 2>/dev/null) # -f is better than multiple calls to readlink, but it is not supported on Mac
    if [[ $link ]]; then
        path="$link"
        debug "  link-->$link"
    else
        while link=$(readlink "$path"); do
            path="$link"
            debug "  link->$path"
        done
    fi
    # path = 
    #    mac:  /System/Library/Frameworks/JavaVM.framework/Versions/Current/Commands/java
    #    lin:  /usr/lib/jvm/java-8-oracle/jre/bin/java
    #    win:  /cygdrive/c/Program Files/Java/jre1.8.0_121/bin/java.exe
    path=$(dirname "$path")
    # path = 
    #    mac:  /System/Library/Frameworks/JavaVM.framework/Versions/Current/Commands
    #    lin:  /usr/lib/jvm/java-8-oracle/jre/bin
    #    win:  /cygdrive/c/Program Files/Java/jre1.8.0_121/bin
    path="${path%[\\/]bin}"
    path="${path%[\\/]jre}"
    # path = 
    #    mac:  /System/Library/Frameworks/JavaVM.framework/Versions/Current
    #    lin:  /usr/lib/jvm/java-8-oracle
    #    win:  /cygdrive/c/Program Files/Java/jre1.8.0_121
    
    # On mac, using readlink does not find JAVA_HOME, but it does result in the directory
    # that contains the 'java_home' command, which outputs the true JAVA_HOME when invoked.
    [[ $OS == mac && -x $path/java_home ]] && path=$($path/java_home)
    verifyHome "$path" && return 0

    debug "  FAILED to find JAVA_HOME from binary $1"
}

# Find JAVA_HOME by running `java -showversion -verbose` and parsing
# the output to get the path to the main Java runtime library (rt.jar)
function getHomeFromJava() {
    local java rtlib jre home
    
    getHomeFromBinary javac && return 0
    getHomeFromBinary java  && return 0

    java=$(type -p java) || { verboseMsg "Java is not installed or is not available in the PATH."; return 1; }
    debug -f "querying $java for JAVA_HOME"

    # If the above failed, try to use the output from `java -showversion -verbose:class` to find the home directory.
    local rtlib=$(java -showversion -verbose:class 2>&1 | head -1)
    # rtlib = [Opened <PATH>]
    rtlib="${rtlib#\[Opened\ }"
    rtlib="${rtlib%\]}"
    # rtlib =
    #    mac: /Library/Java/JavaVirtualMachines/jdk1.8.0_77.jdk/Contents/Home/jre/lib/rt.jar
    #    lin: /usr/lib/jvm/java-8-oracle/jre/lib/rt.jar
    #    lin: /usr/lib/jvm/java-8-openjdk-amd64/jre/lib/rt.jar
    #    win: C:\Program Files\Java\jre1.8.0_51\lib\rt.jar
    #    win: C:\Program Files\Java\jdk1.8.0_101\jre\lib\rt.jar
    local jre="${rtlib%[\\/]lib[\\/]rt.jar}"
    # jre =
    #    mac: /Library/Java/JavaVirtualMachines/jdk1.8.0_77.jdk/Contents/Home/jre
    #    lin: /usr/lib/jvm/java-8-oracle/jre
    #    lin: /usr/lib/jvm/java-8-openjdk-amd64/jre
    #    win: C:\Program Files\Java\jre1.8.0_51
    #    win: C:\Program Files\Java\jdk1.8.0_101\jre
    local home="${jre%[\\/]jre}"
    # home = 
    #    mac: /Library/Java/JavaVirtualMachines/jdk1.8.0_77.jdk/Contents/Home
    #    lin: /usr/lib/jvm/java-8-oracle
    #    lin: /usr/lib/jvm/java-8-openjdk-amd64
    #    win: C:\Program Files\Java\jre1.8.0_51
    #    win: C:\Program Files\Java\jdk1.8.0_101
    verifyHome "$home" && return 0

    # failed
    debug "  query FAILED"
    return 1
}

#######################################################################################
# Windows - specific
#######################################################################################
# Sets a number of variables to point to windows directories.
# Each directory has a unix version (*_U) and a Windows version (*_W) suitable for 
#   windows programs. BOTH versions should contain only forward slashes and should NOT
#   contain spaces (which could mess with `make`).   
#   In addition there can be a *_RAW version which is a standard Windows fornat that 
#   contains spaces and forward slashes.
function winGetDirs() {
  # on Windows, there are many variants of the root drive, depending on the OS/shell (e.g. cygwin, MSys2)
  # $SYSTEMROOT is a native environment variable--something like C:\Windows
  # SYSROOT_U will be something like  "/c" or "/cygdrive/c"  
  SYSROOT_U=$(cd "$SYSTEMROOT/.." && pwd)
  # SYSROOT_W will be something like  "c:" (it must be used for windows programs that do not understand /c or /cygdrive/c )
  SYSROOT_W=$(dirname $SYSTEMROOT)
  SYSROOT_W="${SYSROOT_W//\\/\/}" # it will probably not contain slashes, but just in case.

  # Program Files directory
  PF_RAW="$PROGRAMFILES" # Windows form that includes spaces
  PF_W=$SYSROOT_W/Progra~1  # Windows form that avoids spaces
  PF_U=$SYSROOT_U/Progra~1  # use the ~1 form to avoid spaces.

  # Program Files (x86) directory (for 32-bit programs)
  PF86_RAW=$(printenv 'ProgramFiles(x86)')
  PF86_W=$SYSROOT_W/Progra~2
  PF86_U=$SYSROOT_U/Progra~2
}

# Query the Windows registry to get the current JDK version. 
# The perform a second query using this version number to retrieve 
# the location of JAVA_HOME from the Windows registry.
function winGetJavaHomeRegistry() {
    local REG_KEY='HKLM\SOFTWARE\JavaSoft\Java Development Kit'
    local regval=$(reg query "$REG_KEY" -v "CurrentVersion" 2>/dev/null)
    [[ -z $regval ]] && { debug -f "JDK not registered."; return 1; }
    # regval is something like:
    #
    #  HKEY_LOCAL_MACHINE\SOFTWARE\JavaSoft\Java Development Kit
    #     CurrentVersion    REG_SZ    1.8
    #
    regval="${regval##*REG_SZ}" #trim off everything up to and including REG_SZ
    local version=$(trim "$regval") # remove spaces
    debug -f "Current JDK Version: $version"
    
    # query for path to JAVA_HOME
    regval=$(reg query "$REG_KEY\\$version" -v "JavaHome" 2>/dev/null)
    [[ -z $regval ]] && { debug "  JDK home not registered."; return 1; }
    # regval is something like:
    #
    #  HKEY_LOCAL_MACHINE\SOFTWARE\JavaSoft\Java Development Kit\1.8
    #      JavaHome    REG_SZ    C:\Program Files\Java\jdk1.8.0_101
    regval="${regval##*REG_SZ}" #trim off everything up to and including REG_SZ
    local path=$(trim "$regval") # remove spaces
    debug "  JDK path from registry: $path"
    path=$(resolve -v "$path")

    # Now path should be something like "C:\Program Files\Java\jdk1.8.0_101"
    # Unfotunately this is NOT usuable by make (because done of the make functions support spaces in path names)
    # winGetDirs
    # if [[ $path == "$PF_RAW"* ]]; then
    # 	path=${path/#"$PF_RAW"/"$PF_U"}
    # elif [[ $path == "$PF86_RAW"* ]]; then
    # 	path=${path/#"$PF86_RAW"/"$PF86_U"}
    # fi
    # path=$(fixPath "$path")  # get path with \ replaced by /
      
      verifyHome "$path" && return 0
      #failed
      debug '  FAILED'
      return 1
  }

function winGetJavaHomeSearch() {
    # This function works by finding all the paths of the form "C:\Program Files\Java\*\include\jni.h"
    # the folders matching the * wildcard represent JAVA_HOME. 
    # This is applied to all dirs in $dirs, so both 32-bit and 64-bit program files are checked.
    local dir dirs=(
        #"$PF_U"
        "$PF86_U"
    )
    local file found=()
    for dir in "${dirs[@]}"; do
      # This will add the results of glob expansion to found.
      # However, if there are no matches, an invalid path (containing *) will be added.
      found+=($dir/Java/*/include/jni.h)
    done
    debug -f "Found ${#found[@]} possible JDK paths."
    for file in "${found[@]}"; do
        if [[ -f $file ]]; then
              dir="${file%'/include/jni.h'}"
            debug "  possible home: $dir"
              verifyHome "$dir" && return 0
        fi
    done

      debug '  FAILED'
    return 1 #failed
}


##################################################
# Run the script
startup "$@"
##################################################
