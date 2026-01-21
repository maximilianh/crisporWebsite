#!/bin/bash
set -eo pipefail # set options for better error detection
trap 'rc=$?; echo >&2 "ERROR at line ${LINENO} (rc: $rc)"; exit $rc' ERR
trap 'rc=$?; (($rc)) && { echo >&2 "EXITING DUE TO ERROR ($rc)"; exit $rc; } || echo Done.' EXIT

# This script accepts the current test output as the reference "OK" files.
# This should only be run after carefully examining the output to ensure that it is correct.
SELF="$BASH_SOURCE"
SCRIPT_DIR=$(dirname "$SELF")
SCRIPT_NAME=$(basename "$SELF")
cd "$SCRIPT_DIR"/.. # cd to RNAstructure/tests folder

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

EM="$FMT_CYN"
DBG="$FMT_BLD"
F0="$FMT_RST"
ERR="$FMT_YLW"
WARN="$FMT_YLW"
NL=$'\n' #newline

function dbg() {
	echo "${DBG}DEBUG: ${*//"$F0"/"$F0$DBG"}$F0"
}
function warn() {
	echo >&2 "${WARN}WARN: ${*//"$F0"/"$F0$WARN"}$F0"
}

function die() {
	echo >&2 "${ERR}ERROR: ${*//"$F0"/"$F0$ERR"}$F0"
	return 1
}

handled_output=()
accepted_output=()
error_output=()

# read through log of tests commands and pick those involving DIFF
function parseLogs() {
	local output ref diff testname list item
	local prog="$1"
	list=()
	#               diff   (  <FLAGS>  )* (FILE1 )   (FILE2 )   '>&'   (DIFF_FILE) .*
	perl_regex='^\s*diff\s+(?:-[^\s]+|\s)*([^\s]+)\s+([^\s]+)\s*\>\&\s*([^\s]+).*'
	for log in ${prog}_OUTPUT/*_commands.log; do
		#echo "File: $log"
		while read ref output diff; do
			dbg "OUT=$output  REF=$ref  DIFF=$diff"
			testname=${diff#"$prog"_}  # remove program name
			testname=${testname%_diff.txt}  # remove _diff.txt suffix
			dbg "testname=$testname"
			list+=("${prog}_OUTPUT/$output $ref $testname")
		done < <(perl -n -e 's/'"$perl_regex"'/\1 \2 \3/   and print $_' $log)
	done

	for item in "${list[@]}"; do
		handleTest $item || { dbg "handleTest Failed for $item"; return 1; }
	done
}


# loop through *_output.* files in the test_OUTPUT folder
function loopOutputs() {
	local prog=$1 progbase=$2 refdir=$3
	local file output refdir ref testname
	for output in ${prog}_OUTPUT/*_output.*; do
		[[ -f $output ]] || continue # make sure it exists (i.e. not just literal '<prog>_OUTPUT/*_output.*')
		file=${output##*/}      # remove directory

		if [[ -z $progbase ]]; then
			progbase=${prog%-smp}     # get the "base-name" of the program e.g. Fold-smp -> Fold
			progbase=${progbase%-gui}    #    .. continued
		fi

		if [[ -z $refdir ]]; then
			refdir=$progbase
		fi

		ref=${file/#"$prog"_/"$progbase"_} # replace prog with progbase in the name of the OK file
		ref=${file/_output./_OK.} # replace _output with _OK
		ref=${refdir}/$ref # prepend directory (progbase)

		testname=${file#"$prog"_}  # remove program name
		testname=${testname%_output.*}  # remove _output.ext suffix

		handleTest $output $ref $testname || return 1
	done
}

# Choose whether or not to accept a given test output as the reference file
function handleTest() {
	local output=$1 ref=$2 testname=$3
	local diffresult answer

	# See if we've already handled this output. If so, just return. Otherwise add it.
	for handled in "${handled_output[@]}"; do [[ $handled == $output ]] && return; done
	handled_output+=($output)

	#dbg "handleTest $output $ref $testname"

	if [[ -e $ref ]]; then
		#dbg "Found OK-ref: $ref"
		diffresult=`diff --strip-trailing-cr  "$ref"  "$output" |& head -n 10`
	else
		diffresult="(Current OK file missing)"
	fi

	if [[ -z $diffresult ]]; then
		echo "Skipping $EM$testname$F0 -- matches."
		return
	fi

	if [[ ! -e $output ]]; then
		error_output+=("$testname: Missing output: '$output'")
		warn "Expected output for $EM$testname$F0 is missing!${NL}  Path: '$output'"
		return
	fi

	echo "Accept result for $EM$testname$F0 and copy '$output'  to  '$ref'  ?"
	echo "$FMT_PRP$FMT_BLD$diffresult$F0"

	while read -p 'Accept? [y,n,d,o,r,?]: ' answer; do
		case "$answer" in
			[Dd]) 	# view full diff
				  	diff --strip-trailing-cr  $ref  $output | less
				  	;;
			[Oo]) 	# view full output
				  	less $output
				  	;;
			[Rr]) 	# view full OK file
				  	less $ref
				  	;;
			[Yy]) 	# Accept the file
					if ERR_MSG=$(cp $output $ref 2>&1); then
						accepted_output+=($output $ref)
					else
						warn "Failed to copy  $output  to  $ref: $ERR_MSG"
						error_output+=("$testname: Failed to copy '$output' to '$ref'")
						return 0
						#return 1
					fi
					break
				  	;;
			[NnSs])	# Skip (do NOT accept) the file
				  	break
				  	;;
			[Qq]) 	# Quit
				  	echo "Test acceptance canceled by user."
					return 1
				  	;;
		    ''|'?'|*) # show help
					echo $'\t d=show full diff\n\t o=show output\n\t r=show ref OK\n\t y=Accept\n\t n=Skip (default) i.e. do NOT accept\n\t q=Quit'
				 	;;
		esac
	done
}

function acceptProg() {
	# parseLogs -- read through log of tests commands and pick those involving DIFF
	# loopOutputs -- loop through *_output.* files in the test_OUTPUT folder
	# If either returns an error code, then quit.
	{ parseLogs $1 || die "Parse-method FAILED for $1"; } && { loopOutputs $1 || die "Loop-method FAILED for $1"; }
}

function usage() {
	echo >&2 "Usage: $SCRIPT_NAME all | <PROGRAMS>..."
}

###############################################
##### Actual program entry point is here: #####
###############################################
progs=("$@")
if [[ ! $1 ]]; then
	echo >&2 "Missing required PROGRAM name!"; 	usage; exit 1
elif [[ $1 == all ]]; then
	progs=( *_OUTPUT/ )
fi

# If the user passed in a full directory path, remove the final slash, the parent directory, and the '_OUTPUT' suffix.
progs=( "${progs[@]%\/}" ) #remove final slash
progs=( "${progs[@]##*\/}" ) #remove parent directory
progs=( "${progs[@]%_OUTPUT}" ) #remove _OUTPUT suffix

foundDirs=0
for prog in "${progs[@]}"; do
	if [[ -d ${prog}_OUTPUT ]]; then
		echo "Analyzing output for program: $prog"
		acceptProg $prog
		((++foundDirs))
	else
		echo >&2 "Program output directory not found:  ${prog}_OUTPUT"
	fi
done

if [[ $foundDirs -eq 0 ]]; then
	# Show usage information if none of the programs were found.
	usage; exit 1
else
	if (( ${#accepted_output[@]} )); then
		echo "Accepted the Following: "
		printf "\t%s --> %s\n" "${accepted_output[@]}"
	else
		echo "NO output files were accepted."
	fi

	if (( ${#error_output[@]} )); then
		warn "The Following Errors Occurred: "
		printf "\t$WARN%s$F0\n" "${error_output[@]}"
	else
		echo "(No errors occurred.)"
	fi
fi
