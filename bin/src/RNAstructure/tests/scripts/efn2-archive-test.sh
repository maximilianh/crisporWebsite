#!/bin/bash

set -euo pipefail # exit on all errors and return the exit code of any failing program in a pipeline

# Exclude files with known errors.
EXCLUDE=(
    M.tub 
    C.lim 
    lyc-es-c 
    zea-ma-a 
    Esch.coli._AE014075_1-363-brokenPK 
    met-fer 
    zea-ma-h 
    Esch.coli._AE014075_1-363 
    tet-the 
    bac-bre 
    yar-li-a
)

INIT_DIR="$PWD"
SELF="$BASH_SOURCE"
SELF_DIR=$(dirname "$SELF")             # <ROOT>/RNAstructure/tests/scripts
cd "$SELF_DIR"/../..                    # <ROOT>/RNAstructure
ARCHIVE_DIR=../archive                  # <ROOT>/archive
REF_DIR=tests/efn2/archive              # <ROOT>/RNAstructure/tests/efn2/archive
OUTPUT_DIR=tests/efn2_OUTPUT/archive    # <ROOT>/RNAstructure/tests/efn2_OUTPUT/archive
EXE=exe/efn2                            # <ROOT>/RNAstructure/exe/efn2
export DATAPATH=data_tables             # <ROOT>/RNAstructure/data_tables
ZIPNAME=+efn2_archive.tar.lzma         # Name of compressed archive file (the '+' prefix is to keep it at the top of directory listings)

if [[ ! -d $ARCHIVE_DIR ]] ; then
  ARCHIVE_DIR=../rna-archive
  [[ -d $ARCHIVE_DIR ]] || { echo >&2 'RNA Archive Directory not found!'; exit 1; }
fi

#echo SELF=$(realpath "$SELF")
#echo SELF_DIR=$(realpath  "$SELF_DIR")
#echo ARCHIVE_DIR=$(realpath  "$ARCHIVE_DIR")
#echo REF_DIR=$(realpath -m "$REF_DIR")
#echo OUTPUT_DIR=$(realpath -m "$OUTPUT_DIR")
#echo EXE=$(realpath  "$EXE")

suffix=_output.txt
ref_suffix=_OK.txt

isREF= ; skipDone=1 ; keepOutput= ; operation=run
for arg; do 
    case "$arg" in
        -r) isREF=1; echo "Building reference output." ;;
        -w) skipDone=; echo "Overwriting existing output." ;;
        -k) keepOutput=1; echo "Keeping all output files." ;;
        -z)  operation=zip              ;;
        -zx) operation=unzip            ;;
        -rmf) operation=removeOutput    ;;
        *) echo >&2 "Unknown command-line flag: $arg"; exit 1 ;;
    esac
done

#setup
[[ $isREF ]] && { OUTPUT_DIR="$REF_DIR"; suffix=$ref_suffix; }  # output files into the reference directory if user passed -r flag

errexit=0
case $operation in 
    zip)    echo "Building reference archive: $ZIPNAME"
            rm -f $ZIPNAME
            cd $REF_DIR && tar --lzma -cf $ZIPNAME *$ref_suffix || errexit=$?
            ;;
    unzip)  echo 'Extracting reference output.'
            cd $REF_DIR && tar --lzma -xf $ZIPNAME || errexit=$?
            ;;
    removeOutput)
            [[ $isREF ]] && type='reference' || type='program'
            printf 'Removing all %s output from %s.\n' $type "$OUTPUT_DIR"
            rm -f "$OUTPUT_DIR"/*$suffix
            ;;
esac

[[ $errexit -eq 0 ]] || { echo "Failed with error code $errexit"; exit $errexit;  }
[[ $operation == run ]] || { echo 'done.';  exit; }

mkdir -p "$OUTPUT_DIR" # program will exit if this fails.

#clean-up last output
[[ $skipDone ]] || rm -f "$OUTPUT_DIR"/*$suffix

#trap 'trap - INT; kill -INT $$' INT CHLD # if the user presses Ctrl+C, exit the loop. (Otherwise sometimes efn2 will catch the SIGNAL and the loop will NOT exit)

# If we are doing a standard test (NOT generating reference files) and no reference files (*_OK.txt) exist, but the reference zip-archive does exist, then extract it.
# Probably works, but untested: [[ ! $isREF && -f $REF_DIR/$ZIPNAME ]] && ! ls >&/dev/null $REF_DIR/*$ref_suffix && "$INIT_DIR/$SELF" -zx

count=0
list=( "$ARCHIVE_DIR"/*.ct )
total=${#list[@]}
echo "Number of archive structures: $total"
for f in "${list[@]}"; do
    : $((count++))
    name="${f##*/}"; name="${name%.ct}" # remove directory and extension
    
    # skip files excluded due to known errors
    for exc in "${EXCLUDE[@]}"; do [[ $exc == "$name" ]] && { echo "Skipped $name (excluded)"; continue 2; }; done
    
    #Generate names of output and error files.
    output="$OUTPUT_DIR"/"$name"$suffix
    errout="$OUTPUT_DIR"/"$name"_errors.txt

    [[ -s "$output" && $skipDone ]] && { echo "Skipped $name (already exists)"; continue; }

    echo "$name ($count of $total)"
    exit=0
    "$EXE" "$f" "$output" -w  2>"$errout" 1>/dev/null || exit=$?
    [[ $exit -ne 0 ]] && echo >&2 "FAILED $name ($exit)";
    if [[ -s "$errout" ]]; then
        echo >&2 "ERRORS: $name";
        cat >&2 "$errout"
    else
        rm "$errout" # remove empty file
    fi

    # If this is NOT a reference run, then do a diff comparison
    if [[ ! $isREF ]]; then
        ref="$REF_DIR"/"$name"$ref_suffix
        diffout="$OUTPUT_DIR"/"$name"_diff.txt
        diff -b "$ref" "$output" >& "$diffout" || echo "(diff exit $?)";
        if [[ -s "$diffout" ]]; then 
            # if the diff output is not empty, show an error message.
            echo >&2 "CHANGED: $name";
        else
            rm "$diffout" # remove empty file
            [[ $keepOutput ]] || rm -f "$output"
        fi
    fi
done

# function postProcessOutput() {
#     sed -i.tmp -e 's/Stack = //' -e 's/ for \(stack\|closure\) of //' -e 's/^\t//' -e 's/Helix total = /\t==/' "$1"
# }