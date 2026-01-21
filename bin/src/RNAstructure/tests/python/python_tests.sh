# This runs the python tests in both python2 and python3 (if available)

# cause this script to exit with an error code if ANY command fails.
set -eo pipefail

SELF="${BASH_SOURCE[0]}"; SELF_DIR=$(dirname "$SELF") # RNAstructure/tests/python

# Set the current working directory to RNAstructure/tests
cd "$SELF_DIR/.." # RNAstructure/tests

# Include RNAstructure/exe in PYTHONPATH, so it can import RNAstructure.py (and _RNAstrcture_wrap.so)
export PYTHONPATH="$PYTHONPATH":../exe
export DATAPATH="${DATAPATH:-../data_tables}"

# Array to hold versions of python that have been found.
versions=()

# determine whether a specific version was already found
function found_version() {
    local len=${#versions[@]} # number of elements in versions array
    for (( i=0; i<len; i++ )); do
        if [[ "$1" == "${versions[i]}" ]]; then return 0; fi
    done
    return 1 # not found
}

# Test whether a command is available on the PATH
function is_cmd() {
    type -p $1 &>/dev/null
}

PYTHON_LIST=( python python2 python3 ) # list of possible python binaries to test.

for py in ${PYTHON_LIST[@]}; do
    # if this python binary exists, test it
    if is_cmd $py; then
        # Get the version number of this python binary
        ver=$( $py --version 2>&1)
        
        # If we already tested this version, skip it
        # (for example, if `python` is the same as `python2` etc)
        found_version "$ver" && continue
        
        versions+=("$ver")
        echo "Building RNAstructure-Python Interface for $ver"
        make -C .. python_interface PYTHON=$py
        $py ./python/python_tests.py
    else # this python binary does not exist
        echo "python_test: Skipping $py -- command not found."
    fi
done

# If at least one python was found, exit successfully
if ((${#versions[@]}>0)); then exit 0; fi

echo >&2 "No version of python was found on the PATH"
exit 1
