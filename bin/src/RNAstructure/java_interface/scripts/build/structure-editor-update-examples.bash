####################################################
# This script copies example files from RNAstructure/Examples 
# to the resources/examples folder in the StructureEditor 
# program.
####################################################
####### Configuration ##############################

# Path of the root java folder (e.g. RNAstructure/java_interface) 
# relative to THIS script's location
JAVA_ROOT="../.."  

# Path of source directory relative to JAVA_ROOT
SRC="../examples"

# Path of local examples directory relative to JAVA_ROOT.
DST="resources/StructureEditor/examples"

####################################################
function resolve() { (cd "$1"; echo "$PWD") }
SCRIPT=$(realpath "${BASH_SOURCE[0]}")   # RNAstructure/java_interface/scripts/build/script-name.bash
SCRIPT_DIR=$(dirname "$SCRIPT")          # RNAstructure/java_interface/scripts/build

ROOT=$(resolve "$SCRIPT_DIR/$JAVA_ROOT") # RNAstructure/java_interface
SRC=$(resolve  "$ROOT/$SRC")
DST=$(resolve  "$ROOT/$DST")

echo -e "SRC:$SRC\nDST:$DST"

copy_files=(
######## LIST EXAMPLE FILES TO COPY #########
# Copy these files from RNAstructure/examples
bmorivector.ct
bmorivector.dot
ca5s.ct
ca5s.pfs
ec16s.fasta
P546.fasta
P546.shape
########## END OF FILES ############
)
#rm "$ROOT/$DST"/*.*
cd "$SRC"
cp -f "${copy_files[@]}" "$DST"
find "$DST" -type f -printf "%f\n" > "$DST/examples.dir"
echo
echo "Examples Directory List:"
cat "$DST/examples.dir"
echo 
read -p "Done. Press a key to quit." 