# Some of the test input files are generated from other input files.
# In order to improve the speed and reliability of regression tests, 
# These files are NOT generated on-the-fly. Instead this script can be
# run to (re-)generate the files only when necessary.

cd $(dirname $BASH_SOURCE)/../.. # cd to RNAstructure folder
make Fold partition bifold ProbabilityPlot dynalign || exit 1

cd tests # # cd to RNAstructure/tests folder
. scripts/test-env.sh #include definitions of $SINGLESEQ etc
PATH+=:../exe # add RNAstructure/exe to the PATH


bifold $DOUBLESEQ $BIMOLCT --intramolecular
Fold $SINGLESEQ $SINGLECT -s $SINGLESAV
partition $SINGLESEQ $SINGLEPFS
partition $SINGLESEQ2 $SINGLEPFS2
partition $SINGLESEQ3 $SINGLEPFS3
partition $SINGLESEQ4 $SINGLEPFS4

ProbabilityPlot $SINGLEPFS $SINGLEPFS.txt --text
ProbabilityPlot $SINGLEPFS3 $SINGLEPFS3.txt --text



# Run dynalign to prepare the input file for DynalignDotPlot
echo > $DYNDOTPLOT~.conf \
"inseq1 = testFiles/testFile_RD0260.seq
inseq2 = testFiles/testFile_RD0500.seq
outct = $DYNDOTPLOT~1.ct
outct2 = $DYNDOTPLOT~2.ct
aout = $DYNDOTPLOT~.ali
savefile = $DYNDOTPLOT.dsv"

dynalign $DYNDOTPLOT~.conf
rm $DYNDOTPLOT~* # remove temporary output files.

echo "Done rebuilding test input files."