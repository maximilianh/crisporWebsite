###############################################################################
# Set variables for input files. 
# (Ideally these will all be moved into app-specific folders)
# This file is sourced by test-runner.sh
###############################################################################
[[ $INDIR ]] || INDIR=testFiles; export INDIR
PREFIX=testFiles/testFile_ 

export SINGLESEQ=${PREFIX}RA7680.seq
export SINGLESEQ_SHORT=${PREFIX}RA7680_short.fasta
export SINGLESEQ2=${PREFIX}met-vol.seq
export SINGLESEQ2_FASTA=${PREFIX}met-vol.fasta
export SINGLESEQ3=${PREFIX}ivslsu.seq
export SINGLESEQ4=${PREFIX}ca5s.seq
export SINGLESEQ5=${PREFIX}ec5s.seq
export DOUBLESEQ="$SINGLESEQ5 $SINGLESEQ4"
export SINGLEPFS=${PREFIX}RA7680.pfs
export SINGLEPFS2=${PREFIX}met-vol.pfs
export SINGLEPFS3=${PREFIX}ivslsu.pfs
export SINGLEPFS4=${PREFIX}ca5s.pfs
export KNOTSCT=${PREFIX}knotted.ct
export SINGLECT=${PREFIX}RA7680.ct
export ENSEMBLECT=${PREFIX}ensemble.ct
export OLIGOLIST=${PREFIX}oligoscreen_example.lis
export SINGLESAV=${PREFIX}RA7680.fsv

export BIMOLCT=${PREFIX}ec5s_ca5s_bimol.ct
export CIRC_PREDICTED=${PREFIX}CircleCompare_predicted.ct
export CIRC_ACCEPTED=${PREFIX}CircleCompare_accepted.ct
export DOUBLECT="$CIRC_PREDICTED $CIRC_ACCEPTED"
export DOUBLECT2="$KNOTSCT $CIRC_ACCEPTED"
export DYNDOTPLOT=${PREFIX}DynalignDotPlot
export PROBPLOT=${PREFIX}ProbabilityPlot
