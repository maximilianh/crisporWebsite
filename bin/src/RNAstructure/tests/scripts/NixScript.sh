#!/bin/bash

# *****************************************************************************
# Note that the individual scripts sourced by this main script for individual
# test types are not independent; to work correctly the individual scripts
# must be called from within this one. The individual scripts are separated
# from this main script only for maintenance and readability purposes; they
# are essentially part of this script. Note also that the individual scripts
# are not executable on their own.
# *****************************************************************************

###############################################################################
## Error checking function and variables for common file names
###############################################################################
function checkErrors {
	#Usage: checkErrors TEST_NAME [FILES...]
	#Example: checkErrors efn2_without_options efn2_without_options_errors.txt efn2_without_options_diff_output.txt
	local testName=$1
	local errorLevel=0
	local file failureReason msg errNum fileType
	local errorDir=${PROGRAM_NAME}_OUTPUT
	local errorList=()
	shift

	for file in "$@"; do
		if [[ -s $file ]]; then
			case $file in 
				*_errors.txt) 
					errNum=10 #major error
					fileType='Error'
					;;

				*_diff_output.txt)
					errNum=5 #possible error
					fileType='Diff'
					;;
				*)
					errNum=1 #possible error
					fileType='Non-empty'
					;;
			esac
			if ((errNum>errorLevel)); then
				((errorLevel=errNum))
				numLines=$(wc -l $file | { read -r first _; echo $first; } || echo UNKNOWN )
				failureReason=" (see $fileType File $file; $numLines lines)"
			fi
		fi
	done

	if ((errorLevel)); then
		# Remove EMPTY diff and error files  (The "-s $file"  test is true only if the file has a size > 0)
		for file in ${testName}_errors.txt ${testName}_*_errors.txt ${testName}_diff_output.txt ${testName}_*_diff_output.txt; do
			[[ -e $file && ! -s $file ]] && rm -f $file  #i.e. if the file exists, but has size==0, then delete it.
		done
		# Copy error files to {PROGRAM_NAME}_OUTPUT folder
		[[ -d $errorDir ]] || mkdir $errorDir  #Make the output folder if it doesn't already exists.
		cp ${testName}* $errorDir/
		#if ((errorLevel>=10)); then
		#	msg='Test Failure'
		#else
			msg='Possible Test Failure'
		#fi
		echo "      ${msg}: $testName $failureReason" >&2
	else
		echo '      Test Passed.'
	fi
}


OS=$(uname -s 2>/dev/null)
if [[ $OS == CYGWIN_* || $OS == MSYS_* ]]; then
	export DIFF_CMD='diff --strip-trailing-cr'
else
	export DIFF_CMD=diff
fi

PROGRAM_NAME=$1

###############################################################################
# Set variables for input files.
###############################################################################

SINGLESEQ=testFiles/testFile_RA7680.seq
SINGLESEQ_SHORT=testFiles/testFile_RA7680_short.fasta
SINGLESEQ2=testFiles/testFile_met-vol.seq
SINGLESEQ2_FASTA=testFiles/testFile_met-vol.fasta
SINGLESEQ3=testFiles/testFile_ivslsu.seq
SINGLESEQ4=testFiles/testFile_ca5s.seq
DOUBLESEQ='testFiles/testFile_ec5s.seq testFiles/testFile_ca5s.seq'
SINGLEPFS=testFile_RA7680.pfs
SINGLEPFS2=testFile_met-vol.pfs
SINGLEPFS3=testFile_ivslsu.pfs
SINGLEPFS4=testFile_ca5s.pfs
KNOTSCT=testFiles/testFile_knotted.ct
SINGLECT=testFiles/testFile_RA7680.ct
ENSEMBLECT=testFiles/testFile_ensemble.ct
OLIGOLIST=testFiles/testFile_oligoscreen_example.lis
SINGLESAV=Fold_sav_test.sav
DOUBLECT='testFiles/testFile_CircleCompare_predicted.ct testFiles/testFile_CircleCompare_accepted.ct'
DOUBLECT2='testFiles/testFile_knotted.ct testFiles/testFile_CircleCompare_accepted.ct'
BIMOLCT=multi.ct

###############################################################################
# Call the appropriate individual script.
###############################################################################

if [[ $1 == AllSub ]]; then source AllSub/AllSub_Script;
elif [[ $1 == bifold ]]; then source bifold/bifold_Script;
elif [[ $1 == bifold-smp ]]; then source bifold/bifold_Script;
elif [[ $1 == bipartition ]]; then source bipartition/bipartition_Script;
elif [[ $1 == bipartition-smp ]]; then source bipartition/bipartition_Script;
elif [[ $1 == CircleCompare ]]; then source CircleCompare/CircleCompare_Script;
elif [[ $1 == ct2dot ]]; then source ct2dot2ct/ct2dot_Script;
elif [[ $1 == dot2ct ]]; then source ct2dot2ct/dot2ct_Script;
elif [[ $1 == DuplexFold ]]; then source DuplexFold/DuplexFold_Script;
elif [[ $1 == design ]]; then source design/design_Script;
elif [[ $1 == draw ]]; then source draw/draw_Script;
elif [[ $1 == dynalign ]]; then source dynalign/dynalign_Script;
elif [[ $1 == dynalign-smp ]]; then source dynalign/dynalign_Script;
elif [[ $1 == dynalign_ii ]]; then source dynalign_ii/dynalign_ii_Script;
elif [[ $1 == dynalign_ii-smp ]]; then source dynalign_ii/dynalign_ii_Script;
elif [[ $1 == DynalignDotPlot ]]; then source DynalignDotPlot/DynalignDotPlot_Script;
elif [[ $1 == efn2 ]]; then source efn2/efn2_Script;
elif [[ $1 == efn2-smp ]]; then source efn2/efn2_Script;
elif [[ $1 == EnergyPlot ]]; then source EnergyPlot/EnergyPlot_Script;
elif [[ $1 == EnsembleEnergy ]]; then source EnsembleEnergy/EnsembleEnergy_Script;
elif [[ $1 == Fold ]]; then source fold/Fold_Script;
elif [[ $1 == Fold-smp ]]; then source fold/Fold_Script;
elif [[ $1 == fold-cuda ]]; then source fold-cuda/fold-cuda_Script;
elif [[ $1 == MaxExpect ]]; then source MaxExpect/MaxExpect_Script;
elif [[ $1 == multilign ]]; then source multilign/multilign_Script;
elif [[ $1 == multilign-smp ]]; then source multilign/multilign_Script;
elif [[ $1 == NAPSS ]]; then source NAPSS/NAPSS_Script;
elif [[ $1 == oligoscreen ]]; then source oligoscreen/oligoscreen_Script;
elif [[ $1 == oligoscreen-smp ]]; then source oligoscreen/oligoscreen_Script;
elif [[ $1 == partition ]]; then source pfunction/partition_Script;
elif [[ $1 == partition-smp ]]; then source pfunction/partition_Script;
elif [[ $1 == partition-cuda ]]; then source partition-cuda/partition-cuda_Script;
elif [[ $1 == PARTS ]]; then source PARTS/PARTS_Script;
elif [[ $1 == ProbabilityPlot ]]; then source ProbabilityPlot/ProbabilityPlot_Script;
elif [[ $1 == ProbablePair ]]; then source ProbablePair/ProbablePair_Script;
elif [[ $1 == ProbScan ]]; then source ProbScan/ProbScan_Script;
elif [[ $1 == ProbKnot ]]; then source ProbKnot/ProbKnot_Script;
elif [[ $1 == refold ]]; then source refold/refold_Script;
elif [[ $1 == RemovePseudoknots ]]; then source RemovePseudoknots/RemovePseudoknots_Script;
elif [[ $1 == scorer ]]; then source scorer/scorer_Script;
elif [[ $1 == ShapeKnots ]]; then source ShapeKnots/ShapeKnots_Script;
elif [[ $1 == stochastic ]]; then source stochastic/stochastic_Script;
elif [[ $1 == stochastic-smp ]]; then source stochastic/stochastic_Script;
elif [[ $1 == TurboFold ]]; then source TurboFold/TurboFold_Script;
elif [[ $1 == TurboFold-smp ]]; then source TurboFold/TurboFold_Script;
else
    echo "    $1 Testing Protocol Not Determined Yet."
    echo "The testing protocol for $1 has not been determined yet." > "$ERROR_OUTPUT_DIR$1_tests_missing.txt"
    echo "This should mean that tests for $1 are coming soon." >>"$ERROR_OUTPUT_DIR$1_tests_missing.txt"
    cat "$ERROR_OUTPUT_DIR$1_tests_missing.txt" >&2
    echo >&2
fi
