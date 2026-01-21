# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r refold Fold-cuda # run `make` for the main exe and any other listed programs
EXT=.ct #ct file is final output from refold

# This function is run after each `runFullTest` is performed (but before output is validated).
# For Fold-cuda, the purpose of post-processing is to run `refold`` to convert the 
# SAV file (produced by Fold-Cuda) into a CT file for comparison with the reference.
function TEST_POST_PROCESS() {
	runSubTest 'refold' refold @OUT.sav @OUT.ct
}

# Test fold-cuda_without_options.
runFullTest 'without_options' time/ivslsu.seq @OUT.sav # sav file will be post-processed with refold

# Todo: enable DNA option (or determine why DNA option is not available)
#runFullTest 'dna_option' -d time/ivslsu.seq @OUT

EXT=.out
runFullTest 'array_option' -v time/ivslsu.seq @OUT.sav ---stdout  # File is quite large. 3.8 MB

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
