# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs

DESIGN_INPUT1=design/input/hairpin_structure.ct
DESIGN_INPUT2=design/input/bmorivector_structure.ct
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# For tests to be reproducible, we need to specify a random seed.
# So even the test "without options" has a seed option. (-s 1)
# (The number 1 was arbitrarily chosen, and has no special significance.)
# Test design_without_options.
runFullTest 'without_options' $DESIGN_INPUT1 -o @OUT.ct -s 1

# Test design_dna_option.
runFullTest 'dna_option' $DESIGN_INPUT1 -o @OUT.ct -s 1    -d

# Test design_random_seed_option.
# Test another arbitrary seed option (123456789). In general, different random seeds lead to different results.
runFullTest 'seed_option' $DESIGN_INPUT1 -o @OUT.ct        -s 123456789

# Test design_max_defect_option.
runFullTest 'max_defect_option' $DESIGN_INPUT1 -o @OUT.ct -s 1    -e 0.001

# Test design_preselect_option.
runFullTest 'preselect_option' $DESIGN_INPUT1 -o @OUT.ct -s 1    -p

# Test design_maxmutate_option.
#	Set the maximum number of times a nucleotide will be mutated during
#    defect-weighted reoptimization. The default is 4.
runFullTest 'maxmutate_option' $DESIGN_INPUT2 -o @OUT.ct -s 1    -e 0.1 -p -mm 2

# Test design_maxleaf_option.
# The maximum number of times a leaf can be re-optimized at random. The default is 3.
runFullTest 'maxleaf_option' $DESIGN_INPUT2 -o @OUT.ct -s 1    -e 0.1 -p -ml 1

# Test design_maxredesign_option.
# The maximum number of redesigns per parent node.  The default is 10.
runFullTest 'maxredesign_option' $DESIGN_INPUT2 -o @OUT.ct -s 1    -e 0.1 -p -mr 2

# Test design_maxdepth_option.
# Max-depth: The maximum extent to which the structure will be sub-divided in the binary decomposition. The default is 5.
runFullTest 'maxdepth_option' $DESIGN_INPUT2 -o @OUT.ct -s 1     -e 0.1 -p -md 2

# Test design_medium_structure.
runFullTest 'medium_structure' $DESIGN_INPUT2 -o @OUT.ct -s 1  -e 0.1 -p

endTestBlock # End a group of tests. Also cleans up orphaned files etc.
