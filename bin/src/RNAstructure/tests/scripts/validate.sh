# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.err  # we are testing the error output
DIFF_FLAGS+=' -b' # ignore differences in spaces

EXPECT_OK=---exit=0   # expect a success exit code (0)
EXPECT_ERROR=---exit=1  # expect an error exit code (1)
TEST_ERR=---stderr  # we are testing the stderr output.

# Test validate with the general self-contained example given.
DIR=validate/input
runFullTest 'good_seq'				$DIR/good.seq 				$TEST_ERR $EXPECT_OK
runFullTest 'good_seq_plain'		$DIR/plainseq.txt -t seq    $TEST_ERR $EXPECT_OK  # plain text sequence
runFullTest 'good_seq_dna_option'	$DIR/good.seq     -a "dna" 	$TEST_ERR $EXPECT_OK
runFullTest 'good_ct'				$DIR/good.ct    			$TEST_ERR $EXPECT_OK
runFullTest 'good_dbn'				$DIR/good.dbn    			$TEST_ERR $EXPECT_OK
runFullTest 'good_shape'			$DIR/good.shape    			$TEST_ERR $EXPECT_OK
runFullTest 'good_olis'				$DIR/good.olis    			$TEST_ERR $EXPECT_OK
runFullTest 'good_fcon'				$DIR/good.fcon    			$TEST_ERR $EXPECT_OK
runFullTest 'good_mseq'				$DIR/good.multiseq   -t mseq    $TEST_ERR $EXPECT_OK
runFullTest 'good_mfasta'		    $DIR/good.multifasta -t mseq    $TEST_ERR $EXPECT_OK
runFullTest 'good_mseq_plain'		$DIR/plainseq.txt    -t mseq    $TEST_ERR $EXPECT_OK

runFullTest 'bad_seq_baseZ' 		$DIR/base-Z.seq  			$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_olis_baseZ'		$DIR/base-Z.olis    		$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_mseq_baseZ'		$DIR/base-Z.multiseq   -t mseq  $TEST_ERR $EXPECT_ERROR
runFullTest 'bad_mfasta_baseZ'		$DIR/base-Z.multifasta -t mseq  $TEST_ERR $EXPECT_ERROR

runFullTest 'good_seq_baseZ' 	    $DIR/base-Z.seq        -u   $TEST_ERR $EXPECT_OK
runFullTest 'good_olis_baseZ'		$DIR/base-Z.olis       -u   $TEST_ERR $EXPECT_OK
runFullTest 'good_mseq_baseZ'		$DIR/base-Z.multiseq   -u  -t mseq  $TEST_ERR $EXPECT_OK
runFullTest 'good_mfasta_baseZ'		$DIR/base-Z.multifasta -u  -t mseq  $TEST_ERR $EXPECT_OK

runFullTest 'bad_unknown_ext'		$DIR/good.seq.unk 			$TEST_ERR $EXPECT_ERROR
runFullTest 'good_unknown_ext'		$DIR/good.seq.unk  -t seq	$TEST_ERR $EXPECT_OK

runFullTest 'bad_fasta'				$DIR/bad.fasta 				$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_ct'				$DIR/bad.ct 				$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_dbn_char'			$DIR/bad.dbn 				$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_dbn_unmatched'		$DIR/bad-unmatched.dbn		$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_shape_index'		$DIR/bad-index.shape    	$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_shape_value'		$DIR/bad-value.shape    	$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_fcon_section'		$DIR/bad-section.fcon    	$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_fcon_index'		$DIR/bad-index.fcon     	$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_mseq_empty'		$DIR/bad-empty.multiseq   -t mseq  $TEST_ERR $EXPECT_ERROR
runFullTest 'bad_mfasta_empty'		$DIR/bad-empty.multifasta -t mseq  $TEST_ERR $EXPECT_ERROR
runFullTest 'bad_mseq_past1'		$DIR/bad-past1.multiseq   -t mseq  $TEST_ERR $EXPECT_ERROR

runFullTest 'good_seq_length'	    $DIR/good.seq 	-l 75		$TEST_ERR $EXPECT_OK
runFullTest 'good_ct_length'		$DIR/good.ct    -l 26		$TEST_ERR $EXPECT_OK
runFullTest 'good_dbn_length'	    $DIR/good.dbn   -l 26		$TEST_ERR $EXPECT_OK
runFullTest 'good_olis_length'		$DIR/good.olis  -l 26		$TEST_ERR $EXPECT_OK
runFullTest 'good_mseq_length'		$DIR/good.multiseq   -l 75 -t mseq    $TEST_ERR $EXPECT_OK
runFullTest 'good_mfasta_length'	$DIR/good.multifasta -l 75 -t mseq    $TEST_ERR $EXPECT_OK


runFullTest 'bad_seq_length' 	    $DIR/good.seq 	-l 74		$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_ct_length'		    $DIR/good.ct    -l 25		$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_dbn_length' 	    $DIR/good.dbn   -l 25		$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_olis_length'		$DIR/good.olis  -l 25		$TEST_ERR $EXPECT_ERROR
runFullTest 'bad_mseq_length'		$DIR/good.multiseq   -l 74 -t mseq    $TEST_ERR $EXPECT_ERROR
runFullTest 'bad_mfasta_length'     $DIR/good.multifasta -l 74 -t mseq    $TEST_ERR $EXPECT_ERROR

BIN=$DIR/binary.dat
runFullTest 'bad_ct_binary'		    "$BIN"      -t ct           $TEST_ERR $EXPECT_ERROR
runFullTest 'bad_dbn_binary'		"$BIN"      -t dbn          $TEST_ERR $EXPECT_ERROR
runFullTest 'bad_seq_binary'		"$BIN"      -t seq          $TEST_ERR $EXPECT_ERROR
runFullTest 'bad_shape_binary'		"$BIN"      -t chem         $TEST_ERR $EXPECT_ERROR
runFullTest 'bad_olis_binary'		"$BIN"      -t olis         $TEST_ERR $EXPECT_ERROR
runFullTest 'bad_fcon_binary'		"$BIN"      -t fcon         $TEST_ERR $EXPECT_ERROR

endTestBlock