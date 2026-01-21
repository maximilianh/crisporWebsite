# This script is not intended to be executed directly.  It is sourced by test-runner.sh and most functions you see here are defined in test-tools.sh
# The variables EXENAME, EXEBASE, OKDIR, and EXEDIR are defined externally.

beginTestBlock  # Begin a group of tests.
runMake -r @EXE  # run `make` for the main exe and any other listed programs
EXT=.ct  # set extension for @OKFILE and @OUTFILE placeholders (defaults for runDiff)

# Run a TurboHomology test. This includes:
#   1. Writing the TurboHomology configuration file
#   2. running TurboHomology
#   3. Comparing the output (alignment file and CT files) with references ("OK" files)
# Usage: 
#   turboTest <TEST_NAME>  <CONF_TEXT> [--standard | <SEQ_NAMES...>]
# Where SEQ_NAMES is a list of sequence names for which the CT file output 
# should be compared to analogously-named reference ("OK") files.
function turboTest() { 
	# Initilize the test with the testing system
	initTest "$1" || return  # exit if the test is excluded etc.
	
	# Write the conf file
	# Note that whereas runTest and runDiff automatically replace test placeholders (such as @TEST and @OUT etc),
	# we have to do this explicitly with `replaceVarsInText` because we are writing to a file.
	replaceVarsInText "$2" > "$TEST.conf"  # Write output to $TEST.conf e.g. "TurboFold_general_without_options.conf"
	
	# Now run TurboFold with the given conf
	runTest @EXE "$TEST.conf"

	# Now compare the output alignment files and (optionally) CT files
	# We can compare with the "standard" output files (from the general_without_options test)
	# or with a custom set of reference files.
	local okfile sequences
	if [[ $3 == --standard ]]; then
		setTestRef 'general_without_options'  # use the  'general_without_options' OK files as references for this test. This affects @OKBASE, @OK, and @OKFILE
		set -- seq1 seq2 seq3  # set the argument list to these three names
	else
		shift 2   # sequences are listed as all remaining arguments (after the first two)
	fi

	# Run a diff on the alignment file 
	# diff output is written to ${TEST}_aln_diff.txt e.g.  "TurboFold_general_without_options_aln_diff.txt"
	runDiff @OUT.aln @OK.aln 'aln'  

	# Now for each of the listed input sequences, do a diff with the expected output CT
	local seqName 
	for seqName; do  # loop through each sequence name in the remaining arguments (i.e. "$@")
		runDiff "@TEST_${seqName}_output.ct"  "@OKBASE_${seqName}_OK.ct"  "${seqName}"  # diff output written to  ${TEST}_${seqName}_diff.txt
	done

	endTest # perform end-of-test duties
}

# This is the directory where Turbofold input files are stored
INPUT=TurboHomology/input

# Test TurboHomology SRP testset 
turboTest 'TurboHomology_Fasta_output_options_groups_SRP' "
InSeq = {$INPUT/Dani.reri._EH487580.seq;$INPUT/Dros.viri._GSP-7244.seq;$INPUT/Cion.inte._BW230361.seq;$INPUT/Cani.spec._SDB-9616.seq;$INPUT/Taen.gutt._CK235759.seq;$INPUT/Orni.anat._GSP-9258.seq;$INPUT/Stro.purp._GSP-7668.seq;$INPUT/Homo.sapi._CK004881.seq;$INPUT/Lott.giga._GSP-225164.seq;$INPUT/Dros.mela._X00952.seq;$INPUT/Homo.sapi._BF908693.seq;$INPUT/Batr.dend._GSP-403673.seq;$INPUT/Mono.dome._GSP-13616.seq;$INPUT/Homo.sapi._NR_002715.seq;$INPUT/Caen.eleg._Z30973.seq;$INPUT/Caen.eleg._Z34533.seq;$INPUT/Tetr.nigr._GSP-99883.seq;$INPUT/Dros.yaku._GSP-7245.seq;$INPUT/Homo.sapi._X02067.seq;$INPUT/Sus.scro._BP459377.seq;$INPUT/Caen.eleg._AY948606.seq;$INPUT/Aply.cali._GSP-6500.seq;$INPUT/Mus.musc._DQ285765.seq;$INPUT/Dros.will._GSP-7260.seq;$INPUT/Ratt.norv._AC091616.seq;$INPUT/Homo.sapi._AA559987.seq;$INPUT/Cion.inte._GSP-7719.seq;$INPUT/Dani.rrer._CR925739.seq;$INPUT/Dros.erec._GSP-7220.seq;$INPUT/Homo.sapi._BF875598.seq;$INPUT/Dros.mela._X01056.seq;$INPUT/Phak.pach._GSP-170000-1.seq;$INPUT/Maca.mula._GSP-9544.seq;$INPUT/Homo.sapi._AL627171.seq;$INPUT/Onch.volv._AI249203.seq;$INPUT/Homo.sapi._BF909425.seq;$INPUT/Fugu.rubr._U40756.seq;$INPUT/Dani.reri._EH513970.seq;$INPUT/Dros.simu._GSP-7240.seq;$INPUT/Oryc.cuni._AB021174.seq;$INPUT/Bran.flor._GSP-7739.seq;$INPUT/Schm.medi._GSP-79327.seq;$INPUT/Homo.sapi._AL139099.seq;$INPUT/Caen.eleg._AY948608.seq;$INPUT/Dani.rrer._CU179657.seq;$INPUT/Caen.rema._GSP-31234.seq;$INPUT/Dani.rrer._BX469912.seq;$INPUT/Pris.paci._GSP-54126.seq;$INPUT/Caen.brig._AC084446.seq;$INPUT/Ratt.norv._AA923995.seq;$INPUT/Loxo.afri._GSP-9785.seq;$INPUT/Cion.savi._BW551179.seq;$INPUT/Bos.taur._GSP-6282.seq;$INPUT/Call.jacc._GSP-9483.seq;$INPUT/Cion.inte._BW088710.seq;$INPUT/Anop.gamb._CD747410.seq;$INPUT/Cypr.carp._CF662133.seq;$INPUT/Caen.eleg._AY948607.seq;$INPUT/Tric.spir._GSP-6334.seq;$INPUT/Aede.aegy._DW714815.seq;$INPUT/Gall.gall._AB073218.seq;$INPUT/Lemu.catt._GSP-9447.seq;$INPUT/Dasy.nove._GSP-9361.seq;$INPUT/Onch.volv._GSP-6282.seq;$INPUT/Onch.volv._BF727619.seq;$INPUT/Mus.musc._AC126244.seq;$INPUT/Cion.savi._BW520975.seq;$INPUT/Nema.vect._GSP-45351.seq;$INPUT/Homo.sapi._X04249.seq;$INPUT/Homo.sapi._BE083383.seq;$INPUT/Aede.aegy._GSP-7159.seq;$INPUT/Caen.eleg._U41993.seq;$INPUT/Cion.inte._BW223536.seq;$INPUT/Homo.sapi._BU566906.seq;$INPUT/Phak.meib._GSP-169999.seq;$INPUT/Xeno.laev._X01055.seq;$INPUT/Dros.anan._GSP-7217.seq;$INPUT/Aede.aegy._DW714804.seq;$INPUT/Mono.brev._GSP-81824.seq;$INPUT/Brug.mala._GSP-6279.seq;$INPUT/Dros.moja._GSP-7230.seq;$INPUT/Dani.reri._CR854846.seq;$INPUT/Dros.mela._AC002512.seq;$INPUT/Anop.gamb._SDB-377271.seq;$INPUT/Pan.trog._GSP-37010.seq;$INPUT/Dani.reri._BX649540.seq;$INPUT/Dros.anan._CJ969290.seq;$INPUT/Dani.rrer._BX470258.seq;$INPUT/Schi.mans._GSP-6183.seq;$INPUT/Phak.pach._GSP-170000-2.seq;}
RefCT = {$INPUT/Dros.viri._GSP-7244.ct;$INPUT/Cion.inte._BW230361.ct;$INPUT/Cani.spec._SDB-9616.ct;$INPUT/Taen.gutt._CK235759.ct;$INPUT/Orni.anat._GSP-9258.ct;$INPUT/Stro.purp._GSP-7668.ct;$INPUT/Homo.sapi._CK004881.ct;$INPUT/Lott.giga._GSP-225164.ct;$INPUT/Dros.mela._X00952.ct;$INPUT/Homo.sapi._BF908693.ct;$INPUT/Batr.dend._GSP-403673.ct;$INPUT/Mono.dome._GSP-13616.ct;$INPUT/Homo.sapi._NR_002715.ct;$INPUT/Caen.eleg._Z30973.ct;$INPUT/Caen.eleg._Z34533.ct;$INPUT/Tetr.nigr._GSP-99883.ct;$INPUT/Dros.yaku._GSP-7245.ct;$INPUT/Homo.sapi._X02067.ct;$INPUT/Sus.scro._BP459377.ct;$INPUT/Caen.eleg._AY948606.ct;$INPUT/Aply.cali._GSP-6500.ct;$INPUT/Mus.musc._DQ285765.ct;$INPUT/Dros.will._GSP-7260.ct;$INPUT/Ratt.norv._AC091616.ct;$INPUT/Homo.sapi._AA559987.ct;$INPUT/Cion.inte._GSP-7719.ct;$INPUT/Dani.rrer._CR925739.ct;$INPUT/Dros.erec._GSP-7220.ct;$INPUT/Homo.sapi._BF875598.ct;$INPUT/Dros.mela._X01056.ct;$INPUT/Phak.pach._GSP-170000-1.ct;$INPUT/Maca.mula._GSP-9544.ct;$INPUT/Homo.sapi._AL627171.ct;$INPUT/Onch.volv._AI249203.ct;$INPUT/Homo.sapi._BF909425.ct;$INPUT/Fugu.rubr._U40756.ct;$INPUT/Dani.reri._EH513970.ct;$INPUT/Dros.simu._GSP-7240.ct;$INPUT/Oryc.cuni._AB021174.ct;$INPUT/Bran.flor._GSP-7739.ct;$INPUT/Schm.medi._GSP-79327.ct;$INPUT/Homo.sapi._AL139099.ct;$INPUT/Caen.eleg._AY948608.ct;$INPUT/Dani.rrer._CU179657.ct;$INPUT/Caen.rema._GSP-31234.ct;$INPUT/Dani.rrer._BX469912.ct;$INPUT/Pris.paci._GSP-54126.ct;$INPUT/Caen.brig._AC084446.ct;$INPUT/Ratt.norv._AA923995.ct;$INPUT/Loxo.afri._GSP-9785.ct;$INPUT/Cion.savi._BW551179.ct;$INPUT/Bos.taur._GSP-6282.ct;$INPUT/Call.jacc._GSP-9483.ct;$INPUT/Cion.inte._BW088710.ct;$INPUT/Anop.gamb._CD747410.ct;$INPUT/Cypr.carp._CF662133.ct;$INPUT/Caen.eleg._AY948607.ct;$INPUT/Tric.spir._GSP-6334.ct;$INPUT/Aede.aegy._DW714815.ct;$INPUT/Gall.gall._AB073218.ct;$INPUT/Lemu.catt._GSP-9447.ct;$INPUT/Dasy.nove._GSP-9361.ct;$INPUT/Onch.volv._GSP-6282.ct;$INPUT/Onch.volv._BF727619.ct;$INPUT/Mus.musc._AC126244.ct;$INPUT/Cion.savi._BW520975.ct;$INPUT/Nema.vect._GSP-45351.ct;$INPUT/Homo.sapi._X04249.ct;$INPUT/Homo.sapi._BE083383.ct;$INPUT/Aede.aegy._GSP-7159.ct;$INPUT/Caen.eleg._U41993.ct;$INPUT/Cion.inte._BW223536.ct;$INPUT/Homo.sapi._BU566906.ct;$INPUT/Phak.meib._GSP-169999.ct;$INPUT/Xeno.laev._X01055.ct;$INPUT/Dros.anan._GSP-7217.ct;$INPUT/Aede.aegy._DW714804.ct;$INPUT/Mono.brev._GSP-81824.ct;$INPUT/Brug.mala._GSP-6279.ct;$INPUT/Dros.moja._GSP-7230.ct;$INPUT/Dani.reri._CR854846.ct;$INPUT/Dros.mela._AC002512.ct;$INPUT/Anop.gamb._SDB-377271.ct;$INPUT/Pan.trog._GSP-37010.ct;$INPUT/Dani.reri._BX649540.ct;$INPUT/Dros.anan._CJ969290.ct;$INPUT/Dani.rrer._BX470258.ct;$INPUT/Schi.mans._GSP-6183.ct;$INPUT/Phak.pach._GSP-170000-2.ct;}
OutCT = {@TEST_seq1_output.ct;}
ExistingAln = $INPUT/template_SRP.fasta
Mode = MEA
AlnFormat = Fasta
OutAln = @OUT.aln
$SMP_SETTING"   seq1

# Test TurboHomology_without_options_groups.
# Note that the general test has a mode specified of MEA; MEA is the default mode.
turboTest 'general_without_options_groups' "
InSeq = {$INPUT/testFile_E00001.seq;$INPUT/testFile_E00002.seq;$INPUT/testFile_E00003.seq;$INPUT/testFile_E00004.seq;$INPUT/testFile_E00005.seq;$INPUT/testFile_E00006.seq;$INPUT/testFile_E00007.seq;$INPUT/testFile_E00008.seq;$INPUT/testFile_E00009.seq;$INPUT/testFile_E00010.seq;$INPUT/testFile_E00011.seq;$INPUT/testFile_E00012.seq;$INPUT/testFile_E00014.seq;$INPUT/testFile_E00015.seq;}
RefCT = {$INPUT/testFile_E00002.ct;$INPUT/testFile_E00003.ct;$INPUT/testFile_E00004.ct;$INPUT/testFile_E00005.ct;$INPUT/testFile_E00006.ct;$INPUT/testFile_E00007.ct;$INPUT/testFile_E00008.ct;$INPUT/testFile_E00009.ct;$INPUT/testFile_E00010.ct;$INPUT/testFile_E00011.ct;$INPUT/testFile_E00012.ct;$INPUT/testFile_E00014.ct;$INPUT/testFile_E00015.ct;}
OutCT = {@TEST_seq1_output.ct;}
ExistingAln = $INPUT/testFile_TurboHomology.fasta
Mode = MEA
OutAln = @OUT.aln
$SMP_SETTING"   seq1

# Test TurboHomology_Fasta_output_option.
turboTest 'TurboHomology_Fasta_output_option' "
InSeq = {$INPUT/testFile_E00001.seq;$INPUT/testFile_E00002.seq;$INPUT/testFile_E00003.seq;$INPUT/testFile_E00004.seq;$INPUT/testFile_E00005.seq;$INPUT/testFile_E00006.seq;$INPUT/testFile_E00007.seq;$INPUT/testFile_E00008.seq;$INPUT/testFile_E00009.seq;$INPUT/testFile_E00010.seq;$INPUT/testFile_E00011.seq;$INPUT/testFile_E00012.seq;$INPUT/testFile_E00014.seq;$INPUT/testFile_E00015.seq;}
RefCT = {$INPUT/testFile_E00002.ct;$INPUT/testFile_E00003.ct;$INPUT/testFile_E00004.ct;$INPUT/testFile_E00005.ct;$INPUT/testFile_E00006.ct;$INPUT/testFile_E00007.ct;$INPUT/testFile_E00008.ct;$INPUT/testFile_E00009.ct;$INPUT/testFile_E00010.ct;$INPUT/testFile_E00011.ct;$INPUT/testFile_E00012.ct;$INPUT/testFile_E00014.ct;$INPUT/testFile_E00015.ct;}
OutCT = {@TEST_seq1_output.ct;}
ExistingAln = $INPUT/testFile_TurboHomology.fasta
Mode = MEA
OutAln = @OUT.aln
AlnFormat = Fasta
$SMP_SETTING"   seq1

# Test TurboHomology_output_ColumnNumber_option.
turboTest 'TurboHomology_output_ColumnNumber_option' "
InSeq = {$INPUT/testFile_E00001.seq;$INPUT/testFile_E00002.seq;$INPUT/testFile_E00003.seq;$INPUT/testFile_E00004.seq;$INPUT/testFile_E00005.seq;$INPUT/testFile_E00006.seq;$INPUT/testFile_E00007.seq;$INPUT/testFile_E00008.seq;$INPUT/testFile_E00009.seq;$INPUT/testFile_E00010.seq;$INPUT/testFile_E00011.seq;$INPUT/testFile_E00012.seq;$INPUT/testFile_E00014.seq;$INPUT/testFile_E00015.seq;}
RefCT = {$INPUT/testFile_E00002.ct;$INPUT/testFile_E00003.ct;$INPUT/testFile_E00004.ct;$INPUT/testFile_E00005.ct;$INPUT/testFile_E00006.ct;$INPUT/testFile_E00007.ct;$INPUT/testFile_E00008.ct;$INPUT/testFile_E00009.ct;$INPUT/testFile_E00010.ct;$INPUT/testFile_E00011.ct;$INPUT/testFile_E00012.ct;$INPUT/testFile_E00014.ct;$INPUT/testFile_E00015.ct;}
OutCT = {@TEST_seq1_output.ct;}
ExistingAln = $INPUT/testFile_TurboHomology.fasta
Mode = MEA
OutAln = @OUT.aln
ColumnNumber = 500
$SMP_SETTING"   seq1

# Test TurboHomology_output_ColumnNumber_option_2.
turboTest 'TurboHomology_output_ColumnNumber_option_2' "
InSeq = {$INPUT/testFile_E00001.seq;$INPUT/testFile_E00002.seq;$INPUT/testFile_E00003.seq;$INPUT/testFile_E00004.seq;$INPUT/testFile_E00005.seq;$INPUT/testFile_E00006.seq;$INPUT/testFile_E00007.seq;$INPUT/testFile_E00008.seq;$INPUT/testFile_E00009.seq;$INPUT/testFile_E00010.seq;$INPUT/testFile_E00011.seq;$INPUT/testFile_E00012.seq;$INPUT/testFile_E00014.seq;$INPUT/testFile_E00015.seq;}
RefCT = {$INPUT/testFile_E00002.ct;$INPUT/testFile_E00003.ct;$INPUT/testFile_E00004.ct;$INPUT/testFile_E00005.ct;$INPUT/testFile_E00006.ct;$INPUT/testFile_E00007.ct;$INPUT/testFile_E00008.ct;$INPUT/testFile_E00009.ct;$INPUT/testFile_E00010.ct;$INPUT/testFile_E00011.ct;$INPUT/testFile_E00012.ct;$INPUT/testFile_E00014.ct;$INPUT/testFile_E00015.ct;}
OutCT = {@TEST_seq1_output.ct;}
ExistingAln = $INPUT/testFile_TurboHomology.fasta
Mode = MEA
OutAln = @OUT.aln
ColumnNumber = 60
$SMP_SETTING"   seq1

# Test TurboHomology_Fasta_output_ColumnNumber_option.
turboTest 'TurboHomology_Fasta_output_ColumnNumber_option' "
InSeq = {$INPUT/testFile_E00001.seq;$INPUT/testFile_E00002.seq;$INPUT/testFile_E00003.seq;$INPUT/testFile_E00004.seq;$INPUT/testFile_E00005.seq;$INPUT/testFile_E00006.seq;$INPUT/testFile_E00007.seq;$INPUT/testFile_E00008.seq;$INPUT/testFile_E00009.seq;$INPUT/testFile_E00010.seq;$INPUT/testFile_E00011.seq;$INPUT/testFile_E00012.seq;$INPUT/testFile_E00014.seq;$INPUT/testFile_E00015.seq;}
RefCT = {$INPUT/testFile_E00002.ct;$INPUT/testFile_E00003.ct;$INPUT/testFile_E00004.ct;$INPUT/testFile_E00005.ct;$INPUT/testFile_E00006.ct;$INPUT/testFile_E00007.ct;$INPUT/testFile_E00008.ct;$INPUT/testFile_E00009.ct;$INPUT/testFile_E00010.ct;$INPUT/testFile_E00011.ct;$INPUT/testFile_E00012.ct;$INPUT/testFile_E00014.ct;$INPUT/testFile_E00015.ct;}
OutCT = {@TEST_seq1_output.ct;}
ExistingAln = $INPUT/testFile_TurboHomology.fasta
Mode = MEA
OutAln = @OUT.aln
AlnFormat = Fasta
ColumnNumber = 60
$SMP_SETTING"   seq1
