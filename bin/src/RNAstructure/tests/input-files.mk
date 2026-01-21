################################################
# This makefile generates input files for tests.
################################################

ALL=

# List of RNAstructure executables required to build test input files.
TOOLS=partition Fold ProbabilityPlot bifold dynalign

# Add RNAstructure executables to the PATH
PATH+=:../exe
export PATH

# Map simple tool name to the executable in ../exe
$(TOOLS): %: ../exe/%

# Rule to build standard RNAstructure tools (e.g. partition)
../exe/%:
	make -C .. $*

# Rule to build PFS files from fasta files
%.pfs: %.fasta partition
	partition  $< $@

# Rule to build PFS files from seq files
%.pfs: %.seq partition
	partition  $< $@

# Rule to build Folding Save Files (SAV, FSV) from seq files
%.fsv %.sav: %.seq Fold
	Fold  $< /dev/null -s $@

# Rule to build Folding Save Files (SAV, FSV) from fasta files
%.fsv %.sav: %.fasta Fold
	Fold  $< /dev/null -s $@

# Rule to build probabilty text files from PFS files
%.pfs.txt: %.pfs ProbabilityPlot
	ProbabilityPlot  $< --text $@

# Build input for DynalignDotPlot (dsv)
DDPLOT=testFiles/testFile_DynalignDotPlot
define DDPLOT_CONF
 inseq1 = testFiles/testFile_RD0260.seq
 inseq2 = testFiles/testFile_RD0500.seq
 outct  = $(DDPLOT)~1.ct
 outct2 = $(DDPLOT)~2.ct
 aout   = $(DDPLOT)~.ali
 savefile = $(DDPLOT).dsv
endef
export DDPLOT_CONF
$(DDPLOT).dsv: dynalign \
	testFiles/testFile_RD0260.seq \
	testFiles/testFile_RD0500.seq \

	echo "$$DDPLOT_CONF" >$(DDPLOT)~.conf
	dynalign $(DDPLOT)~.conf
	rm $(DDPLOT)~* # remove temporary output files.
