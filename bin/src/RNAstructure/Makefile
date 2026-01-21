################################################################################
#  Main RNAstructure Makefile
################################################################################
#  Contents:
#  Section 1: List of Programs.
#  Section 2: Main Makefile Targets.
#  Section 3: Standard Program Dependencies.
#  Section 4: Standard Linking Recipes.
#  Section 5: Non-Standard Programs and Recipes.
#  Section 6: Special Interfaces (Java, Python)
#  Section 7: RNAstructure regression tests
#  Section 8: Cleanup and Other Maintenance Goals.
#  Section 9: Make Help Message.
################################################################################
ROOTPATH=.# Define the relative path to the RNAstructure root directory.
## Include all macro, dependency, and variable definitions.
include ${ROOTPATH}/compiler.h
include ${ROOTPATH}/library_defines.h
include ${ROOTPATH}/dependencies.h

################################################################################
#  Section 1: List of Programs
#
#  RNAstructure programs are listed below in groups. The variables defined 
#  in this section are important because they are used to:
#    * Display the 'Help' message which lists Makefile targets.
#    * Generate linking & compilation rules for each program.
#    * Perform cleanup (e.g. delete programs during `make realclean`)
################################################################################

####  SERIAL Programs  ####
# List all serial (single-threaded) programs here. 
# If a serial program also has an SMP counterpart, suffix its name with a 
# tilde (~), e.g. "bifold~".
#$(MF) will include Multifind in build if user uncomments the definition in compiler.h.
define SERIAL_PROGS
  AlignmentFold~ AlignmentPartition~ AccessFold  AllSub  bifold~  bipartition~  CircleCompare  ct2dot  CycleFold
  dot2ct design~ DistanceMatrix draw  dynalign~  dynalign_ii~  DuplexFold~  DynalignDotPlot
  EDcalculator~ efn2~ EnergyPlot  EnsembleEnergy  ETEcalculator~ Fold~  MaxExpect~  multilign~  NAPSS 
  oligoscreen~  OligoWalk   orega~
  partition~  PARTS  phmm  ProbabilityPlot
  ProbablePair~  ProbKnot~  ProbScan  refold  RemovePseudoknots~  scorer Shannon
  ShapeKnots~  StructureProb~ stochastic~  SymmetryTester TurboFold~  TurboHomology validate Rsample~ $(MF)
endef

####  SMP Programs  ####
# 'SMP' (Symmetric Multi-Processing) in general means programs that are 
# designed to run tasks simultaneously on multiple threads, ideally making use of 
# multiple processor cores, for parallel computation. 
# Different RNAstructure programs use different methods to achieve this: 
#  (1) OS-specific thread-handling  (e.g. TurboFold-smp via PHMM's ansi_thread etc)
#  (2) POSIX Threads (pthreads)     (e.g. dynalign-smp  using rankmanager etc.)
#  (3) OpenMP (OMP).                (e.g. Fold-smp, partition-smp, using "#pragma omp ...")
# The latter two require dedicated compiler support and system libraries.
# Serial programs suffixed with "~" are already automatically included.
# Only SMP programs that do NOT have serial analogs need to be listed here.
SMP_PROGS= # (Also includes all serial programs ending in ~)

####  CUDA Programs  ####
# CUDA is a parallel computing platform and API created by Nvidia.
# It allows a CUDA-enabled graphics processing unit (GPU) to be used 
# for general purpose processing. GPUs are highly parallel multi-core 
# systems allowing very efficient manipulation of large blocks of data.
CUDA_PROGS=Fold-cuda  partition-cuda  orega-cuda

####  MISC programs  ####
# These are other programs & tools used for testing, etc., that are not typically 
# compiled along with the above groups.
MISC_PROGS=unit-tests CycleFold-test

####  ALL programs  ####
# If new program groups are added, please make sure to add them to the following
# variable, which should contain ALL executable targets available in this Makefile.
ALL_PROGS=$(SERIAL_PROGS) $(SMP_PROGS) $(CUDA_PROGS) $(MISC_PROGS)

#### List Preparation -- (Users should not typically edit this.) ####
# Filter out serial vs smp tests and replace the ~ suffix with -smp for the latter.
SMP_FROM_SERIAL:=$(patsubst %~, %-smp, $(filter %~, $(strip $(SERIAL_PROGS))))
SMP_PROGS    :=$(sort  $(strip $(SMP_PROGS))   $(SMP_FROM_SERIAL)) # Added `strip` because some versions of Make have trouble with the embedded linefeeds.
SERIAL_PROGS :=$(sort $(patsubst %~, %, $(strip $(SERIAL_PROGS))))
CUDA_PROGS   :=$(sort $(strip $(CUDA_PROGS)))
MISC_PROGS   :=$(sort $(strip $(MISC_PROGS)))
################################################################################


################################################################################
#  Section 2: Main Makefile Targets.
################################################################################
.DEFAULT_GOAL=help # Set the default Makefile target

 # '.PHONY' indicates that these targets should always be run when invoked (even if a file of the same name exists)
.PHONY: all ALL serial SERIAL smp SMP CUDA cuda install instructions HELP help

# Default Make Target -- Show a help message (defined at bottom of Makefile)
instructions HELP help:
	@$(eval export MAKE_HELP_MESSAGE) echo -e "$$MAKE_HELP_MESSAGE"

# Make all serial and SMP executables
ALL all: serial smp

# Make all the serial executables.
SERIAL serial: $(SERIAL_PROGS) ; @:

# Make all the SMP (multi-processor) executables.
SMP smp: $(SMP_PROGS) ; @:

# Make all the CUDA (GPU) executables.
CUDA cuda: $(CUDA_PROGS) ; @:

# Copy the executables to the /usr/local directory.
install:
	cp -r exe/ /usr/local/RNAstructure

####  Define Makefile Targets for all programs ####
# Redirect `make PROGRAM` to `make exe/PROGRAM` for all programs listed above.
$(ALL_PROGS): %: exe/%

# Create aliases for programs that are commonly mistyped.
# For example, without the following aliases:
#  `make fold` would NOT build the 'Fold' program (because the case is different). 
#  Instead the command would always succeed with the message "'fold' is up to date." 
#  because the RNAstructure/fold/ directory already exists.
fold: Fold
napss: NAPSS
pfunction: partition
Turbofold: TurboFold
Multilign: multilign
Orega: orega

# Also allow the lower-case name of each program to function as an alias.
#  E.g. `make oligowalk` will make OligoWalk.
#  Note that we do NOT add an alias if the name is already lower-case, because
#  doing so would create a "circular dependency" that would result in a warning.
$(foreach prog,$(ALL_PROGS),$(call fnAddDep,$(call fnToLower,$(prog)),$(prog)))


################################################################################
#  Section 3: Standard Program Dependencies.
#
#  Define the following variables for each RNAstructure program:
#  * OBJ_<PROGRAM>   -- (required) A list of object-file (*.o) dependencies 
#    that should be linked together to create the executable.
#  * FLAGS_<PROGRAM> -- (optional) A list of extra flags that should be added 
#    when compiling/linking the program.
#  * LIBS_<PROGRAM>  -- (optional) A list of extra linker flags that should 
#    be appended *after* all object-files during linking to specify additional 
#    system libraries to be linked with the executable.  (rarely necessary)
#  * LD_<PROGRAM>  -- (optional) A list of extra linker flags that should 
#    be appended *before* object-files during linking. (rarely necessary)
#
################################################################################

#  Note: Many of the object-file-collection variables used below (e.g. 
#  RNA_FILES, HYBRID_FILES, etc) are defined in dependencies.h
CLI=${CMD_LINE_PARSER} # The command-line parser is used for most text-interface programs.

# List text interface programs and their object-file (*.o) dependencies.
OBJ_AccessFold=      AccessFold/AccessFold.o 	${CLI} ${HYBRID_FILES}
OBJ_AlignmentFold=          AlignmentFold/AlignmentFold.o ${CLI} ${RNA_FILES}
OBJ_AlignmentFold-smp=      AlignmentFold/AlignmentFold.o ${CLI} ${RNA_FILES_SMP}
OBJ_AlignmentPartition=     AlignmentPartition/AlignmentPartition.o ${CLI} ${PARTITION_FILES}
OBJ_AlignmentPartition-smp= AlignmentPartition/AlignmentPartition.o ${CLI} ${PARTITION_FILES_SMP}
OBJ_AllSub=          AllSub/AllSub.o 			${CLI} ${RNA_FILES}
OBJ_bifold=          bifold/bifold.o 			${CLI} ${HYBRID_FILES}
OBJ_bifold-smp=      bifold/bifold.o 			${CLI} ${HYBRID_FILES_SMP}
OBJ_bipartition=     bipartition/bipartition.o 	${CLI} ${HYBRID_FILES}
OBJ_bipartition-smp= bipartition/bipartition.o  ${CLI} ${HYBRID_FILES_SMP}
OBJ_CircleCompare=   CircleCompare/CircleCompare_Interface.o ${CLI} ${RNA_FILES} \
                     ${STRUCTURE_SCORER} ${STRUCTURE_IMAGER} ${STRUCTURE_COMPARER} \
                     ${PROBSCAN_FILES}
OBJ_ct2dot=          ct2dot/ct2dot.o ${CLI} ${RNA_FILES}
OBJ_dot2ct=          dot2ct/dot2ct.o ${CLI} ${RNA_FILES}
OBJ_CycleFold=       ${CYCLEFOLD_FILES}
OBJ_CycleFold-test= $(patsubst %/main.o,%/main-test.o,${CYCLEFOLD_FILES}) # Use main-test.o instead of main.o for CycleFold-test
OBJ_validate=        validate/validate.o ${CLI} ${RNA_FILES} ${ROOTPATH}/src/RNAFileType.o
OBJ_design=          design/design_interface.o ${CLI} ${RNA_FILES} # NOTE: RNA_class/design.o is included in RNA_FILES
OBJ_design-smp=          design/design_interface.o ${CLI} ${RNA_FILES_SMP} 
OBJ_EDcalculator=    EDcalculator/EDcalculator_Interface.o ${CLI} ${RNA_FILES}
OBJ_EDcalculator-smp=EDcalculator/EDcalculator_Interface.o ${CLI} ${RNA_FILES_SMP}
OBJ_draw=            draw/DrawStructure.o ${CLI} ${STRUCTURE_IMAGER} ${RNA_FILES} ${PROBSCAN_FILES}
OBJ_DuplexFold=      DuplexFold/DuplexFold.o ${CLI} ${HYBRID_FILES}
OBJ_DistanceMatrix=  DistanceMatrix/DistanceMatrix.o ${CLI} ${RNA_FILES}
OBJ_DuplexFold-smp=  DuplexFold/DuplexFold.o ${CLI} ${HYBRID_FILES_SMP}
OBJ_dynalign=        dynalign/dynaligninterface.o ${DYNALIGN_SERIAL_FILES}
OBJ_dynalign-smp=    dynalign/dynaligninterface-smp.o ${DYNALIGN_SMP_FILES}
OBJ_dynalign_ii=     dynalign/dynaligninterface_ii.o ${DYNALIGN_II_SERIAL_FILES}
OBJ_dynalign_ii-smp= dynalign/dynaligninterface_ii-smp.o ${DYNALIGN_II_SMP_FILES}
OBJ_DynalignDotPlot= DynalignDotPlot/DynalignDotPlot.o ${CLI} ${PLOT_HANDLER} ${DYNALIGN_SERIAL_FILES}
OBJ_efn2=            efn2/efn2.o     ${CLI} ${RNA_COUNT_FILES} ${PSEUDOKNOT_FILES_COUNT}
OBJ_efn2-smp=        efn2/efn2-smp.o ${CLI} ${RNA_FILES} ${PSEUDOKNOT_FILES}
OBJ_EnergyPlot=      EnergyPlot/EnergyPlot.o ${CLI} ${PLOT_HANDLER} ${RNA_FILES}
OBJ_EnsembleEnergy=  EnsembleEnergy/EnsembleEnergy_Interface.o ${CLI} ${RNA_FILES}
OBJ_ETEcalculator=    ETEcalculator/ETEcalculator_Interface.o ${CLI} ${RNA_FILES}
OBJ_ETEcalculator-smp=ETEcalculator/ETEcalculator_Interface.o ${CLI} ${RNA_FILES_SMP}
OBJ_Fold=            fold/Fold.o ${CLI} ${RNA_FILES}
OBJ_Fold-smp=        fold/Fold.o ${CLI} ${RNA_FILES_SMP}
OBJ_MaxExpect=       MaxExpect/MaxExpectInterface.o ${CLI} ${RNA_FILES}
OBJ_MaxExpect-smp=   MaxExpect/MaxExpectInterface.o ${CLI} ${RNA_FILES_SMP}
OBJ_multilign=       multilign/Multilign_Interface.o RNA_class/Multilign_object.o ${CLI} \
                     ${DYNALIGN_SERIAL_FILES}
OBJ_multilign-smp=   multilign/Multilign_Interface-smp.o RNA_class/Multilign_object.o ${CLI}\
                     ${DYNALIGN_SMP_FILES}
OBJ_Multifind=       Multifind/Multifind_Interface.o RNA_class/Multifind_object.o RNA_class/Multilign_object-Multifind.o \
                     ${CLI} ${DYNALIGN_SERIAL_FILES} ${INCLUDESVM}/svm.o
OBJ_Multifind-smp=   Multifind/Multifind_Interface-smp.o RNA_class/Multifind_object.o RNA_class/Multilign_object-Multifind.o \
                     ${CLI} ${DYNALIGN_SMP_FILES} ${INCLUDESVM}/svm.o
OBJ_NAPSS=           napss/NAPSS_Interface.o napss/napss.o src/algorithm_instrumented.o  \
                     ${PSEUDOKNOT_FILES} ${CLI} ${CONF_PARSER} ${RNA_FILES}
OBJ_oligoscreen=     oligoscreen/oligoscreen.o ${CLI} ${OLIGO_FILES}
OBJ_oligoscreen-smp= oligoscreen/oligoscreen.o ${CLI} ${OLIGO_FILES_SMP}
OBJ_OligoWalk=       oligowalk/src/oligowalk.o oligowalk/src/globals.o ${CMD_LINE_PARSER} ${OLIGO_FILES}
OBJ_orega=           orega/orega_interface.o orega/GeneticAlgorithm.o orega/RnaContainerClass.o \
                     ${CLI} ${RNA_FILES}
OBJ_orega-smp=       orega/orega_interface.o orega/GeneticAlgorithm.o orega/RnaContainerClass.o \
                     ${CLI} ${RNA_FILES_SMP}
OBJ_PARTS=           PARTS/src/parts/main.o ${PARTS_FILES} ${RNA_FILES} ${PHMM_FILES}
OBJ_partition=       pfunction/partition.o ${CLI} ${PARTITION_FILES}
OBJ_partition-smp=   pfunction/partition.o ${CLI} ${PARTITION_FILES_SMP}
OBJ_phmm=            phmm/phmm_interface.o ${DYNALIGN_SERIAL_FILES} ${CLI}
OBJ_ProbabilityPlot= ProbabilityPlot/ProbabilityPlot.o ${CLI} ${PLOT_HANDLER} ${RNA_FILES}
OBJ_ProbablePair=    ProbablePair/ProbablePair.o ${CLI} ${RNA_FILES}
OBJ_ProbablePair-smp=ProbablePair/ProbablePair.o ${CLI} ${RNA_FILES_SMP}
OBJ_ProbKnot=        ProbKnot/ProbKnot_Interface.o ${CLI} ${RNA_FILES}
OBJ_ProbKnot-smp=    ProbKnot/ProbKnot_Interface.o ${CLI} ${RNA_FILES_SMP}
OBJ_ProbScan=        ProbScan/ProbScan_Interface.o ${CLI} ${RNA_FILES} ${PROBSCAN_FILES}
OBJ_ProbScan-smp=    ProbScan/ProbScan_Interface.o ${CLI} ${PROBSCAN_FILES_SMP}
OBJ_refold=          refold/refold.o ${CLI} ${RNA_FILES}
OBJ_RemovePseudoknots=     RemovePseudoknots/RemovePseudoknots.o ${CLI} ${RNA_FILES}
OBJ_RemovePseudoknots-smp= RemovePseudoknots/RemovePseudoknots.o ${CLI} ${RNA_FILES_SMP}
OBJ_scorer=          scorer/Scorer_Interface.o ${CLI} ${STRUCTURE_SCORER} ${RNA_FILES}
OBJ_Shannon=            Shannon/Shannon.o     ${CLI} ${RNA_FILES}
OBJ_ShapeKnots=      ${CLI} ${RNA_FILES} ${PSEUDOKNOT_FILES} ${SHAPEKNOTS_FILES}
OBJ_ShapeKnots-smp=  ${CLI} ${RNA_FILES_SMP} ${PSEUDOKNOT_FILES} ${SHAPEKNOTS_FILES}
OBJ_stochastic=      stochastic/stochastic.o ${CLI} ${RNA_FILES}
OBJ_stochastic-smp=  stochastic/stochastic.o ${CLI} ${RNA_FILES_SMP}
OBJ_StructureProb=   StructureProb/StructureProb.o ${CLI} ${RNA_FILES} ${PSEUDOKNOT_FILES}
OBJ_StructureProb-smp=StructureProb/StructureProb.o ${CLI} ${RNA_FILES} ${PSEUDOKNOT_FILES}
OBJ_SymmetryTester= SymmetryTester/SymmetryTester.o ${CLI} ${RNA_FILES}
OBJ_TurboFold=       TurboFold/TurboFold_Interface.o  ${TURBOFOLD_FILES} \
                     ${CLI} ${RNA_FILES} ${CONF_PARSER}
# Note: Turbofold SMP uses the SERIAL versions of the RNA_FILES
OBJ_TurboFold-smp=   TurboFold/TurboFold_Interface-smp.o ${TURBOFOLD_FILES_SMP} \
                     ${CLI} ${RNA_FILES} ${CONF_PARSER}
OBJ_TurboHomology=   TurboHomology/TurboHomology_Interface.o ${TURBOHOMOLOGY_FILES} \
                     ${CONF_PARSER} ${CLI} ${RNA_FILES}
OBJ_Rsample=         Rsample/Rsample.o ${CLI} ${RNA_FILES}
OBJ_Rsample-smp=     Rsample/Rsample.o ${CLI} ${RNA_FILES_SMP}

#### CycleFold ####
# CycleFold uses some non-standard compiler and linking flags.
FLAGS_CycleFold=-std=c++0x -pthread
FLAGS_CycleFold-test=$(FLAGS_CycleFold) -D TEST
LIBS_CycleFold-test=-lboost_unit_test_framework #CycleFold-test must be compiled with the boost unit test library.


#### CUDA PROGRAMS ####

OBJ_partition-cuda-base=  partition-smp/base.o partition-smp/param.o partition-smp/prna.o partition-smp/util.o
OBJ_partition-cuda= partition-smp/main.o ${OBJ_partition-cuda-base}

OBJ_Fold-cuda-base= fold-smp/frna.o  fold-smp/fparam.o fold-smp/fbase.o \
                    partition-smp/util.o # Note: Fold-cuda uses util.o from partition-cuda

OBJ_Fold-cuda=      fold-smp/main.o ${OBJ_Fold-cuda-base}

OBJ_orega-cuda=      ${OBJ_orega} ${OBJ_partition-cuda-base} ${OBJ_Fold-cuda-base}# orega-cuda uses the same list of object files, but it is built with CUDA, while orega is not.

#### MISC PROGRAMS ####
OBJ_unit-tests= unit-tests/unit-tests.o ${RNA_FILES}

####  Non-Standard Programs  ####
# Most programs can be linked using the standard recipes in Section 4.
# Any programs that should NOT be compiled with the standard linking 
# recipes should be listed in NON_STANDARD_PROGRAMS below.
define NON_STANDARD_PROGRAMS
  
endef

################################################################################
#  Section 4: Standard Linking Recipes.
################################################################################

# Add object-file dependencies and program-specific flags for all programs listed in Section 1.
 $(foreach prog,$(ALL_PROGS),$(eval exe/$(prog): $$(OBJ_$(prog)))) # For all programs listed above, add object-file dependencies
 $(foreach prog,$(ALL_PROGS),$(eval exe/$(prog): CXXFLAGS+=$$(FLAGS_$(prog)))) # For all programs listed above, add flags defined in FLAGS_<PROGRAM> to CXXFLAGS
 $(foreach prog,$(ALL_PROGS),$(eval exe/$(prog): LDLIBS+=$$(LIBS_$(prog)))) # For all programs listed above, add flags defined in LIBS_<PROGRAM> to LDLIBS
 $(foreach prog,$(ALL_PROGS),$(eval exe/$(prog): LDFLAGS+=$$(LD_$(prog)))) # For all programs listed above, add flags defined in LD_<PROGRAM> to LDFLAGS

# This defines the default linking rule that is applied to all SERIAL and MISC 
# programs (unless they are listed in NON_STANDARD_PROGRAMS).
STANDARD_PROGS=$(addprefix exe/,$(filter-out $(NON_STANDARD_PROGRAMS), $(SERIAL_PROGS) $(MISC_PROGS)))
$(STANDARD_PROGS):
	@# LINK_DEPS Uses CXX to link all object-file dependencies. (See compiler.h)
	${LINK_DEPS}


# This defines the default linking rule that is applied to 
# ALL serial and SMP and misc programs (unless they are listed in NON_STANDARD_PROGRAMS).
STANDARD_SMP_PROGS=$(addprefix exe/,$(filter-out $(NON_STANDARD_PROGRAMS), $(SMP_PROGS)))
$(STANDARD_SMP_PROGS):
	@# LINK_SMP_DEPS Uses CXX to link all object-file dependencies using SMP and OpenMP flags (See compiler.h)
	${LINK_SMP_DEPS}

# The linking rule that applies to all CUDA programs (unless they are listed in NON_STANDARD_PROGRAMS).
STANDARD_CUDA_PROGS=$(addprefix exe/,$(filter-out $(NON_STANDARD_PROGRAMS), $(CUDA_PROGS)))
$(STANDARD_CUDA_PROGS): export RNA_INCLUDE_CUDA=1 # this is only needed for combined CPU/GPU programs, but it should also be compatible with pure-CUDA programs.
$(STANDARD_CUDA_PROGS):
	@# Uses CUDA_COMPILER to link all object-file dependencies.
	${LINK_CUDA}    $^

################################################################################
#  Section 5: Non-Standard Programs and Recipes.
################################################################################

# -----PYTHON-----
PY_DIR=python_interface
# These goals are aliases that all build some variant of the python interface.
#   they are differentiated by Make variables that are set below.
PYTHON_INTERFACES  = python_interface python python_setup python_dist      # standard python build (attempts to build both python2 and python3 interfaces)
PYTHON3_INTERFACES = python3_interface python3 python3_setup python3_dist # build for python3 only
PYTHON2_INTERFACES = python2_interface python2 python2_setup python2_dist # build for python2 only

ALL_PYTHON_INTERFACES = $(PYTHON_INTERFACES) $(PYTHON2_INTERFACES) $(PYTHON3_INTERFACES)
.PHONY: python_interface_sources  $(ALL_PYTHON_INTERFACES)

IS_CMD=type -p &>/dev/null
IS_PYDIST=$(if $(PYDIST),(using distutils),(using makefile))# message to show distutils is being used
NO_PY=$(if $(PYTHON),$(C_HASH)) # If PYTHON is set, this will evaluate to "#" causing the command to be commented out (ignored) 
YES_PY=$(if $(PYTHON),,$(C_HASH)) # If PYTHON is NOT set, this will evaluate to "#" causing the command to be commented out (ignored) 

$(info PY: $(YES_PY), $(NO_PY))

# If the PYTHON variable is NOT set, this attemps to build the python
#   interface for both python2 and python3 (and for 'python' if neither 
#   python2 nor python3 are valid commands).
# If PYTHON is set, then this builds the python interface ONLY for the specified 
#   python binary.
$(PYTHON_INTERFACES):
	@# if python2 exists, build the python2 interface
	@$(NO_PY) ! $(IS_CMD) python2 || { echo "Building python2 Interface $(IS_PYDIST)"; $(MAKE_CD) $(PY_DIR) interface PYTHON=python2; }
	@# if python3 exists, build the python2 interface
	@$(NO_PY) ! $(IS_CMD) python3 || { echo "Building python3 Interface $(IS_PYDIST)"; $(MAKE_CD) $(PY_DIR) interface PYTHON=python3; }
	@# if neither python3 nor python2 exists, try simply 'python'
	@$(NO_PY) $(IS_CMD) python3 || $(IS_CMD) python2 || { echo "Building Python Interface $(IS_PYDIST)"; $(MAKE_CD) $(PY_DIR) interface PYTHON=python; }
	@# if the user has specified PYTHON=... then use their value.
	@$(YES_PY) echo "Building Python Interface (for $(PYTHON)) $(IS_PYDIST)"; $(MAKE_CD) $(PY_DIR) interface

# Build the Python interface library (_RNAstructure_wrap.so).
$(PYTHON2_INTERFACES) $(PYTHON3_INTERFACES):
	@$(MAKE_CD) $(PY_DIR) interface

# Export a variable for goals that build the Python interface using python3
python3_interface python3 python3_setup python3_dist: export PYTHON=python3
# Export a variable for goals that build the Python interface using python2
python2_interface python2 python2_setup python2_dist: export PYTHON=python2
# Export a variable for goals that build the Python interface using distutils
python_setup python_dist python3_setup python3_dist python2_setup python2_dist:  export PYDIST=1

# Delegate other Python-related goals to the python Makefile. For example:
#    python-all     --  Run SWIG, update python sources and finally build the interface.
#    python-debug   --  Show some configuration information.
#    python-source-list --  Update list of source files used by distutils (setup.py)
#    python-scripts --  Copy RNAstructure.py and other scripts to exe/
#    python-clean   --  Clean the distutils build directory
python-%:
	@$(MAKE_CD) $(PY_DIR) $(subst python-,,$@)

python_interface_sources: $(PYTHON_LIBRARY_FILES)



################################################################################
#  Section 6: Special Interfaces (Java, Python)
################################################################################

#######################################
################ Java #################
#######################################
#### Build the RNAstructure Java GUI ####
.PHONY: gui GUI editor StructureEditor gui-tester
gui GUI:
	@$(MAKE_CD) java_interface gui

#### Build the Structure Drawing Program ####
editor StructureEditor:
	@$(MAKE_CD) java_interface editor

#### Build the Regression-Test Framework for the Java GUI ####
gui-tester:
	@$(MAKE_CD) java_interface gui tester

java-%:
	@$(MAKE_CD) java_interface $(subst java-,,$@)


################################################################################
#  Section 7: RNAstructure regression tests
################################################################################
.PHONY: tests
tests:
	@echo "Running RNAstructure SERIAL Tests" && $(MAKE_CD) tests/ serial
	@[ '$(OPSYSTEM)' != Mac ] && echo "Running RNAstructure SMP Tests" && $(MAKE_CD) tests/ smp

################################################################################
#  Section 8: Cleanup and Other Maintenance Goals.
#
# `make clean`     removes all temporary build objects (e.g. object and 
#                  dependency files)
# `make realclean` performs `clean` and also removes final build products 
#                  (i.e. executables).
################################################################################

# Remove object files and any temporary files from building.
.PHONY: clean
CLEAN_FILES= *~  *.o  *.class *.pyc  *_make.log *_make-errors.log
clean:
	@echo 'Deleting all object, java-class, pyc, and temporary files.'
	@find . \( -false $(patsubst %,-o -name '%',$(CLEAN_FILES)) \) -type f -delete

# Remove Visual-Studio-specific object files
.PHONY: cleanVS
CLEAN_FILES_VS= *~  *.obj  *.pdb  *.tlog  *.VC.db  *.ilk 
cleanVS:
	@echo 'Deleting all Visual Studio object and temp files.'
	@find . \( -false $(patsubst %,-o -name '%',$(CLEAN_FILES_VS)) \) -type f  -print  -delete

.PHONY: clean-empty-logs
clean-empty-logs:
	@find . \( -name '*_make.log' -o -name '*_make-errors.log' \) -empty -delete

# Read the file ".gitignore" in the exe/ directory to get a list of files and folders to keep
#KEEP_EXE_ITEMS:=$(shell  sed -nE 's/^[ \t]*!(.*)/\1/p' $(ROOTPATH)/exe/.gitignore)
#List any files below that should NOT be deleted from the exe folder.
#  Directories MUST end in a slash (/)
#KEEP_EXE_ITEMS+=
#KEEP_EXE_FILES:=$(filter-out %/,$(KEEP_EXE_ITEMS)) # Remove items that end in slash
#KEEP_EXE_DIRS:= $(patsubst %/,%,$(filter %/,$(KEEP_EXE_ITEMS))) # remove items that do not end in slash

# Remove object files and executables.
.PHONY: realclean
realclean: clean cleanVS
	@rm -f RNA_class/*_class
	@echo 'Removing executables in exe/'
	@cd exe && rm -f $(ALL_PROGS)
	-@$(MAKE_CD) java_interface realclean
	-@$(MAKE_CD) $(PY_DIR)      realclean

#	@find exe -maxdepth 1 -type f  $(KEEP_EXE_FILES:%=! -name "%") -print  -delete
#	@echo "Keeping the following items (as listed in 'exe/.keep')"
#	@echo "Files: $(KEEP_EXE_FILES)"
#	@echo "Folders: $(KEEP_EXE_DIRS)"
#	find exe -maxdepth 1 -type f  $(patsubst %,! -name "%",$(KEEP_EXE_FILES))   -delete
#	find exe -mindepth 1 -maxdepth 1 -type d $(patsubst %, ! -name "%",$(KEEP_EXE_DIRS)) | xargs rm -rf

####  Run Regression Tests with `make test+<PROGRAM>` ####
# This is simply a shortcut for running Make commands from within the tests/ folder.
# Examples:	`make test+Fold` will run tests for the Fold program.
#          	`make test+all` will run all tests
test+%:
	@$(MAKE_CD) tests $(subst test+,,$@)

####  Get source-file dependencies of a Make target ####
# `make src+<PROGRAM>` will output a list of all cpp, cxx, and h source files 
# that <PROGRAM> depends on (even recursively).
# This is a simplified version of deps+<PROGRAM>
#  (Remove '@' to echo the command before executing it, e.g. for debugging.)
src+%:
	@target=$(subst src+,,$@); result='$$1'; $(~SHELL_GET_DEPS_CODE) | sort -u

####  Get object-file and source-file dependencies of a Make target ####
# `make deps+<PROGRAM>` will output all object file (*.o) dependencies of <PROGRAM>
#  along with the the object file's source file dependencies (cpp, cxx, and h).
#  (Remove '@' to echo the command before executing it, e.g. for debugging.)
deps+%:
	@target=$(subst deps+,,$@); result='$$2: $$1'; $(~SHELL_GET_DEPS_CODE) | sort -u


#shell code to list all Makefile dependencies for a given target
define ~SHELL_GET_DEPS_CODE
missing=$$(make -dn $$target | perl -ne "print if s/File '(.*\.o)' does not exist\./-W \$$1/"); \
regex="\s*Prerequisite '([^']*\.(?:cpp|cxx|h))' is (?:newer|older) than target '([^']*\.o)'\."; \
make -dn $$missing $$target | perl -ne "print if s/$$regex/$$result/"
endef

# Show some information about the OS and environment, for debugging.
show-env:
	@echo "uname: $$(uname -a)"
	@type -p &>/dev/null lsb_release && echo "OS Release: $$(lsb_release -a)" || true
	@echo "BASH_VERSION=$$BASH_VERSION"
	-@echo "PWD=$$PWD"
	@-type -p &>/dev/null cygpath && { echo -n "Cygpath: "; cygpath -wa .; } || true
	-@echo "SHELL=$$SHELL"
	-@echo "PATH=$$PATH"

################################################################################
#  Section 7: Make Help Message
################################################################################
define MAKE_HELP_MESSAGE
Make Goals: type `make <PROGRAM>` for any of the following:
  Serial (single-threaded) Programs:
    $(call fnToColumns,$(SERIAL_PROGS))

  SMP (multi-processor) Programs:
    $(call fnToColumns,$(SMP_PROGS))

  CUDA (GPU) Programs: 
    $(call fnToColumns,$(CUDA_PROGS))

  Note: Lower-Case variants of program names are also accepted
        E.g. `make oligowalk` will make OligoWalk.

  Other Programs:
    gui       Build the Java Graphical User Interface (GUI)
              and the corresponding native RNAstructure_GUI library.
    editor    Build the RNA Structure Drawing Editor (and dependencies).

    python-interface  -- Build the python interface library (used to 
              access RNAstructure tools from python scripts)

    unit-tests  --  Build and run some simple RNAstructure unit tests.


Other Commands and Tools:
  Build Multiple Programs:
    serial    Compile all the serial (single-threaded) executables.
    SMP       Compile all SMP (multi-processor) executables.
    CUDA      Compile all the CUDA (GPU) executables.
    all       Compile all serial and SMP executables.

  Do Cleanup:
    clean     Delete all object files (*.o) and Java class files (*.class)
    cleanVS   Delete Visual Studio object and temp files.
    realclean Runs clean, cleanVS, and deletes all files inside exe/

  Get Information:
    show-vars List Makefile variables (for debugging)
    help      Show this list of goals.

endef
