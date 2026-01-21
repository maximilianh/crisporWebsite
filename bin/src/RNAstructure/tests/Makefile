SHELL=/bin/bash
include ../common.mk

# Define the standard test-runner script that does the actual script execution 
# for most tests:
runtest=scripts/test-runner.sh

############################## SERIAL TESTS ####################################
# Add new SERIAL tests to the list below. 
# If a serial test also has an SMP counterpart, suffix its name with a tilde (~).
define SERIAL_TESTS
  AccessFold   AllSub
  bifold~   bipartition~
  CircleCompare   ct2dot   CycleFold
  design   dot2ct   draw   DuplexFold~   DynalignDotPlot   dynalign_ii~   dynalign~
  EDcalculator~ efn2~   EnergyPlot   EnsembleEnergy
  Fold~
  MaxExpect~   multilign~  NAPSS 
  oligoscreen~ orega~
  partition~   PARTS   ProbabilityPlot   ProbablePair~   ProbKnot~   ProbScan
  refold   RemovePseudoknots~  Rsample~
  scorer   ShapeKnots~   stochastic~
  TurboFold~  TurboHomology
  unit-tests   validate
  OligoWalk
  ETEcalculator~
endef

################################ SMP TESTS #####################################
# All serial tests listed above with a "~" suffix will ALSO be added as SMP tests.
# But if there are any SMP tests that do NOT have serial analogs, they can be listed below:
SMP_TESTS= # (Also includes "~" suffixed SERIAL tests from above.)

############################### CUDA TESTS #####################################
# List CUDA tests here:
CUDA_TESTS=Fold-cuda partition-cuda orega-cuda # New Cuda tests can be added here.

###################### MISCELLANEOUS TESTS #####################################
define GUI_TESTS
  AllSub-gui  bifold-gui  bipartition-gui  
  draw-gui  dynalign-gui  efn2-gui Fold-gui
  MaxExpect-gui OligoWalk-gui oligoscreen-gui
  ProbKnot-gui partition-gui
  stochastic-gui
endef
# draw-gui       -- no corresponding java functionality (java prints single pages)
# DuplexFold-gui -- no corresponding java functionality

###################### MISCELLANEOUS TESTS #####################################
MISC_TESTS=python

########################## LIST PREPARATION ####################################
# Filter out serial vs smp tests and replace the ~ suffix with -smp for the latter.
SMP_FROM_SERIAL=$(patsubst %~, %-smp, $(filter %~, $(strip $(SERIAL_TESTS))))  # Note: $(patsubst pattern,replacement,text)  $(filter pattern…,text)
SMP_TESTS    :=$(sort  $(strip $(SMP_TESTS))    $(SMP_FROM_SERIAL)) # Added `strip` because some versions of Make have trouble with the embedded linefeeds.
SERIAL_TESTS :=$(sort $(patsubst %~, %, $(strip $(SERIAL_TESTS))))
CUDA_TESTS   :=$(sort $(strip $(CUDA_TESTS)))
MISC_TESTS   :=$(sort $(strip $(MISC_TESTS)))
GUI_TESTS    :=$(sort $(strip $(GUI_TESTS)))
################################################################################

.DEFAULT_GOAL := help

.PHONY: all
all: serial smp cuda misc 

############ SERIAL TESTS ##############
.PHONY: SERIAL serial  $(SERIAL_TESTS)
SERIAL serial: $(SERIAL_TESTS)

$(SERIAL_TESTS):        # run the same recipe for all serial files
	@$(runtest) --serial $@

+%:  # allows running serial scripts (during development) that are not yet listed in SERIAL_TESTS by prepending '+' to the name
	$(runtest) --serial $(patsubst +%,%,$@)

############ SMP TESTS ##############
.PHONY: smp SMP $(SMP_TESTS)
SMP smp: $(SMP_TESTS)

$(SMP_TESTS): 
	@$(runtest) --smp $@

+%-smp:  # allows running smp scripts (during development) that are not yet listed in SMP_TESTS by prepending '+' to the name
	$(runtest) --smp $(patsubst +%,%,$@)


############ CUDA TESTS ##############
.PHONY: cuda CUDA $(CUDA_TESTS)
CUDA cuda: $(CUDA_TESTS)

$(CUDA_TESTS):         # run the same recipe for all cuda files
	@$(runtest) --cuda $@

############ GUI TESTS ##############
.PHONY: GUI gui $(GUI_TESTS)
GUI gui: $(GUI_TESTS)

$(GUI_TESTS): 
	@$(runtest) --gui $@

+%-gui:  # allows running gui scripts (during development) that are not yet listed in GUI_TESTS  by prepending '+' to the name
	$(runtest) --gui $(patsubst +%,%,$@)

############ MISC TESTS ##############
.PHONY: MISC misc $(MISC_TESTS)
MISC misc:  $(MISC_TESTS)

## Note:  There are no standard recipes for MISC tests. Each misc test should define its own recipe.

.PHONY: python
python: 
	$(SHELL) ./python/python_tests.sh

####### Tests for the testing system itself ####
system-tests:
	@$(runtest) --serial $@

# Simply show the results of the tests (i.e. pass/failed how many tests) and
# return an error code if !FAILED_TESTS.log exists and is not empty.
# In the code below, PL and FL are the logs of passed and failed tests.
# NP and NF are the number of lines in each respectively.
show-results:
	@lines() { local ct=0; [[ -s $$1 ]] && ct=$$(wc -l<"$$1" 2>/dev/null); echo $$ct; };\
	PL='!PASSED_TESTS.log'; FL='!FAILED_TESTS.log'; NF=$$(lines $$FL); NP=$$(lines $$PL);\
	if [[ -s $$FL ]]; then \
	  echo "$R------ FAILED $$NF TESTS (of $$((NF+NP))) ------$;";\
	  cat $$FL; echo "$R-------------------------------------$;"; exit 1;\
	else echo "$G------ PASSED ALL $$NP TESTS ------$;"; fi

############ Test-Related Environtment Variables #########

# Update a setting if the user passed a command-line flag to override it.
# $(call fnUpdateSetting,SETTING_VAR,FLAG_NAME)
fnUpdateSetting=$(if $(call fnGetOption,$2),$(eval export $1 := $($2)))

DEFAULT_TEST_OPTS = DumpErrors CleanExit
export TEST_OPTS ?= $(DEFAULT_TEST_OPTS)
DEFAULT_TEST_PROG_MAKE_ARGS= debug=on
export TEST_PROG_MAKE_ARGS ?= $(DEFAULT_TEST_PROG_MAKE_ARGS)

$(call fnUpdateSetting,INCLUDE_TESTS,include) # Update the INCLUDE_TESTS setting if include='...' is on the command-line
$(call fnUpdateSetting,EXCLUDE_TESTS,exclude) # Update the EXCLUDE_TESTS setting if exclude='...' is on the command-line
$(call fnUpdateSetting,TEST_OPTS,opts)        # Update the TEST_OPTS setting if opts='...' is on the command-line
$(call fnUpdateSetting,TEST_CONFIG,conf)      # Update the TEST_CONFIG setting if conf='...' is on the command-line

# By default, test programs are always rebuilt by Make if their source code 
# is newer than the existing executable. But if TEST_SKIP_REMAKE is non-empty
# and the program already exists, it will NOT be remade.
# On the command-line, this can be set with `remake=OFF`
ifneq ($(call fnIsOptionOFF,remake),)
  export TEST_SKIP_REMAKE=ALL
endif

# valgrind=ON turns on memory analysis using valgrind for all tests that support it.
ifneq ($(call fnIsOptionON,valgrind),)
  export TEST_PROG_ANALYZER=$(VALGRIND_TEST_CMD) # See below.
endif

## Common Test Analyzers ##
VALGRIND_TEST_CMD = valgrind --log-file=@TEST_valgrind.log --error-exitcode=199

############ Utility functions ##############
CLEAN_KEEP_FILES=Makefile .gitignore input-files.mk
FIND_CLEAN_KEEP_FILES:=$(patsubst %,! -name '%',$(CLEAN_KEEP_FILES))
.PHONY: clean
clean:
	@ # Move unknown files in the test directory to TRASH/
	@ rm \!*_TESTS.log 2>/dev/null && echo "Deleted test logs." || echo "No test logs to clean."
	@ list=$$(find . -maxdepth 1 -type f $(FIND_CLEAN_KEEP_FILES));\
	if [ -n "$$list" ]; then\
	    mkdir -p TRASH;\
	    find . -maxdepth 1 -type f $(FIND_CLEAN_KEEP_FILES) -exec mv -f '{}' TRASH \;  ;\
	    printf "$(R)The following files have been moved to ./TRASH/$;\n  %s\n"  "$${list//$$'\n'/$$'\n  '}";\
	else \
	    printf "No unknown files to clean.\n";\
	fi
	@ # Delete test output directories (*_OUTPUT)
	@ list=$$(find . -maxdepth 1 -type d -name '*_OUTPUT');\
	  if [ -n "$$list" ]; then\
	    rm -rf ./*_OUTPUT/;\
	    printf "$(R)The following directories and their contents have been deleted:$;\n  %s\n"  "$${list//$$'\n'/$$'\n  '}";\
	  else\
	    printf "No output directories to delete.\n";\
	  fi

.PHONY: clean-all
clean-all: clean+ clean-obj

.PHONY: clean+
clean+: clean
	@ if [[ -d TRASH ]]; then rm -rf TRASH && echo "Emptied TRASH."; else echo "No TRASH to empty."; fi

.PHONY: clean-obj
clean-obj:
	cd .. && make realclean
	cd ../java_interface && make realclean

$(C_BAR)=$(EMPTY)# Use $| to represent 'Empty', e.g. for the definition of LEFT_MARGIN below.

#### Definitions for formatted terminal output (bold, underline, colors etc) ###
R=$(TFMT_BR)#   Bold/Bright Red
G=$(TFMT_BG)#   Bold/Bright Green
Y=$(TFMT_Y)#    Dark Yellow
B=$(TFMT_BB)#   Bold/Bright Blue
D=$(TFMT_DIM)#   Dim text
_=$(TFMT_UND)#  Underline with "$_"
$(C_EXCL)=$(TFMT_BLD)#  Bold Text with "$!"
$(C_SEMI)=$(TFMT_RST)#  Reset format with "$;"


.PHONY: help
help:
	@$(eval export HELP_MESSAGE) echo "$$HELP_MESSAGE"
	@make -s list

#################### <-- 80 Characters --> #####################################
define HELP_MESSAGE
 $!RNAstucture: Regression Test Help$;
   $_Special `make` Targets:$;
   $G  list    $;  - Show a list of all tests.
   $G  all     $;  - Run all tests.
   $G  options $;  - Show a list of environment variables that influence tests.
                 (aka $Gopts$;)

   $G  clean   $;  - Cleanup output from previous tests (in *_OUTPUT/ dirs). Also
                 moves any misc/unknown files in the test directory to TRASH/
   $G  clean+  $;  - Do standard clean then also delete files in TRASH/ (which are
                 normally preserved for caution).
   $G  clean-obj$; - Delete all object and binary files (i.e. 'realclean').
   $G  clean-all$; - Do all cleaning (i.e. `$!make clean clean+ clean-obj$;`).

  $_Test-Category Targets:$;
   $G  serial  smp  cuda  misc $;

  $_Targets can be combined:$;
     e.g. $Y$$ make clean serial smp $;

endef

.PHONY: list lsit 
list lsit:  # List all tests
	@$(eval export MAKE_TEST_LIST) echo -e "$$MAKE_TEST_LIST"

define MAKE_TEST_LIST
$!RNAstucture: Individual Test Targets$;
  $_Serial Tests:$;
    $(call fnToColumns,$(SERIAL_TESTS))
  $_SMP Tests:$;
    $(call fnToColumns,$(SMP_TESTS))
  $_Cuda Tests:$;
    $(call fnToColumns,$(CUDA_TESTS))
  $_Java GUI Tests:$;
    $(call fnToColumns,$(GUI_TESTS))
  $_Miscellaneous Tests:$;
    $(call fnToColumns,$(MISC_TESTS))
endef


.PHONY: options opts
options opts:  # List environment vars
	@$(eval export SPECIAL_TEST_VARS) echo "$$SPECIAL_TEST_VARS" | less -R

define SPECIAL_TEST_VARS
$!RNAstucture: Test-Related Environment Variables$;
  $_DATAPATH$;
    The path to the data_tables folder that contains thermodynamic parameter 
    files. This is set automatically if it is empty or invalid.
  $_EXCLUDE_TESTS$;  (aka $(Y)exclude=...$;)
    A list of test names (which can include wildcards) that should be skipped.
      e.g. $(Y)export EXCLUDE_TESTS="Fold_without_options  *_dna_option  partition_shape_*"$;
      e.g. $(Y)make Fold $!exclude=$;$(Y)'without_options dna_option'$;  $D(as `make` argument)$;
  $_INCLUDE_TESTS$;  (aka $(Y)include=...$;)
    A list of test names (which can include wildcards) that should be included in the test run.
    If this is unset or empty, all tests are run. Otherwise, only tests that match a
    name in the list are run.  Importantly, EXCLUDE_TESTS has higher precedence than
    INCLUDE_TESTS, so a file that matches in both will be excluded.
      e.g. $(Y)export INCLUDE_TESTS="*_without_options  *_dna_option"$;
      e.g. $(Y)make Fold $!include=$;$(Y)'without_options dna_option'$;  $D(as `make` argument)$;
  $_TEST_OPTS$;  (aka $(Y)opts=...$;)
    This environment variable contains a list of flags that control the test
    system, especially regarding verbosity of output. These flags are described
    below. The flags are NOT case sensitive, and 1-letter forms can be used and
    combined with or without spaces.
      e.g. $(Y)export TEST_OPTS="KeepAll DumpErrors ShowOutput DebugGui CleanExit"$;
      e.g. $(Y)$Yexport TEST_OPTS="kdsgc"$; $D(same as above)$;
      e.g. $(Y)make Fold $!opts=$;$(Y)KQX$;      $D(as `make` argument)$;
    The DEFAULT value of TEST_OPTS (if not set explicitly) is:
        "$(Y)$(DEFAULT_TEST_OPTS)$;"
    $!KeepAll (K)$; - Keep all test output (including stdout, error, and diff 
      files etc.) for ALL tests. (Usually this is only done when a test fails.)
      These files are stored in $!<TESTNAME>$;_OUTPUT/
    $!ListErrors (L)$; - List all "error" and "diff" files that indicate a failed
      test (as detected by checkErrorFiles). This can help to quickly identify 
      the source of each error.
    $!DumpErrors (D)$; - Display the first few lines of each of the error files 
      listed by $!ListErrors$;. This can help quickly identify the specific
       details of a potentially failed test. (ListErrors is also implied).
    $!ShowOutput (S)$; - Display the stdout from commands invoked with "runTest".
      (It is also logged to file whether or not this flag is present.)
    $!ExitOnError (X)$; - Immediately exit a test script after a failed test. 
      Note that the status of the $!CleanExit$; option affects whether or not
      $!Make$; will continue to run susbequent tests.
    $!CleanExit (C)$; - Exit a test script with a clean exit code (0), even if
      errors occurred (allowing $!Make$; to run subsequent test scripts).
      This suppresses the error code that is usually returned from a test script
      when errors have been detected.
    $!DebugGui (G)$; - Show additional information during GUI tests, including the
      full gui-script as well as internal tracing of the gui-script as it runs.
    $!Quiet (Q)$; - Suppress all output except test results and warnings.
      This setting will not be completely effective if options that INCREASE 
      verbosity (e.g. ListErrors, DebugGui, etc) are also set.
    $!NoWarnings (W)$; - Suppress warnings. This is discouraged because
      Warnings can indicate test-authoring errors and skipped tests.
  $_DIFF_FLAGS$;
    Special flags to pass to the $!diff$; command when comparing test output.
    On Windows, this is set to $!--strip-trailing-cr$; to ignore differences 
    that are due to differing line endings.
  $_DIFF_CMD$;
    The command to use to compare test output with reference files. By default,
    this is "$!diff$;", but the user may have an alternate preference.
  $_RUNTEST_CMD$;
    The runFullTest command internally invokes $(Y)runTest$; to run the executable
    under test. Setting RUNTEST_CMD to the name of another command or function
    will cause it to be used instead. For example, this is used to run Java GUI
    test scripts using the same command as the Text interfaces (so runFullTest
    internally calls $(Y)runGuiTest$; instead of runTest).
  $_TEST_SKIP_REMAKE$;  (aka $(Y)remake=off$;)
    By default, test programs are always rebuilt by Make if their source code 
    is newer than the existing executable. But if TEST_SKIP_REMAKE is non-empty
    and a program already exists, it will NOT be remade.
    On the command-line, this can be set with `remake=OFF`
      e.g. $(Y)export TEST_SKIP_REMAKE=1$;
      e.g. $(Y)make Fold $!remake=$;$(Y)off$;      $D(as `make` argument)$;
  $_valgrind=on$;
    Passing this on the command-line to Make will setup valgrind to monitor each
    test program. Internally, this is done by setting the $_TEST_PROG_ANALYZER$;
    environment variable to a value similar to the following:
      $!$(VALGRIND_TEST_CMD)$;  $D(*program and args...*)$;
  $_TEST_PROG_MAKE_ARGS$;
    These arguments are passed to `make` when rebuilding test programs.
    The default is: $!$(DEFAULT_TEST_PROG_MAKE_ARGS)$;
endef

####### Define default message/action for files with no matching rule. #######
.DEFAULT:
	@echo "Unknown target \"$@\".  For a list of goals, type \"make help\"."; \
	exit 1
##############################################################################
