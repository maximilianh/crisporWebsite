##########################################################################################
# Define the relative path to the RNAstructure root directory
# and include all macro, dependency, and variable definitions.
##########################################################################################
ROOTPATH=..
JAVA=yes

include ${ROOTPATH}/compiler.h
include ${ROOTPATH}/library_defines.h
include ${ROOTPATH}/dependencies.h
include compiler-java.h

##########################################################################################
#  Define the help message to list known Makefile goals
##########################################################################################
.DEFAULT_GOAL=help

#############################<--80-chars-->#####################################
help:
	@$(eval export MAKEFILE_HELP_MSG) echo "$$MAKEFILE_HELP_MSG"
define MAKEFILE_HELP_MSG
 Type "make <goal>..." to perform one of the actions listed below.

   native - Compile the native (C++) class library ($(RNASTRUCTURE_LIBRARY)).
   gui    - Build the Java user interface ($(JARFILE_GUI)).
   editor - Build the RNA Structure Drawing Editor ($(JARFILE_EDITOR)).
   tester - Build GUI testing platform ($(JARFILE_TESTER)).
   swig   - Generate SWIG glue code for JNI (only required when breaking changes
            have been made to the public c++ interface.)
 
   all    - Perform all actions listed above (swig, native, gui, tester, editor)
   help   - Show this list of goals.
 
 Note: `gui`, `editor` and `tester` each compile the native library if it is not
   already built. You can build the JAR file without the native lib by appending
   appending the suffix "-jar" (e.g. `make gui-jar`) or you can compile the
   Java classes (without the JAR) by apending "-java", e.g. `make gui-java`.
endef

################################################################################
# RNAstructure Java Modules
################################################################################
# These are logical units of code that correspond either to shared libraries or 
# top-level programs. Each is stored in a separate subfolder of the source tree.
# Each has a short makefile name (e.g. UTIL) and a "real" name (e.g. Utilities)
# that corresponds to a folder and/or package name.

# The Utilities module contains code common to many of the other modules
MODULE_UTIL 	= Utilities
# RNAstructureUI is the main RNAstructure GUI application.
MODULE_GUI 		= RNAstructureUI
# GUITester runs unit tests on the RNAstructure GUI
MODULE_TESTER 	= GUITester
# StructureEditor is the RNA Structure Drawing Editor standalone application.
MODULE_EDITOR 	= StructureEditor
# RNAstructure is the C++ backend RNA class
MODULE_RNACLS 	= RNAstructure
# List of all modules. This is used by fn_EvalForeachModule to perform tasks and setup variables for all listed modules.
ALL_MODULES = UTIL GUI TESTER EDITOR RNACLS

################################################################################################
#                           Utility Functions
################################################################################################
# Define some function MACROS used by $(call fn_...)

# fn_FindJavaSrc -- Recursively list java source files in the specified folder.
# Usage: $(call fn_FindJavaSrc,<SEARCH_DIR>)
# This will call `find <SEARCH_DIR> -type f -name '*.java'`  and will return the results of the search.
fn_FindJavaSrc = $(shell find '$1' $(SHELL_FIND_EXCLUDE) -type f -name '*.java')

# fn_EvalForeachModule -- Evaluate an expression once for each module in 
#    ALL_MODULES with the text "{?}" replaced with the short-name of the module
#    (e.g. GUI)
# Usage: $(call fn_EvalForeachModule,<EXPRESSION_WITH{?}>)
# See notes in `fn_EvalForeach` about double-escaping ($$) variables.
fn_EvalForeachModule=$(foreach ~~tmp_var,$(ALL_MODULES),$(eval $(subst {?},$(~~tmp_var),$1)))

################################################################################################
#                  Source and Output folders for Java compilation
################################################################################################
# Java organizes code into a tree-like structure of packages. BASE_PKG defines the root package prefix for all RNAstructure code.
BASE_PKG=ur_rna
# Java requires that compiled Java class files are stored in subdirectories correspond to their packages.
# Some Java IDEs also prefer or require source files to be stored in a directory  hierarchy that matches each file's package (with dot (.) replaced by slash (/))
BASE_PKG_SUBDIR=$(subst /,.,$(BASE_PKG))#### same as  BASE_PKG except that all dots (.) are replaced with commas.

# The root directory for all java-related stuff, including SWIG etc.
JAVA_ROOT_DIR =.# (CURRENT DIRECTORY)
# The root directory for source code
JAVA_SRC_DIR  =$(JAVA_ROOT_DIR)/src
# A subdirectory of JAVA_SRC_DIR where the real source code begins. (See notes above for BASE_PKG_SUBDIR). Since the root package is "ur_rna", all code has to be in src/ur_rna/
BASE_PKG_SRC_DIR =$(JAVA_ROOT_DIR)/src/$(BASE_PKG_SUBDIR)
# Set the source directory for each specific module listed in ALL_MODULES
$(call fn_EvalForeachModule,SRC_DIR_{?}=$$(BASE_PKG_SRC_DIR)/$$(MODULE_{?}))  ### e.g. SRC_DIR_GUI=<java_root>/src/ur_rna/RNAstructureUI
# Set the package subdirectory for each specific module listed in ALL_MODULES. This is the same as SRC_DIR_XXX, but is relative to JAVA_SRC_DIR
$(call fn_EvalForeachModule,PKG_SUBDIR_{?}=$$(BASE_PKG_SUBDIR)/$$(MODULE_{?}))### e.g. PKG_SUBDIR_GUI=ur_rna/RNAstructureUI

# Directory where resources (icons, images, etc) are stored.
RES_SRC_DIR  =$(JAVA_ROOT_DIR)/resources
# Set the module-specific subdirectory where resources (icons, images, etc) are stored (for each module).
$(call fn_EvalForeachModule,RES_SRC_DIR_{?}=$$(RES_SRC_DIR)/$$(MODULE_{?}))  ### e.g. RES_SRC_DIR_GUI=<bin>/java/ur_rna/RNAstructureUI

# Root output directory
BIN_OUT_ROOT      =${ROOTPATH}/exe
# Subdirectory of BIN_OUT_ROOT for all java output
JAVA_OUT_DIR      =${BIN_OUT_ROOT}/java
# Subdirectory of BIN_OUT_ROOT for java unit testing output
JAVA_TEST_OUT_DIR =${BIN_OUT_ROOT}/java_test
# Subdirectory where all real output goes. (See note for BASE_PKG_SRC_DIR about Java package requirements)
BASE_PKG_OUT_DIR   =${BIN_OUT_ROOT}/java/$(BASE_PKG_SUBDIR)
# Set the output directory for each specific module listed in ALL_MODULES
$(call fn_EvalForeachModule,OUT_DIR_{?}=$$(BASE_PKG_OUT_DIR)/$$(MODULE_{?}))  ### e.g. MODULE_OUT_GUI=<bin>/java/ur_rna/RNAstructureUI
# Set the directory where resources (icons, images, etc) are output, for each different module.
$(call fn_EvalForeachModule,RES_OUT_DIR_{?}=$$(OUT_DIR_{?})/resources)    ### e.g. RES_OUT_DIR_GUI=<bin>/java/ur_rna/RNAstructureUI

# Source directory for SWIG-relatd files (*.i etc)
SWIG_SRC    = $(JAVA_ROOT_DIR)/SWIG
# Source directory for scripts. Some scripts are for Makefile use or maintenance. Others are for distribution in releases.
SCRIPTS_SRC = $(JAVA_ROOT_DIR)/scripts

############## List ABBOT libs required for running GUITester #############################################
JAVA_LIB_DIR=$(JAVA_ROOT_DIR)/lib
#specify the required dependency libraries for java compilation
JAVA_CLASSPATH:=.
ABBOT_LIBS:=$(JAVA_LIB_DIR)/abbot.jar $(wildcard $(JAVA_LIB_DIR)/gnu-regexp-*.jar $(JAVA_LIB_DIR)/jdom-*.jar $(JAVA_LIB_DIR)/junit-*.jar)
ABBOT_LIBS_OUT:=$(patsubst  $(JAVA_LIB_DIR)/%,$(BIN_OUT_ROOT)/lib/%,$(ABBOT_LIBS))

space:=$(C_SPACE)
JAVA_CLASSPATH:=$(subst $(space),$(OS_PATH_SEP),$(JAVA_CLASSPATH))#Separate paths with ":" on linux/mac and with ";" on Windows.
TESTER_CLASSPATH:=$(subst $(space),$(OS_PATH_SEP),$(ABBOT_LIBS))#Separate paths with ":" on linux/mac and with ";" on Windows.

#################################################################
# Scripts to include in Software Releases for use by end-users
#################################################################
GUI_SCRIPTS=RNAstructureScript
TESTER_SCRIPTS=run-gui-test
EDITOR_SCRIPTS=StructureEditor
ifeq (${OPSYSTEM},Windows)
  GUI_SCRIPTS+=RNAstructure.bat RNAstructure_CLI.bat
  EDITOR_SCRIPTS+=StructureEditor.bat
endif
# Prefix the output directory (../exe/) to the scripts listed above.
GUI_SCRIPTS:=	$(addprefix $(BIN_OUT_ROOT)/,$(GUI_SCRIPTS))
TESTER_SCRIPTS:=$(addprefix $(BIN_OUT_ROOT)/,$(TESTER_SCRIPTS))
EDITOR_SCRIPTS:=$(addprefix $(BIN_OUT_ROOT)/,$(EDITOR_SCRIPTS))

#################################################################
#  Definition of source files and dependencies for each Java Module
#################################################################
# Recursively list *.java files in each module's subdirectory.
$(call fn_EvalForeachModule,SRC_FILES_{?}=$$(call fn_FindJavaSrc,$$(SRC_DIR_{?})))#  e.g. SRC_FILES_UTIL=./src/ur_rna/Utilities/**/*.java
# Convert each source file path into a class file, rooted in the output directory.
$(call fn_EvalForeachModule,CLASS_FILES_{?}=$$(patsubst $$(JAVA_SRC_DIR)/%.java,$$(JAVA_OUT_DIR)/%.class,$$(SRC_FILES_{?})))#  e.g. CLASS_FILES_UTIL=../exe/java/ur_rna/Utilities/**/*.class

# Resource files (images, sounds etc) for RNAstructure GUI
RES_SRC_FILES_GUI:=$(call fnFind,'$(RES_SRC_DIR_GUI)',-type f)
RES_OUT_FILES_GUI:=$(patsubst $(RES_SRC_DIR_GUI)/%,$(RES_OUT_DIR_GUI)/%,$(RES_SRC_FILES_GUI))
# Resource files for StructureEditor
~RES_EXCLUDE = -path '$(RES_SRC_DIR_EDITOR)/images/icon-src' # remove SVG source files for icons from the resources
RES_SRC_FILES_EDITOR:=$(call fnFind,'$(RES_SRC_DIR_EDITOR)',-type f,$(~RES_EXCLUDE))
RES_OUT_FILES_EDITOR:=$(patsubst $(RES_SRC_DIR_EDITOR)/%,$(RES_OUT_DIR_EDITOR)/%,$(RES_SRC_FILES_EDITOR))

##########################################################################################
## PHONY Goals
##########################################################################################
# Always run these goals, even if a file of the same name already exists.
.PHONY: help all clean realclean debug swig swig-rna depend #native java gui GUI tester

##########################################################################################
## Top-Level Goals
##########################################################################################

JARFILE_GUI=RNAstructure.jar
JARFILE_TESTER=GUITester.jar
JARFILE_EDITOR=StructureEditor.jar
JARFILE_OUT_GUI=$(BIN_OUT_ROOT)/$(JARFILE_GUI)
JARFILE_OUT_TESTER=$(BIN_OUT_ROOT)/$(JARFILE_TESTER)
JARFILE_OUT_EDITOR=$(BIN_OUT_ROOT)/$(JARFILE_EDITOR)

ALL_JARFILE_NAMES=$(JARFILE_GUI) $(JARFILE_TESTER) $(JARFILE_EDITOR)
ALL_JARFILE_PATHS=$(patsubst %,$(BIN_OUT_ROOT)/%,$(ALL_JARFILE_NAMES))

# Rule to build the C++ part of the RNAstructure Java interface only.
native: $(BIN_OUT_ROOT)/${RNASTRUCTURE_LIBRARY}
# Rule to build the Java interface, C++ and Java code.
GUI gui: gui-jar     native
tester:  tester-jar  native
editor:  editor-jar  native
all:    swig native gui tester editor

#Specify that each JAR file should be deleted if there is an error making it. (Otherwise the jar is created but may not have all required files inside.)
.DELETE_ON_ERROR: $(ALL_JARFILE_PATHS)

##########################################################################################
## Compilation Rules
##########################################################################################

# IMPORTANT:  The '$(fn_VerifyJDK)' call should be added to any recipe that requires 
#   the JDK (either javac, jar, or jni.h/jni_md.h etc). When it is expanded, the JDK is 
#   auto-detected and verified. 
#   The result is cached, so auto-detection and verification will happen at most once per
#   makefile invocation.

# Test conditions for SMP usage.
ifeq ($(OPSYSTEM),Windows)
  # if using windows, do not use SMP for 32-bit arch. (Causes UnsatisfiedlinkError. Not sure exactly what isn't compatible.)
  UseSMP=$(findstring 64,$(TARGET_ARCH))
else
  UseSMP=$(filter-out Mac,$(OPSYSTEM))# do not use SMP on mac.
endif

#$(info UseSMP=$(UseSMP) TARGET_ARCH=$(TARGET_ARCH) OPSYSTEM=$(OPSYSTEM))
ifneq ($(UseSMP),)
  #$(info USING SMP)
  # When compiling the RNAstructure_GUI library, use OpenMP (SMP) and pThreads, but only if using Linux or Windows (Not Mac)
  JAVA_LIBRARY_FILES=$(JAVA_LIBRARY_FILES_SMP)
  # COMPILE_JNI  - The default recipe to compile the c++ files (.cpp and .cxx) 
  COMPILE_JNI  = ${COMPILE_SMP} ${JNI_CXXFLAGS}
  # LINK_JNI_LIB - rule for linking object files to make the shared library that interfaces with Java. 
  LINK_JNI_LIB = ${LINK_SMP_DEPS}
endif


### Rule for compiling *.java source files into *.class files.
# Note: the -sourcepath flag is usually not required if the classpath (-cp) includes the source directory. 
# However if it is NOT specified, and javac finds source files (*.java) inside a jar file in the CP, 
# then it may try to recompile them, leading to compilation errors related to third-party pre-compiled code
$(JAVA_OUT_DIR)/%.class: $(JAVA_SRC_DIR)/%.java
	@$(fn_VerifyJDK)#Verify the JDK.
	$(!MAKE_OUTDIR) # make the output directory for the class file
	"${JXX}" ${JXXFLAGS} -cp "$(JAVA_SRC_DIR)$(OS_PATH_SEP)$(JAVA_CLASSPATH)" -d "$(JAVA_OUT_DIR)" -sourcepath $(JAVA_SRC_DIR)  $<

# Rule to compile the native library (RNAstructure.dll, libRNAstructure.so, etc)
$(BIN_OUT_ROOT)/${RNASTRUCTURE_LIBRARY}: LDFLAGS+=${JNI_CXXFLAGS} ${LIBFLAGS}
$(BIN_OUT_ROOT)/${RNASTRUCTURE_LIBRARY}: $(JAVA_LIBRARY_FILES)
	@$(fn_VerifyJDK)#Verify the JDK.
	$(!MAKE_OUTDIR)
	@#Note that CXXFLAGS are modified to include ${JNI_CXXFLAGS} and ${LIBFLAGS}
	${LINK_JNI_LIB}

# Rule to compile all SWIG-generated c++ files. (COMPILE_JNI is defined in compiler-java.h)
$(SWIG_OUT)/%.o: $(SWIG_SRC)/%.cpp $(SWIG_SRC)/%.h
	@$(fn_VerifyJDK)#Verify the JDK.
	$(!MAKE_OUTDIR)
	${COMPILE_JNI} $<

# Rule to compile all SWIG-generated cxx files. (COMPILE_JNI is defined in compiler-java.h)
$(SWIG_OUT)/%.o: $(SWIG_SRC)/%.cxx
	@$(fn_VerifyJDK)#Verify the JDK.
	$(!MAKE_OUTDIR)
	${COMPILE_JNI} $<

$(CLASS_FILES_TESTER): JAVA_CLASSPATH:=$(JAVA_CLASSPATH)$(OS_PATH_SEP)$(TESTER_CLASSPATH)# Temporarily set the JAVA_CLASSPATH to include abbot libraries.
$(CLASS_FILES_TESTER): $(ABBOT_LIBS)

gui-java: 		$(CLASS_FILES_GUI)
util-java: 		$(CLASS_FILES_UTIL)
editor-java: 	$(CLASS_FILES_EDITOR)
rna-java: 		$(CLASS_FILES_RNACLS)
tester-java: 	$(CLASS_FILES_TESTER)

gui-jar:   	$(JARFILE_OUT_GUI)   	$(GUI_SCRIPTS)
tester-jar:	$(JARFILE_OUT_TESTER)	$(TESTER_SCRIPTS)	$(ABBOT_LIBS_OUT)
editor-jar:	$(JARFILE_OUT_EDITOR)	$(EDITOR_SCRIPTS)

# Rule to build the RNAstructure.jar
$(JARFILE_OUT_GUI): $(CLASS_FILES_UTIL)  $(CLASS_FILES_RNACLS)  $(CLASS_FILES_GUI)  \
                    $(RES_OUT_FILES_GUI)  $(JAVA_ROOT_DIR)/Manifest.txt
	@$(fn_VerifyJDK)#Verify the JDK.
	$(!MAKE_OUTDIR) #Create the output directory.
	@# Build the JAR file
	"$(JAR)" cvmf $(JAVA_ROOT_DIR)/Manifest.txt $@ \
		-C $(JAVA_OUT_DIR)  $(PKG_SUBDIR_UTIL)  \
		-C $(JAVA_OUT_DIR)  $(PKG_SUBDIR_RNACLS) \
		-C $(JAVA_OUT_DIR)  $(PKG_SUBDIR_GUI)   \
		 > $@.log

$(JARFILE_OUT_TESTER): $(CLASS_FILES_UTIL)  $(CLASS_FILES_RNACLS)  $(CLASS_FILES_GUI)  $(CLASS_FILES_TESTER) \
                       $(RES_OUT_FILES_GUI)  $(JAVA_ROOT_DIR)/GUITester_Manifest.txt
	@$(fn_VerifyJDK)#Verify the JDK.
	$(!MAKE_OUTDIR) #Create the output directory.
	@# Build the JAR file
	"$(JAR)" cvmf $(JAVA_ROOT_DIR)/GUITester_Manifest.txt $@ \
		-C $(JAVA_OUT_DIR)  $(PKG_SUBDIR_UTIL)   \
		-C $(JAVA_OUT_DIR)  $(PKG_SUBDIR_TESTER) \
		 > $@.log

# Rule to build the StructureEditor java program
$(JARFILE_OUT_EDITOR): $(CLASS_FILES_UTIL)  $(CLASS_FILES_RNACLS)  $(CLASS_FILES_EDITOR)  \
                       $(JAVA_ROOT_DIR)/StructureEditor_Manifest.txt $(RES_OUT_FILES_EDITOR)
	@$(fn_VerifyJDK)#Verify the JDK.
	$(!MAKE_OUTDIR) #Create the output directory.
	@# Build the JAR file
	"$(JAR)" cvmf $(JAVA_ROOT_DIR)/StructureEditor_Manifest.txt $@ \
		-C $(JAVA_OUT_DIR)  $(PKG_SUBDIR_UTIL)   \
		-C $(JAVA_OUT_DIR)  $(PKG_SUBDIR_RNACLS) \
		-C $(JAVA_OUT_DIR)  $(PKG_SUBDIR_EDITOR) \
		 > $@.log

ifeq ($(OPSYSTEM),Mac)
  # if this is a mac, add an additional '' to the sed command
  $(GUI_SCRIPTS) $(TESTER_SCRIPTS) $(EDITOR_SCRIPTS): $(BIN_OUT_ROOT)/%: $(SCRIPTS_SRC)/%
	@$(fn_VerifyJDK)#Verify the JDK.
	cp -fp '$<' '$@'
	sed -i '' "s,__JAVA_BIN_DIR__,$(JDK_BIN),g" '$@'
	chmod u+x '$@'
else
  #for Windows and Linux:
  $(GUI_SCRIPTS) $(TESTER_SCRIPTS) $(EDITOR_SCRIPTS): $(BIN_OUT_ROOT)/%: $(SCRIPTS_SRC)/%
	@$(fn_VerifyJDK)#Verify the JDK.
	cp -fp '$<' '$@'
	sed -i "s,__JAVA_BIN_DIR__,$(JDK_BIN),g" '$@'
	chmod u+x '$@'
endif


#	cp -fp '$<' '$@' && chmod u+x '$@'

!CPRES=$(!MAKE_OUTDIR); cp -fp '$<' '$@' && echo "cp resource $@"
$(RES_OUT_FILES_GUI):    $(RES_OUT_DIR_GUI)/%:    $(RES_SRC_DIR_GUI)/%
	$(!CPRES)

$(RES_OUT_FILES_EDITOR): $(RES_OUT_DIR_EDITOR)/%: $(RES_SRC_DIR_EDITOR)/%
	$(!CPRES)

# Rule to build SWIG files.
# This rule is almost never necessary for the end user.
EDITOR_SWIG_OUT_DIR=$(SRC_DIR_RNACLS)/backend
swig:
	# Removing old invalid swig files.
	find ${EDITOR_SWIG_OUT_DIR} $(SRC_DIR_GUI) -name 'SWIGTYPE_p_*.java' -delete
	@# Remove old generated swig files.
	@#grep -l 'This file was automatically generated by SWIG' ${EDITOR_SWIG_OUT_DIR}/*.java | tr '\n' '\0' | xargs -0 -I FILE mv -f FILE FILE.old
	# Generating SWIG files for RNAstructure GUI
	swig -java -c++ -package $(BASE_PKG).$(MODULE_GUI)               -outdir $(SRC_DIR_GUI)               ${SWIG_SRC}/RNAstructureBackendCalculator.i
	swig -java -c++ -package $(BASE_PKG).$(MODULE_GUI).drawing.proxy -outdir $(SRC_DIR_GUI)/drawing/proxy ${SWIG_SRC}/drawing/StructureBackend.i
	swig -java -c++ -package $(BASE_PKG).$(MODULE_GUI).drawing.proxy -outdir $(SRC_DIR_GUI)/drawing/proxy ${SWIG_SRC}/drawing/DotPlotBackend.i
	@rm -f $(SRC_DIR_GUI)/*Proxy.java
	@rm -f $(SRC_DIR_GUI)/drawing/proxy/*Proxy.java
	@#bash ${SWIG_SRC}/remove-default-args.sh "${SWIG_SRC}/RNA.i" > "${SWIG_SRC}/remove-default-args.i"
	@echo > "${SWIG_SRC}/remove-default-args.i"
	# Generating SWIG files for Structure Editor
	swig -java -c++ -package ur_rna.RNAstructure.backend -outdir ${EDITOR_SWIG_OUT_DIR} ${SWIG_SRC}/RNA.i # this is for the structure drawing editor and alignment-editor	
	@#rm -f ${EDITOR_SWIG_OUT_DIR}/RNABackend.java
	# Searching for SWIGTYPE_p_XXX files (which may indicate that an internal type was incorrectly exposed).
	@found=($$(find ${EDITOR_SWIG_OUT_DIR} $(SRC_DIR_GUI) -name 'SWIGTYPE_p_*.java')); \
	    if [[ $$found ]]; then printf >&2 "$(TFMT_BR)Error: internal SWIGTYPE file generated:$(TFMT_BY)%s$(TFMT0)\n" "$${found[@]}"; exit 1; fi

##########################################################################################
# Special Goals (debug, cleanup etc)
##########################################################################################
DEBUG_VARS=BASE_PKG_SUBDIR MODULE_RES_SRC MODULE_OUT MODULE_RES_OUT #SHOW debugging information
DEBUG_VARS_MULTILINE=SRC_FILES CLASS_FILES RES_OUT_FILES_SRC RESOURCE_FILES
debug: ; @:
	$(foreach MODULE_NAME,$(ALL_MODULES), \
	  $(info ------- DEBUG INFO for Module: $(MODULE_NAME) -------) \
	  $(foreach var,$(DEBUG_VARS),$(info $(var) = $($(var)_$(MODULE_NAME)))) \
	  $(foreach var,$(DEBUG_VARS_MULTILINE),\
	    $(info $n---$(var)_$(MODULE_NAME)---$n$(call fnSpacesToLines,$($(var)_$(MODULE_NAME)))) \
	  ) \
	  $(info ------- END of INFO for Module: $(MODULE_NAME) -------$n) \
	)

# Cleanup.
# Object cleanup removes all temporary build objects.
# Executable cleanup removes all possible executables.

# Rule to clean all possible object files and compiled Java class files.
clean: clean-java clean-native

clean-java:
	@[ -d $(JAVA_OUT_DIR) ]       && find $(JAVA_OUT_DIR) -name '*.class' -delete       || echo 'No java classes to clean.'
	@[ -d $(JAVA_TEST_OUT_DIR) ] && find $(JAVA_TEST_OUT_DIR) -name '*.class' -delete || echo 'No java test classes to clean.'
	rm -f $(BIN_OUT_ROOT)/*.jar.log

clean-native:
	@echo 'Deleting all Java native library object files.'
	@rm -f $(sort ${JAVA_LIBRARY_FILES} $(JAVA_LIBRARY_FILES_SMP))

# Rule to clean all possible objects and executables.
realclean: clean
	rm -f $(addprefix $(BIN_OUT_ROOT)/,\
        ${RNASTRUCTURE_LIBRARY} \
        RNAstructure.jar Utilities.jar GUITester.jar StructureEditor.jar \
        $(GUI_SCRIPTS) $(TESTER_SCRIPTS) $(EDITOR_SCRIPTS) \
      )
	rm -rf $(JAVA_OUT_DIR)  $(BIN_OUT_ROOT)/java_test  $(BIN_OUT_ROOT)/lib

.PHONY: abbot
abbot:  $(ABBOT_LIBS)

$(ABBOT_LIBS_OUT): $(BIN_OUT_ROOT)/lib/%: $(JAVA_LIB_DIR)/%
	@$(!MAKE_OUTDIR)
	cp -fp '$^' '$@'

$(JAVA_LIB_DIR)/%.jar:
	bash $(JAVA_LIB_DIR)/get-abbot
