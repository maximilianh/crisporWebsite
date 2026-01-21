################################################################################
# This file includes Makefile variables and function definitions to 
# detect and configure the Java compiler and the c++ compiler as it pertains 
# to building the Java Native Interface component of RNAstructure.
#
# The environment variables JAVA_HOME or JDK_BIN and JDK_INCLUDE can be set to
# speed up or circumvent auto-detection.
################################################################################

# Include 'java-config.mk', which is a user-editable file that sets 
#   JAVA_HOME or JDK_BIN and JDK_INCLUDE
JAVA_CONFIG_FILE=java-config.make
-include $(JAVA_CONFIG_FILE)
# List the include-file here, with no recipe to avoid an error if not found.
$(JAVA_CONFIG_FILE):

################################################################################
# Java Compiler Configuration
################################################################################
# Set java compiler and jar-archiver (usually included in the PATH, except possibly on Windows)
JXX=$(JDK_BIN)/javac
JAR=$(JDK_BIN)/jar
# JXXFLAGS - Java compiler flags
JXXFLAGS = -encoding UTF-8 -Xlint:-options 

# Modern installations of Java typically have the binary commands (e.g. javac) 
# as well as the JNI header files (e.g. jni.h) all housed under a single parent 
# directory which is *sometimes* specified in the environement variable JAVA_HOME.
# 
# In this scheme, the binary commands are located in $JAVA_HOME/bin and the headers
# are in $JAVA_HOME/include  and  $JAVA_HOME/include/$JNI_OS_NAME is an OS-dependent 
# subfolder that contains jni_md.h, another required C++ header file.
# (JNI_OS_NAME is defined later in this file and is "darwin" on Mac, 
# "linux" on GNU/Linux, and "win32" on Microsoft Windows)
#
# In addition to the JAVA_HOME location, the binaries and/or headers can be found
# in other locations--usually via symlinks to the files under JAVA_HOME.
# Thus, if JAVA_HOME has been defined or can be determined by searching etc, then the 
# paths to the binary files and header files can usually be determined from it.
# 
# We use a different variable JDK_HOME so that it doesn't conflict with the 
# JAVA_HOME environment variable, which may be used by sub-processes.
JDK_HOME ?= $(JAVA_HOME)
# JDK_BIN -- The path to the JDK binary commands (e.g. javac, jar etc) This variable can 
#   usually be left alone, as long as the typical java installation has been performed. 
JDK_BIN  ?= $(JDK_HOME)/bin

################################################################################
# C++ Configuration related to Java (compilation of RNAstructure_GUI)
################################################################################
# CXXFLAGS: C++ compiler flags for compiling and linking the native library 
#           (RNAstructure_GUI) that interfaces with Java.
CXXFLAGS+= # Add JAVA-specific CXX flags here.

# Uncomment below to enable multithreading (SMP: OpenMP and pThreads) in the 
# RNAstructure_GUI library, but only if using Linux or Windows (Not Mac which 
# doesn't support OpenMP)
# TODO: This is not possible yet. Changes to the src/rank*.cpp files are necessary 
#       before SMP can be used.
# ifneq ($(OPSYSTEM),Mac)
#   CXXFLAGS+=$(CXXFLAGS_OPENMP) -pthread
# endif

# JNI-specific flags: 
# ('JNI' stands for Java Native Interface. It is the 'glue' that connects
# the Java Runtime with 'native' libraries, such as the RNAstructure class library.)
# JNI_CXXFLAGS - Additional flags that must be passed to the C++ compiler when
#    compiling the SWIG-generated JNI *.cxx source files.
JNI_CXXFLAGS = $(JNI_INCLUDE_FLAGS) -O ### (This is usually different for Mac, so it is redefined below.)

# COMPILE_JNI  - The default recipe to compile the c++ files (.cpp and .cxx) 
#    required for the Java native library.
COMPILE_JNI  = ${COMPILE_CPP} ${JNI_CXXFLAGS}

# LINK_JNI_LIB - rule for linking object files to make the shared library that interfaces with Java. 
LINK_JNI_LIB = ${LINK_DEPS} #LINK_DEPS is defined in compiler.h

# JDK_INCLUDE - The path to the Java JDK headers (including jni.h and other headers for JNI). 
#               See the section at the end of this file called "About JDK_INCLUDE" for details.
JDK_INCLUDE ?= $(JDK_HOME)/include

# JNI_INCLUDE_FLAGS - Flags listing the folders to include for JNI compilation
#    (e.g. -I<directory>).  Specifying them separately from JNI_CXXFLAGS allows
#    better portability between operating systems.
#    See the section called "About JNI_INCLUDE_FLAGS" for details.
JNI_INCLUDE_FLAGS=-I"${JDK_INCLUDE}" -I"${JDK_INCLUDE}"/$(JNI_OS_NAME)

################################################################################
#  Operating-System-Specific Configuration
################################################################################
# OPSYSTEM is defined in ../compiler.h
ifeq (${OPSYSTEM},Linux)
  ############# LINUX #############
  # LINK_JNI_LIB is fine as defined above.
  JNI_OS_NAME=linux
else ifeq (${OPSYSTEM},Mac)
  ############# MAC ###############
  LINK_JNI_LIB += -framework JavaVM  # Required object files will be appended here.
  JNI_OS_NAME=darwin
else ifeq (${OPSYSTEM},Windows) 
  ############# WINDOWS ###########
  LINK_JNI_LIB += -Wl,--add-stdcall-alias # Required object files will be appended here.
  JNI_OS_NAME=win32
endif

################################################################################
#  Verification and Auto-detection of Java installation
################################################################################
# Define some macros for auto-detecting and verifying and Java.
# These are not run on each make invocation. Instead *each* goal defined in 
# Makefile should include $(fn_VerifyJDK) as the first instruction. 

# This will invoke auto-detection (if necessary) and verification
# of javac and jni.h
# Usage: $(call fn_FindJavaDir,<TARGET>,<GUESS>) where TARGET is one of 
# (home, jni, bin, or jre) and GUESS is the best guess for the location.
fn_FindJavaDir=$(shell bash scripts/build/find-java.sh $(FIND_JAVA_FLAGS) '$1' '$2')
FIND_JAVA_FLAGS=-v # -v is verbose (useful messages for end-users) 
                   # -d is debug   (not useful for users, but good for debugging.)

# Usually we can just locate the JDK_HOME dir and define JDK_BIN and JDK_INCLUDE
# implicitly as subdirectores of JDK_HOME.
fn_FindJDK_Home=$(call fn_FindJavaDir,home,$(JDK_HOME))
# Sometimes these folders are in separate locations. In that case, we must verify 
# them individually.
fn_FindJDK_Bin=$(call fn_FindJavaDir,bin,$(JDK_BIN))
fn_FindJDK_Inc=$(call fn_FindJavaDir,jni,$(JDK_INCLUDE))

# The function 'fn_VerifyJDK_Home' is appropriate when the Java installed folder 
# contains both 'bin' and 'include' directories. (As opposed to having them in 
# separate locations). The verify method should expand to a non-empty value on 
# success and empty on failure. Here we override JDK_HOME to ensure it receives 
# the value from the find-java script.
# Then we return $(JDK_HOME) which will be empty on failure.
fn_VerifyJDK_Home=$(eval override JDK_HOME:=$(fn_FindJDK_Home))$(JDK_HOME)


# The function 'fn_VerifyJDK_Dirs' is appropriate when the JDK bin' and 'include' folders
# are not contained in a common parent directory.
# This function overrides JDK_BIN and JDK_INCLUDE to ensure they receive the 
# value from the find-java.sh script.
# Returns $(and $(JDK_BIN),$(JDK_INCLUDE)) which will be empty if EITHER 
# variable could not be auto-detected by the script.
fn_VerifyJDK_Dirs=$(eval override JDK_BIN:=$(fn_FindJDK_Bin))$(eval override JDK_INCLUDE:=$(fn_FindJDK_Inc))$(and $(JDK_BIN),$(JDK_INCLUDE))

# These functions are run depending on the outcome of fn_VerifyJDK (regardless of which verification method is used -- fn_VerifyJDK_Dirs or fn_VerifyJDK_Home).
fn_OnVerifySuccess=$(eval export PATH:=$(JDK_BIN):$(PATH))# Other variables are (JDK_HOME etc) are always exported (below)
fn_OnVerifyFailed=$(error $(MSG_JDK_MISSING))

# tmpBinUsesHome:=$(findstring JDK_HOME,$(value $(JDK_BIN))) # will be non-empty only if the definition of $(JDK_BIN) includes the variable JDK_HOME
# tmpIncUsesHome:=$(findstring JDK_HOME,$(value $(JDK_INCLUDE))) # will be non-empty only if the definition of $(JDK_INCLUDE) includes the variable JDK_HOME

# fn_UseJDKHome -- determines which verification method is appropriate. 
# (*) fn_VerifyJDK_Home should be used if JDK_BIN and JDK_INCLUDE both depend 
#     on JDK_HOME (which they do by default).
# (*) fn_VerifyJDK_Dirs should be used otherwise, because if either JDK_BIN 
#     or JDK_INCLUDE do not contain JDK_HOME, they were likely customized by 
#     the user and therefore might not be in the same parent directory. So each
#     directory should be located individually.
fn_ChooseJDKVerifyBy=$(if $(and $(findstring JDK_HOME,$(value JDK_BIN)),$(findstring JDK_HOME,$(value JDK_INCLUDE))),fn_VerifyJDK_Home,fn_VerifyJDK_Dirs)

# fn_VerifyJDK -- Call this at the top of any recipe that requires Java.
# This function performs Java detection and verification and uses the variable
# `JDK_VERIFIED` to cache the result. If JDK_VERIFIED is empty, a function is
# called that invokes a shell script that returns the auto-detected, verified 
# paths to JDK_BIN and JDK_INCLUDE. 
#
# This function always returns an empty string so it is safe to call it at the 
# top of each recipe that requires Java. The results are cached so repeated 
# calls are fast (and will not result in additional calls to the shell script.
#
# The variable DISABLE_JAVA_AUTODETECT can be set to non-empty to completely 
# disable Java auto-detection and verification.
fn_VerifyJDK=$(if $(JDK_VERIFIED)$(DISABLE_JAVA_AUTODETECT),,$(if $($(fn_ChooseJDKVerifyBy)),$(eval JDK_VERIFIED=YES)$(fn_OnVerifySuccess),$(fn_OnVerifyFailed)))

define MSG_JDK_MISSING

 *** The Java JDK could not be auto-detected.
 *** Please make sure you have the JDK installed (not just the JRE).
 *** Java version 8 (aka 1.8) or later is required.

endef
ifeq (${OPSYSTEM},Linux)
  MSG_JDK_MISSING+=*** ("OpenJDK" on Linux has both a JRE-only version and full JDK version.)$(newline)
endif

#JDK_VERIFIED=#  This must be left blank for auto-detection

# Export variables so they will be accessible to recursive calls to make (preventing additional auto-detection etc)
export JDK_HOME JDK_INCLUDE JDK_BIN PATH JDK_VERIFIED

# Auto-detect and verify the JDK, showing verbose *and* debugging
# information.  This is just for debugging purposes.
# (Do not use find-java in recipes or as a dependency for other goals because it
# is PHONY so it is ALWAYS invoked and forces dependtant goals to be rebuilt).
# Instead, run $(fn_VerifyJDK) in recipes that require Java.
.PHONY: find-java
find-java: FIND_JAVA_FLAGS= -v -d # turn on both verbose and debugging flags.
find-java:
	@$(fn_VerifyJDK)

# Clear the default goal if it was set to one of the targets in this file. 
# (It should be set in the Makefile itself.)
ifeq ($(.DEFAULT_GOAL),$(JAVA_CONFIG_FILE))
   .DEFAULT_GOAL=
endif

define !jdk-show-vars
	printf 'Make vars:\n  JDK_BIN=`$(value JDK_BIN)`=$(JDK_BIN)\n  JDK_INCLUDE=`$(value JDK_INCLUDE)`=$(JDK_INCLUDE)\n  JDK_HOME=`$(value JDK_HOME)`=$(JDK_HOME)\n'
	printf 'Bash vars:\n  JDK_BIN=%s\n  JDK_INCLUDE=%s\n  JDK_HOME=%s\n  PATH=%s\n' "$$JDK_BIN" "$$JDK_INCLUDE" "$$JDK_HOME" "$$PATH"
	javac -version || true
endef

# Test auto-detection and cacheing of results.
ALL_TESTS=test1 test2 test3 test4
test-verify: $(ALL_TESTS)
$(ALL_TESTS): FIND_JAVA_FLAGS= -v -d # turn on both verbose and debugging flags.
$(ALL_TESTS): TEST_HEADER=$(info $n----------$n$!$@ -- $(DESC)$;$n----------)

test1: DESC=Force testing by fn_VerifyJDK_Home
test1:
	$(TEST_HEADER)
	$(eval JDK_BIN=$$(JDK_HOME)/bin)
	$(eval JDK_INCLUDE=$$(JDK_HOME)/include)
	$(eval JDK_VERIFIED=)
	$(info "Verifing By $(fn_ChooseJDKVerifyBy)")
	$(fn_VerifyJDK)
	@$(!jdk-show-vars)

test2: DESC=Force testing by fn_VerifyJDK_Dirs
test2:
	$(TEST_HEADER)
	$(eval JDK_BIN=bananas)
	$(eval JDK_INCLUDE=$$(JDK_HOME)/include)
	$(eval JDK_VERIFIED=)
	$(info "Verifing By $(fn_ChooseJDKVerifyBy)")
	$(fn_VerifyJDK)
	@$(!jdk-show-vars)

test3: DESC=Verification should NOT occur (already cached).
test3:
	$(TEST_HEADER)
	$(info JDK_VERIFIED=$(JDK_VERIFIED))
	$(fn_VerifyJDK)
	@$(!jdk-show-vars)

test4: DESC=Running Sub-Make for test3. Variables should be inherited.
test4:
	$(TEST_HEADER)
	$(info JDK_VERIFIED=$(JDK_VERIFIED))
	$(MAKE) test3

################################################################################
######### Older Method of Java detection without Bash Script ###################
################################################################################
# The makefile OR function expands each argument UNTIL it reaches one that evaluates to NON-empty.
#FOUND_JDK_INCLUDE=$(if $(wildcard $(JDK_INCLUDE)/jni.h),$(call ))
#FOUND_JDK_BIN=$(wildcard $(JDK_BIN)/javac $(JDK_BIN)/javac.exe)

# Get the list of all 'goals' passed to make. e.g. make "GUI" etc. If there are no goals (i.e. just "make"), use the term "default"
#GOALS:=$(or $(MAKECMDGOALS),default)
#SKIP_RESOLVE -- A list of goals that DO NOT require resolution or verification of of paths (like JAVA_HOME, JXX etc)
#SKIP_RESOLVE:=help instructions clean realclean default swig
#SKIP_VERIFY -- A list of goals that DO NOT require verification of paths (like JAVA_HOME, JXX etc), but do require simple resolution
#SKIP_VERIFY:=$(SKIP_RESOLVE) javaconfig 
# Get the list of all goals that DO require path resolution.
#MUST_RESOLVE:=$(filter-out $(SKIP_RESOLVE),$(GOALS))
# Get the list of all goals that DO require verification.
#MUST_VERIFY:=$(filter-out $(SKIP_VERIFY),$(GOALS))
# If the user has specified "javaconfig" as a goal, set the following flag
#DO_CONFIG:=$(filter javaconfig,$(GOALS))

# # DO NOT perform these additional steps if we are only doing 'make clean' or 'make help' etc (hence the test for goals in the 'MUST_RESOLVE' category)
# ifneq ($(strip $(MUST_RESOLVE)),)
#   ifeq ($(FOUND_JDK_INCLUDE),)
#     $(info Detecting JDK 'include' Location... (set JDK_INCLUDE as described in compiler-java.h to avoid this step.))
#     JDK_INCLUDE:=$(call fn_FindJava,jni)
#     $(info JDK_INCLUDE=$(JDK_INCLUDE))
#   endif
#   ifeq ($(FOUND_JDK_BIN),)
#     $(info Detecting JDK 'bin' Location... (set JDK_BIN as described in compiler-java.h to avoid this step.))
#     JDK_BIN:=$(call fn_FindJava,jdk-bin)
#     $(info JDK_BIN=$(JDK_BIN))
#   endif
# endif

# # Run verification of folder and command configuration (e.g. JXX, JDK_INCLUDE).
# JDK_MSG=Please make sure you have the JDK installed and not just the JRE. (Note that "OpenJDK-<VERSION>-jre" is NOT the JDK; It is the JRE.)
# ifeq ($(FOUND_JDK_INCLUDE),)#Still not found!!
#   #Just show a warning because JAVA_HOME is not necessary (as long as JDK_INCLUDE and either JDK_BIN or JXX have been set)
#   $(error The Java JDK 'include' directory could not be located. $(JDK_MSG)$n)
# endif

# ifeq ($(FOUND_JDK_BIN)),)
#   #Just show a warning because JDK_BIN is not necessary (as long as JAR and JXX have been set)
#   $(error The Java JDK 'bin' directory could not be located. $(JDK_MSG)$n)
# endif

  # ifeq ($(strip $(FOUND_JXX)),)
  #   $(error The Java Compiler was not defined or could not be located. $(JAVA_HOME_MSG)$n)
  # endif

  # ifeq ($(strip $(FOUND_JAR)),)
  #   $(error The Java Jar program could not be located.  $(JAVA_HOME_MSG)$n)
  # endif

  # ifeq ($(strip $(FOUND_JDK_INCLUDE)),)
  #   $(error The Java JDK header file 'jni.h' could not be located.$n$t--JDK_INCLUDE: $(JDK_INCLUDE).$n$t--$(JAVA_HOME_MSG)$n)
  # endif

  # ifeq ($(strip $(FOUND_JDK_INCLUDE2)),)
  #   $(error The Java JDK header file 'jni_md.h' could not be located.$n$t--JDK_INCLUDE_OS: $(JDK_INCLUDE)/$(JNI_OS_NAME).$n$t--$(JAVA_HOME_MSG)$n))
  # endif

  #$(info Verification Complete.)
#endif

# Perform the following operations ONLY if the user is running 'make javaconfig'
# ifeq ($(filter javaconfig,$(GOALS)),javaconfig)
#   FAILED_MSG=(Please consider setting the JAVA_HOME or JDK_BIN environment variables manually.)
#   $(info Makefile Java Auto-Config)
#   ifeq ($(FOUND_JAVA_HOME),)
#     $(info Searching for JAVA_HOME...)
#     JAVA_HOME:=$(fn_FIND_JAVA_HOME)
#     $(info JAVA_HOME: $(if $(FOUND_JAVA_HOME),$(JAVA_HOME),Not Found. $(FAILED_MSG)$n))
#   endif

#   ifeq ($(FOUND_JDK_INCLUDE),)
#     JDK_INCLUDE=$(JAVA_HOME)/include
#     ifeq ($(FOUND_JDK_INCLUDE),)
#       $(info Searching for JDK_INCLUDE...)
#       JDK_INCLUDE:=$(fn_FIND_JDK_INCLUDE)
#       $(info JDK_INCLUDE: $(if $(FOUND_JDK_INCLUDE),$(JDK_INCLUDE),Not Found. $(FAILED_MSG)$n))
#     endif
#   endif

#   ifeq ($(FOUND_JDK_BIN),)
#     JDK_BIN=$(JAVA_HOME)/bin
#     ifeq ($(FOUND_JDK_BIN),)
#       $(info Searching for JDK_BIN...)
#       JDK_BIN:=$(fn_FIND_JDK_BIN)
#       $(info JDK_BIN: $(if $(FOUND_JDK_BIN),$(JDK_BIN),Not Found. $(FAILED_MSG)$n))
#     endif
#   endif

#   #$(info JXX==$(FOUND_JXX))
#   ifeq ($(FOUND_JXX),)
#     #Try this alternate first, before showing error
#     JXX=$(JAVA_HOME)/bin/javac
#     #$(info JXX==$(JXX))
#     $(info cmdv=$(shell command -v $(JXX)))
#     #If still not found, show an error
#     $(if $(FOUND_JXX),,$(error The Java Compiler could not be located. $(FAILED_MSG)$n))
#   endif

#   ifeq ($(FOUND_JAR),)
#     #Try this alternate first, before showing error
#     JAR=$(JAVA_HOME)/bin/jar
#     #If still not found, show an error
#     $(if $(FOUND_JAR),,$(error The Java Jar program could not be located. $(FAILED_MSG)$n))
#   endif
# endif



# #fn_REPLACE_WITH_HOME=$(subst $(JAVA_HOME),$$(JAVA_HOME),$1)

# ############# Save Configuration ##################################################
# # Run some tests and save some settings.
# .PHONY: javaconfig
# javaconfig:
# 	$(call fn_FindJava,output "$(JAVA_CONFIG_FILE)")

# 	$(info Writing Java Auto-Config File: $(JAVA_CONFIG_FILE))
# 	$(if $(FOUND_JAVA_HOME),     $(file >  $(JAVA_CONFIG_FILE),JAVA_HOME   = $(JAVA_HOME)))
# 	$(if $(FOUND_JDK_INCLUDE),   $(file >> $(JAVA_CONFIG_FILE),JDK_INCLUDE = $(call fn_REPLACE_WITH_HOME,$(JDK_INCLUDE))))
# 	$(if $(FOUND_JDK_BIN),       $(file >> $(JAVA_CONFIG_FILE),JDK_BIN     = $(call fn_REPLACE_WITH_HOME,$(JDK_BIN))))
# 	$(if $(FOUND_JXX),           $(file >> $(JAVA_CONFIG_FILE),JXX         = $(call fn_REPLACE_WITH_HOME,$(JXX))))
# 	$(if $(FOUND_JAR),           $(file >> $(JAVA_CONFIG_FILE),JAR         = $(call fn_REPLACE_WITH_HOME,$(JAR))))

# 	$(info JAVA_HOME   = "$(JAVA_HOME)")
# 	$(info JDK_INCLUDE = "$(JDK_INCLUDE)")
# 	$(info JDK_BIN     = "$(JDK_BIN)")
# 	$(info JXX         = "$(JXX)")
# 	$(info JAR         = "$(JAR)")

# 	$(info Makefile Java Auto-Config Complete)

##############################################################################
# About JNI_INCLUDE_FLAGS
##############################################################################
#    The JNI_INCLUDE_FLAGS variable contains library-include flags (e.g. -I...) 
# for JDK folders that contain required linking headers for compiling c++ code that 
# interfaces with Java applications via Java Native Interface (JNI).
#
# NOTE: The Java Development Kit (JDK) is required!
#      The Java runtime (JRE) is not sufficient.
#
#    The directories that are included must sometimes be specified manually to 
# coincide with the location where Java native headers are installed. For best results, 
# it is recommended that you base the flags on the JDK_INCLUDE variable.
# This approach is easier to maintain between different machines due to the fact 
# that the location of the JDK can vary, even between machines with the same OS.
#    When looking for these directories on your system, you can search for the files
# "jni.h" and "jni_md.h". On linux and windows, these are usually found
# in two separate folders, but both are subdirectories of the JDK folder.
##############################################################################
##################### EXAMPLE: Microsoft WINDOWS #############################
#  On WINDOWS, the JDK is usually installed in:
#          "C:\Program Files\Java\jdk1.x.x_x"   (1.x.x_x is the java version)
#  If the 32-bit version of java is installed on a 64-bit Windows OS, 
#    it would be in "C:\Program Files (x86)\Java\jdk1.x.x_x"
#  You should set the JDK_INCLUDE variable to the 'include' subfolder under the 
#  JDK folder. I.e. JDK_INCLUDE='C:\Program Files (x86)\Java\jdk1.x.x_x\include'
#  Then the required files "jni.h" and "jni_md.h" can be found in
#      $(JDK_INCLUDE)   and    $(JDK_INCLUDE)\win32    respectively.
#  The JNI_INCLUDE_FLAGS would then be simply:
#      -I"${JDK_INCLUDE}" -I"${JDK_INCLUDE}"/win32 
#     The above "Windows-style" path format ("C:\...") is required for typical 
#  Windows-centric compilers (e.g. Intel, Visual Studio). However when using 
#  a POSIX-aware g++ variabt (e.g. MinGW etc) it may be necessary to specify the 
#  paths in POSIX format: "/c/Program Files/Java/jdk1.x.x_x/include" (or /cygdrive/c/... etc)
#  *** Importantly, quotes are required if the JDK path contains spaces.
#  However, this can be avoided by creating the following **highly-recommended**
#  symlinks to the "Program Files" folder:
#      "C:\Program Files"        ==>  C:\bin64   (or name it Apps64, Programs64, Progs etc.)
#      "C:\Program Files (x86)"  ==>  C:\bin32   
#  This can be done (as an administrator) with the command:  
#       mklink /D C:\bin64 "C:\Program Files"
#  Then JDK_INCLUDE becomes e.g.  C:\bin64\Java\jdk1.x.x_x\include  which does 
#  not require quotes. Alternatively you can refer to these folders using 
#  their DOS 8-character abbreviations: C:\Progra~1 and C:\Progra~2
#  C:\Progra~1 = "C:\Program Files" and C:\Progra~2 = "C:\Program Files (x86)"
##############################################################################
##################### EXAMPLE: Linux #########################################
#  The path to the JDK folder differs between versions of linux, but
# the required subdirectories are usually the same. For example:
#    Fedora Linux:   JDK_INCLUDE=/usr/java/jdk1.x.x_x/include
#    Ubuntu Linux:   JDK_INCLUDE=/usr/lib/jvm/java-x-sun-1.x.x.x/include
# But on both Fedora and Ubuntu, the directories that must be included are:
#    JNI_INCLUDE_FLAGS= -I$(JDK_INCLUDE) -I$(JDK_INCLUDE)/linux
##############################################################################
##################### EXAMPLE: Mac OSX #######################################
#    For modern versions of Java, the include path setup on Mac OSX is 
#  usually similar to:
#  	JDK_INCLUDE=/Library/Java/JavaVirtualMachines/jdk1.x.x_x.jdk/Contents/Home/include
#  	JNI_INCLUDE_FLAGS=-I$(JDK_INCLUDE) -I$(JDK_INCLUDE)/darwin
#  The folder can be found, for example using a command such as:
#     find /Library/Java/JavaVirtualMachines -name jni.h 
#            or perhaps 
#     find /System/Library/Frameworks -name jni.h
#  These commands will print one line for each JDK version found. You only
#  need to use the most recent version, and ONLY the *directory* part
#  (NOT the full path to jni.h).
#  Also, please note that this setup may differ between versions of OSX 
#  and/or versions of Java. For example with Java 1.6 there used to be just 
#  one folder to include:
#     JNI_INCLUDE_FLAGS=-I/System/Library/Frameworks/JavaVM.framework/Headers
##############################################################################
# !!! Note that the example paths given above are only meant to be a general
#  guide to the path forms on different systems. The actual paths on any given
#  system may differ from these !!!
##############################################################################
