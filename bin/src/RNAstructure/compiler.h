SHELL=/bin/bash # Use bash as the shell
ROOTPATH:=$(or $(ROOTPATH),$(patsubst %/,%,$(dir $(lastword $(MAKEFILE_LIST)))))# (See ROOTPATH in Section 1.)
~~PRESERVE_DEFAULT_GOAL:=$(.DEFAULT_GOAL) # Store default goal (to be restored later)
include $(ROOTPATH)/common.mk # general-purpose utility functions and variable definitions. 
################################################################################
#  This file is user-editable and allows configuration of the RNAstructure 
#  build system.
#
#  Contents:
#  Section 1: Description of High-level Configuration Options.
#  Section 2: Standard Compiler Settings, Flags, and Recipes.
#  Section 3: Configuration Customization (Command-Line Options & OS detection).
#  Section 4: Operating System-Specific Configuration.
#  Section 5: Common Linking and Compiling Rules.
#  Section 6: Program-specific Linking and Compiling Recipes.
#  Section 7: CUDA and Combined CPU/CUDA Builds.
#  Section 8: Debugging and Utility Macros/Variables.
#
################################################################################


################################################################################
#  Section 1: Description of High-level Configuration Options.
################################################################################
# The following variables can be set to change the behavior of this makefile.
#
#### Environment variables: ####
#   These are unset by default. But they can be set here or defined in the 
#   environment as desired.
#
#   RNA_MAKE_OS   -  If set, RNA_MAKE_OS will override the default
#                    operating system auto-detection.  Values: Mac, Windows or Linux
#   RNA_MAKE_ARCH -  If set, RNA_MAKE_ARCH will override the target machine 
#                    architecture (e.g. x86 or x86_64), e.g. for cross-compiling.
#
#### Command-line variables: ####
#   These can be passed on the command-line to make in the form "name=value"
#     Example:  ` make Fold debug=on arch=x86 `
#   Note that these must be set on the command-line to have any effect -- simply
#      setting them in the environment or in this file will NOT work.
#
#   static=on|off   - Force static or dynamic linking (default is OS-dependent)
#   dynamic=on|off  - Opposite of static option above.
#   debug=on|off    - Add debug flags (i.e. -g) to CXXFLAGS.
#   optimize=off|on - Turn off optimization flags (i.e. -O3). Default is "on".
#   warnings=on|off - Turn on additional compiler warnings (e.g. -Wall)
#   os=Mac|Linux|Windows - Compile for the specified operating system
#   arch=x86|x86_64 - Compile for the specified machine architecture
#   display=...     - Influences what information is displayed during compilation.
#                     This option affects the VERBOSE flag in compilation and linking rules.
#                     Valid values are described below. Single-character or full word forms are acceptable.
#   display=Q|quiet   - Do not echo recipe commands.  Redirect all command output to files (stdout to $@_make.log and stderr to $@_make-errors.log).
#   display=R|results - Do not echo recipe commands.  But do display all output (stderr and stdout) from the commands.
#   display=V|verbose - Echo recipe commands and show all output from the commands (default Make behavior).
#   display=+T|+target - In addition to the above values for 'display', the text "+target" (or +T) can also be appended, 
#                     which causes the name of the target to be displayed. 
#                     Example:  display=quiet+target  -- Only display the name of the target, but do not show any other output.
#   PYTHON            - The PYTHON variable can be set on the Make command-line to 
#                     set the path to the python executable. (Used when building the python-interface)
#                     For example: `make python-interface PYTHON=/usr/bin/python2.7`

#### Make Variables: ####
#   ROOTPATH  - Usually set in the client Makefile. If not, it is auto-detected at the top of this file.
#               It is the relative path to THIS file's directory FROM the current directory (i.e. the client Makefile's directory)
#               # Note -- we could verify the ROOTPATH passed by the client using $(wildcard ${ROOTPATH}/src/version.h) or $(wildcard ${ROOTPATH}/$(basename $(lastword $(MAKEFILE_LIST)))))
#

################################################################################
#  Section 2: Standard Compiler Settings, Flags, and Recipes.
################################################################################
# These variables are common among several operating systems and compilers, but
# they can also be modified in the OS-specific setup in "Section 4".

# CXX           - Define the c++ compiler. 
#                 Examples: g++ (GNU C++ compiler), icc (Intel C compiler),
#                 CC (SGI and Sun compiler), xlC (xlc compiler, e.g. on AIX)
# CXXFLAGS      - Compiler flags that are common to all platforms and all builds.
#                 (Additional platform-specific or program-specific flags are appended in a later section.)
#                 GCC Options can be found here:  https://gcc.gnu.org/onlinedocs/gcc-4.9.2/gcc/Invoking-GCC.html#Invoking-GCC
#                 Do not directly include optimization flags (-O3) debugging symbols (-g), or warning flags here. 
#                 Those flags are added dynamically, depending on configuration settings and command-line options.
#                 See DEBUG_FLAGS and OPTIMIZE_FLAGS below.
# COMPILE_CPP   - The default recipe to compile a .cpp source file into an object file.
#                 (Do not include $< because COMPILE_CPP is the base for other recipes that include it.)
# VERBOSE       - A macro which influences what information is displayed during target compilation.
#                 See the 'display' option in Section 3 for more information.
  CXX         = g++

  CXXFLAGS    = -DNDEBUG ${ARCHITECTURE} -std=c++11 $(CPP_FLAGS) # C++11 requirement added for log_double.  CPP_FLAGS can be specified by the user (make ... CPP_FLAGS=-DBANANA)
  CXXFLAGS    +=  ${??DEBUG_FLAGS} ${??OPTIMIZE_FLAGS} ${??WARNINGS_FLAGS} #  '??' macros are calculated based on configuration settings described in Section 1.
  COMPILE_CPP = ${VERBOSE}${CXX} -c  -o $@    $(strip ${CXXFLAGS})

#  CXXFLAGS    = -DNDEBUG ${ARCHITECTURE} ${??DEBUG_FLAGS} ${??OPTIMIZE_FLAGS} ${??WARNINGS_FLAGS} # Note: '??' macros are calculated based on configuration settings described in Section 1.
  # _CXX_ is identical to CXX and _CXXFLAGS_ is identical to CXXFLAGS unless RNA_INCLUDE_CUDA is set (see "Combined CPU/CUDA Builds" below)
  COMPILE_CPP = ${VERBOSE}${_CXX_} -c  -o $@    $(strip ${_CXXFLAGS_})

# LINK       - Recipe for linking object files into an executable.
# LDFLAGS    - Additional flags that should be included in the linking command (before object files)
# LDLIBS     - List of additional libraries to pass to the linker. These are listed AFTER all object files (order is important).
# LINK_DEPS  - Full recipe for linking all listed dependencies with the standard recipe and libraries
  LINK      = ${VERBOSE}${_CXX_}  -o $@    $(strip ${_CXXFLAGS_} ${LDFLAGS} ${??STATIC_FLAGS})
  LDLIBS    = -lstdc++ # added to the LINK_DEPS recipe (used for linking standard executables).
  LINK_DEPS = ${LINK}    $^    ${LDLIBS}

# DEBUG_FLAGS - compiler flag(s) to add debugging symbols. These are only added
#   to CXXFLAGS if Make is invoked with the option "debug=on" or 
#   if DEFAULT_DEBUG is ON.
  DEBUG_FLAGS= -g
  DEFAULT_DEBUG=OFF

# OPTIMIZE_FLAGS - compiler flag(s) to enable optimization. These are added 
#   to CXXFLAGS unless make is invoked with the option "optimize=off" 
#   or DEFAULT_OPTIMIZE is OFF.
  OPTIMIZE_FLAGS= -O3
  DEFAULT_OPTIMIZE=ON

# WARNINGS_FLAGS   - Compiler flag(s) to enable additional warnings. 
# DISABLE_WARNINGS - Compiler flag(s) to disable specific warnings. 
#   If make is invoked with the option "warnings=on" WARNINGS_FLAGS are added to CXXFLAGS.
#   otherwise, DISABLE_WARNINGS are added to CXXFLAGS.
  WARNINGS_FLAGS=-Wall
  DISABLE_WARNINGS=-w # More can be added in OS-dependent section.
  DEFAULT_SHOW_WARNINGS=OFF

# STATIC_FLAGS  - Flags for compiling with static libraries. (Only used if DO_STATIC_LINK is "ON").
# DYNAMIC_FLAGS - Flags for compiling with dynamic/shared libraries. (Only used if DO_STATIC_LINK is "OFF").
# DO_STATIC_LINK- Set to 'ON' to enable static linking (by including STATIC_FLAGS in the LINK command) 
#                 or set to 'OFF' to use dynamic linking.  
#                 Calling `make [TARGETS] static=on|off` will override this setting.
  DO_STATIC_LINK=# Off by default (but OS-dependent section turns it on by default for some OSes)
  STATIC_FLAGS  =-static -static-libgcc -static-libstdc++ -Wl,-Bstatic # Seems redundant but necessary for fully static GCC build.
  DYNAMIC_FLAGS =# can be modified in OS-dependent section

# SMP Compiling and Linking:
  COMPILE_SMP = ${COMPILE_CPP} -D COMPILE_SMP #(for SMP without OpenMP -- e.g. Dynalign, Turbofold).
  LINK_SMP    = ${LINK} ${CXXFLAGS_OPENMP}
  LDLIBS_SMP = -lpthread -lgomp # added to the LINK_SMP_DEPS recipe in addition to LDLIBS (used for linking multi-processor executables).
  LINK_SMP_DEPS = ${LINK_SMP}    $^   ${??STATIC_FLAGS} ${LDLIBS} ${LDLIBS_SMP}
  # OpenMP compiling and linking:
  COMPILE_OMP=${COMPILE_CPP} ${CXXFLAGS_OPENMP}
  CXXFLAGS_OPENMP = -fopenmp -D SMP # g++ compiler
  #CXXFLAGS_OPENMP = -openmp -D SMP # intel compiler

  COMPILE_COUNT = ${COMPILE_CPP}  -DCOUNTING

# LIBFLAGS   - compiler flags when producing a shared library (.dll, .so, etc)
  LIBFLAGS = -shared

################################################################################
#  Section 3: Configuration Customization (Command-Line Options & OS detection).
################################################################################
# The config file can be created to contain additional, directory or user-specific 
# variable definitions.
CONFIG_FILE=config.mk
-include $(CONFIG_FILE)

#### Detection of command-line flags passed to make: ####
# Note: fnHandleOption is defined in common.mk. It allows certain flags
# to be turned on or off depending on a default setting which can be overriden
# by a command-line flag of the form 'setting=ON|OFF'

# detect whether the "debug" flag is set to "on"
??DEBUG_FLAGS=$(call fnHandleOption,$(DEFAULT_DEBUG),debug,,$(DEBUG_FLAGS))
# detect whether the "optimize" flag is set to "on" or "off"
??OPTIMIZE_FLAGS=$(call fnHandleOption,$(DEFAULT_OPTIMIZE),optimize,,$(OPTIMIZE_FLAGS))
# detect whether the "warnings" flag is set to "off"
??WARNINGS_FLAGS=$(call fnHandleOption,$(DEFAULT_SHOW_WARNINGS),warnings,,$(WARNINGS_FLAGS),$(DISABLE_WARNINGS))
# detect whether the "static" or "dynamic" flags were set on the command-line to make.
??STATIC_FLAGS=$(call fnHandleOption,$(DO_STATIC_LINK),static,dynamic,$(STATIC_FLAGS),$(DYNAMIC_FLAGS))

#### The 'display' option ####
# Set the default display option (if not specified on the command-line)
DEFAULT_DISPLAY_MODE ?= verbose+target

VERBOSE=$(~DISPLAY_ECHO)$(~DISPLAY_TARGET)$(~DISPLAY_REDIRECT)
  # Get the value of the 'display' option, if set on the command-line.
  # Convert it to lower-case and replace '+'' with ' ' to separate it into words e.g. 'V+T' --> 'V T'
  ~DISPLAY_MODE := $(subst +,$s,$(call fnToLower,$(call fnGetOption,display,$(DEFAULT_DISPLAY_MODE))))
  ~DISPLAY_ECHO=$(if $(filter verbose v,$(~DISPLAY_MODE)),,@)# Start the command with '@' to suppress echo, unless the 'VERBOSE' option is specified
  ~DISPLAY_TARGET=$(if $(filter target t,$(~DISPLAY_MODE)),$(call fnInfo,$(TFMT_BC)$(call fnGetAction,$@) $@$;))# Use $(info) to display the target name if the 'TARGET' option is specified.
  ~DISPLAY_REDIRECT=$(if $(filter quiet q,$(~DISPLAY_MODE)), >'$(call fnMkFileName,$@)_make.log' 2>'$(call fnMkFileName,$@)_make-errors.log' )# Redirect the output if the 'QUIET' option is specified.
  fnMkFileName=$(subst .,_,$(subst /,_,$1))# Convert a path to a flat file name.
  fnGetAction=$(if $(filter %.o,$1),Compiling,$(if $(filter exe/% $(ROOTPATH)/exe%,$1),Linking,Making))# Guess the type of recipe (Compiling vs Linking) from the target name.
  #DEBUG: $(info DISPLAY_MODE=$(~DISPLAY_MODE))$(info DISPLAY_ECHO=$(~DISPLAY_ECHO))$(info DISPLAY_TARGET=$(~DISPLAY_TARGET))$(info VERBOSE=$(VERBOSE))

# OPSYSTEM -- The desired operating system. 
#   Can be set by: (1) Command-line (e.g. os=Windows)   (2) environment RNA_MAKE_OS   (3) auto-detected (in common.mk)
  OPSYSTEM := $(or $(call fnGetOption,os),${RNA_MAKE_OS},$(OS_NAME),$(OPSYSTEM))

# TARGET_ARCH -- Target processor architecture (x86 or x86_64)
#   Can be used to cross-compile for a different architecture if supported by the compiler etc.
#   Can be set by: (1) Command-line (e.g. arch=x86)   (2) environment RNA_MAKE_ARCH   (3) auto-detected (in common.mk)
  TARGET_ARCH := $(or $(call fnGetOption,arch),${RNA_MAKE_ARCH},$(OS_ARCH))

# ARCH_BITS -- The number of bits (32 or 64) of the target architecture.
#   This is used to name some files (see library_defines.h for details.)
  ARCH_BITS := $(if $(findstring 64,$(TARGET_ARCH)),64,32)

# NATIVE_ARCH_BITS -- The number of bits (32 or 64) of the actual current operating system.
#   (May differ from target arch if cross-compiling.)
  NATIVE_ARCH_BITS := $(OS_BITS)# defined in common.mk, based on kernel architecture given by uname

################################################################################
#  Section 4: Operating-System Specific Configuration.
################################################################################
    ifeq (${OPSYSTEM},Linux)
      ############ LINUX ##########################################################
      CXXFLAGS+= -fPIC
      #DYNAMIC_FLAGS=-Wl,--as-needed # "as-needed" forces linker to ignore listed dlls that are not really necessary.
      # Add -m32 if compiling for 32-bit on a 64-bit machine.
      ifneq ($(NATIVE_ARCH_BITS),$(ARCH_BITS))
        ARCHITECTURE = -m$(ARCH_BITS)
      endif
    else ifeq (${OPSYSTEM},Mac)
      ############ MAC ############################################################
      DISABLE_WARNINGS+=-Wno-write-strings -Wno-unused-value -Wno-c++11-compat-deprecated-writable-strings -Wno-logical-op-parentheses -Wno-parentheses-equality
      CXXFLAGS+= -fPIC -arch $(TARGET_ARCH)
      LIBFLAGS = -dynamiclib -arch $(TARGET_ARCH)
      DO_STATIC_LINK=OFF # Do not use static linking on mac. Apple strongly discourages it.
      #DYNAMIC_FLAGS= # Note: '--as-needed' not supported by clang
    else ifeq (${OPSYSTEM},Windows)
      ############ WINDOWS ########################################################
      #  The following settings are for compilation on Windows using the MinGW-w64 
      #  compiler from within the Cygwin POSIX environment. Other choices of shell
      #  or compiler may work as well. 
      #  Compiler executable names may be dependent on the flavor of MinGW 
      #  you have installed. The names below correspond to the "MinGW-w64" build 
      #  obtained directly from the Cygwin repository. Note that 32-bit and 64-bit 
      #  targeted versions of the compiler can be installed side by side.
      WIN64-g++ = x86_64-w64-mingw32-g++
      WIN32-g++ = i686-w64-mingw32-g++
      CXX=$(if $(findstring 64,${TARGET_ARCH}),$(WIN64-g++),$(WIN32-g++))
      CXXFLAGS+= -fsched-spec-load -D__USE_MINGW_ANSI_STDIO=1 ## -fPIC is not supported on windows  ## removed -pthread  for testing
      # On Windows, statically link libgcc, libstdc++, winpthreads, libgomp (OpenMP). 
      # Read about static and shared linking on Windows: https://access.redhat.com/documentation/en-US/Red_Hat_Enterprise_Linux/4/html/Using_ld_the_GNU_Linker/win32.html
      DO_STATIC_LINK=ON # This can be overriden by the `static=ON|OFF` flag to make.
      # INCLUDE_GAMMA: When compiling on Windows with certain compilers (Intel, VS etc), 
      #   we need to include a file to redefine the 'tgamma' function that is missing from cmath. 
      INCLUDE_GAMMA = $(if $(findstring g++,$(CXX)),,${ROOTPATH}/src/gamma.o)
      #DYNAMIC_FLAGS=-Wl,--as-needed # "as-needed" forces linker to ignore listed dlls that are not really necessary.

      #CUDAFLAGS+= -ccbin $$(cygpath -aw $$(which ${CXX}))
    else ifeq (${OPSYSTEM},)
      ############ NO OS ########################################################
      $(error No Operating system defined!!!)
    else
    ############ UNKNOWN OS ###################################################
      $(error Unknown Operating system defined: $(OPSYSTEM))
    endif

################################################################################
#  Section 5: Common Linking and Compiling Rules.
################################################################################
# Default Compilation Rules - generic rules for compiling c++ object files (.o) from
# source files (cpp). These rules will be overridden if a more specific rule 
# is defined (see dependencies.h) 
# Add a rule for cpp files with corresponding header files:
%.o: %.cpp %.h
	$(!MAKE_OUTDIR)
	${COMPILE_CPP} $<

# If no header file is found, this second rule will apply:
%.o: %.cpp
	$(!MAKE_OUTDIR)
	${COMPILE_CPP} $<

################################################################################
#  Section 6: Program-Specific Linking and Compiling Flags.
################################################################################
# This section contains compilation and linking recipes and definitions that
# are used in specific circumstances or with specific programs.
# Usually nothing below this point needs to be changed by end-users, unless certain 
# standard compiler flags need to be modified for a given compiler (e.g. -c or -o )
# It may also be useful to look at the commands below for reference.
  # Define a Dynalign II compiling rule.
  COMPILE_DYNALIGN_II = ${COMPILE_CPP} -D DYNALIGN_II

  # Define a Dynalign II SMP compiling rule.
  COMPILE_DYNALIGN_II_SMP = ${COMPILE_CPP} -D DYNALIGN_II -D COMPILE_SMP

  # Define the MULTIFIND compiling rule.
  COMPILE_MULTIFIND = ${COMPILE_CPP} -D MULTIFIND

  # Define the TURBOHOMOLOGY compiling rule.
  COMPILE_TURBOHOMOLOGY = ${COMPILE_CPP} -D TURBOHOMOLOGY

  # Define CUDA compiling and linking rules (for use with GPU programs).
  COMPILE_CUDA = nvcc ${CUDAFLAGS} -c -o $@
  LINK_CUDA = nvcc ${CUDAFLAGS} -o $@

  # Define an SVM compiling rule.
  COMPILE_SVM = ${COMPILE_CPP} -I${INCLUDESVM}
  # Define an SVM/SMP compiling rule.
  COMPILE_SVM_SMP = ${COMPILE_CPP} -I${INCLUDESVM} -D COMPILE_SMP
  

  # If making Multifind, INCLUDESVM must indicate the location of svm.o from the libsvm package. Also, uncomment the following line, MF = Multifind~.  This will include Multifind in the list of programs to build.
  INCLUDESVM = /usr/local/libsvm-3.24
  #MF = Multifind~

  # Define an INSTRUMENTED compiling rule.
  COMPILE_INSTRUMENTED = ${COMPILE_CPP} -D INSTRUMENTED

######################### GPU/CUDA ############################
#  Section 7: CUDA and Combined CPU/CUDA Builds.
###############################################################
# CUDA_COMPILER - NVIDIA CUDA compiler, for GPU calculations
  CUDA_COMPILER = nvcc
# CUDAFLAGS - flags for the CUDA compiler. These will be merged with CXXFLAGS using -Xcompiler.
  CUDAFLAGS = -D FLOAT -D SHORT -D _FORCE_INLINES -D _CUDA_CALC_ -D NDEBUG \
              ${??DEBUG_FLAGS} ${??OPTIMIZE_FLAGS} \
              -use_fast_math
              # -Wno-deprecated-gpu-targets

# Combined CPU/CUDA Builds: 
# The RNA_INCLUDE_CUDA variable indicates that RNAstructure programs should use
# the CUDA versions of partition and fold. This requires building the "normal"
# classes with the CUDA compiler (nvcc) and CUDA flags as well as compiling the
# CUDA-versions of partition and Fold.
# The RNA_INCLUDE_CUDA variable can be set as a dependency of a Make goal, e.g.:
#          exe/partition: RNA_INCLUDE_CUDA=1
#          exe/partition: pfunction/partition.o ...
#                  ${LINK_DEPS}

# _CXX_ -- This is normally the same as ${CXX} (the cpp compiler), but if
# RNA_INCLUDE_CUDA is set, _CXX_ will be set to ${CUDA_COMPILER}.
# This will cause the normal RNAstructure C++ code to be compiled by nvcc so it
# can be linked with CUDA versions of partition and/or fold.
_CXX_      = $(if ${RNA_INCLUDE_CUDA},${CUDA_COMPILER},$(CXX))
# _CXXFLAGS_ -- This is normally the same as ${CXXFLAGS}, but if
# RNA_INCLUDE_CUDA is set, _CXXFLAGS_ will be set to ${CUDAFLAGS} along with
# multiple '-Xcompiler ...' flags that cause nvcc to pass the real C++ flags 
# through to the compiler it uses internally (which is gcc on Linux or cl.exe
# on Windows).
_CXXFLAGS_ = $(if ${RNA_INCLUDE_CUDA},$(CUDAFLAGS) $(patsubst %,-Xcompiler %,$(CXXFLAGS)),$(CXXFLAGS))

# This preserves the previous method of turning on combined CPU/CUDA builds by
# setting the CUDA environment variable.
ifneq (${CUDA},)
  $(info CUDA environment variable detected: compiling for GPU)
  RNA_INCLUDE_CUDA=1
endif

# For compilation of programs using CUDA, the cuda flags need to be set
  #  CUDAFLAGS = -DFLOAT -DSHORT -O3 -use_fast_math
#    CUDAFLAGS = -DFLOAT -DSHORT -g -use_fast_math

##############################################################################
#  Section 8: Debugging and Utility Macros/Variables.
##############################################################################

# List variables we are interested in when debugging (with `make show-config`)
# TODO: This section is out-of-date. Update with new variable names.
CPP_MAKE_VARS=COMPILE_CPP CXX CXXDEFINES CXXFLAGS LIBFLAGS LINK DO_STATIC_LINK STATIC_FLAGS OPSYSTEM TARGET_ARCH 
APP_MAKE_VARS=COMPILE_CUDA COMPILE_DYNALIGN_II COMPILE_DYNALIGN_II_SMP COMPILE_INSTRUMENTED  COMPILE_MULTIFIND COMPILE_TURBOHOMOLOGY COMPILE_SMP COMPILE_SVM COMPILE_SVM_SMP CUDA_COMPILER CUDAFLAGS CXXFLAGS_OPENMP INCLUDE_GAMMA INCLUDESVM LINK_CUDA LINK_SMP
JAVA_MAKE_VARS=COMPILE_JNI JNI_CXXFLAGS JNI_DIR JNI_INCLUDE_FLAGS JNI_INCLUDE JXX JXXFLAGS LINK_JNI_LIB

.PHONY: show-config
show-config: ; @:
	$(info \
$n---------------------------- Makefile Configuration ----------------------------\
$n  Building on OS: $(OPSYSTEM) (uname: $(SHELL_UNAME))\
$n  Target Architecture: $(TARGET_ARCH) ($(ARCH_BITS)-bits)\
$n$n-----------C++ Configuration------------$(call fnShowVarList,$(CPP_MAKE_VARS))\
$n$n-----------Java Configuration-----------$(call fnShowVarList,$(JAVA_MAKE_VARS))\
$n$n-----------App-specific Configuration-----------$(call fShowVarList,$(APP_MAKE_VARS))\
$n$n       (For a list of ALL variables, use the target 'show-all-vars')\
$n--------------------------------------------------------------------------------)


##############################################################################
#  Final Modifications
##############################################################################
# Define the config file here so make doesn't complain if it is missing.
$(CONFIG_FILE):

# Disable built-in rules and implicit suffix rules. 
# All the rules we need are specified explicitly.
MAKEFLAGS += --no-builtin-rules  --no-builtin-variables
.SUFFIXES:

# Define default message/action for files with no matching rule.
.DEFAULT:
	$(error Unknown target "$@".  For a list of goals, type "make help".)

# Restore the default goal in case it was set to one of the targets in this file. (It should be set in the Makefile itself.)
.DEFAULT_GOAL:=$(~~PRESERVE_DEFAULT_GOAL)
