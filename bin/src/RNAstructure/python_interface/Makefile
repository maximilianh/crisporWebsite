#### The PYTHON interface library ####

ROOTPATH=..# Define the relative path to the RNAstructure root directory.
## Include all macro, dependency, and variable definitions.
include ${ROOTPATH}/compiler.h
include ${ROOTPATH}/library_defines.h
include ${ROOTPATH}/dependencies.h

# The subdirectory *this* Makefile is in, relative to the RNAstructure root.
PY_DIR = python_interface

# The PYTHON variable can be set on the Make command-line to set the path to the python executable.
# For example: `make python-interface PYTHON=/usr/bin/python2.7`
#          or  `make python-interface PYTHON=python3`
PYTHON ?= python

# The output library filename e.g. _RNAstructure_wrap.dll
RNA_WRAP_NAME=_RNAstructure_wrap$(PY_LIB_EXT)
# In python 3 the RNA_WRAP_NAME is longer and includes platform info:
#  e.g.  _RNAstructure_wrap.cpython-36m-x86_64-linux-gnu.so
#    or  _RNAstructure_wrap.cpython-36m-x86_64-cygwin.dll

# The output directory, relative to this Makefile
BIN_DIR=../exe
# Python scripts and other files that need to be in the same directory as the library
SCRIPTS = RNAstructure.py  RNAstructure_wrap.py  Error_handling.py
# A temporary directory for object file output
BUILD_DIR=build/$(OS_NAME)_py$(PY_VERSION)/

RNA_WRAP_PATH=$(BIN_DIR)/$(RNA_WRAP_NAME)

# fnGetPythonConfig -- calls a python script to get python-specific build configuration settings.
fnGetPythonConfig=$(shell $(PYTHON) -c 'import distutils.sysconfig as S; print(S.get_config_var("$1"));')

# Alternatively, you can modify C:\cygwin64\lib\python2.7\_sysconfigdata.py to point BLDLIBRARY to /lib/libpython2.6.dll
ifeq (${OPSYSTEM},Windows)
	fnGetPythonBLDLIBRARY=$(shell $(PYTHON) -c 'import sys, distutils.sysconfig as S; LD = S.get_config_var("$1") if sys.version[0]=="3" else "/lib/"+S.get_config_var("$1"); print(LD)')
else
	fnGetPythonBLDLIBRARY=$(shell $(PYTHON) -c 'import distutils.sysconfig as S; print(S.get_config_var("$1"));')
endif

######### Get some python configuration settings ##########
# Note: `fnCacheVar` allows these variables to be evaluated on-demand
# and then cached so they will not be re-evaluated in the same Make run.
#  PY_INCLUDE_DIR -- directory containing headers needed for source code compilation
#                    of the wrapper library, e.g. "python.h".
$(call fnCacheVar,PY_INCLUDE_DIR,$$(call fnGetPythonConfig,INCLUDEPY))
#  PY_VERSION -- The "short" version of python. Used to determine the name of the python library 
#                to use in the link command.
$(call fnCacheVar,PY_VERSION,$$(call fnGetPythonConfig,VERSION))
#  PY_LIB_DIR -- Directory containing python shared libraries (needed for linking).
$(call fnCacheVar,PY_LIB_DIR,$$(call fnGetPythonConfig,LIBPL))
#  PY_BLDLIBRARY -- Flag required for linking python shared library, eg. "-lpython3.6m"
$(call fnCacheVar,PY_BLDLIBRARY,$$(call fnGetPythonBLDLIBRARY,BLDLIBRARY))
#  PY_LIB_EXT -- Suffix of python extension libraries. (.dll .so, etc)
#                In python3, the extension will be longer, e.g. .cpython-36m-x86_64-linux-gnu.so
$(call fnCacheVar,PY_LIB_EXT,$$(call fnGetPythonConfig,SO))


PYTHON_CXXFLAGS= -I$(PY_INCLUDE_DIR)
LD_DYNAMIC=-Wl,-Bdynamic
FORCE_DYNAMIC=$(call fnIf,$(DO_STATIC_LINK),$(LD_DYNAMIC))# If we are statically linking (C-runtime and stdc++), make sure to DYNAMICALLY link to the python library.
ifeq ($(OPSYSTEM),Mac)
  #DHM-- On mac, I found that -lpython is needed...
  PYTHON_LIBS= -lpthread $(FORCE_DYNAMIC) -L$(PY_LIB_DIR)  $(PY_BLDLIBRARY) -lpython$(PY_VERSION)
else
  PYTHON_LIBS= -lpthread $(FORCE_DYNAMIC) -L$(PY_LIB_DIR)  $(PY_BLDLIBRARY) #-lpython$(PY_VERSION)
endif
COPY_SCRIPTS = $(!CP) $(SCRIPTS) $(BIN_DIR)/
PY_FILEVER=$(subst .,,$(PY_VERSION))# Python version with dots removed (for making a filename)
RNA_WRAP_OBJ=RNAstructure_wrap.$(PY_FILEVER).o

# The PYDIST variable determines whether distutils is used or not.
ifeq ($(PYDIST),)
	# Use the RNAstructure makefile system to build the RNAstructure python interface (default option)
  PY_BUILD_TYPE=makefile
else
	# Use python's internal distutils system to build the RNAstructure python interface
  PY_BUILD_TYPE=distutils
endif

interface: interface-from-$(PY_BUILD_TYPE) # 'interface-from-distutils' or 'interface-from-makefile'

# This uses the RNAstructure Makefile build system to compile the 
# the Python Interface Native Libarary, _RNAstructure_wrap.{dll,so,dylib}

ifeq ($(OPSYSTEM),Mac)
  #On mac $(RNA_WRAP_PATH) is flagged as an error, expand the name explicitly.
  #DHM -- I have no idea why..
  interface-from-makefile: ../exe/_RNAstructure_wrap.so
  ../exe/_RNAstructure_wrap.so: LDLIBS+=${PYTHON_LIBS}
  ../exe/_RNAstructure_wrap.so: LDFLAGS+=${LIBFLAGS} ${PYTHON_CXXFLAGS}
  ../exe/_RNAstructure_wrap.so: $(PYTHON_LIBRARY_FILES) $(RNA_WRAP_OBJ)
	${LINK_DEPS}
	$(COPY_SCRIPTS)
else
  interface-from-makefile: $(RNA_WRAP_PATH)
  $(RNA_WRAP_PATH): LDLIBS+=${PYTHON_LIBS}
  $(RNA_WRAP_PATH): LDFLAGS+=${LIBFLAGS} ${PYTHON_CXXFLAGS}
  $(RNA_WRAP_PATH): $(PYTHON_LIBRARY_FILES) $(RNA_WRAP_OBJ)
	${LINK_DEPS}
	$(COPY_SCRIPTS)
endif

#RNAstructure_wrap.o
$(RNA_WRAP_OBJ): CXXFLAGS+=$(PYTHON_CXXFLAGS)
$(RNA_WRAP_OBJ): RNAstructure_wrap.cxx
	${COMPILE_CPP}   $<

# This uses python's internal distutils to build the RNAstructure interface library.
# The name of the output file can sometimes contain platform and architecture 
# information, e.g.:   _RNAstructure_wrap.cpython-36m-x86_64-linux-gnu.so
interface-from-distutils: # build the python interface library
	mkdir -p $(BUILD_DIR)
	$(PYTHON) setup.py build_ext --inplace --build-temp $(PY_DIR)/$(BUILD_DIR)
	@# The wrapper is placed in the RNAstructure directory. Move it to exe/
	mv ../$(RNA_WRAP_NAME) $(BIN_DIR)
	$(COPY_SCRIPTS)

scripts: # copy RNAstructure.py and dependencies to ../exe/
	$(COPY_SCRIPTS)

source-list:
	#Generate list of sources (used by setup.py in `interface-distutils`)
	@echo '#This file is auto-generated by `make $@`. Manual edits may be overwritten.' > rna_sources.h
	@$(MAKE_CD) ../  src+python_interface_sources | grep -E '\.(cpp|cxx)$$' >> rna_sources.h
	@echo $(PY_DIR)/RNAstructure_wrap.cxx >> rna_sources.h

# The lines in the 'debug' target below look commented-out, 
# but they are valid bash commands that will be displayed as comments.
debug: 
	# This displays some of the configuration variables obtained from python.
	# PYTHON=$(PYTHON)  ($(shell which $(PYTHON)))
	# PY_VERSION=$(PY_VERSION)
	# PY_INCLUDE_DIR=$(PY_INCLUDE_DIR)
	# PY_LIB_DIR=$(PY_LIB_DIR)
	# PY_LIB_EXT=$(PY_LIB_EXT)
	# PY_BLDLIBRARY=$(PY_BLDLIBRARY)
	# RNA_WRAP_NAME=$(RNA_WRAP_NAME)

#to generate your own wrapper files using SWIG
#should not normally be needed
swig:
	swig -python -c++ -keyword -w511 RNAstructure.i

all: all-from-$(PY_BUILD_TYPE) #PY_BUILD_TYPE is either 'distutils' or 'makefile'
all-from-distutils: source-list swig interface-from-distutils
all-from-makefile:  swig interface-from-makefile

clean:
	rm -rf ./build/
	rm -f RNAstructure_wrap.*.o

realclean: clean
	cd $(BIN_DIR) && rm -f $(SCRIPTS) # Remove python interface scripts
	@rm -rf $(BIN_DIR)/__pycache__  # Remove python cache
	rm -f $(RNA_WRAP_NAME) _RNAstructure_wrap.*.so _RNAstructure_wrap.*.dll # Remove python interface library binary
