##########
## Set Library Names
## Options for Windows, Mac, and Linux (default)
##########

#### fnLibName ####
# This function generates shared library file names that are appropriate
# for each OS.
# For example, a library named "Drawing" should have the following file
# names on the respective OS:
#  Windows: Drawing.dll
#  Linux:   libDrawing.so
#  Mac:     libDrawing.dylib
#
# Usage: $(call fnLibName,NAME_WITHOUT_EXT)
fnLibName=$(subst ?,$1,$(LIB_FILE_${OPSYSTEM}))

# Define the library file name format for each OPSYSTEM.
# (this is used by fnLibName)
LIB_FILE_Windows = ?.dll
LIB_FILE_Linux   = lib?.so
LIB_FILE_Mac     = lib?.dylib

# Set shared library names for Windows OS
DRAWING_LIBRARY       = $(call fnLibName,Drawing)
DRAWING_PLOTS_LIBRARY = $(call fnLibName,DrawingPlots)
DRAWING_STRUCTURES_LIBRARY = $(call fnLibName,DrawingStructures)
DYNALIGN_LIBRARY      = $(call fnLibName,Dynalign)
DYNALIGN_SMP_LIBRARY  = $(call fnLibName,Dynalign_SMP)
HYBRID_RNA_LIBRARY    = $(call fnLibName,HybridRNA)
MULTILIGN_LIBRARY     = $(call fnLibName,Multilign)
MULTILIGN_SMP_LIBRARY = $(call fnLibName,Multilign_SMP)
OLIGO_LIBRARY         = $(call fnLibName,Oligo)
PARTS_LIBRARY         = $(call fnLibName,PARTS)
RNA_LIBRARY           = $(call fnLibName,RNA)
TURBOFOLD_LIBRARY     = $(call fnLibName,TurboFold)
TURBOFOLD_SMP_LIBRARY = $(call fnLibName,TurboFold_SMP)

ifeq ($(OPSYSTEM),Mac)
  RNASTRUCTURE_LIBRARY  = $(call fnLibName,RNAstructure_GUI)
else
  # See note at bottom about the RNAstructure_GUI library name
  RNASTRUCTURE_LIBRARY  = $(call fnLibName,RNAstructure_GUI_${ARCH_BITS})
endif

PYTHON_INTERFACE_LIB = $(call fnLibName,_RNAstructure_wrap)

# Run `make DEBUG=1` to show the calculated library names for debugging.
ifneq ($(DEBUG),)
  $(info -----------------$(C_NEWLINE)Library Names:)
  $(foreach var,DRAWING_LIBRARY DRAWING_PLOTS_LIBRARY DRAWING_STRUCTURES_LIBRARY DYNALIGN_LIBRARY DYNALIGN_SMP_LIBRARY HYBRID_RNA_LIBRARY MULTILIGN_LIBRARY MULTILIGN_SMP_LIBRARY OLIGO_LIBRARY PARTS_LIBRARY RNA_LIBRARY RNASTRUCTURE_LIBRARY TURBOFOLD_LIBRARY TURBOFOLD_SMP_LIBRARY,\
  	 $(info $(C_TAB)$(var)$(C_TAB)= $($(var))))
  $(info -----------------$(C_NEWLINE))
endif

############# Note about the RNAstructure_GUI library name ####################
# 32- and 64-bit Java can be installed side-by side on a 64-bit machine, so the
# architecture of the Java VM is not necessarily the same as that of the 
# operating system. So two versions of the native library need to be produced--
# one for each potential VM architecture. The ARCH_BITS variable is set in 
# compiler.h and is either "32" or "64" depending on the target architecture.
# (See also TARGET_ARCH in compiler.h)
###############################################################################