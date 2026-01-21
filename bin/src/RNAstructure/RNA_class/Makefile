##########
## Describe which targets the Makefile can use.
##########

instructions:
	@echo 'All libraries are placed in the exe directory after building.'
	@echo 'Use "make Dynalign_Library" to make the Dynalign library.'
	@echo 'Use "make DynalignSMP_Library" to make the Dynalign SMP library.'
	@echo 'Use "make HybridRNA_Library" to make the HybridRNA library.'
	@echo 'Use "make Multilign_Library" to make the Multilign Library.'
	@echo 'Use "make MultilignSMP_Library" to make the Mulitlign SMP Library.'
	@echo 'Use "make Oligo_Library" to make the Oligo library.'
	@echo 'Use "make RNA_Library" to make the RNA library.'
	@echo 'All class example executables remain in the RNA_class directory after building.'
	@echo 'Use "make Dynalign_class" to make "Dynalign_class."'
	@echo 'Use "make DynalignSMP_class" to make "DynalignSMP_class."'
	@echo 'Use "make HybridRNA_class" to make "HybridRNA_class."'
	@echo 'Use "make OligoWalk_class" to make "OligoWalk_class."'
	@echo 'Use "make RNA_class" to make "RNA_class."'

##########
## Define the relative path to the RNAstructure root directory.
## Include all macro, dependency, and variable definitions.
##########

ROOTPATH=..
include ${ROOTPATH}/compiler.h
include ${ROOTPATH}/library_defines.h
include ${ROOTPATH}/dependencies.h

##########
## Define targets.
##########

# Build the Dynalign shared library for serial calculations.
Dynalign_Library: ${OUTDIR}/${DYNALIGN_LIBRARY}
${OUTDIR}/${DYNALIGN_LIBRARY}: ${DYNALIGN_SERIAL_FILES}
	ar rc ${OUTDIR}/${DYNALIGN_LIBRARY} ${DYNALIGN_SERIAL_FILES}

# Build the Dynalign shared library for SMP calculations.
DynalignSMP_Library: ../exe/${DYNALIGN_SMP_LIBRARY}
../exe/${DYNALIGN_SMP_LIBRARY}: ${DYNALIGN_SMP_FILES}
	ar rc ${ROOTPATH}/exe/${DYNALIGN_SMP_LIBRARY} ${DYNALIGN_SMP_FILES}

# Build the HybridRNA class shared library.
HybridRNA_Library: ${OUTDIR}/${HYBRID_RNA_LIBRARY}
${OUTDIR}/${HYBRID_RNA_LIBRARY}: ${HYBRID_FILES}
	ar rc ${OUTDIR}/${HYBRID_RNA_LIBRARY} ${HYBRID_FILES}

# Build the Multilign library for serial calculations.
Multilign_Library: ${OUTDIR}/${MULTILIGN_LIBRARY}
${OUTDIR}/${MULTILIGN_LIBRARY}: Multilign_object.o ${DYNALIGN_SERIAL_FILES}
	ar rc ${OUTDIR}/${MULTILIGN_LIBRARY} Multilign_object.o ${DYNALIGN_SERIAL_FILES}

# Build the Multilign library for SMP calculations.
MultilignSMP_Library: ../exe/$(MULTILIGN_SMP_LIBRARY)
../exe/${MULTILIGN_SMP_LIBRARY}: Multilign_object.o ${DYNALIGN_SMP_FILES}
	ar rc ${ROOTPATH}/exe/${MULTILIGN_SMP_LIBRARY} Multilign_object.o ${DYNALIGN_SMP_FILES}

# Build the Oligo shared library.
Oligo_Library: ${OUTDIR}/${OLIGO_LIBRARY}
${OUTDIR}/${OLIGO_LIBRARY}: $(OLIGO_FILES)
	ar rc ${OUTDIR}/${OLIGO_LIBRARY} ${OLIGO_FILES}

# Build the RNA class shared library.
RNA_Library: ${OUTDIR}/${RNA_LIBRARY}
${OUTDIR}/${RNA_LIBRARY}: ${RNA_FILES}
	ar rc ${OUTDIR}/${RNA_LIBRARY} ${RNA_FILES}

# Build the Dynalign class example program for Dynalign_object.
Dynalign_class: ${DYNALIGN_SERIAL_FILES} Dynalign_class.o
	${LINK} ${DYNALIGN_SERIAL_FILES} Dynalign_class.o

# Build the Dynalign class example SMP program for Dynalign_object.
DynalignSMP_class: ${DYNALIGN_SMP_FILES} DynalignSMP_class.o
	${LINK} ${DYNALIGN_SMP_FILES} DynalignSMP_class.o -lpthread

# Build the HybridRNA class example program for HybridRNA.
HybridRNA_class: ${HYBRID_FILES} HybridRNA_class.o
	${LINK} ${HYBRID_FILES} HybridRNA_class.o

# Build the OligoWalk_class example program for Oligowalk_object.
OligoWalk_class: ${OLIGO_FILES} OligoWalk_class.o
	${LINK} ${OLIGO_FILES} OligoWalk_class.o

# Build the RNA class example program for RNA.
RNA_class: ${RNA_FILES} RNA_class.o 
	$(LINK) ${RNA_FILES} RNA_class.o

##########
## Cleanup.
## Object cleanup removes all temporary build objects.
## Executable cleanup removes all possible executables.
##########

clean:
	rm -f *.o
	rm -f ${DYNALIGN_SERIAL_FILES}
	rm -f ${DYNALIGN_SMP_FILES}
	rm -f ${HYBRID_FILES}
	rm -f ${OLIGO_FILES}
	rm -f ${RNA_FILES}

realclean: clean
	rm -f ${ROOTPATH}/exe/${DYNALIGN_LIBRARY}
	rm -f ${ROOTPATH}/exe/$(DYNALIGN_SMP_LIBRARY)
	rm -f ${ROOTPATH}/exe/${HYBRID_RNA_LIBRARY}
	rm -f ${ROOTPATH}/exe/${MULTILIGN_LIBRARY}
	rm -f ${ROOTPATH}/exe/${MULTILIGN_SMP_LIBRARY}
	rm -f ${ROOTPATH}/exe/${OLIGO_LIBRARY}
	rm -f ${ROOTPATH}/exe/${RNA_LIBRARY}
	rm -f Dynalign_class
	rm -f DynalignSMP_class
	rm -f HybridRNA_class
	rm -f OligoWalk_class
	rm -f RNA_class

