##########
## Describe which targets the Makefile can use.
##########

instructions:
	@echo 'Use "make PARTS" to create executable "PARTS."'
	@echo 'Use "make PARTS_Library" to make the PARTS library.'

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

PARTS: $(ROOTPATH)/exe/PARTS
$(ROOTPATH)/exe/PARTS: src/parts/main.o ${PARTS_FILES} ${RNA_FILES} ${PHMM_FILES}
	${LINK} src/parts/main.o ${PARTS_FILES} ${RNA_FILES} ${PHMM_FILES}

PARTS_Library: ${PARTS_FILES} ${RNA_FILES} ${PHMM_FILES}
	ar rc ${ROOTPATH}/exe/${PARTS_LIBRARY} ${PARTS_FILES} ${RNA_FILES} ${PHMM_FILES}

##########
## Cleanup.
## Object cleanup removes all temporary build objects.
## Executable cleanup removes all possible executables.
##########

clean:
	rm -f src/parts/main.o
	rm -f ${PARTS_FILES}
	rm -f ${RNA_FILES}
	rm -f ${PHMM_FILES}

realclean: clean
	rm -f ${ROOTPATH}/exe/PARTS
	rm -f ${ROOTPATH}/exe/${PARTS_LIBRARY}

