# define the C compiler to use
CC = gcc

# define any compile-time flags
CFLAGS = -Wall -g

# define any directories containing header files other than /usr/include
#
INCLUDES = -I./include

# define the C source files
APIS = ./src/performance.c ./src/words.c 
MAIN1 = ./src/SSC.c
MAIN2 = ./src/Fasta2Spacer.c 

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
API_OBJS = $(APIS:.c=.o)
MAIN1_OBJS = $(MAIN1:.c=.o)
MAIN2_OBJS = $(MAIN2:.c=.o)

# define the executable file 
MAIN1_APP = ./bin/SSC
MAIN2_APP = ./bin/Fasta2Spacer

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

all:    $(MAIN1_APP) $(MAIN2_APP) 

$(MAIN1_APP): $(API_OBJS) $(MAIN1_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN1_APP) $(API_OBJS) $(MAIN1_OBJS) -lm 


$(MAIN2_APP): $(API_OBJS) $(MAIN2_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN2_APP) $(API_OBJS) $(MAIN2_OBJS) -lm 

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) $(API_OBJS) $(MAIN1_OBJS) $(MAIN2_OBJS)

depend: $(SRCS)
	makedepend $(INCLUDES) $^
