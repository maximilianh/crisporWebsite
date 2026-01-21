CC=g++
CFLAGS=-std=c++98
SRC=../src
RNA_class=../RNA_class

all: main

clean:
	rm *.o
	rm main

main: rna_library.o main.cpp nn_model.o structure.o nn_model_linear.o nn_model_efn2.o energy.o thermodynamics.o
	$(CC) $(CFLAGS) -o main main.cpp *.o

nn_model_efn2.o: nn_model_efn2.hpp nn_model_efn2.cpp nn_model_linear.o
	$(CC) $(CFLAGS) -c -o nn_model_efn2.o nn_model_efn2.cpp

nn_model_linear.o: nn_model_linear.hpp nn_model_linear.cpp nn_model.o
	$(CC) $(CFLAGS) -c -o nn_model_linear.o nn_model_linear.cpp

nn_model.o: nn_model.hpp nn_model.cpp primary_structure.o secondary_structure.o rna_library.o energy.o
	$(CC) $(CFLAGS) -c -o nn_model.o nn_model.cpp

secondary_structure.o: secondary_structure.hpp secondary_structure.cpp primary_structure.o
	$(CC) $(CFLAGS) -c -o secondary_structure.o secondary_structure.cpp

energy.o: energy.hpp energy.cpp
	 $(CC) $(CFLAGS) -c -o energy.o energy.cpp

primary_structure.o: primary_structure.hpp primary_structure.cpp
	$(CC) $(CFLAGS) -c -o primary_structure.o primary_structure.cpp

thermodynamics.o: $(RNA_class)/thermodynamics.cpp $(RNA_class)/thermodynamics.h rna_library.o
	$(CC) $(CFLAGS) -c -o thermodynamics.o $(RNA_class)/thermodynamics.cpp

rna_library.o: $(SRC)/rna_library.cpp $(SRC)/rna_library.h stackclass.o stackstruct.o structure.o
	$(CC) $(CFLAGS) -c -o rna_library.o $(SRC)/rna_library.cpp

stackclass.o: $(SRC)/stackclass.cpp $(SRC)/stackclass.h
	$(CC) $(CFLAGS) -c -o stackclass.o $(SRC)/stackclass.cpp

stackstruct.o: $(SRC)/stackstruct.cpp $(SRC)/stackstruct.h
	$(CC) $(CFLAGS) -c -o stackstruct.o $(SRC)/stackstruct.cpp

structure.o: $(SRC)/structure.cpp $(SRC)/structure.h
	$(CC) $(CFLAGS) -c -o structure.o $(SRC)/structure.cpp
