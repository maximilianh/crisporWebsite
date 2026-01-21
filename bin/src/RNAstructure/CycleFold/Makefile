CC=g++
CFLAGS=-std=c++0x -c -O3 -Wall -pedantic -pthread
LIBS=-lstdc++ -lpthread
SOURCES=NCM_parameters.cpp options.cpp arrays.cpp io.cpp sequence.cpp mainloop.cpp main.cpp constraints.cpp maxexpect.cpp ../src/ParseCommandLine.cpp extrinsic.cpp alignment.cpp $(shell find ../src/phmm/ -name "*.cpp" -type f) turbo_calculation.cpp
TEST_SOURCES=$(SOURCES)
TESTS=test.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../exe/CycleFold

all: $(SOURCES) $(EXECUTABLE)

clean:
	rm -rf *.o

realclean:
	rm -rf *.o $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ ${LIBS}

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

test:
	$(MAKE) CFLAGS="$(CFLAGS) -D TEST" LIBS="$(LIBS) -lboost_unit_test_framework"

main.o: test.h $(shell find ./Tests -type f) mainloop.h
#logdouble.o: logdouble.h
arrays.o: arrays.h constants.h
NCM_parameters.o: NCM_parameters.h constants.h
io.o: io.h sequence.h
sequence.o: sequence.h
#logdouble.o: logdouble.h
constraints.o: constraints.h
mainloop.o: mainloop.h NCM_parameters.h extrinsic.h
alignment.o: alignment.h
extrinsic.o: extrinsic.h
