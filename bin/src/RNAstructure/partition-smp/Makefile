include config.h

all: $(BIN)/partition-lowmem $(BIN)/partition-cuda $(BIN)/partition-save-params
$(BIN)/partition-cuda: cuda_base.o cuda_main.o cuda_param.o cuda_prna.o cuda_util.o
	$(CUDA_CC) -o $@ $(CUDA_CFLAGS) $^
$(BIN)/partition-lowmem: base.o main.o param.o prna.o util.o
	$(CC) -o $@ $(CFLAGS) $^ -lm
$(BIN)/partition-save-params: base.o param.o save-params.o util.o
	$(CC) -o $@ $(CFLAGS) $^ -lm
clean:
	rm -f *.o *.cu

# symlinks for cuda source files
cuda_base.cu: base.c
	ln -s $< $@
cuda_main.cu: main.c
	ln -s $< $@
cuda_param.cu: param.c
	ln -s $< $@
cuda_prna.cu: prna.c
	ln -s $< $@
cuda_util.cu: util.c
	ln -s $< $@

# suffix rules
.SUFFIXES: .c .o .cu
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
.cu.o:
	$(CUDA_CC) $(CUDA_CFLAGS) -o $@ -c $<

# dependencies on header files
cuda_base.o base.o: base.c base.h cu.h util.h
cuda_main.o main.o: main.c prna.h param.h base.h cu.h real.h util.h
cuda_param.o param.o: param.c param.h base.h cu.h real.h util.h
cuda_prna.o prna.o: prna.c base.h cu.h prna.h param.h real.h util.h
save-params.o: save-params.c param.h base.h cu.h real.h util.h
cuda_util.o util.o: util.c util.h
