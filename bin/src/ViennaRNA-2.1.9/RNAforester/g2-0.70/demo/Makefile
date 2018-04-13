#
#
# Makefile for g2_library (demo)
#
#


CC = gcc
CXX = g++
##Linux
F77 = g77
CPP = /lib/cpp

##OSF
#F77 = f77
#F77_FLAGS=-i8

##AIX
#F77 = xlf
#CPP = /usr/ccs/lib/cpp

##IRIX
#F77 = f77 -n32
#CPP = /lib/cpp -P -DFORTRAN

CFLAGS = -I../src -I../src/X11 -I../src/PS -I../src/Win32 -I../src/GD -I../src/FIG -g -O2  -I/usr/local/include -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DLINUX=1 -DX_DISPLAY_MISSING=1 -DDO_PS=1 -DDO_FIG=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1
CXXFLAGS = -I../src -I../src/X11 -I../src/PS -I../src/Win32 -I../src/GD/ -I../src/FIG -g -O2 -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DLINUX=1 -DX_DISPLAY_MISSING=1 -DDO_PS=1 -DDO_FIG=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1
LDFLAGS = -L.. -lg2   -L/usr/local/lib -lm

DEMO_PS = simple_PS
DEMO_FIG = simple_FIG
#DEMO_X11 = simple_X11 pointer handles
#DEMO_GD = simple_gd
#DEMO_WIN32 = simple_win32


DEMO = g2_test g2_arc g2_splines_demo demo_cpp $(DEMO_PS) $(DEMO_FIG) $(DEMO_X11) $(DEMO_GD) $(DEMO_WIN32)

all:	$(DEMO)

demo:	$(DEMO)


g2_test:  g2_test.c
	$(CC) $(CFLAGS) g2_test.c $(LDFLAGS) -o $@

g2_arc: g2_arc.c
	$(CC) $(CFLAGS) g2_arc.c $(LDFLAGS) -o $@

g2_splines_demo:  g2_splines_demo.c
	$(CC) $(CFLAGS) g2_splines_demo.c $(LDFLAGS) -o $@

simple_PS:  simple_PS.c
	$(CC) $(CFLAGS) simple_PS.c $(LDFLAGS) -o $@

simple_FIG: simple_FIG.c
	$(CC) $(CFLAGS) simple_FIG.c $(LDFLAGS) -o $@

simple_X11:  simple_X11.c
	$(CC) $(CFLAGS) simple_X11.c $(LDFLAGS) -o $@

simple_win32:  simple_win32.c
	$(CC) $(CFLAGS) simple_win32.c $(LDFLAGS) -lkernel32 -lgdi32 -luser32 -lcomdlg32 -o $@

pointer:  pointer.c
	$(CC) $(CFLAGS) pointer.c $(LDFLAGS) -o $@

handles: handles.c
	$(CC) $(CFLAGS) handles.c $(LDFLAGS) -o $@

simple_gd:  simple_gd.c
	$(CC) $(CFLAGS) simple_gd.c $(LDFLAGS) -o $@

demo_cpp: demo_cpp.cpp
	$(CXX) $(CXXFLAGS) demo_cpp.cpp $(LDFLAGS) -o $@	

demo_f: demo_f.f
	$(F77) $(F77_FLAGS) demo_f.f $(LDFLAGS) -o $@

g2testf: g2_test.f
	$(CPP) -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DLINUX=1 -DX_DISPLAY_MISSING=1 -DDO_PS=1 -DDO_FIG=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1 -E g2_test.f > g2testf.f
	$(F77) g2testf.f $(LDFLAGS) -o g2testf

clean:
	-rm $(DEMO) g2testf g2testf.f demo_f *~ core *.ps *.fig *.eps *.png *.jpg *.fig.bak


