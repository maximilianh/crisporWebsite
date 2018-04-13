bin_PROGRAMS = Kinfold
SUBDIRS = Example

Kinfold_LDFLAGS=$(OPENMP_CFLAGS)

Kinfold_SOURCES = baum.c cache.c globals.c main.c nachbar.c \
		  baum.h cache_util.h globals.h   nachbar.h \
		  cmdline.c cmdline.h
EXTRA_DIST = Laplace/extract_data.pl \
             Laplace/laplace.sh \
             Laplace/README \
             Laplace/seq1 \
             Laplace/seq1_ref.ps \
             Laplace/to_boxplot.R \
             kinfold.ggo

globals.o: cmdline.h

cmdline.h cmdline.c: kinfold.ggo
	gengetopt -i kinfold.ggo

dist_man_MANS = Kinfold.1
