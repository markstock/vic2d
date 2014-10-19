# makefile for vic-moc program

# build options
#DEBUG=1
OPENMP=1

# build targets
all: vic2d
#all: vic3d
#all: vic2d vic3d

# no user-serviceable parts below

UNAME := $(shell uname)
CC=gcc
LINKER=gcc
FC=gfortran
CFLAGS=


ifdef DEBUG
  CFLAGS+=-g -p -ggdb -fbounds-check
else
  CFLAGS+=-O2 -funroll-loops -ffast-math -fomit-frame-pointer
  CFLAGS+=-mtune=native
endif
ifdef OPENMP
  CFLAGS+=-fopenmp
endif
CFLAGS+=-std=c99

ifeq ($(UNAME), Linux)
  LDFLAGS=
endif
ifeq ($(UNAME), Darwin)
  LDFLAGS=-L/Developer/SDKs/MacOSX10.5.sdk/usr/X11/lib
  CFLAGS=-I/Developer/SDKs/MacOSX10.5.sdk/usr/X11/include
endif
LDFLAGS+=-lm -lgfortran -lpng

gr2.o : gr2.f Makefile
	$(FC) $(CFLAGS) $(MACH) -c $<

gr3.o : gr3.f Makefile
	$(FC) $(CFLAGS) $(MACH) -c $<

mud2sp_full.o : mud2sp_full.f Makefile
	$(FC) $(CFLAGS) $(MACH) -c $<

mud3sp.o : mud3sp.f Makefile
	$(FC) $(CFLAGS) $(MACH) -c $<

%.o : %.c vicmoc.h Makefile
	$(CC) $(CFLAGS) $(MACH) -c $<

LIB2D=libvicmoc2d.o utility.o maskops.o gr2.o mud2sp_full.o
LIB3D=libvicmoc3d.o utility.o gr3.o mud3sp.o mud2sp_full.o

libvicmoc2d.a: $(LIB2D) vicmoc.h Makefile
	ar -rcs libvicmoc2d.a $(LIB2D)

libvicmoc3d.a: $(LIB3D) vicmoc.h Makefile
	ar -rcs libvicmoc3d.a $(LIB3D)

vic2d: libvicmoc2d.a inout.o vic2d.o vicmoc.h Makefile
	$(LINKER) -o $@ $(CFLAGS) $(MACH) inout.o vic2d.o libvicmoc2d.a $(LDFLAGS)

vic3d: libvicmoc3d.a inout.o vic3d.o vicmoc.h Makefile
	$(LINKER) -o $@ $(CFLAGS) $(MACH) $(LDFLAGS) inout.o vic3d.o libvicmoc3d.a $(LDFLAGS)

clean :
	rm -f *.o *.a gmon* a.out
