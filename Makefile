# makefile for vic-moc program

# build options
#DEBUG=1
OPENMP=1
#MINGW=1

# no user-serviceable parts below

UNAME := $(shell uname)
CC=gcc
LINKER=gcc
FC=gfortran
CFLAGS=
LDFLAGS=
EXE=vic2d
MACH=-march=native

ifdef DEBUG
  CFLAGS+=-g -p -ggdb -fbounds-check
else
  CFLAGS+=-O2 -funroll-loops -ffast-math -fomit-frame-pointer
  #CFLAGS+=-mtune=native
endif
ifdef OPENMP
  CFLAGS+=-fopenmp
endif
CFLAGS+=-std=c99

ifeq ($(UNAME), Linux)
  LDFLAGS+=
endif
ifeq ($(UNAME), Darwin)
  #LDFLAGS+=-L/opt/X11/lib
  LDFLAGS+=-L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk/usr/X11/lib
  LDFLAGS+=-L/usr/local/gfortran/lib
  #CFLAGS+=-I /opt/X11/include
  CFLAGS+=-I /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk/usr/X11/include
endif
ifdef MINGW
  CC=x86_64-w64-mingw32-gcc
  LINKER=x86_64-w64-mingw32-gcc
  FC=x86_64-w64-mingw32-gfortran
  EXE=vic2d.exe
  LDFLAGS+=-static
endif
LDFLAGS+=-lgfortran -lpng -lm -lz

# build targets
all: $(EXE)
#all: vic3d
#all: vic2d vic3d

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

LIB2D=libvicmoc2d.o utility.o maskops.o gr2.o mud2sp_full.o mud2sp_extern.o
LIB3D=libvicmoc3d.o utility.o gr3.o mud3sp.o mud2sp_full.o

libvicmoc2d.a: $(LIB2D) vicmoc.h Makefile
	ar -rcs libvicmoc2d.a $(LIB2D)

libvicmoc3d.a: $(LIB3D) vicmoc.h Makefile
	ar -rcs libvicmoc3d.a $(LIB3D)

$(EXE): libvicmoc2d.a inout.o vic2d.o vicmoc.h Makefile
	$(LINKER) -o $@ $(CFLAGS) $(MACH) inout.o vic2d.o libvicmoc2d.a $(LDFLAGS)

vic3d: libvicmoc3d.a inout.o vic3d.o vicmoc.h Makefile
	$(LINKER) -o $@ $(CFLAGS) $(MACH) $(LDFLAGS) inout.o vic3d.o libvicmoc3d.a $(LDFLAGS)

clean :
	rm -f *.o *.a gmon* a.out

