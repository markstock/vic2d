# makefile for vic-moc program

CC = gcc
LINKER = gcc
FC = gfortran

#MACH = -march=pentium4 -malign-double -mno-ieee-fp
#MACH = -march=athlon-xp -malign-double -mno-ieee-fp
#MACH = -malign-double -mno-ieee-fp
#MACH = -malign-double
#MACH =-m64

#OPT = -O2
#OPT = -O3 -funroll-loops -ffast-math -fomit-frame-pointer
OPT = -O2 -funroll-loops -ffast-math -fomit-frame-pointer -fopenmp
#OPT = -O2 -funroll-loops -ffast-math -fomit-frame-pointer
#OPT = -O3 -funroll-loops -ffast-math -fbounds-check
#OPT =-g -p -ggdb
#OPT =-I/sw/include
#OPT =-g -p -ggdb -I/sw/include
#OPT =-g -p -ggdb -I/usr/local/include -I/sw/include
#OPT =-I/Developer/SDKs/MacOSX10.5.sdk/usr/include -I/Developer/SDKs/MacOSX10.5.sdk/usr/X11/include
#OPT =-I/Developer/SDKs/MacOSX10.5.sdk/usr/X11/include
#OPT =-I/Developer/SDKs/MacOSX10.5.sdk/usr/X11/include -O3 -funroll-loops -ffast-math -fomit-frame-pointer

LFLAGS =-lm -lgfortran -lpng
#LFLAGS =-L/sw/lib -lm -lgfortran -lpng
#LFLAGS =-L/usr/local/lib -lm -lgfortran -L/sw/lib -lpng
#LFLAGS =-L/Developer/SDKs/MacOSX10.5.sdk/usr/lib -L/Developer/SDKs/MacOSX10.5.sdk/usr/X11/lib -lm -lgfortran -lpng
#LFLAGS =-L/Developer/SDKs/MacOSX10.5.sdk/usr/X11/lib -lm -lgfortran -lpng

#all: diffconv
all: vic2d
#all: vic2d vic3d

gr2.o : gr2.f Makefile
	$(FC) $(OPT) $(MACH) -c $<

gr3.o : gr3.f Makefile
	$(FC) $(OPT) $(MACH) -c $<

mud2sp_full.o : mud2sp_full.f Makefile
	$(FC) $(OPT) $(MACH) -c $<

mud3sp.o : mud3sp.f Makefile
	$(FC) $(OPT) $(MACH) -c $<

%.o : %.c vicmoc.h Makefile
	$(CC) $(OPT) $(MACH) -c $<

#libvicmoc.a: libvicmoc.o gr23.o mud2sp_full.o vicmoc.h Makefile
#	ar -rcs libvicmoc.a libvicmoc.o gr23.o mud2sp_full.o

libvicmoc.a: libvicmoc.o gr2.o mud2sp_full.o vicmoc.h Makefile
	ar -rcs libvicmoc.a libvicmoc.o gr2.o mud2sp_full.o

#libvicmoc.a: libvicmoc.o gr2.o gr3.o mud2sp_full.o mud3sp.o vicmoc.h Makefile
#	ar -rcs libvicmoc.a libvicmoc.o gr2.o gr3.o mud2sp_full.o mud3sp.o

diffconv: libvicmoc.a inout.o diffconv.o vicmoc.h Makefile
	$(LINKER) -o $@ $(OPT) $(MACH) $(LFLAGS) inout.o diffconv.o libvicmoc.a
#vic2d: libvicmoc.a inout.o vic2d.o vicmoc.h Makefile
#	$(LINKER) -o $@ $(OPT) $(MACH) $(LFLAGS) inout.o vic2d.o libvicmoc.a
vic2d: libvicmoc.o gr2.o mud2sp_full.o inout.o vic2d.o vicmoc.h Makefile
	$(LINKER) -o $@ $(OPT) $(MACH) $(LFLAGS) inout.o vic2d.o libvicmoc.o gr2.o mud2sp_full.o
vic3d: libvicmoc.a inout.o vic3d.o vicmoc.h Makefile
	$(LINKER) -o $@ $(OPT) $(MACH) $(LFLAGS) inout.o vic3d.o libvicmoc.a

clean :
	rm -f *.o *.a gmon* a.out
