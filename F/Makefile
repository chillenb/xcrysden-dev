TOPDIR = $(PWD)/..
include ../Make.sys

.SUFFIXES : .o .f .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

OBJ =   gengeom.o \
	multhexa.o \
	multatom.o \
	box.o \
	ttcplane.o \
	readf1.o \
	auxil.o \
	wigner.o \
	wigner12d.o \
	reduce.o \
	convexhull.o \
	p2c.o

OBJ2 =  calplane.o 

OBJ3 =  atomlab.o \
	readXYZ.o \
	readf1.o \
	auxil.o

OBJ5 =  str2xcr.o \
	auxil.o

OBJ6 =  multislab.o \
	readf1.o \
	auxil.o \
	writeftn34.o

OBJ7 =  savestruct.o \
	readf1.o \
	auxil.o 

OBJ8 =  kPath.o \
	getintcoor.o \
	auxil.o \
	readf1.o 

OBJ88=  pwKPath.o \
	auxil.o \
	readf1.o 

OBJ9 = wn_readbands.o

OBJ10 =	SRC_spaghetti/get_k.o SRC_spaghetti/get_ei.o

OBJ11 = pwi2xsf_old.o
OBJ111= pwi2xsf.o latgen.o auxil.o pw_util.o

OBJ12 = fsReadBXSF.o auxil.o wigner.o convexhull.o reduce.o

OBJ13 = pwo_xsf2xsf.o latgen.o auxil.o pw_util.o

all: compile install

compile: spaghetti gengeom calplane atomlab recvec str2xcr \
	fracCoor multislab ftnunit savestruct kPath pwKPath \
	fhi_inpini2ftn34 fhi_coord2xcr \
	wn_readbands wn_readbakgen pw fsReadBXSF \
	cube2xsf xsf2cube

spaghetti:
	cd SRC_spaghetti; $(MAKE)

gengeom: $(OBJ)
	$(FC) $(FFLAGS) -o gengeom $(OBJ)

calplane: $(OBJ2)
	$(FC) $(FFLAGS) -o calplane $(OBJ2)

atomlab: $(OBJ3)
	$(FC) $(FFLAGS) -o atomlab $(OBJ3)

recvec: recvec.o auxil.o
	$(FC) $(FFLAGS) -o recvec recvec.o auxil.o 

str2xcr: $(OBJ5)
	$(FC) $(FFLAGS) -o str2xcr $(OBJ5)

fracCoor: fracCoor.o auxil.o  readf1.o
	$(FC) $(FFLAGS) -o fracCoor fracCoor.o auxil.o  readf1.o

multislab: $(OBJ6)
	$(FC) $(FFLAGS) -o multislab $(OBJ6)

savestruct: $(OBJ7)
	$(FC) $(FFLAGS) -o savestruct $(OBJ7)

ftnunit: ftnunit.o
	$(FC) $(FFLAGS) -o ftnunit ftnunit.o

kPath: $(OBJ8)
	$(FC) $(FFLAGS) -o kPath $(OBJ8)

pwKPath: $(OBJ88)
	$(FC) $(FFLAGS) -o pwKPath $(OBJ88)

fhi_inpini2ftn34: fhi_inpini2ftn34.o
	$(FC) $(FFLAGS) -o fhi_inpini2ftn34 fhi_inpini2ftn34.o

fhi_coord2xcr: fhi_coord2xcr.o
	$(FC) $(FFLAGS) -o fhi_coord2xcr fhi_coord2xcr.o

wn_readbands: $(OBJ9)
	$(FC) $(FFLAGS) -o wn_readbands $(OBJ9)  $(OBJ10)

wn_readbakgen: wn_readbakgen.o
	$(FC) $(FFLAGS) -o wn_readbakgen wn_readbakgen.o

pw: pwi2xsf pwi2xsf_old pwo_xsf2xsf


pwi2xsf_old: $(OBJ11) 
	$(FC) $(FFLAGS) -o pwi2xsf_old $(OBJ11)

pwi2xsf: $(OBJ111)
	$(FC) $(FFLAGS) -o pwi2xsf $(OBJ111)

pwo_xsf2xsf: $(OBJ13)
	$(FC) $(FFLAGS) -o pwo_xsf2xsf $(OBJ13)

fsReadBXSF: $(OBJ12)
	$(FC) $(FFLAGS) -o fsReadBXSF $(OBJ12)

cube2xsf: cube2xsf.o
	$(FC) $(FFLAGS) -o cube2xsf cube2xsf.o

xsf2cube: xsf2cube.o
	$(FC) $(FFLAGS) -o xsf2cube xsf2cube.o

PROGRAMS = \
	gengeom calplane atomlab recvec str2xcr fracCoor \
	multislab savestruct ftnunit kPath pwKPath \
	fhi_inpini2ftn34 fhi_coord2xcr \
	wn_readbands wn_readbakgen pwi2xsf pwi2xsf_old pwo_xsf2xsf fsReadBXSF \
	cube2xsf xsf2cube

install:
	for prog in $(PROGRAMS); do \
		if test -x $$prog.exe ; then cp $$prog.exe ../bin ; \
		elif test -x $$prog ; then cp $$prog ../bin ; fi; \
	done

clean:
	-rm -f *.o *.mod gdb.com
	for prog in $(PROGRAMS); do \
		if test -f $$prog.exe ; then rm $$prog.exe; \
		elif test -f $$prog ; then rm $$prog; fi; \
	done
	-cd SRC_spaghetti; rm -f *.o
	-cd SRC_nn; rm -f nn *.o

$(OBJ):  mode.inc
