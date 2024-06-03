include ../../Make.sys

OBJ =   nn.o \
	dirlat.o \
	ord2.o

all: nn install

nn:     $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o nn 

install:
	if test -f nn.exe ; then cp nn.exe ../../bin ; \
	elif test -f nn     ; then cp nn     ../../bin ; fi

clean: 
	-rm *.o
	if test -f nn.exe; then \
		rm nn.exe; \
	elif test -f nn; then \
		rm nn; \
	fi

