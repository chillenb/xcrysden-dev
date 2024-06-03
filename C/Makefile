TOPDIR = $(PWD)/..
include ../make.include
include ../Make.sys

#
# don't touch below this
#

INCS      = $(FFTW3_INCDIR) $(MESCHACH_INCDIR) $(TCL_INCDIR) $(TK_INCDIR) $(GL_INCDIR) $(X_INCDIR) $(EXTERNAL_INCDIR)

include make-objects

.SUFFIXES: .o .c
.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCS) -c $<

objects = $(XC_OBJS) $(EXTERNAL_OBJS)

default: xcrys xctclsh xsf2xsf install

compile: 
	$(MAKE) $(CYGWIN) TOPDIR=$(TOPDIR)

#objects: $(objects)
#
#$(objects): %.o: %.c
#	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCS) -c $<

xcrys: $(XC_OBJS) $(EXTERNAL_OBJS)
	$(CC)  $(CFLAGS) $(XC_OBJS) $(EXTERNAL_OBJS) \
	-o xcrys $(FFTW3_LIB) $(MESCHACH_LIB) $(TK_LIB) $(TCL_LIB) $(TOGL_LIB) $(GLU_LIB) $(GL_LIB) $(X_LIB) $(MATH) $(LDLIB) $(EXT_LIB)

xctclsh: xcTclsh.o
	$(CC) $(CFLAGS) xcTclsh.o -o xctclsh $(TCL_LIB) $(MATH) $(LDLIB)

XFS_OBJS = xsf2xsf.o
xsf2xsf: $(XFS_OBJS) 
	$(CC) $(CFLAGS) $(XFS_OBJS) -o xsf2xsf $(MATH) $(LDLIB)


install:
	-cp xcrys   ../bin
	-cp xctclsh ../bin
	-cp xsf2xsf ../bin

# Note: this target is executed via the cproto.sh script !!!
cproto:
	cproto $(INCS) $(CFILE)

clean:
	-rm -f *.o *~ xcrys xctclsh xsf2xsf

header:
	$(TOPDIR)/sys_utils/c_header.sh

depend:
	-if test "no`type makedepend`" != "no" ; then makedepend $(INCS) $(XC_OBJS:.o=.c); fi	

