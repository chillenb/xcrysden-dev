TOPDIR = $(PWD)/..
TCL_INDEX = $(TOPDIR)/util/tcl_index

include ../Make.sys

.tcl:

tclIndex: *.tcl
	XCRYSDEN_TOPDIR=$(TOPDIR); \
	export XCRYSDEN_TOPDIR; $(TCL_INDEX)

clean:
	-rm -f *~
