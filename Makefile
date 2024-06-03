#------------------------------------------------------------------------
#
# XCRYSDEN Top Makefile
#
#------------------------------------------------------------------------

SHELL     = /bin/sh
TOPDIR    = $(CURDIR)
TCL_INDEX = $(TOPDIR)/util/tcl_index

include make.include

what:
	@clear;	$(SHELL) make-usage 


# some varibales ...

all: tcl tk mesa togl fftw xcrysden
#  meschach 
xcrysden: usage bwidget bindir src-C src-F src-Tcl

usage: $(TOPDIR)/docs/xcrysden.1
	man $(TOPDIR)/docs/xcrysden.1 | awk 'BEGIN {lprint=0; print "## do not edit changes will be lost (file automatically generated)\n"; } /SYNOPSIS/ { lprint=1; } /SEE ALSO/ { lprint=0; } /a*/ { if (lprint) print; }' > usage

check-make-sys:
	@if test ! -f Make.sys ; then \
	echo ""; \
	echo "   First copy an appropriate system/Make.* file to ./Make.sys "; \
	echo "   and edit it to suit your needs."; \
	echo ""; exit 1; fi

help:
	@. make-usage


bindir:
	if test ! -d bin ; then mkdir bin; fi

tcl: check-make-sys
	if test -d external/src ; then \
		cd external/src; $(MAKE)  "TOPDIR=$(TOPDIR)" tcl; \
	fi

tk: check-make-sys
	if test -d external/src ; then \
		cd external/src; $(MAKE)  "TOPDIR=$(TOPDIR)" tk; \
	fi

togl: check-make-sys
	if test -d external/src ; then \
		cd external/src; $(MAKE)  "TOPDIR=$(TOPDIR)" togl; \
	fi

mesa: check-make-sys
	if test -d external/src ; then \
		cd external/src; $(MAKE)  "TOPDIR=$(TOPDIR)" mesa; \
	fi

meschach: check-make-sys
	if test -d external/src ; then \
		cd external/src; $(MAKE)  "TOPDIR=$(TOPDIR)" meschach; \
	fi

fftw: check-make-sys
	if test -d external/src ; then \
		cd external/src; $(MAKE)  "TOPDIR=$(TOPDIR)" fftw; \
	fi

bwidget: check-make-sys
	if test -d external/src ; then \
		cd external/src; $(MAKE)  "TOPDIR=$(TOPDIR)" bwidget; \
	fi

src-C: check-make-sys
	@echo
	@echo "#------------------------------#"
	@echo "#                              #"
	@echo "#   Compiling XCRYSDEN C-code  #"
	@echo "#                              #"
	@echo "#------------------------------#"
	@echo
	cd C; $(MAKE) "TOPDIR=$(TOPDIR)" compile


src-F: check-make-sys
	@echo
	@echo "#------------------------------#"
	@echo "#                              #"
	@echo "#   Compiling XCRYSDEN F-code  #"
	@echo "#                              #"
	@echo "#------------------------------#"
	@echo
	cd F/SRC_nn; $(MAKE)  "TOPDIR=$(TOPDIR)"
	cd F/SRC_spaghetti; $(MAKE)  "TOPDIR=$(TOPDIR)"
	cd F; $(MAKE)  "TOPDIR=$(TOPDIR)"

src-Tcl:
	@echo
	@echo "#------------------------------#"
	@echo "#                              #"
	@echo "#  Managing XCRYSDEN Tcl-code  #"
	@echo "#                              #"
	@echo "#------------------------------#"
	@echo
	cd Tcl; $(MAKE)
	cd Tcl/fs; $(MAKE)

tests:
	if test -x bin/xcrys ; then cd tests; ./make_all_tests.sh; fi

# ------------------------------------------------------------------------
#
# clean-targets
#
# ------------------------------------------------------------------------

clean: clean-C clean-F clean-Tcl wrappers_clean	
	cd examples; $(MAKE) clean

veryclean: clean clean-bin clean-docs clean-external 

distclean: veryclean clean-bck
	if test -f usage;    then rm -f usage; fi
	if test -d external; then cd external; $(MAKE) distclean; fi

clean-bck:
	-rm -f *~
	-rm -f */*~
	-rm -f */*/*~
	-rm -f */*/*/*~

clean-C: check-make-sys
	cd C; $(MAKE) "TOPDIR=$(TOPDIR)" clean

clean-F: check-make-sys
	-cd F/SRC_nn; $(MAKE) "TOPDIR=$(TOPDIR)" clean
	-cd F/SRC_spaghetti; $(MAKE) "TOPDIR=$(TOPDIR)" clean
	-cd F; $(MAKE) "TOPDIR=$(TOPDIR)" clean

clean-Tcl:
	cd Tcl/fs; $(MAKE) clean
	cd Tcl; $(MAKE) clean

clean-bin:
	if test -d bin/; then \
		rm -f bin/*; \
	fi

clean-docs:
	if test -d docs/; then \
		cd docs; $(MAKE) "TOPDIR=$(TOPDIR)" veryclean; \
	fi

clean-external:
	if test -d external/src ; then \
	   cd external/src ; $(MAKE) "TOPDIR=$(TOPDIR)" clean; \
	fi

wrappers_clean:
	for file in $(addsuffix .wrapper,$(PROGS)); do \
		if test -f $$file; then rm -f $$file; fi; \
	done 

# ========================================================================
#
# here are targets for making various distributions
#
# ========================================================================

PROGS = xcrysden pwi2xsf pwo2xsf ptable unitconv

README_FILES = \
	AUTHORS \
	COPYING \
	COPYRIGHT \
	ChangeLog \
	NEWS \
	README README.cygwin \
	THANKS \
	otherLICENSES/ 

IRON_ITEMS = usage version xcrysden
IRON_DIRS = \
	Awk/ \
	Tcl/*.tcl Tcl/tclIndex Tcl/fs/*.tcl Tcl/fs/tclIndex \
	Tcl/Xcrysden_resources Tcl/custom-definitions \
	contrib/ \
	examples/ \
	images/ \
	scripts/ \
	tests/ \
	util/ 

MAN_PAGES  = $(addsuffix .1,$(PROGS))
MAN_FILES  = $(addprefix docs/,$(MAN_PAGES))
BAT_FILES  = xcrysden.bat 
IRON_FILES = $(IRON_ITEMS) $(IRON_DIRS) $(MAN_FILES)

SRC_ONLY_FILES = \
	Makefile make-usage make.include \
	C/*.c C/*.h C/Makefile C/make-objects C/*.cygwin \
	F/*.f F/*.f90 F/*.inc F/Makefile \
	F/*/*.f F/*/*.inc F/*/Makefile \
	Tcl/Makefile Tcl/fs/Makefile \
	docs/Makefile \
	sys_utils/*.sh \
	system/ \

EXTERNAL_SRC_FILES    = external/Makefile external/src/Makefile
EXTERNAL_LIB_BWIDGET  = external/lib/bwidget-$(BWIDGET_VER)
EXTERNAL_LIB_TCLTK    = external/lib/tcl$(TCL_VER2) external/lib/tk$(TCL_VER2) external/lib/Togl$(TOGL_VER)
EXTERNAL_LIB          = external/lib

EXTERNAL_SHAREDLIB_FILES = $(wildcard \
	external/lib/libtcl$(TCL_VER2).so* external/lib/libtcl$(TCL_VER2).dylib \
	external/lib/libtk$(TCL_VER2).so* external/lib/libtk$(TCL_VER2).dylib \
	external/lib/libTogl$(TOGL_VER).so* external/lib/libTogl$(TOGL_VER).dylib \
	external/lib/libGL.so* external/lib/libGL.*.dylib* \
	external/lib/libGLU.so* external/lib/libGLU.*.dylib* \
	external/lib/libmeschach.so* external/lib/libmeschach*.dylib* \
	external/lib/libfftw3.so* external/lib/libfftw3*.dylib* \
	external/lib/libquadmath.so* external/lib/libquadmath.dylib*)

EXTERNAL_GFORTRAN_LIB =
#EXTERNAL_GFORTRAN_LIB = external/lib/libgfortran.so.$(GFORTRAN_MINOR_VERSION)

BIN_FILES = bin/

# ------------------------------------------------------------------------
#
# Target for installing system wide: default localtion is /usr/local
# To install to different location run as: prefix=XXXXX make iunstall
#
# ------------------------------------------------------------------------

prefix   ?= /usr/local
version  := $(shell cat version)
xcrysden  = xcrysden-$(version)

install: xcrysden
	@echo
	@echo "#-----"
	@echo "#                          "
	@echo "#   Installing XCRYSDEN to: $(prefix)"
	@echo "#                          "
	@echo "#-----"
	@echo	
	install -m755 -d $(prefix)/share/$(xcrysden)
	cp -a $(IRON_ITEMS) $(prefix)/share/$(xcrysden)
	\
	for subdir in Awk $(EXTERNAL_LIB_BWIDGET) images scripts Tcl util; do \
		if test -d $$subdir; then \
			install -m755 -d $(prefix)/share/$(xcrysden)/$$subdir; \
			cp -a   $$subdir/* $(prefix)/share/$(xcrysden)/$$subdir; \
		fi; \
	done
	\
	if test -d examples; then \
		install -m755 -d $(prefix)/share/doc/$(xcrysden)/examples; \
		cp -a   examples/* $(prefix)/share/doc/$(xcrysden)/examples; \
		ln -sf $(prefix)/share/doc/$(xcrysden)/examples $(prefix)/share/$(xcrysden)/examples; \
	fi; \
	\
	install -m755 -d $(prefix)/share/man/man1
	install -m644 $(MAN_FILES) $(prefix)/share/man/man1/
	gzip -f $(addprefix $(prefix)/share/man/man1/,$(MAN_PAGES))
	\
	install -m755 -d $(prefix)/lib/$(xcrysden)
	install -m755 bin/*  $(prefix)/lib/$(xcrysden)/
	\
	prefix=$(prefix) xcrysden=$(xcrysden) sh sys_utils/wrappers.sh
	if test ! -d $(prefix)/bin; then install -m755 -d $(prefix)/bin; fi
	for prog in $(PROGS); do \
		install -m755 $$prog.wrapper   $(prefix)/bin/$$prog; \
	done



#------------------------------------------------------------------------
#
# SOURCE-distributions
#
#------------------------------------------------------------------------

srcdist:      clean-bck src-Tcl usage _src-dist    

#------------------------------------------------------------------------
#
# BINARY (i.e. precompiled) distributions 
#
#------------------------------------------------------------------------

bindist: bindist-shared; # the default bindist is shared

bindist-static:      clean-bck xcrysden _bin-dist-static
bindist-shared:      clean-bck xcrysden _bin-dist-fully-shared
bindist-semishared:  clean-bck xcrysden _bin-dist-semishared


#
# sorce distribution 
#
_src-dist: 
	tar -cvf xcrysden.tar $(README_FILES) $(IRON_FILES) $(BAT_FILES) $(SRC_ONLY_FILES) $(EXTERNAL_SRC_FILES)
	sys_utils/xcRepackage.sh $(TOPDIR) xcrysden.tar


#
# statically linked binary distribution (with external/lib/tcl/ and external/lib/tk/)
#
_bin-dist-static: 
	tar -cvf xcrysden.tar $(README_FILES) $(IRON_FILES) $(EXTERNAL_LIB_BWIDGET) $(EXTERNAL_LIB_TCLTK) $(BIN_FILES)
	sys_utils/xcRepackage.sh $(TOPDIR) xcrysden.tar bin-static


# shared linked distribution without shared libs
_bin-dist-fully-shared: 
	if test -f $(EXTERNAL_GFORTRAN_LIB); then \
		tar -cvf xcrysden.tar $(README_FILES) $(IRON_FILES) $(EXTERNAL_LIB_BWIDGET) $(BIN_FILES) $(EXTERNAL_GFORTRAN_LIB); \
	else \
		tar -cvf xcrysden.tar $(README_FILES) $(IRON_FILES) $(EXTERNAL_LIB_BWIDGET) $(BIN_FILES); fi
	sys_utils/xcRepackage.sh $(TOPDIR) xcrysden.tar bin-shared


# semi-shared linked distribution with shared libs in external/lib
_bin-dist-semishared:
	-tar -cvf xcrysden.tar \
		$(README_FILES) $(IRON_FILES) $(EXTERNAL_LIB_BWIDGET) $(EXTERNAL_LIB_TCLTK) $(EXTERNAL_SHAREDLIB_FILES) $(BIN_FILES); \
		status=$$?; \
		if test $$status -gt 0; then \
			tar -cvf xcrysden.tar \
			$(README_FILES) $(IRON_FILES) $(EXTERNAL_LIB) $(BIN_FILES); \
		fi
	sys_utils/xcRepackage.sh $(TOPDIR) xcrysden.tar bin-semishared

