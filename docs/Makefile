# $Id: Makefile,v 1.4 2019/08/07 15:29:12 tone Exp $

# ------------------------------------------------------------------------
# variables
# ------------------------------------------------------------------------

ROBODOC = robodoc  --html --multidoc --index
#ROBODOC = /temp/tone/building_home/tone/bin/robodoc4.0  --html --multidoc --index

EXAMPLES_DIR = ../examples/Scripting
XC_TCL_DIR   = ../Tcl
XC_TCL_FILES = \
	scripting.tcl \
	scriptingFilter.tcl \
	scriptingMakeMovie.tcl \
	scriptingScalarField2D.tcl \
	scriptingScalarField3D.tcl \
	scriptingAtomicLabels.tcl \
	auxil.tcl

XC_TCL_DIRFILES   = $(addprefix $(XC_TCL_DIR)/,$(XC_TCL_FILES))
EXAMPLES_DIRFILES = $(shell ls $(EXAMPLES_DIR)/*.tcl) $(shell ls $(EXAMPLES_DIR)/*.sh)

all: scripting_docs

scripting_docs:
	if test ! -d scripting_all/examples; then mkdir -p scripting_all/examples; fi
	cp $(XC_TCL_DIRFILES) scripting_all/
	cp $(EXAMPLES_DIRFILES)  scripting_all/examples/
	$(ROBODOC) --src ./scripting_all/ --doc ./html_docs	
	rm -rf scripting_all/

clean:
	-rm -rf html_docs/
	-rm *~

veryclean: clean	
distclean: clean
