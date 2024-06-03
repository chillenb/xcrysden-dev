TOPDIR = $(PWD)/../..
include ../../make.include
include ../../Make.sys

default:
	@echo; echo "*** Possible targets: tcl tk togl mesa fftw clean"; echo

all: tcl tk mesa togl fftw meschach bwidget

tcl: 
	@if test "$(COMPILE_TCLTK)" = yes  -o  "$(COMPILE_ALL)" = yes; then \
		\
		if test ! -f $(TCL_TGZ); then \
			echo ;\
			echo "#------------------------------------#";\
			echo "#                                    #";\
			echo "#   Downloading external Tcl-library #";\
			echo "#                                    #";\
			echo "#------------------------------------#";\
			echo ;\
			$(DOWNLOAD) $(TCL_DOWNLOAD); \
		fi; \
		\
		echo ;\
		echo "#----------------------------------#";\
		echo "#                                  #";\
		echo "#   Compiling external Tcl-library #";\
		echo "#                                  #";\
		echo "#----------------------------------#";\
		echo ;\
		\
		if [ \( ! -d tcl$(TCL_VER3) \)  -a \( -f $(TCL_TGZ) \) ]; then \
			gunzip -c $(TCL_TGZ) | tar xvf -; \
		fi; \
		\
		cd tcl$(TCL_VER3)/unix; \
		./configure $(TCLTK_OPTIONS) --prefix=$(TOPDIR)/external; \
		$(MAKE); $(MAKE) install; \
	fi


tk:
	@if test "$(COMPILE_TCLTK)" = yes  -o  "$(COMPILE_ALL)" = yes; then \
		\
		if test ! -f $(TK_TGZ); then \
			echo ;\
			echo "#-----------------------------------#";\
			echo "#                                   #";\
			echo "#   Downloading external Tk-library #";\
			echo "#                                   #";\
			echo "#-----------------------------------#";\
			echo ;\
			$(DOWNLOAD) $(TK_DOWNLOAD); \
		fi; \
		\
		echo ;\
		echo "#---------------------------------#";\
		echo "#                                 #";\
		echo "#   Compiling external Tk-library #";\
		echo "#                                 #";\
		echo "#---------------------------------#";\
		echo ;\
		\
		if [ \( ! -d tk$(TCL_VER3) \)  -a \( -f $(TK_TGZ) \) ]; then \
			gunzip -c $(TK_TGZ) | tar xvf -; \
		fi; \
		\
		cd tk$(TCL_VER3)/unix; \
		CPPFLAGS=-I$(TOPDIR)/external/include LDFLAGS=-L$(TOPDIR)/external/lib ./configure $(TCLTK_OPTIONS) --prefix=$(TOPDIR)/external; \
		$(MAKE); $(MAKE) install; \
	fi


togl:
	@if test "$(COMPILE_TOGL)" = yes  -o  "$(COMPILE_ALL)" = yes; then \
		\
		if test ! -f $(TOGL_TGZ); then \
			echo ;\
			echo "#-------------------------------------#";\
			echo "#                                     #";\
			echo "#   Downloading external Togl-library #";\
			echo "#                                     #";\
			echo "#-------------------------------------#";\
			echo ;\
			$(DOWNLOAD) $(TOGL_DOWNLOAD); \
		fi; \
		\
		echo ;\
		echo "#-----------------------------------#";\
		echo "#                                   #";\
		echo "#   Compiling external Togl-library #";\
		echo "#                                   #";\
		echo "#-----------------------------------#";\
		echo ;\
		\
		if [ \( ! -d Togl$(TOGL_VER) \)  -a \( -f $(TOGL_TGZ) \) ]; then \
			gunzip -c $(TOGL_TGZ) | tar xvf -; \
		fi; \
		\
		(cd Togl$(TOGL_VER); \
		CPPFLAGS=-I$(TOPDIR)/external/include LDFLAGS=-L$(TOPDIR)/external/lib ./configure $(TOGL_OPTIONS) --prefix=$(TOPDIR)/external; \
		$(MAKE); $(MAKE) install); \
		(cd ../lib; togl=`ls libTogl$(TOGL_VER).*`; if test x$$togl = x; then ln -sf Togl$(TOGL_VER)/libTogl$(TOGL_VER).* .; fi); \
	fi


#mesa: mesa-traditional
mesa: mesa-new glu-new

mesa-new: external_dirs
	@if test "$(COMPILE_MESA)" = yes  -o  "$(COMPILE_ALL)" = yes; then \
		\
		if test ! -f $(MESA_TGZ); then \
			echo ;\
			echo "#-------------------------------------#";\
			echo "#                                     #";\
			echo "#   Downloading external Mesa-library #";\
			echo "#                                     #";\
			echo "#-------------------------------------#";\
			echo ;\
			$(DOWNLOAD) $(MESA_DOWNLOAD); \
		fi; \
		\
		echo ;\
		echo "#-----------------------------------#";\
		echo "#                                   #";\
		echo "#   Compiling external Mesa-library #";\
		echo "#                                   #";\
		echo "#-----------------------------------#";\
		echo ;\
		\
		if [ \( ! -d mesa-$(MESA_VER) \)  -a \( -f $(MESA_TGZ) \) ]; then \
			gunzip -c $(MESA_TGZ) | tar -xvf -; \
		fi; \
		cd mesa-$(MESA_VER); \
		./configure $(MESA_OPTIONS) --prefix=$(TOPDIR)/external; \
		$(MAKE); $(MAKE) install; \
	fi

glu-new: external_dirs
	@if test "$(COMPILE_MESA)" = yes  -o  "$(COMPILE_ALL)" = yes; then \
		\
		if test ! -f $(GLU_TGZ); then \
			echo ;\
			echo "#------------------------------------#";\
			echo "#                                    #";\
			echo "#   Downloading external GLU-library #";\
			echo "#                                    #";\
			echo "#------------------------------------#";\
			echo ;\
			$(DOWNLOAD) $(GLU_DOWNLOAD); \
		fi; \
		\
		echo ;\
		echo "#----------------------------------#";\
		echo "#                                  #";\
		echo "#   Compiling external GLU-library #";\
		echo "#                                  #";\
		echo "#----------------------------------#";\
		echo ;\
		\
		if [ \( ! -d glu-$(GLU_VER) \)  -a \( -f $(GLU_TGZ) \) ]; then \
			gunzip -c $(GLU_TGZ) | tar -xvf -; \
		fi; \
		cd glu-$(GLU_VER); \
		CPPFLAGS=-I$(TOPDIR)/external/include LDFLAGS=-L$(TOPDIR)/external/lib ./configure $(GLU_OPTIONS) --prefix=$(TOPDIR)/external; \
		$(MAKE); $(MAKE) install; \
	fi


mesa-traditional: external_dirs
	@if test "$(COMPILE_MESA)" = yes  -o  "$(COMPILE_ALL)" = yes; then \
		\
		if test ! -f $(MESA_TGZ); then \
			echo ;\
			echo "#-------------------------------------#";\
			echo "#                                     #";\
			echo "#   Downloading external Mesa-library #";\
			echo "#                                     #";\
			echo "#-------------------------------------#";\
			echo ;\
			$(DOWNLOAD) $(MESA_DOWNLOAD); \
		fi; \
		\
		echo ;\
		echo "#-----------------------------------#";\
		echo "#                                   #";\
		echo "#   Compiling external Mesa-library #";\
		echo "#                                   #";\
		echo "#-----------------------------------#";\
		echo ;\
		\
		if [ \( ! -d Mesa-$(MESA_VER) \)  -a \( -f $(MESA_TGZ) \) ]; then \
			gunzip -c $(MESA_TGZ) | tar -xvf -; \
		fi; \
		cd Mesa-$(MESA_VER); \
		$(MAKE) $(MESA_TARGET); \
		$(MAKE) install DESTDIR=$(TOPDIR)/external; \
		for dir in usr/local usr/X11; do \
			if test -d $(TOPDIR)/external/$$dir; then \
                                mv $(TOPDIR)/external/$$dir/lib/*     $(TOPDIR)/external/lib/; \
                                mv $(TOPDIR)/external/$$dir/include/* $(TOPDIR)/external/include/; \
                                rm -rf $(TOPDIR)/external/$$dir; \
                        fi; \
                done; \
	fi


fftw:
	@if test "$(COMPILE_FFTW)" = yes   -o  "$(COMPILE_ALL)" = yes; then \
		\
		if test ! -f $(FFTW_TGZ); then \
			echo ;\
			echo "#-------------------------------------#";\
			echo "#                                     #";\
			echo "#   Downloading external FFTW-library #";\
			echo "#                                     #";\
			echo "#-------------------------------------#";\
			echo ;\
			$(DOWNLOAD) $(FFTW_DOWNLOAD); \
		fi; \
		\
		echo ;\
		echo "#----------------------------------#";\
		echo "#                                  #";\
		echo "#   Compiling external FFTW-library #";\
		echo "#                                  #";\
		echo "#----------------------------------#";\
		echo ;\
		\
		if [ \( ! -d fftw-$(FFTW_VER) \)  -a \( -f $(FFTW_TGZ) \) ]; then \
			gunzip -c $(FFTW_TGZ) | tar xvf -; \
		fi; \
		\
		cd fftw-$(FFTW_VER); \
		./configure $(FFTW_OPTIONS) --prefix=$(TOPDIR)/external; \
		$(MAKE); $(MAKE) install; \
	fi


meschach: external_dirs
	@if test "$(COMPILE_MESCHACH)" = yes  -o  "$(COMPILE_ALL)" = yes; then \
		\
		if test ! -f $(MESCHACH_TGZ); then \
			echo ;\
			echo "#-----------------------------------------#";\
			echo "#                                         #";\
			echo "#   Downloading external MESCHACH-library #";\
			echo "#                                         #";\
			echo "#-----------------------------------------#";\
			echo ;\
			$(DOWNLOAD) $(MESCHACH_DOWNLOAD); \
		fi; \
		\
		echo ;\
		echo "#---------------------------------------#";\
		echo "#                                       #";\
		echo "#   Compiling external Meschach-library #";\
		echo "#                                       #";\
		echo "#---------------------------------------#";\
		echo ;\
		\
		if [ \( ! -d meschach$(MESCHACH_VER) \)  -a \( -f $(MESCHACH_TGZ) \) ]; then \
			mkdir meschach$(MESCHACH_VER); \
			cd meschach$(MESCHACH_VER); \
			gunzip -c ../$(MESCHACH_TGZ) | tar xvf -; \
			cd ..; \
		fi; \
		\
		cd meschach$(MESCHACH_VER); \
		$(MAKE) clean; $(MAKE); \
		if test ! -d $(TOPDIR)/external/include/meschach; then \
			mkdir $(TOPDIR)/external/include/meschach; \
		fi; \
		cp *.h $(TOPDIR)/external/include/meschach/; \
		mv meschach.a $(TOPDIR)/external/lib/libmeschach.a; \
	fi


bwidget:
	@if test ! -f $(BWIDGET_TGZ); then \
		echo ;\
		echo "#----------------------------------------#";\
		echo "#                                        #";\
		echo "#   Downloading external BWidget-library #";\
		echo "#                                        #";\
		echo "#----------------------------------------#";\
		echo ;\
		$(DOWNLOAD) $(BWIDGET_DOWNLOAD); \
	fi
	if test ! -d $(TOPDIR)/external/lib; then \
		mkdir $(TOPDIR)/external/lib; \
	fi
	if test -f $(BWIDGET_TGZ); then \
		if test ! -d $(TOPDIR)/external/lib/bwidget-$(BWIDGET_VER); then \
			cd $(TOPDIR)/external/lib; tar zxvf $(TOPDIR)/external/src/bwidget-$(BWIDGET_VER).tar.gz; \
		fi; \
	fi


external_dirs:
	for dir in lib include bin; do \
		if test ! -d $(TOPDIR)/external/$$dir; then \
			mkdir $(TOPDIR)/external/$$dir; \
		fi; \
	done

clean:
	for file in tcl*gz* tk*gz* Mesa*gz* mesa*gz* glu*gz* Togl*gz* bwidget*gz* fftw*gz* mesch*gz*; do \
		if test -f $$file; then rm -f $$file; fi; \
	done; \
	for dir in tcl$(TCL_VER3) tk$(TCL_VER3) Togl$(TOGL_VER) mesa-$(MESA_VER) glu-$(GLU_VER) Mesa-$(MESA_VER) bwidget-$(BWIDGET_VER) fftw-$(FFTW_VER) meschach$(MESCHACH_VER); do \
		if test -d $$dir ; then rm -rf $$dir ; fi; \
	done
