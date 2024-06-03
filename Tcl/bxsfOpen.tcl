#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/bxsfOpen.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc bxsfOpenMenu {{can {}}} {
    global system

    set file [tk_getOpenFile -defaultextension .bxsf \
		  -filetypes {
		      {{All Files}          {.*}      }
		      {{BXSF XSF Files}     {.bxsf}   }
		      {{GZipped BXSF Files} {.bxsf.gz}}
		      {{XSF Files}          {.xsf}    }
		      {{GZipped XSF Files}  {.xsf.gz} }
		  } \
		  -initialdir $system(PWD) \
		  -title "Open BXSF (i.e. Fermi Surface) File"]
    if { $file == "" } {
	return
    }

    genFSInit $file
}

	