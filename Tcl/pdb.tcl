#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/pdb.tcl
# ------                                                                    #
# Copyright (c) 2019 by Anton Kokalj                                        #
#############################################################################

proc pdfUpdate {filedir {can .mesa}} {
    global select

    if { [info exists select(selection_mode)] && $select(selection_mode) } {
	tk_dialog [WidgetName] "WARNING !!!" "WARNING: structure can not be reloaded while in selection mode. First exit from the selection mode, then reload" warning 0 OK
	return
    }

    set file [gunzipXSF $filedir]

    if { ! [UpdateStruct $can $file pdb] } {
        return
    }
    xcUpdateState
    return
}
