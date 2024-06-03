#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/sInfo.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# Tools menu path is is .menu.vmdat
proc Get_sInfoArray {} {
    global sInfo periodic

    xcDebug "######################################"
    if { [array exists sInfo] } {
	xcDebug -debug [array get sInfo *]
    } else {
	xcDebug "array s Info does not EXISTS"
    }
    xcDebug "######################################"
    
    set m .menu.vmdat
    if { $sInfo(ldatagrid2D) || $sInfo(ldatagrid3D) } {
	$m.menu entryconfig "Data Grid" -state normal
    } else {
	$m.menu entryconfig "Data Grid" -state disabled
    }

    if { $sInfo(lprimvec) || $sInfo(lconvvec) } {
	set periodic(igroup) $sInfo(groupn)
    }

    set periodic(dim) $sInfo(dim)
    
    if { !$sInfo(lprimvec) && !$sInfo(lconvvec) } {
	set periodic(dim) 0
    }

    if { $periodic(dim) == 0 } {
	$m.menu entryconfig {k-path Selection} -state disabled
    } else {
	$m.menu entryconfig {k-path Selection} -state normal
    }
    xcDebug -stderr "sInfo(dim) = $sInfo(dim); periodic(dim) = $periodic(dim)"
    
    if { $sInfo(lforce) } {
	.menu.vmdis.menu entryconfig {Forces} -state normal
	.menu.vmmod.menu entryconfig {Force Settings} -state normal
    } else {
	.menu.vmdis.menu entryconfig {Forces} -state disabled
	.menu.vmmod.menu entryconfig {Force Settings} -state disabled
    }

    # initialize WignerSeitz
    WignerSeitzInit
}
