############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/scroll.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
# ------                                                                    #
#############################################################################


proc mouseWheelScroll {w scrollCmd_B4 scrollCmd_B5 scrollCmd_Wheel} {
    global scrollWin tcl_platform

    #if { [info exist scrollWin($w)] } {
    #	# binding already set for $w
    #	return 1
    #}
    #
    #set scrollWin($w) 1

    set scrWin [getAllDescendantWid $w]

    foreach wid $scrWin {

	if { [winfo exists $wid] } {
	    if { $tcl_platform(platform) == "unix" } {
		bind $wid <Button-4> $scrollCmd_B4
		bind $wid <Button-5> $scrollCmd_B5
	    } else {
		# TODO: please tune the %D on windows
		bind $wid <MouseWheel>  $scrollCmd_Wheel
	    }    
	}
    }
}