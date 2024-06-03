#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/popupMenu.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc popupMenu {W x y} {
    if { [winfo exists $W.menu] } {
	destroy $W.menu
    }
    set m [menu $W.menu -tearoff 0]
    tk_popup $m $x $y
    popupMenu:popup $m $W
}

proc popupMenu:popup {m w} {
    $m add cascade  -image colors       -menu $m.vmcolor
    $m add cascade  -label "File"       -menu $m.vmfile 
    $m add cascade  -label "Display"    -menu $m.vmdis
    $m add cascade  -label "Modify"     -menu $m.vmmod
    $m add cascade  -label "AdvGeom"    -menu $m.vmadvg 
    $m add cascade  -label "Properties" -menu $m.vmpro 
    $m add cascade  -label "Tools"      -menu $m.vmdat
    $m add cascade  -label "Help"       -menu $m.vmhelp
    $m add separator
    $m add command  -label "Exit"        -command exit_pr

    ########################################################################
    set mcolor [menu $m.vmcolor -tearoff 0]
    set mfile  [menu $m.vmfile  -tearoff 0]
    set mmod   [menu $m.vmmod   -tearoff 0]
    set mdis   [menu $m.vmdis   -tearoff 0]
    set madvg  [menu $m.vmadvg  -tearoff 0]
    set mpro   [menu $m.vmpro   -tearoff 0]
    set mdat   [menu $m.vmdat   -tearoff 0]
    set mhelp  [menu $m.vmhelp  -tearoff 0]

    mainMenu $w $mcolor $mfile $mmod $mdis $madvg $mpro $mdat $mhelp
    update
    xcUpdateState $m {}
}
