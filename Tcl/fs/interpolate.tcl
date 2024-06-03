#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/interpolate.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################


proc FS_Interpolation:submit {togl spin iband} {
    global FS_Interpolation fs
    
    set fs($spin,$iband,interpolationdegree) [list \
						  $FS_Interpolation($togl) \
						  $FS_Interpolation($togl) \
						  $FS_Interpolation($togl)]
    FS_fsConfig $iband $spin
}
    


