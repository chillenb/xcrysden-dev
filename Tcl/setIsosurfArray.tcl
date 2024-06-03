#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/setIsosurfArray.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc SetIsoSurfArray {} {
    global isosurf

    set isosurf(spin)             {}
    set isosurf(res_type)         "angstroms"
    set isosurf(mb_angs/bohr)     "Angstroms"
    set isosurf(space_sel)        "whole_cell"
    set isosurf(Y_Sel)            "centered"
    set isosurf(Z_Sel)            "centered"
    set isosurf(expand)           "none"
    set isosurf(expand_X)         1
    set isosurf(expand_Y)         1
    set isosurf(expand_Z)         1
    set isosurf(type_of_isosurf)  "solid"
    set isosurf(shade_model)      "smooth"
    set isosurf(transparency)     "off" 
    set isosurf(twoside_lighting) [xc_getGLparam lightmodel -get two_side_iso]
    set isosurf(tessellation_type) cubes
    set isosurf(normals_type)      gradient
}
