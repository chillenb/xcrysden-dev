#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/xcAbout.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc xcAbout {} {
    global system

    set version [ReadFile $system(TOPDIR)/version]
    set text {
 =====================================================================
  About XCrySDen (version: $version)
 =====================================================================
 
 Copyright (C) 1996--2019 Anton Kokalj (tone.kokalj@ijs.si)
                          Jozef Stefan Institute, Ljubljana, Slovenia
 
 XCrySDen was written by Anton Kokalj, following a project of Mauro
 Causa. The project was initiated because Mauro Causa and Anton Kokalj
 felt a growing need for a simple visualization tool aimed at
 displaying the crystalline structures. Programming started in 1996
 and the first implementation of the program was made available in
 1999.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
}

    set text [subst $text]
    set t [xcDisplayVarText $text "XCrySDen: About"]
    if { [winfo exists $t.f1.t] } {
	catch {$t.f1.t configure -state disabled}
    }
}
