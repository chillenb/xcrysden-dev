#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/loadAttributes.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc xsfLoadAttributes {file} {
    global mody
    
    set elementcolor  0
    set elementradius 0
    set genericattr   0         
    foreach line [split [ReadFile $file] \n] {
	xcDebug -stderr "---> $elementradius ;;; $line"
	if { [llength $line] == 0 } {
	    # an empty-line --> skip it
	    continue
	} elseif { [regexp -- {^[a-zA-Z]a*} $line] } {
	    set keyword [string trim $line { }]; # get rid of trailing white-spaces
	    switch -glob -- $keyword {
		"ELEMENTCOL*" { 
		    set elementcolor  1 
		    set elementradius 0
		    set genericattr   0         
		}
		"ELEMENTRAD*" {
		    set elementcolor  0
		    set elementradius 1
		    set genericattr   0         
		}
		default {
		    set elementcolor  0
		    set elementradius 0
		    set genericattr   1         

		    set key     [lindex $keyword 0]
		    set values  [lrange $keyword 1 end]
		    xcDebug -stderr "----> generic option: $keyword $values"
		    if { ! [info exists mody($key)] } {
			ErrorDialog "while parsing attribute file \"$file\", while reading string \"$keyword\""
		    }
		    if { [catch {eval {xc_newvalue .mesa $mody($key)} $values} error ] } {
			ErrorDialog "while parsing attribute file \"$file\", while reading string \"$keyword\".\n\nError message: $error"
		    }
		}	
	    }
	} else {
	    if { $elementcolor } {

		if { [llength $line] != 4 } {
		    ErrorDialog "should get 4 numbers (i.e. atomic-number r g b), but got: $line\nwhile reading attribute file: $file"
		    return
		}
		xc_newvalue .mesa $mody(L_ATCOL_ONE) [lindex $line 0] \
		    [lindex $line 1]  [lindex $line 2]  [lindex $line 3]

	    } elseif { $elementradius } {

		if { [llength $line] != 2 } {
		    ErrorDialog "should get 2 numbers (i.e. atomic-number atomic-radius), but got: $line\nwhile reading attribute file: $file"
		    return
		}
		xc_newvalue .mesa $mody(L_ATRAD_ONE) [lindex $line 0] \
		    [lindex $line 1]  
		xc_newvalue .mesa $mody(L_RCOV_ONE)  [lindex $line 0] \
		    [lindex $line 1]  
	    }
	}
    }
}
		

# ------------------------------------------------------------------------
# loads the attributes from a definition file such as 
# ~/.xcrysden/custom-definitions
# ------------------------------------------------------------------------
proc defLoadAttributes {} {
    global mody atmRad atmCol 
    
    #
    # load atmRad
    #
    foreach atmNum [array names atmRad] {
	if { $atmNum < 0 || $atmNum > 100 } {
	    ErrorDialog "atomic number out of range. Check the atmRad array in the deifinition file"
	}
	foreach code [list $mody(L_ATRAD_ONE) $mody(L_RCOV_ONE)] {
	    xcCatchEval \
		[list xc_newvalue .mesa $code $atmNum $atmRad($atmNum)] \
		"Can't change the atomic radius for atomic-number $atmNum. Check the definition of atmRad($atmNum) in the definition file."
	}
    }
    
    #
    # load atmCol
    #
    foreach atmNum [array names atmCol] {
	if { $atmNum < 0 || $atmNum > 100 } {
	    ErrorDialog "atomic number out of range. Check the atmCol array in the deifinition file"
	}
	xcCatchEval \
	    [concat xc_newvalue .mesa $mody(L_ATCOL_ONE) $atmNum $atmCol($atmNum)] \
	    "Can't change the atomic color for atomic-number $atmNum. Check the definition of atmCol($atmNum) in the definition file."
    }

    #
    # load myParam
    #
    load_myParam
}


proc load_myParam {} {
    global myParam mody

    if { [info exists myParam] } {
	foreach _elem [array names myParam] {
	    switch -exact -- $_elem {
		"ATRAD_SCALE" - 
		"TESSELLATION" -
		"UNIBONDCOLOR" -
		"PERSPECTIVEFOVY" - 
		"PERSPECTIVEFRONT" -
		"PERSPECTIVEBACK" -
		"BALLF" - 
		"RODF" - 
		"WFLINEWIDTH" - 
		"PLLINEWIDTH" - 
		"FRAMELINEWIDTH" -
		"OUTLINEWIDTH" -
		"WF3DLINEWIDTH" -
		"PLRADIUS" - 
		"COV_SCALE" - 
		"FRAMECOL" - 
		"FRAMELINEWIDTH" - 
		"FRAMERODF" - 
		"BACKGROUND" {
		    
		    # myParam array element is OK !!!

		    set elem L_$_elem
		    xcCatchEval \
			[list eval "xc_newvalue .mesa $mody($elem)" $myParam($_elem)] \
			"Can't change the $_elem attribute. Check the definition of myParam($_elem) in the definition file."
		}

		"CRYSTAL_MAXCELL" -
		"SLAB_MAXCELL" -
		"POLYMER_MAXCELL" -
		"MPEG_ENCODE_PARAM_FILE" {
		    continue
		}

		"ATOMIC_LABEL_FONT" -
		"ATOMIC_LABEL_BRIGHTCOLOR" -
		"ATOMIC_LABEL_DARKCOLOR" {
		    global atomLabel
		    if { ! [info exists myParam(ATOMIC_LABEL_BRIGHTCOLOR)] } {
			set myParam(ATOMIC_LABEL_BRIGHTCOLOR) [xc_getvalue .mesa $mody(GET_GLOBALATOMLABEL_BRIGHTCOLOR)]
		    }
		    if { ! [info exists myParam(ATOMIC_LABEL_DARKCOLOR)] } {
			set myParam(ATOMIC_LABEL_DARKCOLOR) [xc_getvalue .mesa $mody(GET_GLOBALATOMLABEL_DARKCOLOR)]
		    }
		    if { ! [info exists myParam(ATOMIC_LABEL_FONT)] } {
			set myParam(ATOMIC_LABEL_FONT) raster
		    }

		    set atomLabel(globalFont.brightColor) $myParam(ATOMIC_LABEL_BRIGHTCOLOR)
		    set atomLabel(globalFont.darkColor)   $myParam(ATOMIC_LABEL_DARKCOLOR)
		    
		    puts stderr "DEBUG> xc_setfont .mesa $myParam(ATOMIC_LABEL_FONT) \
			$myParam(ATOMIC_LABEL_BRIGHTCOLOR) \
			$myParam(ATOMIC_LABEL_DARKCOLOR)"

		    xc_setfont .mesa $myParam(ATOMIC_LABEL_FONT) \
			$myParam(ATOMIC_LABEL_BRIGHTCOLOR) \
			$myParam(ATOMIC_LABEL_DARKCOLOR)
		}		

		"FS_BACKGROUND" - 
		"FS_CELLTYPE" - 
		"FS_CROPBZ" - 
		"FS_DISPLAYCELL" - 
		"FS_CELLDISPLAYTYPE" - 
		"FS_DRAWSTYLE" - 
		"FS_TRANSPARENT" - 
		"FS_SHADEMODEL" - 
		"FS_INTERPOLATIONDEGREE" - 
		"FS_FRONTFACE" - 
		"FS_REVERTNORMALS" - 
		"FS_WIRECELLCOLOR" - 
		"FS_SOLIDCELLCOLOR" -
		"FS_ANTIALIAS" -     
		"FS_DEPTHCUING" {
		    #
		    # Fermi-Surface customizations (do nothing)
		    #
		    DummyProc
		}
		default {
		
		    # myParam array element is not recognized !!!

		    ErrorDialog "syntax error in custom-definition file: myParam array element \"$_elem\" is not recognized. Check the deifinition file."
		}
	    }
	}
    }
}
