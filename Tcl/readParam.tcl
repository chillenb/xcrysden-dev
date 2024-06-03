#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/readParam.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc read_param { w } {
    global groupsel n_groupsel species inp class 
    # THIS PROC DETERMINE A CLASS OF SPECIES AND READ MINIMAL 
    # SET OF LATTICE VECTORS & ANGLES

    # $w         - name of parent's widget
    # $species   - type of species (polymer, slab, ..)
    # n_groupsel - N. of selected group

	
    # ----------------------------------------------
    # DEFINITIONS of needed SET's of LATTICE VECTORS
    # ----------------------------------------------

    # CLASSES for whom only A is NEEDED: "A" class
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # Polymers: -all
    # Slabs:    -square
    #           -hexagonal
    # Crystals: -cubic

    # CLASSES for whom only A,B is NEEDED: "AB" class
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # Slabs:    -rectangular
    
    # CLASSES for whom only A,B,GAMMA is NEEDED: "AB_GAMMA" class
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # Slabs:    -triclinic
	
    # CLASSES for whom only A,C is NEEDED: "AC" class
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # Crystals: -hexagonal
    #           -rhombohedral-hexagonal
    #           -tetragonal
    
    # CLASSES for whom only A,ALFA is NEEDED: "A_ALFA" class
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # Crystals: -rhombohedral-rhombohedral

    # REMARK: the following classes can be handled either RHOMBO.,
    # ------  either HEXA.:   146    R 3
    #                         148    R -3
    #                         155    R 3 2
    #                         160    R 3 M
    #                         161    R 3 C
    #                         166    R -3 M
    #                         167    R -3 C
    # that is ALL groups that are R-centered

    # CLASSES for whom only A,B,C is NEEDED: "ABC" class
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # Crystals: -orthorhombic
    
    # CLASSES for whom only A,B,C,"angle" is NEEDED: "ABC_angle" class
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # Crystals: -monoclinic
    
    # CLASSES for whom A,B,C,ALFA,BETA,GAMMA is NEEDED: "ABC_all"class
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # Crystals: -triclinic

    frame $w.1
    frame $w.2
    pack $w.1 $w.2 -expand 1 -fill both
    # all polymer are in "A" class
    if { $species == "polymer" } { 
	Entries $w.1 A: inp(A) 10 
	set class A
    }
    
    # we must determine class's number 
    # if $species == (slab or crystal); 
    
    if { $species == "slab" } {
	if { $n_groupsel < 8 } {
	    # TRICLINIC LATTICE
	    Entries $w.1 {A: B:} {inp(A) inp(B)} 10
	    Entries $w.2 GAMMA: inp(GAMMA) 10
	    set class {{A B} {GAMMA}}
	}
	if { $n_groupsel >= 8 && $n_groupsel < 49 } { 
	    # RECTANGULAR LATTICE
	    Entries $w.1 {A: B:} {inp(A) inp(B)} 10
	    set class {{A B}}
	}
	if { $n_groupsel >= 49 && $n_groupsel < 65 } { 
	    # SQUARE LATTICE
	    Entries $w.1 A: inp(A) 10
	    set class A
	}
	if { $n_groupsel >= 65 } { 
	    # HEXAGONAL LATTICE
	    Entries $w.1 A: inp(A) 10
	    set class A
	}
    }
    
    if { $species == "crystal" } {
	if { $n_groupsel < 3 } {
	    # TRICLINIC LATTICE
	    Entries $w.1 {A: B: C:} {inp(A) inp(B) inp(C)} 10
	    Entries $w.2 {ALPHA: BETA: GAMMA:} \
		    {inp(ALFA) inp(BETA) inp(GAMMA)} 10
	    set class {{A B C} {ALFA BETA GAMMA}}
	}
	if { $n_groupsel >= 3 && $n_groupsel < 16 } {
	    # MONOCLINIC LATTICE
	    # WARNING: here are three options !!!!
	    # Entries $w.2 {ALFA: BETA: GAMMA:} {ALFA BETA GAMMA} 10
	    Entries $w.1 {A: B: C:} {inp(A) inp(B) inp(C)} 10
	    Entries $w.2 GAMMA: inp(GAMMA) 10
	    set class {{A B C} {GAMMA}}
	}
	if { $n_groupsel >= 16 && $n_groupsel < 75 } {
	    # ORTHOROMBIC LATTICE
	    Entries $w.1 {A: B: C:} {inp(A) inp(B) inp(C)} 10
	    set class {{A B C}}
	}
	if { $n_groupsel >= 75 && $n_groupsel < 143 } {
	    # TETRAGONAL LATTICE
	    Entries $w.1 {A: C:} {inp(A) inp(C)} 10
	    set class {{A C}}
	}
	# 143 -> 167 some TRIGONAL->HEXA, some TRIGONAL->RHOMBO.
	# 168 -> 194 HEXAGONAL
	if [info exist inp(IFHR)] {
	    if { $inp(IFHR) == 1 } {
		# RHOMBOHEDRAL AXES
		Entries $w.1 {A: ALPHA:} {inp(A) inp(ALFA)} 10
		set class {{A} {ALFA}}
		return
	    }
	}
	if { $n_groupsel >= 143 && $n_groupsel < 195 } {
	    # HEXAGONAL LATTICE & TRIGONAL WITH HEXA. AXES
	    Entries $w.1 {A: C:} {inp(A) inp(C)} 10
	    set class {{A C}}
	}
	if { $n_groupsel >= 195 } {
	    # CUBIC LATTICE
	    Entries $w.1 A: inp(A) 10
	    set class A
	}
    }
    # set focus to the first entry
    focus $w.1.frame.entry1
}
	
    
proc WhichPar2Read {input} {
    global groupsel n_groupsel species \
	    inp class distext
    
    #this procs READ&PRINT parameters in OPEN_FILE mode
    set param [gets $input]
    if { $species == "slab" } {
	append distext "> UNIT CELL PARAMETERS - "
	if { $n_groupsel < 8 } {
	    # TRICLINIC LATTICE
	    set inp(A) [lindex $param 0]
	    set inp(B) [lindex $param 1]
	    set inp(GAMMA) [lindex $param 2]
	    append distext "TRICLINIC LATTICE::\n"
	    append distext "A = $inp(A),    B = $inp(B),\
		    \n    GAMMA = $inp(GAMMA)\n\n"
	}
	if { $n_groupsel >= 8 && $n_groupsel < 49 } { 
	    # RECTANGULAR LATTICE
	    set inp(A) [lindex $param 0]
	    set inp(B) [lindex $param 1]
	    append distext "RECTANGULAR LATTICE::\n"
	    append distext "A = $inp(A),    B = $inp(B)\n\n"
	}
	if { $n_groupsel >= 49 && $n_groupsel < 65 } { 
	    # SQUARE LATTICE
	    set inp(A) [lindex $param 0]
	    append distext "SQUARE LATTICE::\n"
	    append distext "A = $inp(A)\n\n"
	}
	if { $n_groupsel >= 65 } { 
	    # HEXAGONAL LATTICE
	    set inp(A) [lindex $param 0]
	    append distext "HEXAGONAL LATTICE::\n"
	    append distext "A =  $inp(A)\n\n"
	}
    }
    
    if { $species == "crystal" } {
	IsNotStandardGroup $groupsel
	append distext "> CRYSTAL FAMILY               : $inp(CRY_FAM)\n"
	append distext "> CRYSTAL CLASS  (GROTH - 1921): $inp(CRY_CLASS)\n\n"
	append distext "> UNIT CELL PARAMETERS::\n"
	if { $n_groupsel <= 230 && [xcNumber n_groupsel int] } {
	    if { $n_groupsel < 3 } {
		# TRICLINIC LATTICE	
		set inp(A) [lindex $param 0]
		set inp(B) [lindex $param 1]
		set inp(C) [lindex $param 2]
		set inp(ALFA) [lindex $param 3]
		set inp(BETA) [lindex $param 4]
		set inp(GAMMA) [lindex $param 5]
		append distext "A = $inp(A),    B = $inp(B),    C = $inp(C),\n"
		append distext "ALPHA = $inp(ALFA),    BETA = $inp(BETA),    \
			GAMMA = $inp(GAMMA)\n\n"
	    }
	    if { $n_groupsel >= 3 && $n_groupsel < 16 } {
		# MONOCLINIC LATTICE
		# WARNING: here are three options !!!!
		set inp(A) [lindex $param 0]
		set inp(B) [lindex $param 1]
		set inp(C) [lindex $param 2]
		set inp(GAMMA) [lindex $param 3]
		append distext "A = $inp(A),    B = $inp(B),    C = $inp(C),\n"
		append distext "GAMMA = $inp(GAMMA)\n\n"
	    }
	    if { $n_groupsel >= 16 && $n_groupsel < 75 } {
		# ORTHOROMBIC LATTICE
		set inp(A) [lindex $param 0]
		set inp(B) [lindex $param 1]
		set inp(C) [lindex $param 2]		
		append distext "A = $inp(A),    B = $inp(B),    C = $inp(C),\n\n"
	    }
	    if { $n_groupsel >= 75 && $n_groupsel < 143 } {
		# TETRAGONAL LATTICE
		set inp(A) [lindex $param 0]
		set inp(B) [lindex $param 1]
		append distext "A = $inp(A),    B = $inp(B)\n\n"
	    }
	    # 143 -> 167 some TRIGONAL->HEXA, some TRIGONAL->RHOMBO.
	    # 168 -> 194 HEXAGONAL
	    if { $n_groupsel >= 143 && $n_groupsel < 168 } {
		# TRIGONAL LATTICE (either rhombo. either hexa. axes)
		set inp(A) [lindex $param 0]			
		set rhombohexa "146 148 155 160 161 166 167"
		foreach grp $rhombohexa {
		    if { $n_groupsel == $grp } {
			# query for rhombo./hexa. axes
			if { $inp(IFHR) == 1 } {
			    # RHOBOHEDRAL AXES!!!
			    set inp(ALFA) [lindex $param 1]
			    append distext "A = $inp(A),    ALPHA = $inp(ALFA)\n\n"
			    return
			}
		    }
		}
		# HEXAGONAL AXES!!!
		set inp(C) [lindex $param 1]
		append distext "A = $inp(A),    C = $inp(C)\n\n"
	    }	
	    if { $n_groupsel >= 168 && $n_groupsel < 195 } {
		# HEXAGONAL LATTICE
		set inp(A) [lindex $param 0]
		set inp(C) [lindex $param 1]
		append distext "A = $inp(A),    C = $inp(C)\n\n"
	    }
	    if { $n_groupsel >= 195 } {
		# CUBIC LATTICE
		set inp(A) [lindex $param 0]
		append distext "A = $inp(A)\n\n"
	    }
	} else {
	    # do it according to class
	    set i 0
	    foreach grp $class {
		foreach item $grp {
		    if { $item == "A" } {
			set inp(A) [lindex $param $i]
			append distext "A = $inp(A)"
			incr i
		    }
		    if { $item == "B" } {
			set inp(B) [lindex $param $i]
			append distext ",    B = $inp(B)"
			incr i
		    }
		    if { $item == "C" } {
			set inp(C) [lindex $param $i]
			append distext ",    C = $inp(C)"
			incr i
		    }
		    if { $item == "ALFA" } {
			set inp(ALFA) [lindex $param $i]
			append distext ",\nALPHA = $inp(ALFA)"
			incr i
		    }
		    if { $item == "BETA" } {
			set inp(BETA) [lindex $param $i]
			append distext ",    BETA = $inp(BETA)"
			incr i
		    }
		    if { $item == "GAMMA" } {
			set inp(GAMMA) [lindex $param $i]
			if { [info exists inp(ALFA)] || \
				[info exists inp(BETA)] } {
			    append distext ",    GAMMA = $inp(GAMMA)"
			} else {
			    append distext "GAMMA = $inp(GAMMA)"
			}
			incr i
		    }			
		}
	    }
	    append distext "\n\n"
	}
    } 
    append distext "--------------------------------------------------\n\n"
}


proc WhichPar2Print {} {
    global groupsel n_groupsel species \
	    inp class distext
    
    #this procs PRINTS parameters in OPEN_FILE mode
    if { $species == "slab" } {
	append distext "> UNIT CELL PARAMETERS - "
	if { $n_groupsel < 8 } {
	    # TRICLINIC LATTICE
	    append distext "TRICLINIC LATTICE::\n"
	    append distext "A = $inp(A),    B = $inp(B),\
		    \n    GAMMA = $inp(GAMMA)\n\n"
	}
	if { $n_groupsel >= 8 && $n_groupsel < 49 } { 
	    # RECTANGULAR LATTICE
	    append distext "RECTANGULAR LATTICE::\n"
	    append distext "A = $inp(A),    B = $inp(B)\n\n"
	}
	if { $n_groupsel >= 49 && $n_groupsel < 65 } { 
	    # SQUARE LATTICE
	    append distext "SQUARE LATTICE::\n"
	    append distext "A = $inp(A)\n\n"
	}
	if { $n_groupsel >= 65 } { 
	    # HEXAGONAL LATTICE
	    append distext "HEXAGONAL LATTICE::\n"
	    append distext "A =  $inp(A)\n\n"
	}
    }
    
    if { $species == "crystal" } {
	IsNotStandardGroup $groupsel
	append distext "> CRYSTAL FAMILY               : $inp(CRY_FAM)\n"
	append distext "> CRYSTAL CLASS  (GROTH - 1921): $inp(CRY_CLASS)\n\n"
	append distext "> UNIT CELL PARAMETERS::\n"	
	puts stdout "IGR:: $n_groupsel"
	if { $n_groupsel <= 230 && [xcNumber n_groupsel int] } {
	    if { $n_groupsel < 3 } {
		# TRICLINIC LATTICE	
		append distext "A = $inp(A),    B = $inp(B),    C = $inp(C),\n"
		append distext "ALPHA = $inp(ALFA),    BETA = $inp(BETA),    \
			GAMMA = $inp(GAMMA)\n\n"
	    }
	    if { $n_groupsel >= 3 && $n_groupsel < 16 } {
		# MONOCLINIC LATTICE
		# WARNING: here are three options !!!!
		append distext "A = $inp(A),    B = $inp(B),    C = $inp(C),\n"
		append distext "GAMMA = $inp(GAMMA)\n\n"
	    }
	    if { $n_groupsel >= 16 && $n_groupsel < 75 } {
		# ORTHOROMBIC LATTICE
		append distext "A = $inp(A),    B = $inp(B),    C = $inp(C),\n\n"
	    }
	    if { $n_groupsel >= 75 && $n_groupsel < 143 } {
		# TETRAGONAL LATTICE
		append distext "A = $inp(A),    B = $inp(B)\n\n"
	    }
	    if { $n_groupsel >= 143 && $n_groupsel < 168 } {
		# TRIGONAL LATTICE (either rhombo. either hexa. axes)
		set rhombohexa "146 148 155 160 161 166 167"
		foreach grp $rhombohexa {
		    if { $n_groupsel == $grp } {
			# query for rhombo./hexa. axes
			if { $inp(IFHR) == 1 } {
			    # RHOBOHEDRAL AXES!!!
			    append distext "A = $inp(A),    ALPHA = $inp(ALFA)\n\n"
			    return
			}
		    }
		}
		# HEXAGONAL AXES!!!
		append distext "A = $inp(A),    C = $inp(C)\n\n"
	    }
	    if { $n_groupsel >= 168 && $n_groupsel < 195 } {
		# HEXAGONAL LATTICE
		append distext "A = $inp(A),    C = $inp(C)\n\n"
	    }
	    if { $n_groupsel >= 195 } {
		# CUBIC LATTICE
		append distext "A = $inp(A)\n\n"
	    }
	} else {
	    # do it according to class
	    set i 0
	    foreach grp $class {
		foreach item $grp {
		    puts stdout "WhichPar2Print:: ITEM> $item"
		    if { $item == "A" } {		
			append distext "A = $inp(A)"
			incr i
		    }
		    if { $item == "B" } {
			append distext ",    B = $inp(B)"
			incr i
		    }
		    if { $item == "C" } {
			append distext ",    C = $inp(C)"
			incr i
		    }
		    if { $item == "ALFA" } {
			append distext ",\nALPHA = $inp(ALFA)"
			incr i
		    }
		    if { $item == "BETA" } {
			append distext ",    BETA = $inp(BETA)"
			incr i
		    }
		    if { $item == "GAMMA" } {
			if { [info exists inp(ALFA)] || \
				[info exists inp(BETA)] } {
			    append distext ",    GAMMA = $inp(GAMMA)"
			} else {
			    append distext "GAMMA = $inp(GAMMA)"
			}
			incr i
		    }			
		}
	    }
	    append distext "\n\n"
	}	    
    }
    append distext "--------------------------------------------------\n\n"
}

