#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/kLabels.tcl                                  #     
# ------                                                                    #
# Copyright (c) 2005 by Anton Kokalj                                        #
#                                                                           #
# ---------------------------------------------------------------------------
# The labelling of k-labels is based on the idea and lookup table of Peter Blaha.
# ---------------------------------------------------------------------------
# Peter BLAHA, Inst.f.Techn.Elektrochemie, TU Vienna, A-1060 Vienna
# Phone: +43-1-58801-5187             FAX: +43-1-5868937
# Email: pblaha@email.tuwien.ac.at    WWW: http://www.tuwien.ac.at/theochem/
# ---------------------------------------------------------------------------


#
# This proc determines the Bravais lattice type on the basis of XSF's
# igroup and primitive and convetional lattice vectors. The following
# lattices are currently supported for k-point labelling:
# 
# P-cubic, F-cubic, I-cubic
# hexagonal
# P-tetragonal, I-tetragonal
# P-orthorhombic
# 
# The proc returns one among above string or the "not-supported" string.
#
proc igroup2BravaisLattice {igroup} {
    global kLabels
  
    # dummy defs, so that these vars exists
    set vec(0,0)   0.0
    set vec_len(0) 0.0
    set dot(0,0)   0.0    
    
    switch -- $igroup {
	1 {
	    # P-lattice, possibilities: cubic, tetragonal, orthorhombic
	    getLatticeVec_ primvec vec vec_len dot

	    # check if lattice is orthogonal
	    if { [IsEqual 1e-6 0.0 $dot(0,1) $dot(0,2) $dot(1,2)] } {

		# check the lengths of lattice vectors
		if { [IsEqual 1e-6 $vec_len(0) $vec_len(1) $vec_len(2)] } {
		    set lattice "P-cubic"
		} elseif { [IsEqual 1e-6 $vec_len(0) $vec_len(1)] } {
		    set lattice "P-tetragonal"
		} else {
		    set lattice "not-supported"
		}
	    } else {
		set lattice "not-supported"
	    }
	}
	2 { set lattice "not-supported" }
	3 { set lattice "not-supported" }
	4 { set lattice "not-supported" }
	5 { set lattice "F-cubic" }
	6 {
	    # I-lattice, possibilities: cubic, tetragonal
	    # made a similar check as for P-lattice, but with conventional vectors 
	    getLatticeVec_ convvec vec vec_len dot

	    # check if lattice is orthogonal
	    if { [IsEqual 1e-6 0.0 $dot(0,1) $dot(0,2) $dot(1,2)] } {

		# check the lengths of lattice vectors
		if { [IsEqual 1e-6 $vec_len(0) $vec_len(1) $vec_len(2)] } {
		    set lattice "I-cubic"
		} elseif { [IsEqual 1e-6 $vec_len(0) $vec_len(1)] } {
		    set lattice "I-tetragonal"
		} else {
		    xcDebug -stderr "*** igroup2BravaisLattice: impossible Bravais lattice (I); check code"
		    set lattice "not-supported"
		}
	    } else {
		set lattice "not-supported"
	    }	    
	}
	7 { set lattice "not-supported" }
	8 { set lattice "hexagonal" }
	9 { set lattice "not-supported" }
	default { set lattice "not-supported" }
    }
    return $lattice
}

#
# Private proc used by igroup2BravaisLattice, to get the attributes
# of primitive or conventional.
#
proc getLatticeVec_ {type vec_ vec_len_ dot_} {
    global sInfo
    upvar $vec_     vec
    upvar $vec_len_ vec_len
    upvar $dot_     dot
        
    #xcDebug -stderr "tk: SINFO(primvec) == $sInfo(primvec)"
    #xcDebug -stderr "tk: SINFO(convvec) == $sInfo(convvec)"

    # get the primitive lattice vectors and its lengths
    for {set i 0} {$i < 3} {incr i} {
	for {set j 0} {$j < 3} {incr j} {
	    set ind [expr $i*3 + $j] 
	    set vec($i,$j) [lindex $sInfo($type) $ind]
	}
	set vec_len($i) [expr sqrt($vec($i,0)*$vec($i,0) + $vec($i,1)*$vec($i,1) + $vec($i,2)*$vec($i,2))]
	#xcDebug -stderr "tk: VEC_LEN($i) == $vec_len($i)"
    }
    set dot(0,1) [expr $vec(0,0)*$vec(1,0)  +  $vec(0,1)*$vec(1,1)  +  $vec(0,2)*$vec(1,2)]
    set dot(0,2) [expr $vec(0,0)*$vec(2,0)  +  $vec(0,1)*$vec(2,1)  +  $vec(0,2)*$vec(2,2)]
    set dot(1,2) [expr $vec(1,0)*$vec(2,0)  +  $vec(1,1)*$vec(2,1)  +  $vec(1,2)*$vec(2,2)]		   		 
    #xcDebug -stderr "tk: $type DOT:: $dot(0,1) $dot(0,2) $dot(1,2)"
}


#
# This proc tries to return the label of selected k-point based on
# table-lookup.
#
proc getKLabel {latticeType kx ky kz} {
    global  kLabels
    
    if { $latticeType == "not-supported" } { 
	return "" 
    }
    if { ! [info exists  kLabels($latticeType)] } { 
	return "" 
    }
    
    foreach {kxx kyy kzz label} $kLabels($latticeType) {
	#xcDebug -stderr "tk: PREDIFINED: $kxx $kyy $kzz  <--- SELECTED: $kx $ky $kz"

	if { [IsEqual 1e-5 $kx $kxx]  &&  [IsEqual 1e-5 $ky $kyy]  &&  [IsEqual 1e-5 $kz $kzz] } {
	    #xcDebug -stderr "tk: LABEL == $label"
	    return $label
	}
    }
    
    return ""
}


proc kLabels_Note {} {
    global Bz periodic kLabels

    #        WARNING: since labeling of the k-points is a new feature, read carefully 
    #                 below warnings:

    set msg {
	NEW FEATURE: automatic labeling of the k-points
	-----------------------------------------------

	1. For a few "supported" Bravais lattices several k-points will be 
           labeled automatically 
           (Ref: http://www.cryst.ehu.es/cryst/get_kvec.html)

	2. The information within the XSF file are sometimes insufficient
	   to determine the Bravais lattice type. The labelling of the k-points 
           will be hopefully correct only if correct Bravais lattice type was 
           determined

	*** CHECK THIS DATA:
	
        - the guessed BRAVAIS LATTICE TYPE : $Bz(lattice_type)
	 ((the XSF's group number is $periodic(igroup)))
    }
    set Msg [subst -nocommands $msg]

    set t [xcDisplayVarText $Msg "Automatic labeling of k-points"]

    tkwait window $t

    #if { ! [info exists kLabels(warning_window)] } {
    #	set t [xcDisplayVarText $Msg "Automatic labeling of k-points"]
    #} elseif { ! [winfo exists $kLabels(warning_window)] } {
    #	set kLabels(warning_window) [xcDisplayVarText $Msg "Automatic labeling of k-points"]
    #}
}

#
# this proc loads the k-label's lookup table
#
proc load_kLabels {} {
    global kLabels

    # not supported lattices
    set kLabels(not-supported) {}

    # lattice type: P cubic (eg. SG 221 Pm-3m)
    set kLabels(P-cubic) {
	0.0 0.0 0.0     GAMMA
	0.5 0.5 0.5     R
	-.5 0.5 0.5     R
	0.5 -.5 0.5     R
	0.5 0.5 -.5     R
	0.5 -.5 -.5     R
	-.5 -.5 0.5     R
	-.5 0.5 -.5     R
	-.5 -.5 -.5     R
	0.5 0.0 0.0     X
	0.0 0.5 0.0     X
	0.0 0.0 0.5     X
	-.5 0.0 0.0     X
	0.0 -.5 0.0     X
	0.0 0.0 -.5     X
	0.5 0.5 0.0     M
	-.5 0.5 0.0     M
	0.5 -.5 0.0     M
	-.5 -.5 0.0     M
	0.5 0.0 0.5     M
	-.5 0.0 0.5     M
	0.5 0.0 -.5     M
	-.5 0.0 -.5     M
	0.0 0.5 0.5     M
	0.0 -.5 0.5     M
	0.0 0.5 -.5     M
	0.0 -.5 -.5     M
    }

    # lattice type: F cubic (eg. SG 225 Fm-3m)
    set kLabels(F-cubic) {
	0.0 0.0 0.0     GAMMA
	0.5 0.5 0.5     L   
	-0.5 -0.5 -0.5  L   
	0.5 0.0 0.0     L   
	0.0 0.5 0.0     L   
	0.0 0.0 0.5     L   
	-.5 0.0 0.0     L   
	0.0 -.5 0.0     L   
	0.0 0.0 -.5     L   
	0.5 0.5 0.0     X   
	0.5 0.0 0.5     X   
	0.0 0.5 0.5     X   
	-0.5 -0.5 0.0     X   
	-0.5 0.0 -0.5     X   
	0.0 -0.5 -0.5     X   
	0.25000  0.50000  0.75000      W
	-0.25000  0.25000  0.50000      W
	-0.25000  0.50000  0.25000      W
	0.25000  0.75000  0.50000      W
	-0.50000  0.25000 -0.25000      W
	-0.75000 -0.25000 -0.50000      W
	-0.50000 -0.25000 -0.75000      W
	-0.25000  0.25000 -0.50000      W
	-0.25000 -0.50000 -0.75000      W
	-0.25000 -0.75000 -0.50000      W
	0.25000 -0.50000 -0.25000      W
	0.25000 -0.25000 -0.50000      W
	0.50000 -0.25000  0.25000      W
	0.25000 -0.25000  0.50000      W
	0.50000  0.25000  0.75000      W
	0.75000  0.25000  0.50000      W
	0.75000  0.50000  0.25000      W
	0.50000  0.75000  0.25000      W
	0.25000  0.50000 -0.25000      W
	0.50000  0.25000 -0.25000      W
	-0.50000 -0.75000 -0.25000      W
	-0.25000 -0.50000  0.25000      W
	-0.50000 -0.25000  0.25000      W
	-0.75000 -0.50000 -0.25000      W
	0.75000  0.37500  0.37500      K
	0.37500  0.00000 -0.37500      K
	0.37500 -0.37500  0.00000      K
	-0.37500 -0.75000 -0.37500      K
	0.00000 -0.37500  0.37500      K
	-0.75000 -0.37500 -0.37500      K
	-0.37500 -0.37500 -0.75000      K
	0.00000  0.37500 -0.37500      K
	-0.37500  0.37500  0.00000      K
	0.37500  0.75000  0.37500      K
	0.37500  0.37500  0.75000      K
	-0.37500  0.00000  0.37500      K
    }


    # lattice type: I cubic (eg. SG 230 Ia-3d) (I == body-centered)
    set kLabels(I-cubic) {
	0.0 0.0 0.0     GAMMA
	.5 .5 -.5        H  
	.5 -.5 .5        H  
	-.5 .5 .5        H  
	-.5 -.5 .5       H  
	.5 -.5 -.5       H  
	-.5 .5 -.5       H  
	-0.25000  0.75000 -0.25000  P
	-0.75000  0.25000  0.25000  P
	-0.25000 -0.25000 -0.25000  P
	-0.25000 -0.25000  0.75000  P
	0.25000 -0.75000  0.25000  P
	0.75000 -0.25000 -0.25000  P
	0.25000  0.25000  0.25000  P
	0.25000  0.25000 -0.75000  P
	0.5 0.0 0.0      N
	0.0 0.5 0.0      N
	0.0 0.0 0.5      N
	-.5 0.0 0.0      N
	0.0 -.5 0.0      N
	0.0 0.0 -.5      N
	0.50000  0.00000 -0.50000  N
	0.00000  0.50000 -0.50000  N
	-0.50000  0.50000  0.00000  N
	-0.50000  0.00000  0.50000  N
	0.00000 -0.50000  0.50000  N
	0.50000 -0.50000  0.00000  N
    }



    # Hexagonal lattice  (eg. 194  P63/mMc)
    set kLabels(hexagonal) {
	0.0 0.0 0.0     GAMMA
	0.0 0.0 0.5   A  
	0.0 0.0 -.5   A  
	0.5 0.0 0.0   M  
	-.5 0.0 0.0   M  
	0.0 0.5 0.0   M  
	0.0 -.5 0.0   M  
	0.5 -.5 0.0   M  
	-.5 0.5 0.0   M  
	0.5 0.0 0.5   L  
	-.5 0.0 0.5   L  
	0.0 0.5 0.5   L  
	0.0 -.5 0.5   L  
	0.5 -.5 0.5   L  
	-.5 0.5 0.5   L  
	0.5 0.0 -.5   L  
	-.5 0.0 -.5   L  
	0.0 0.5 -.5   L  
	0.0 -.5 -.5   L  
	0.5 -.5 -.5   L  
	-.5 0.5 -.5   L  
	0.333333 0.333333 0.0  K
	-.333333 -.333333 0.0  K
	0.333333 -.666667 0.0  K
	0.666667 -.333333 0.0  K
	-.333333 0.666667 0.0  K
	-.666667 0.333333 0.0  K
	0.333333 0.333333 0.5  H
	-.333333 -.333333 0.5  H
	0.333333 -.666667 0.5  H
	0.666667 -.333333 0.5  H
	-.333333 0.666667 0.5  H
	-.666667 0.333333 0.5  H
	0.333333 0.333333 -.5  H
	-.333333 -.333333 -.5  H
	0.333333 -.666667 -.5  H
	0.666667 -.333333 -.5  H
	-.333333 0.666667 -.5  H
	-.666667 0.333333 -.5  H
    }


    # P tetragonal   (123   P4/mmm)  Coordinates require a=b  ne c ! (not guaranteed)
    set kLabels(P-tetragonal) {
	0.0 0.0 0.0     GAMMA
	0.0 0.0 0.5  Z
	0.0 0.0 -.5  Z
	0.5 0.0 0.0  X
	0.0 0.5 0.0  X
	-.5 0.0 0.0  X
	0.0 -.5 0.0  X
	0.5 0.5 0.0  M
	-.5 0.5 0.0  M
	0.5 -.5 0.0  M
	-.5 -.5 0.0  M
	0.5 0.0 0.5  R
	0.0 0.5 0.5  R
	-.5 0.0 0.5  R
	0.0 -.5 0.5  R
	0.5 0.0 -.5  R
	0.0 0.5 -.5  R
	-.5 0.0 -.5  R
	0.0 -.5 -.5  R
	0.5 0.5 0.5  A
	-.5 0.5 0.5  A
	0.5 -.5 0.5  A
	-.5 -.5 0.5  A
	0.5 0.5 -.5  A
	-.5 0.5 -.5  A
	0.5 -.5 -.5  A
	-.5 -.5 -.5  A
    }


    # I-tetragonal,  (c is "tetragonal")  (139 I4/mmm)  (I == body-centered)
    set kLabels(I-tetragonal) {
	0.0 0.0 0.0     GAMMA
	0.0 0.0 0.5  X
	0.0 0.0 -.5  X
	0.5 -.5 0.0  X
	-.5 0.5 0.0  X
	0.5 0.0 0.0  N
	0.0 0.5 0.0  N
	-.5 0.0 0.0  N
	0.0 -.5 0.0  N
	0.5 0.0 -.5  N
	0.0 0.5 -.5  N
	-.5 0.0 0.5  N
	0.0 -.5 0.5  N
	0.25 0.25 0.25  P
	-.25 0.75 -.25  P
	0.75 -.25 -.25  P
	0.25 0.25 -.75  P
	-.25 -.25 -.25  P
	0.25 -.75 0.25  P
	-.75 0.25 0.25  P
	-.25 -.25 0.75  P
	-.5 0.5 -.5  M  
	-.5 0.5 0.5  M  
	0.5 -.5 0.5  M  
	0.5 -.5 -.5  M  
	0.5 0.5 -.5  Z  
	-.5 -.5 0.5  Z  
    }
    # Note: for above M and Z points "acessibility depends on c/a"



    # P orthorhombic  (47  Pmmm)
    set kLabels(P-orthorhombic) {
	0.0 0.0 0.0     GAMMA
	0.0 0.0 0.5  Z
	0.0 0.0 -.5  Z
	0.5 0.0 0.0  X
	0.0 0.5 0.0  Y
	-.5 0.0 0.0  X
	0.0 -.5 0.0  Y
	0.5 0.5 0.0  S
	-.5 0.5 0.0  S
	0.5 -.5 0.0  S
	-.5 -.5 0.0  S
	0.5 0.0 0.5  U
	0.0 0.5 0.5  T
	-.5 0.0 0.5  U
	0.0 -.5 0.5  T
	0.5 0.0 -.5  U
	0.0 0.5 -.5  T
	-.5 0.0 -.5  U
	0.0 -.5 -.5  T
	0.5 0.5 0.5  R
	-.5 0.5 0.5  R
	0.5 -.5 0.5  R
	-.5 -.5 0.5  R
	0.5 0.5 -.5  R
	-.5 0.5 -.5  R
	0.5 -.5 -.5  R
	-.5 -.5 -.5  R
    }

    # additional comments (PB):
    # --------------------------
    # Trigonal (Rhombohedral) case (2 cases, a,c)
    # 
    # B orthorhombic case (2 cases, depending on a,b,c)
    # 
    # F orthorhombic cases (3 different cases, depending on a,b,c)
    # 
    # C orthorhombic lattice (CXY, and a<>b,  CXZ, CYZ, ...)
    # 
    # C monoclinic lattice (with monoclinic angle gamma, depends on a,b,c and gamma)
    # 
    # P monoclinic lattice (with monoclinic angle gamma, depends on a,b, and gamma)
    # 
    # triclinic
}



