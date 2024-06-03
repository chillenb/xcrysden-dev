#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/genFSInit.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc genFSInit {file} {

    #
    # logo window
    #
    if { [winfo exists .title] } { 
	destroy .title 
    }
    XCRYSDEN_Logo $file

    #
    # read the BXSF file
    #
    set spin [FS_readBXSF $file]

    #
    # select bands for Fermi surface plotting ...
    #
    FS_bandSelect $spin

    #
    # render Fermi surface ...
    #
    FS_GoFermi $spin
}

proc FS_readBXSF {file} {
    global wn fs xcMisc system

    set file [gunzipXSF $file]

    #
    # lets read (band)XSF file; xc_readband return info structure which look
    # like:
    #
    # 1 { grid_index 3D grid_ident grid_nband 1 {subgrid0_ident}}
    set spin {}
    
    #
    # we must construct normal XSF file out of BXSF file and
    # run fsReadBXSF program to construct lattice vectors and Brillouin zone
    #
    if { [catch {exec $system(BINDIR)/fsReadBXSF $file \
		     $system(SCRDIR)/bxsf2xsf.xsf}] } {
	tk_dialog [WidgetName] ERROR "ERROR: an error occured while executing fsReadBXSF program" error 0 OK
	exit
    }
    xc_readXSF $system(SCRDIR)/bxsf2xsf.xsf

    #
    # now read the BXSF file
    #
    set fs(titlefile) $file
    set sinfo [xc_readbandXSF $file]

    if { $sinfo == "" } {
	ErrorDialog "error reading the BXSF file: $file"
	return -code return
    }

    # parse $sinfo
    set slist [lindex $sinfo 1]
    set fs($spin,nbands)        [lindex $slist 3]
    set fs($spin,grid_index)    [lindex $slist 0]
    set fs($spin,grid_subindex) [expr [lindex $slist 4] - 1]
 
    if { ! [info exists fs(Efermi)] } {
	# read the Fermi level from BXSF file if it exists,
	# otherwise set it to zero
	set fs(Efermi) [_FSGetFermiEnergy $file]
    }
    
    #
    # 1. set MIN & MAX values of the band-grids
    # 2. write a band-widths file
    #
    set fs($spin,bandwidthfile) [file join $system(SCRDIR) band_widths.dat]
    set  bwID  [open $fs($spin,bandwidthfile) w]
    puts $bwID $fs($spin,nbands)

    for {set i 1} {$i <= $fs($spin,nbands)} {incr i} {

	
	
	set fs($spin,$i,band_selected) 1
	set fs($spin,$i,minE)     [xc_gridvalue min $fs($spin,grid_index) [expr $i - 1]]
	set fs($spin,$i,maxE)     [xc_gridvalue max $fs($spin,grid_index) [expr $i - 1]]
	set fs($spin,$i,isolevel) $fs(Efermi)    
	
	putsFlush stderr "Reading band: $i    Min-value: $fs($spin,$i,minE)   Max-value: $fs($spin,$i,maxE)"
	puts $bwID [format "%3d    %13.6e %13.6e"  $i $fs($spin,$i,minE) $fs($spin,$i,maxE)]
    }
    close $bwID

    return $spin
}

proc _FSGetFermiEnergy {file} {

    putsFlush stderr "Querying for fermi energy ..."

    set fID [open $file r]
    while { 1 } {
	set line [gets $fID]
	#putsFlush stderr "Line: $line"

	if { [string match -nocase *FERMI*ENERGY* $line] } {
	    set fermi [lindex $line end]
	    break
	}
	if { [eof $fID] } {
	    set fermi     0.0
	    break
	}
    }
    close $fID

    putsFlush stderr "   FERMI ENERGY set to $fermi"
    if { ! [string is double $fermi] } {
	set fermi 0.0
    }
    return $fermi
}




