#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/isoControl.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc IsoControl_InitVar {} {
    global isoControl prop isosurf

    #---
    # default font for thermometer
    set def_font [font create]
    eval {font configure $def_font} [font actual fixed]
    font configure $def_font -size [expr int([font actual fixed -size] * 1.5)]
    #---

    if { ![info exists isoControl(cpl_basis)] } {
	# this "if" is used in isoRender.tcl as an indication if 
	# IsoControl_InitVar was already called
	set isoControl(cpl_basis) MONOCHROME
    }
    if { ![info exists isoControl(cpl_function)] } {       
	set isoControl(cpl_function) LINEAR
    }

    if { ![info exists isoControl(colorplane)] } {set isoControl(colorplane) 0}
    if { ![info exists isoControl(isoline)] }    {set isoControl(isoline)    0}
    if { ![info exists isoControl(colorplane_lighting)] } {set isoControl(colorplane_lighting) 0}
    if { ![info exists isoControl(cpl_transparency)] } {set isoControl(cpl_transparency) 0}
    if { ![info exists isoControl(cpl_thermometer)] } {set isoControl(cpl_thermometer) 0}
    if { ![info exists isoControl(cpl_thermoTplw)] } {set isoControl(cpl_thermoTplw) 0}

    if { ![info exists isoControl(cpl_thermoFmt)] } {set isoControl(cpl_thermoFmt) %+8.4f}
    if { ![info exists isoControl(cpl_thermoLabel)] } {set isoControl(cpl_thermoLabel) " Scale:   [encoding convertfrom symbol D] n(r)"}
    if { ![info exists isoControl(cpl_thermoNTics)] } {set isoControl(cpl_thermoNTics) 6}
    if { ![info exists isoControl(cpl_thermoFont)] } {set isoControl(cpl_thermoFont) $def_font}
    if { ![info exists isosurf(2Dexpand)] }    { set isosurf(2Dexpand) none }
    if { ![info exists isosurf(2Dexpand_X)] }  { set isosurf(2Dexpand_X) 1 }
    if { ![info exists isosurf(2Dexpand_Y)] }  { set isosurf(2Dexpand_Y) 1 }
    if { ![info exists isosurf(2Dexpand_Z)] }  { set isosurf(2Dexpand_Z) 1 }

    if { ![info exists isosurf(tessellation_type)] }  { set isosurf(tessellation_type) cubes }
    if { ![info exists isosurf(normals_type)] }       { set isosurf(normals_type) gradient }

    if { ![info exists isoControl(anim_step)] } {
	set isoControl(anim_step) 1
    }
    if { ![info exists isoControl(time_delay)] } {
	set isoControl(time_delay) 100
    }
    if { ![info exists isoControl(cbfn_apply_to_all)] } {
	set isoControl(cbfn_apply_to_all) 0
    }
    if { ![info exists isoControl(disp_apply_to_all)] } {
	set isoControl(disp_apply_to_all) 0
    }
    if { ![info exists isoControl(anim_apply_to_all)] } {
	set isoControl(anim_apply_to_all) 0
    }
    #
    # just in case
    #
    if { ![info exists isoControl(3Dinterpl_degree)] } {
	set isoControl(3Dinterpl_degree) 1
    }

    set n $isoControl(3Dinterpl_degree)    
    if { ![info exists isoControl(1,nslide)] } {
	set isoControl(1,nslide) \
		[expr ([lindex [xc_iso grid] 2] - 1) * $n + 1]
    }
    if { ![info exists isoControl(2,nslide)] } {
	set isoControl(2,nslide) \
		[expr ([lindex [xc_iso grid] 1] - 1) * $n + 1]
    }
    if { ![info exists isoControl(3,nslide)] } {
	set isoControl(3,nslide) \
		[expr ([lindex [xc_iso grid] 0] - 1) * $n + 1]
    }

    if { ![info exists isoControl(current_slide)] } {
	set isoControl(current_slide) 1
	set isoControl(current_text_slide) "Current slide:  $isoControl(current_slide) / $isoControl(1,nslide)"
    }

    if { ![info exists isoControl(2Dlowvalue)] } {
	set isoControl(2Dlowvalue) $isosurf(minvalue)
    }
    if { ![info exists isoControl(2Dhighvalue)] } {
	set isoControl(2Dhighvalue) $isosurf(maxvalue)
    }
    if { ![info exists isoControl(2Dnisoline)] } {
	set isoControl(2Dnisoline) 15
    }
    if { ![info exists isoControl(isoline_color)] } {
	set isoControl(isoline_color) monocolor 
    }
    if { ![info exists isoControl(isoline_width)] } {
	set isoControl(isoline_width) 2 
    }
    if { ![info exists isoControl(isoline_monocolor)] } {
	set isoControl(isoline_monocolor) #000000
    }
    if { ![info exists isoControl(isoline_stipple)] } {
	set isoControl(isoline_stipple) {no stipple}
    }
   
    foreach i {1 2 3} {	
	if { ![info exists isoControl($i,cpl_basis)] } {
	    set isoControl($i,cpl_basis)    $isoControl(cpl_basis) 
	}
	if { ![info exists isoControl($i,cpl_function)] } {       
	    set isoControl($i,cpl_function) $isoControl(cpl_function)
	}

	if { ![info exists isoControl($i,colorplane)] } {
	    set isoControl($i,colorplane) $isoControl(colorplane)
	}
	if { ![info exists isoControl($i,isoline)] } {
	    set isoControl($i,isoline)    $isoControl(isoline)
	}
	if { ![info exists isoControl($i,colorplane_lighting)] } {
	    set isoControl($i,colorplane_lighting) $isoControl(colorplane_lighting)
	}
	if { ![info exists isoControl($i,cpl_transparency)] } {
	    set isoControl($i,cpl_transparency) $isoControl(cpl_transparency)
	}

	if { ![info exists isoControl($i,cpl_thermometer)] } {
	    set isoControl($i,cpl_thermometer) $isoControl(cpl_thermometer)
	}
	if { ![info exists isoControl($i,cpl_thermoTplw)] } {
	    set isoControl($i,cpl_thermoTplw) $isoControl(cpl_thermoTplw)
	}
	if { ![info exists isoControl($i,cpl_thermoFmt)] } {
	    set isoControl($i,cpl_thermoFmt) $isoControl(cpl_thermoFmt)
	}
	if { ![info exists isoControl($i,cpl_thermoLabel)] } {
	    set isoControl($i,cpl_thermoLabel) $isoControl(cpl_thermoLabel) 
	}
	if { ![info exists isoControl($i,cpl_thermoNTics)] } {
	    set isoControl($i,cpl_thermoNTics) $isoControl(cpl_thermoNTics)
	}
	if { ![info exists isoControl($i,cpl_thermoFont)] } {
	    set isoControl($i,cpl_thermoFont) $isoControl(cpl_thermoFont)
	}

	if { ![info exists isosurf($i,2Dexpand)] }    { 
	    set isosurf($i,2Dexpand)   $isosurf(2Dexpand) 
	}
	if { ![info exists isosurf($i,2Dexpand_X)] }  { 
	    set isosurf($i,2Dexpand_X) $isosurf(2Dexpand_X) 
	}
	if { ![info exists isosurf($i,2Dexpand_Y)] }  { 
	    set isosurf($i,2Dexpand_Y) $isosurf(2Dexpand_Y)
	}
	if { ![info exists isosurf($i,2Dexpand_Z)] }  { 
	    set isosurf($i,2Dexpand_Z) $isosurf(2Dexpand_Z)
	}

	if { ![info exists isoControl($i,anim_step)] } {
	    set isoControl($i,anim_step) $isoControl(anim_step)
	}
	if { ![info exists isoControl($i,time_delay)] } {
	    set isoControl($i,time_delay) $isoControl(time_delay)
	}
	if { ![info exists isoControl($i,current_slide)] } {
	    set isoControl($i,current_slide) $isoControl(current_slide)
	}

	if { ![info exists isoControl($i,2Dlowvalue)] } {
	    set isoControl($i,2Dlowvalue) $isoControl(2Dlowvalue) 
	}
	if { ![info exists isoControl($i,2Dhighvalue)] } {
	    set isoControl($i,2Dhighvalue) $isoControl(2Dhighvalue)
	}
	if { ![info exists isoControl($i,2Dnisoline)] } {
	    set isoControl($i,2Dnisoline) $isoControl(2Dnisoline) 
	}
	if { ![info exists isoControl($i,isoline_color)] } {
	    set isoControl($i,isoline_color) $isoControl(isoline_color)
	}
	if { ![info exists isoControl($i,isoline_width)] } {
	    set isoControl($i,isoline_width) $isoControl(isoline_width)
	}
	if { ![info exists isoControl($i,isoline_monocolor)] } {
	    set isoControl($i,isoline_monocolor) $isoControl(isoline_monocolor)
	}
	if { ![info exists isoControl($i,isoline_stipple)] } {
	    set isoControl($i,isoline_stipple) $isoControl(isoline_stipple)
	}
    }

    # check if close-isocontrol button was pressed previously

    if { [info exists isoControl(close,isosurf) ]    } { set isoControl(isosurf)    $isoControl(close,isosurf) }
    if { [info exists isoControl(close,colorplane) ] } { set isoControl(colorplane) $isoControl(close,colorplane) }
    if { [info exists isoControl(close,isoline) ]    } { set isoControl(isoline)    $isoControl(close,isoline) }
    foreach i {1 2 3} {
	if { [info exists isoControl(close,$i,colorplane)] } { set isoControl($i,colorplane) $isoControl(close,$i,colorplane) }
	if { [info exists isoControl(close,$i,isoline)]    } { set isoControl($i,isoline)    $isoControl(close,$i,isoline) }
    }
}



proc IsoControl_SetColorPlaneVar {type i} {
    global isoControl prop isosurf

    xcDebug "#1 Animation: def = $isoControl(current_slide)"
    xcDebug "#1 Animation:   1 = $isoControl(1,current_slide)"
    xcDebug "#1 Animation:   2 = $isoControl(2,current_slide)"
    xcDebug "#1 Animation:   3 = $isoControl(3,current_slide)"
    if { $type == "1st" } {
	set isoControl($i,cpl_basis)           $isoControl(cpl_basis) 
	set isoControl($i,cpl_function)        $isoControl(cpl_function)
	set isoControl($i,colorplane)          $isoControl(colorplane)
	set isoControl($i,isoline)             $isoControl(isoline)
	set isoControl($i,colorplane_lighting) $isoControl(colorplane_lighting)
	set isoControl($i,cpl_transparency)    $isoControl(cpl_transparency)
	set isoControl($i,cpl_thermometer)     $isoControl(cpl_thermometer)
	set isoControl($i,cpl_thermoTplw)      $isoControl(cpl_thermoTplw)
	set isoControl($i,cpl_thermoFmt)       $isoControl(cpl_thermoFmt)
	set isoControl($i,cpl_thermoLabel)     $isoControl(cpl_thermoLabel) 
	set isoControl($i,cpl_thermoNTics)     $isoControl(cpl_thermoNTics)
	set isoControl($i,cpl_thermoFont)      $isoControl(cpl_thermoFont)
	set isosurf($i,2Dexpand)               $isosurf(2Dexpand) 
	set isosurf($i,2Dexpand_X)             $isosurf(2Dexpand_X) 
	set isosurf($i,2Dexpand_Y)             $isosurf(2Dexpand_Y)
	set isosurf($i,2Dexpand_Z)             $isosurf(2Dexpand_Z)
	set isoControl($i,anim_step)           $isoControl(anim_step)
	#set isoControl($i,current_slide)       $isoControl(current_slide)
	set isoControl($i,2Dlowvalue)          $isoControl(2Dlowvalue)       
	set isoControl($i,2Dhighvalue)         $isoControl(2Dhighvalue)
	set isoControl($i,2Dnisoline)          $isoControl(2Dnisoline)      
	set isoControl($i,isoline_color)       $isoControl(isoline_color)      
	set isoControl($i,isoline_width)       $isoControl(isoline_width)      
	set isoControl($i,isoline_monocolor)   $isoControl(isoline_monocolor) 
	set isoControl($i,isoline_stipple)     $isoControl(isoline_stipple)   
    } else {
	set isoControl(cpl_basis)           $isoControl($i,cpl_basis) 
	set isoControl(cpl_function)        $isoControl($i,cpl_function)
	set isoControl(colorplane)          $isoControl($i,colorplane)
	set isoControl(isoline)             $isoControl($i,isoline)
	set isoControl(colorplane_lighting) $isoControl($i,colorplane_lighting)
	set isoControl(cpl_transparency)    $isoControl($i,cpl_transparency)
	set isoControl(cpl_thermometer)     $isoControl($i,cpl_thermometer) 
	set isoControl(cpl_thermoTplw)      $isoControl($i,cpl_thermoTplw) 
	set isoControl(cpl_thermoFmt)       $isoControl($i,cpl_thermoFmt)   
	set isoControl(cpl_thermoLabel)     $isoControl($i,cpl_thermoLabel) 
	set isoControl(cpl_thermoNTics)     $isoControl($i,cpl_thermoNTics) 
	set isoControl(cpl_thermoFont)      $isoControl($i,cpl_thermoFont)  
	set isosurf(2Dexpand)               $isosurf($i,2Dexpand) 
	set isosurf(2Dexpand_X)             $isosurf($i,2Dexpand_X) 
	set isosurf(2Dexpand_Y)             $isosurf($i,2Dexpand_Y)
	set isosurf(2Dexpand_Z)             $isosurf($i,2Dexpand_Z)
	set isoControl(anim_step)           $isoControl($i,anim_step)
	set isoControl(time_delay)          $isoControl($i,time_delay)

	#set isoControl(current_slide)       $isoControl($i,current_slide)
	set isoControl(current_text_slide)  "Current slide:  $isoControl($i,current_slide) / $isoControl($i,nslide)"
	set isoControl(2Dlowvalue)          $isoControl($i,2Dlowvalue)       
	set isoControl(2Dhighvalue)         $isoControl($i,2Dhighvalue)
	set isoControl(2Dnisoline)          $isoControl($i,2Dnisoline)      
	set isoControl(isoline_color)       $isoControl($i,isoline_color)      
	set isoControl(isoline_width)       $isoControl($i,isoline_width)      
	set isoControl(isoline_monocolor)   $isoControl($i,isoline_monocolor) 
	set isoControl(isoline_stipple)     $isoControl($i,isoline_stipple)   
    }
    xcDebug "#2 Animation: def = $isoControl(current_slide)"
    xcDebug "#2 Animation:   1 = $isoControl(1,current_slide)"
    xcDebug "#2 Animation:   2 = $isoControl(2,current_slide)"
    xcDebug "#2 Animation:   3 = $isoControl(3,current_slide)"
    xcDebug "#2 isoControl(plane) = $isoControl(plane)"
}



proc IsoControl_UpdateColorplane {} {
    global isoControl isosurf

    # COLOR_BASIS & SCALE_FUNCTION
    set item $isoControl(plane)
    if $isoControl(cbfn_apply_to_all) {
	set item {1 2 3}
    }
    foreach it $item {
	set isoControl($it,cpl_basis)          $isoControl(cpl_basis)
	set isoControl($it,cpl_function)       $isoControl(cpl_function)
    }

    # DISPLAY/RANGE/EXPAND/ISOLINE
    set item $isoControl(plane)
    if $isoControl(disp_apply_to_all) {
	set item {1 2 3}
    }
    foreach it $item {
	# check the following three vaiables
	if { ![check_var {
	    {isoControl(2Dlowvalue)    real} 
	    {isoControl(2Dhighvalue)   real} 
	    {isoControl(2Dnisoline)    posint}
	    {isoControl(isoline_width) posint}
	} [list \
	       $isoControl(2Dlowvalue_entry) \
	       $isoControl(2Dhighvalue_entry) \
	       $isoControl(2Dnisoline_entry) \
	       $isoControl(2Disolinewidth_entry)]] } {
	    return	    
	}
	set isoControl($it,colorplane)	        $isoControl(colorplane)
	set isoControl($it,isoline)	        $isoControl(isoline)
	set isoControl($it,colorplane_lighting)	$isoControl(colorplane_lighting)    
	set isoControl($it,cpl_transparency)    $isoControl(cpl_transparency)
	set isoControl($it,cpl_thermometer)     $isoControl(cpl_thermometer)
	set isoControl($it,cpl_thermoTplw)      $isoControl(cpl_thermoTplw)
	set isoControl($it,cpl_thermoFmt)       $isoControl(cpl_thermoFmt)
	set isoControl($it,cpl_thermoLabel)     $isoControl(cpl_thermoLabel) 
	set isoControl($it,cpl_thermoNTics)     $isoControl(cpl_thermoNTics)
	set isoControl($it,cpl_thermoFont)      $isoControl(cpl_thermoFont)

	if { $isoControl(2Dnisoline) > $isoControl(max_allowed_2Dnisoline) } {
	    tk_dialog [WidgetName] WARNING "WARNING: more then $isoControl(max_allowed_2Dnisoline) isolines was requested; Maximum number is $isoControl(max_allowed_2Dnisoline) !!! Setting number of isolines to $isoControl(max_allowed_2Dnisoline)" \
		    warning 0 OK
	    set isoControl(2Dnisoline) 50
	}
	set isoControl($it,2Dlowvalue)         $isoControl(2Dlowvalue)       
	set isoControl($it,2Dhighvalue)        $isoControl(2Dhighvalue)
	set isoControl($it,2Dnisoline)         $isoControl(2Dnisoline)      
	set isoControl($it,isoline_color)      $isoControl(isoline_color)      
	set isoControl($it,isoline_width)      $isoControl(isoline_width)      
	set isoControl($it,isoline_monocolor)  $isoControl(isoline_monocolor) 
	set isoControl($it,isoline_stipple)    $isoControl(isoline_stipple)   
	set isosurf($it,2Dexpand)       $isosurf(2Dexpand)    
	set isosurf($it,2Dexpand_X)	$isosurf(2Dexpand_X)
	set isosurf($it,2Dexpand_Y)	$isosurf(2Dexpand_Y)
	set isosurf($it,2Dexpand_Z)	$isosurf(2Dexpand_Z)
    }
    UpdateIsosurf
}



proc IsoControl_Show {type fiso fplane biso bplane1 bplane2 bplane3} {
    global isoControl
    
    if { $isoControl(plane) != {} } {
	IsoControl_SetColorPlaneVar 1st $isoControl(plane)
    }

    if { $type == "isosurf" } {	
	set isoControl(plane) {}
	$biso   config -bd 3
	$bplane1 config -bd 1
	$bplane2 config -bd 1
	$bplane3 config -bd 1
	pack forget $fplane
	pack $fiso -side top -fill both -expand 1
    } else {
	if { $type == "colorplane1" } {
	    set isoControl(plane) 1
	    $biso   config -bd 1
	    $bplane1 config -bd 3
	    $bplane2 config -bd 1
	    $bplane3 config -bd 1
	    pack forget $fiso
	    pack $fplane -side top -fill both -expand 1
	} elseif { $type == "colorplane2" } {
	    set isoControl(plane) 2
	    $biso    config -bd 1
	    $bplane1 config -bd 1
	    $bplane2 config -bd 3
	    $bplane3 config -bd 1
	    pack forget $fiso
	    pack $fplane -side top -fill both -expand 1
	} elseif { $type == "colorplane3" } {
	    set isoControl(plane) 3
	    $biso    config -bd 1
	    $bplane1 config -bd 1
	    $bplane2 config -bd 1
	    $bplane3 config -bd 3
	    pack forget $fiso
	    pack $fplane -side top -fill both -expand 1
	}
	IsoControl_SetColorPlaneVar 2nd $isoControl(plane)
    }
}



proc IsoControl_IsoLineShow {type d r e i df rf ef ifr} {
    global isoControl

    set dr 1
    set rr 1
    set er 1
    set ir 1
    if { $type == "display" } {
	set dr 3
	pack $df -padx 3 -pady 3 -ipady 3 -side top -fill x -expand 1
	pack forget $rf $ef $ifr
    } elseif { $type == "ranges" } {
	set rr 3
	pack $rf -padx 3 -pady 3 -ipady 3 -side top -fill x -expand 1
	pack forget $df $ef $ifr
    } elseif { $type == "expand" } {
	set er 3
	pack $ef -padx 3 -pady 3 -ipady 3 -side top -fill x -expand 1
	pack forget $df $rf $ifr
    } elseif { $type == "isoline" } {
	set ir 3
	pack $ifr -side top -fill x -expand 1
	pack forget $df $rf $ef
    }
    $d config -bd $dr
    $r config -bd $rr
    $e config -bd $er
    $i config -bd $ir
}


proc IsoControl_Hide {{t .iso}} {
    global unmapWin
    wm withdraw $t
    xcUnmapWindow unmap $t $t $unmapWin(frame,main) isosurf_control
}



# isoControl(color_button)
# isoControl(blend_button)
proc IsoControl {} {
    global isosurf nxdir nydir nzdir periodic xcFonts XCTrace prop \
	    isoControl unmapWin

    if {[winfo exists .iso] } { return }
    
    #
    # initializations
    #
    set isoControl(plane) {}
    set isoControl(datagridDim) 3

    IsoControl_InitVar
    reloadButtonDisable
    
    set t [xcToplevel .iso "Isosurface/Property-plane Controls" "IsoControls" . -0 0 1]
    pack propagate .iso

    xcRegisterUnmapWindow $t $unmapWin(frame,main) \
	    isosurf_control -image unmap-isosurf
    bind $t <Unmap> [list xcUnmapWindow unmap %W $t \
	    $unmapWin(frame,main) isosurf_control]
    bind $t <Map>   [list xcUnmapWindow map %W $t \
	    $unmapWin(frame,main) isosurf_control]
    set ft  [frame $t.ft]
    set fb1 [frame $t.fb1]
    set fb2 [frame $t.fb2]
    if { $prop(type_of_run) == "UHF" } {    
	set fm [frame $t.fm]
	pack $ft $fm -fill x
    } else {
	pack $ft -fill x
    }
    pack $fb1 -fill both -expand 1
    
    set b1 [button $ft.b1 \
	    -text "Isosurface" \
	    -bd 3 \
	    -highlightthickness 0 \
	    -command [list IsoControl_Show isosurf $fb1 $fb2 \
	    $ft.b1 $ft.b2 $ft.b3 $ft.b4]]
    set b2 [button $ft.b2 \
	    -text "Plane #1" \
	    -bd 1 \
	    -highlightthickness 0 \
	    -command [list IsoControl_Show colorplane1 $fb1 $fb2 \
	    $ft.b1 $ft.b2 $ft.b3 $ft.b4]]
    set b3 [button $ft.b3 \
	    -text "Plane #2" \
	    -bd 1 \
	    -highlightthickness 0 \
	    -command [list IsoControl_Show colorplane2 $fb1 $fb2 \
	    $ft.b1 $ft.b2 $ft.b3 $ft.b4]]
    set b4 [button $ft.b4 \
	    -text "Plane #3" \
	    -bd 1 \
	    -highlightthickness 0 \
	    -command [list IsoControl_Show colorplane3 $fb1 $fb2 \
	    $ft.b1 $ft.b2 $ft.b3 $ft.b4]]
    pack $b1 $b2 $b3 $b4 -side left -fill both -expand 1
    
    ########################################
    # if UHF
    if { $prop(type_of_run) == "UHF" } {
	set f0 [frame $fm.f0 -relief raised -bd 2]
	frame $f0.1 -relief groove -bd 2
	xcMenuEntry $f0.1 "What SPIN to take:" 30  \
		isosurf(spin) {ALPHA BETA ALPHA+BETA ALPHA-BETA} \
		-labelwidth 17 \
		-labelfont $xcFonts(small) \
		-entryfont $xcFonts(small_entry)
	pack $f0 -side top -fill both
	pack $f0.1 -padx 2 -pady 5 -ipady 3 -fill x
    }

    ########################
    ### IsoSURFACE frame ###
    ########################
    set f     [frame $fb1.f]
    set left  [frame $f.left  -relief raised -bd 2]
    set right [frame $f.right -relief raised -bd 2]
    set bot   [frame $fb1.bot   -relief raised -bd 2]
    pack $f -side top -expand 1 -fill both
    pack $left $right -side left -fill both -expand 1
    pack $bot -side bottom -fill both -expand 1  

    ########################################
    # LEFT FRAME
    if { ![info exists isoControl(isosurf)] } {
	set isoControl(isosurf) 1
    }
    set ckb [checkbutton $left.cb \
	    -text "Display Isosurface" \
	    -command "UpdateIsosurf" \
	    -relief raised -bd 2 \
	    -anchor w \
	    -variable isoControl(isosurf)]
    pack $ckb -side top -padx 5 -pady 5 -fill x

    scale $left.sc -from 1 -to 4 -length 100 \
	    -variable isosurf(3Dinterpl_degree) -orient horizontal \
	    -label "Degree of triCubic Spline:" \
	    -tickinterval 1 -resolution 1 \
	    -width 7 -sliderlength 20 \
	    -showvalue true -relief groove -bd 2 \
	    -font $xcFonts(small)
    pack $left.sc -side top -padx 2 -pady 5 -fill both -expand 1
    
    set f1 [frame $left.f1 -relief groove -bd 2]
    set f2 [frame $left.f2 -relief groove -bd 2]
    set f3 [frame $left.f3 -relief groove -bd 2]
    set f4 [frame $left.f4 -relief groove -bd 2]
    pack $f3 $f4 $f1 $f2 -side top -padx 2 -pady 5 -fill both -expand 1

    set frame_color [lindex \
	    [GetWidgetConfig frame -background] end]
    FillEntries $f1 {"Minimum grid value:"} \
	    isosurf(minvalue) \
	    18 9 top left \
	    -e_relief flat -e_state disabled -e_bg $frame_color \
	    -l_font $xcFonts(small) -e_font $xcFonts(small_entry)
    FillEntries $f1 {"Maximum grid value:"} \
	    isosurf(maxvalue) \
	    18 9 top left \
	    -e_relief flat -e_state disabled -e_bg $frame_color \
	    -l_font $xcFonts(small) -e_font $xcFonts(small_entry)
    set isosurf(isovalue_entry) [FillEntries $f1 {"Isovalue:"} \
	    prop(isolevel) \
	    18 9 top left \
	    -l_font $xcFonts(small) -e_font $xcFonts(small_entry)]
    focus $isosurf(isovalue_entry)
    #
    # make checkbutton for specifying +/- option
    #
    set ck [checkbutton $f1.ckb \
		-text "Render +/- isovalue" \
		-variable prop(pm_isolevel) \
		-command IsoControlCommand \
		-anchor w]
    pack $ck -side bottom -expand 1 -fill both -padx 5 -pady 2

    #
    # isosurface triangulation algorithm
    #
    set wlist_ [RadioButCmd $f3 "Isosurface tessellation type:" \
		    isosurf(tessellation_type) UpdateIsosurf \
		    top left 0 1 2 "cubes" "tetrahedrons"]
    foreach w $wlist_ {
	$w configure -font $xcFonts(small)
    }

    #
    # isosurface normals computation algorithm
    #
    set wlist_ [RadioButCmd $f4 "Isosurface normals type:" \
		    isosurf(normals_type) UpdateIsosurf \
		    top left 0 1 2 "gradient" "triangles"]
    foreach w $wlist_ {
	$w configure -font $xcFonts(small)
    }

    ###############
    # EXPAND option; just for periodic systems    

    set disabled_color [lindex \
	    [GetWidgetConfig button -disabledforeground] end]
    set enabled_color  [lindex \
	    [GetWidgetConfig scale -foreground] end]
    set enable "\
	    $f2.b.scX configure -foreground $enabled_color; \ 
            $f2.b.scY configure -foreground $enabled_color; \ 
            $f2.b.scZ configure -foreground $enabled_color"
    set disable " \
	    $f2.b.scX configure -foreground $disabled_color; \ 
            $f2.b.scY configure -foreground $disabled_color; \ 
            $f2.b.scZ configure -foreground $disabled_color"

    if { $periodic(dim) > 0 } {
	label $f2.l -text "Expand Isosurface:" -relief flat
	set lfont [ModifyFont [$f2.l cget -font] $f2.l -underline 1]
	$f2.l configure -font $lfont
	set r1 [radiobutton $f2.r1 \
		-text "do not expand" \
		-variable isosurf(expand) \
		-value "none" \
		-anchor w \
		-font $xcFonts(small) \
		-command "xcDisableAll $f2.b; catch {eval $disable}"]
	set r2 [radiobutton $f2.r2 \
		-text "to whole structure" \
		-variable isosurf(expand) \
		-value "whole" \
		-anchor w \
		-font $xcFonts(small) \
		-command "xcDisableAll $f2.b; catch {eval $disable}"]
	set r3 [radiobutton $f2.r3 \
		-text "separately in each direction" \
		-variable isosurf(expand) \
		-value "specify" \
		-anchor w \
		-font $xcFonts(small) \
		-command "xcEnableAll $f2.b; catch {eval $enable}"]
	pack $f2.l -side top -expand 1
	pack $r1 $r2 $r3 -side top -fill x

	set f2a [frame $f2.a -relief flat -width 40]
	set f2b [frame $f2.b -relief flat]
	pack $f2a -side left
	pack $f2b -side left -expand 1 -fill x

	set XCTrace(scX) [scale $f2b.scX -from 1 -to $nxdir -length 100 \
		-variable isosurf(expand_X) -orient horizontal \
		-label "repeat in X-dir:" -tickinterval 1 -resolution 1 \
		-showvalue true \
		-font $xcFonts(small) \
		-width 7 -sliderlength 20]
	# TRACE -- nxdir
	trace variable nxdir w xcTrace
	pack $XCTrace(scX) -side top
	if { $periodic(dim) > 1} {
	    set XCTrace(scY) [scale $f2b.scY -from 1 -to $nydir -length 100 \
		    -variable isosurf(expand_Y) -orient horizontal \
		    -label "repeat in Y-dir:" -tickinterval 1 -resolution 1 \
		    -showvalue true \
		    -font $xcFonts(small) \
		    -width 7 -sliderlength 20]
	    # TRACE -- nydir
	    trace variable nydir w xcTrace
	    pack $XCTrace(scY) -side top	
	}
	if { $periodic(dim) > 2 } {
	    set XCTrace(scZ) [scale $f2b.scZ -from 1 -to $nzdir -length 100 \
		    -variable isosurf(expand_Z) -orient horizontal \
		    -label "repeat in Z-dir:" -tickinterval 1 -resolution 1 \
		    -showvalue true \
		    -font $xcFonts(small) \
		    -width 7 -sliderlength 20]
	    # TRACE -- nzdir
	    trace variable nzdir w xcTrace
	    pack $XCTrace(scZ) -side top
	}
    }    
	
    ########################################
    # RIGHT FRAME
    set f1 [frame $right.f1 -relief groove -bd 2]
    set f2 [frame $right.f2 -relief groove -bd 2]
    set f3 [frame $right.f3 -relief groove -bd 2]
    set f4 [frame $right.f4 -relief groove -bd 2]
    set f5 [frame $right.f5]
    pack $f1 $f2 $f3 $f4 -side top -padx 2 -pady 5 -fill x
    pack $f5  -side top -padx 2 -pady 5 -fill both -expand 1
    
    set wlist1 [RadioButCmd $f1 "Render isosurface as:" \
	    isosurf(type_of_isosurf) UpdateIsosurf \
	    top left 0 1 2 "solid" "wire" "dot"]
    
    set wlist2 [RadioButCmd $f2 "Isosurface's ShadeModel:" \
	    isosurf(shade_model) UpdateIsosurf \
	    top left 0 1 2 "smooth" "flat"]
    
    set wlist3 [RadioButCmd $f3 "Two-sided lighting:" \
	    isosurf(twoside_lighting) ConvertTwoSideVar \
	    top left 0 1 2 "off" "on"]

    set wlist4 [RadioButCmd $f4 "Transparency of isosurface:" \
	    isosurf(transparency) UpdateIsosurf \
	    top left 0 1 2 "off" "on"]

    foreach w [concat $wlist1 $wlist2 $wlist3 $wlist4] {
	$w configure -font $xcFonts(small)
    }

    # REVERT FRONT&BACK SIDE    
    button $f5.b10 \
	    -font $xcFonts(small) \
	    -text "Revert (+) sides" \
	    -command [list RevertIsoSides pos]    
    set isoControl(revert_button2) [button $f5.b10a \
	    -font $xcFonts(small) \
	    -text "Revert (-) sides" \
	    -command [list RevertIsoSides neg]]

    # REVERT NORMALS
    button $f5.b11 \
	    -font $xcFonts(small) \
	    -text "Revert (+) normals" \
	    -command [list RevertIsoNormals pos]    
    set isoControl(revert_button1) [button $f5.b12 \
	    -font $xcFonts(small) \
	    -text "Revert (-) normals" \
	    -command [list RevertIsoNormals neg]]

    set isoControl(smooth_button) [button $f5.b2a \
	    -text "Surface smoothing" \
	    -font $xcFonts(small) \
	    -command SurfaceSmoothing]
    $isoControl(smooth_button) config -state disabled

    set isoControl(color_button) [button $f5.b2 \
	    -text "Set COLOR parameters" \
	    -font $xcFonts(small)]	    

    set isoControl(blend_button) [button $f5.b3 \
	    -text "Set TRANSPARENCY\nparameters" \
	    -font $xcFonts(small)]

    pack $f5.b10 $f5.b10a $f5.b11 $f5.b12 $f5.b2a $f5.b2 $f5.b3 \
	    -fill x -side top -padx 2 -pady 5

    set hid [button $bot.hid -text "Hide" \
	    -command [list IsoControl_Hide $t]]
    set can [button $bot.can -text "Close" \
	    -command [list IsoControlCan $t]]
    set sav [button $bot.sav -text "Save Grid" \
	    -command IsoControlSave]
    set sub [button $bot.sub -text "Submit" \
	    -command UpdateIsosurf]
    pack $hid $can $sav $sub -side left -pady 5 -expand 1
    
    #
    # get the state according to prop(pm_isolevel) variable
    # (look the "Render +/- isovalue" checkbutton)
    #
    IsoControlCommand


    ########################
    ### ColorPLANE frame ###
    ########################
    set f0 [frame $fb2.f0 -relief raised -bd 2]
    frame $f0.1 -relief groove -bd 2
    set ckb [checkbutton $f0.1.cb \
	    -text "Apply to Planes #1/#2/#3" \
	    -relief raised -bd 2 \
	    -anchor w \
	    -variable isoControl(cbfn_apply_to_all)]
    pack $ckb -side top -padx 5 -pady 5 -fill x

    set f1 [frame $f0.f1]
    set f2 [frame $f0.f2]
    xcMenuEntry $f1 "Select color basis:" 30  \
	    isoControl(cpl_basis) {MONOCHROME RAINBOW RGB GEOGRAPHIC BLUE-WHITE-RED BLACK-BROWN-WHITE} \
	    -labelwidth 17 \
	    -labelfont $xcFonts(small) \
	    -entryfont $xcFonts(small_entry) \
	    -labelanchor w
    xcMenuEntry $f2 "Select scale function:" 30  \
	    isoControl(cpl_function) {LINEAR LOG10 SQRT 3th-ROOT EXP(x) EXP(x^2)} \
	    -labelwidth 17 \
	    -labelfont $xcFonts(small) \
	    -entryfont $xcFonts(small_entry) \
	    -labelanchor w
    pack $f0 -side top -fill both -expand 1
    pack $f0.1 -padx 2 -pady 5 -ipady 3 -fill x
    pack $f1 $f2 -side top -in $f0.1 -fill x
    
    #
    # DISPLAY/RANGES/EXPAND/ISOLINE
    #    
    set f1  [frame $fb2.f1 -relief raised -bd 2]
    set ckb [checkbutton $f1.cb \
	    -text "Apply to Planes #1/#2/#3" \
	    -relief raised -bd 2 \
	    -anchor w \
	    -variable isoControl(disp_apply_to_all)]
    pack $ckb -side top -padx 5 -pady 5 -fill x

    set f11 [frame $f1.1]
    pack $f1 -side top -fill both -expand 1 
    pack $f11 -side top -fill both -expand 1 -padx 4 -pady 8
    set mfb [frame $f11.mfb]
    set mf  [frame $f11.mf -relief raised -bd 2]
    pack $mfb $mf -side top -padx 2 -fill both -expand 1
    
    set df  [frame $mf.d -relief groove -bd 2]
    set rf  [frame $mf.r -relief groove -bd 2]
    set ef  [frame $mf.e -relief groove -bd 2]
    set ifr [frame $mf.i]

    set d [button $mfb.d -text "Display" \
	    -highlightthickness 0 -bd 1 \
	    -command [list IsoControl_IsoLineShow display \
	    $mfb.d $mfb.r $mfb.e $mfb.i \
	    $df $rf $ef $ifr]]
    set r [button $mfb.r -text "Ranges" \
	    -highlightthickness 0 -bd 1 \
	    -command [list IsoControl_IsoLineShow ranges \
	    $mfb.d $mfb.r $mfb.e $mfb.i \
	    $df $rf $ef $ifr]]
    set e [button $mfb.e -text "Expand" \
	    -highlightthickness 0 -bd 1 \
	    -command [list IsoControl_IsoLineShow expand \
	    $mfb.d $mfb.r $mfb.e $mfb.i\
	    $df $rf $ef $ifr]]
    set i [button $mfb.i -text "Isoline" \
	    -highlightthickness 0 -bd 1 \
	    -command [list IsoControl_IsoLineShow isoline \
	    $mfb.d $mfb.r $mfb.e $mfb.i\
	    $df $rf $ef $ifr]]

    if { $periodic(dim) == 0 } {
	$e config -state disabled
    }
    pack $d $r $e $i -side left -padx 0 -pady 0 -fill both -expand 1

    IsoControl_IsoLineShow display $mfb.d $mfb.r $mfb.e $mfb.i \
	    $df $rf $ef $ifr

    
    #
    # DISPLAY
    #
    set l     [label $df.l -text "Property-plane display option:"]
    set lfont [ModifyFont [$df.l cget -font] $df.l -underline 1]
    $l configure -font $lfont

    set ck1 [checkbutton $df.c1 \
	    -text "display color-plane" \
	    -variable isoControl(colorplane) \
	    -width 21 \
	    -anchor w]
    set ck2 [checkbutton $df.c2 \
	    -text "display isolines" \
	    -variable isoControl(isoline) \
	    -width 21 \
	    -anchor w]
    set ck3 [checkbutton $df.c3 \
		 -text "lighting of color-plane" \
		 -variable isoControl(colorplane_lighting) \
		 -width 21 \
		 -onvalue 1 \
		 -offvalue 0 \
		 -anchor w]
    set ck4 [checkbutton $df.c4 \
	    -text "transparent color-plane" \
	    -variable isoControl(cpl_transparency) \
	    -width 21 \
	    -onvalue 1 -offvalue 0 \
	    -anchor w]
    set ck5 [checkbutton $df.c5 \
		 -text "display thermometer" \
		 -variable isoControl(cpl_thermometer) \
		 -width 21 \
		 -onvalue 1 -offvalue 0 \
		 -anchor w]
    set ck6 [checkbutton $df.c6 \
		 -text "thermometer in toplevel" \
		 -variable isoControl(cpl_thermoTplw) \
		 -width 21 \
		 -onvalue 1 -offvalue 0 \
		 -anchor w]
    set thermo [frame $df.f]
    grid configure $l      -column 0 -row 0 -columnspan 2 
    grid configure $ck1    -column 0 -row 1
    grid configure $ck4    -column 0 -row 2
    grid configure $ck2    -column 1 -row 1
    grid configure $ck3    -column 1 -row 2
    grid configure $ck5    -column 0 -row 3
    grid configure $ck6    -column 1 -row 3
    grid configure $thermo -column 0 -row 4 -columnspan 2

    # thermometer-widgets
    set tf1  [frame $thermo.1 -relief groove -bd 2]
    set tf11 [frame $tf1.1]
    set tf12 [frame $tf1.2]
    label $tf11.__l -text "Thermometer settings:" -anchor w
    pack $tf11.__l -side top -fill x
    FillEntries $tf11 {
	"Format string:" 
	"Label:" 
	"No. of tics:"
    } {
	isoControl(cpl_thermoFmt)
	isoControl(cpl_thermoLabel)
	isoControl(cpl_thermoNTics)
    } 14 20
    button $tf12.font -text "Set Font" -command isoControl_thermoFont
    pack $tf1 -side top -fill both -expand 1 -padx 5 -pady 5
    pack $tf11 $tf12 -side left -fill x -expand 1 -padx 5 -pady 5
    pack $tf12.font -side top -fill x -expand 1

    #
    # RANGES
    #
    set frame_color [lindex \
	     [GetWidgetConfig frame -background] end]
    FillEntries $rf {"Minimum 3D grid value:"} \
	     isosurf(minvalue) \
	     25 9 top left \
	     -e_relief flat -e_state disabled -e_bg $frame_color \
	     -l_font $xcFonts(small) -e_font $xcFonts(small_entry)
    FillEntries $rf {"Maximum 3D grid value:"} \
	     isosurf(maxvalue) \
	     25 9 top left \
	     -e_relief flat -e_state disabled -e_bg $frame_color \
	     -l_font $xcFonts(small) -e_font $xcFonts(small_entry)

    set isoControl(2Dlowvalue_entry) [FillEntries $rf \
	    {"Lowest rendered value:"} isoControl(2Dlowvalue) \
	     25 9 top left \
	     -l_font $xcFonts(small) -e_font $xcFonts(small_entry)]
    set isoControl(2Dhighvalue_entry) [FillEntries $rf \
	    {"Highest rendered value:"} isoControl(2Dhighvalue) \
	    25 9 top left \
	    -l_font $xcFonts(small) -e_font $xcFonts(small_entry)]
    set isoControl(2Dnisoline_entry) [FillEntries $rf \
	    {"Number of isolines:"} \
	    isoControl(2Dnisoline) \
	    25 9 top left \
	    -l_font $xcFonts(small) -e_font $xcFonts(small_entry)]


    #
    # EXPAND
    #
    set f21  [frame $ef.top]
    set f22  [frame $ef.bot]
    pack $f21 $f22 -side top -fill both -expand 1

    set enable "\
	    $f22.scX configure -foreground $enabled_color; \
            $f22.scY configure -foreground $enabled_color; \
            $f22.scZ configure -foreground $enabled_color"
    set disable " \
	    $f22.scX configure -foreground $disabled_color; \
            $f22.scY configure -foreground $disabled_color; \
            $f22.scZ configure -foreground $disabled_color"
    if { $periodic(dim) > 0 } {
	label $f21.l -text "Expand Property-plane:" -relief flat
	$f21.l configure -font $lfont
	set r1 [radiobutton $f21.r1 \
		-text "do not expand" \
		-variable isosurf(2Dexpand) \
		-value "none" \
		-anchor w \
		-font $xcFonts(small) \
		-command "xcDisableAll $f22; catch {eval $disable}"]
	set r2 [radiobutton $f21.r2 \
		-text "to whole structure" \
		-variable isosurf(2Dexpand) \
		-value "whole" \
		-anchor w \
		-font $xcFonts(small) \
		-command "xcDisableAll $f22; catch {eval $disable}"]
	set r3 [radiobutton $f21.r3 \
		-text "separately in each direction" \
		-variable isosurf(2Dexpand) \
		-value "specify" \
		-anchor w \
		-font $xcFonts(small) \
		-command "xcEnableAll $f22; catch {eval $enable}"]
	pack $f21.l -side top -expand 1
	pack $r1 $r2 $r3 -side top -fill x
	
	set XCTrace(2DscX) [scale $f22.scX -from 1 -to $nxdir -length 100 \
		-variable isosurf(2Dexpand_X) -orient horizontal \
		-label "repeat in X-dir:" -tickinterval 1 -resolution 1 \
		-showvalue true \
		-font $xcFonts(small) \
		-width 7 -sliderlength 20]
	# TRACE -- nxdir
	trace variable nxdir w xcTrace
	pack $XCTrace(2DscX) -side left -pady 5 -expand 1 -fill x
	if { $periodic(dim) > 1} {
	    set XCTrace(2DscY) [scale $f22.scY -from 1 -to $nydir \
		    -length 100 \
		    -variable isosurf(2Dexpand_Y) -orient horizontal \
		    -label "repeat in Y-dir:" -tickinterval 1 -resolution 1 \
		    -showvalue true \
		    -font $xcFonts(small) \
		    -width 7 -sliderlength 20]
	    # TRACE -- nydir
	    trace variable nydir w xcTrace
	    pack $XCTrace(2DscY) -side left -pady 5 -expand 1 -fill x
	}
	if { $periodic(dim) > 2 } {
	    set XCTrace(2DscZ) [scale $f22.scZ -from 1 -to $nzdir \
		    -length 100 \
		    -variable isosurf(2Dexpand_Z) -orient horizontal \
		    -label "repeat in Z-dir:" -tickinterval 1 -resolution 1 \
		    -showvalue true \
		    -font $xcFonts(small) \
		    -width 7 -sliderlength 20]
	    # TRACE -- nzdir
	    trace variable nzdir w xcTrace
	    pack $XCTrace(2DscZ) -side left -pady 5 -expand 1 -fill x 
	}
    }

    #
    # Isoline
    #
    set mf1 [frame $ifr.1 -relief groove -bd 2]
    set mf2 [frame $ifr.2 -relief groove -bd 2]
    set mf3 [frame $ifr.3 -relief groove -bd 2]
    pack $mf1 $mf2 $mf3 -side top -expand 1 -fill both -padx 3 -pady 3
    
    set isoControl(bmc) [button $mf1.b -text "set color" \
	    -command IsoControl_SetIsolineColor]

    RadioButVarCmd $mf1 "Isoline Color:" isoControl(isoline_color) \
	    IsoControl_IsolineColor left top 0 0 \
	    {monocolor} {property color}
    pack $isoControl(bmc) -side left -padx 2

    RadioBut $mf2 "Isoline Stipple:" isoControl(isoline_stipple) left top 0 1 \
	    {no stipple} {stipple negative} {full stipple}    
    # this is temporal, since "Isoline Stipple" is not yet working
    xcDisableAll $mf2

    # here goes isoline_width entry
    set isoControl(2Disolinewidth_entry) [FillEntries \
					      $mf3 {"Isoline width:"} \
					      isoControl(isoline_width) \
					      14 10 left left]
    
    ###################
    # ANIMATION FRAME #
    ###################
    set f4 [frame $fb2.f4 -relief raised -bd 2]
    frame $f4.1 -relief groove -bd 2
    pack $f4 -side top -fill x -expand 1
    pack $f4.1 -padx 2 -pady 5 -ipady 3 -side top -fill x -expand 1

    set ckb [checkbutton $f4.1.cb \
	    -text "Apply to Planes #1/#2/#3" \
	    -relief raised -bd 2 \
	    -anchor w \
	    -variable isoControl(anim_apply_to_all)]
    pack $ckb -side top -padx 5 -pady 5 -fill x

    # slide label
    set f41 [frame $f4.1.1]
    pack $f41 -side top -expand 1 -padx 2 
    label $f41.l1 -textvariable isoControl(current_text_slide) -anchor c
    set lfont [ModifyFont [$f41.l1 cget -font] $f41.l1 -size 16 -underline 1]
    $f41.l1 configure -font $lfont
    pack $f41.l1 -side left -padx 5 -fill x -expand 1
    
    # scales
    set f41a [frame $f4.1.1a]
    pack $f41a -side top -expand 1 -fill x -padx 2 
    scale $f41a.s1 -from 1 -to 10 -length 170 \
	-variable isoControl(anim_step) -orient horizontal \
	-label "Animation Step:" -tickinterval 3 -resolution 1 \
	-showvalue true \
	-font $xcFonts(small) \
	-width 7 -sliderlength 20    
    scale $f41a.s2 -from 0 -to 1000 -length 170 \
	-variable isoControl(time_delay) -orient horizontal \
	-label "Delay between slides (in msec):" -tickinterval 200 -resolution 10 \
	-showvalue true \
	-font $xcFonts(small) \
	-width 7 -sliderlength 20
    pack $f41a.s1 $f41a.s2 -side left -padx 5 -expand 1

    set f42 [frame $f4.1.2]
    pack $f42 -side bottom -expand 1 -padx 2 -pady 5
    set first [button $f42.1st -image first    -anchor center \
	    -command [list IsoControl_Animate first]]
    set backw [button $f42.bck -image backward -anchor center \
	    -command [list IsoControl_Animate backward]]
    set previ [button $f42.prv -image previous -anchor center \
	    -command [list IsoControl_Animate previous]]
    set stop [button $f42.sto -image stop -anchor center \
	      -command [list IsoControl_Animate stop]]
    set next  [button $f42.nxt -image next     -anchor center \
	    -command [list IsoControl_Animate next]]
    set forw  [button $f42.frw -image forward  -anchor center \
	    -command [list IsoControl_Animate forward]]
    set last  [button $f42.lst -image last     -anchor center \
	    -command [list IsoControl_Animate last]]
    pack $first $backw $previ $stop $next $forw $last -side left

    foreach {wid text} [list \
			    $first  "First slide" \
			    $backw  "Play backward" \
			    $previ  "Previous slide" \
			    $stop   "Stop playing" \
			    $next   "Next slide" \
			    $forw   "Play forward" \
			    $last   "Last slide"] {
	DynamicHelp::register $wid balloon $text
    }

    ########################################
    # BOTTOM FRAME
    set f5  [frame $fb2.f5 -relief raised -bd 2]
    pack $f5 -side top -fill both -expand 1

    set hid [button $f5.hid -text "Hide" \
	    -command [list IsoControl_Hide $t]]
    set can [button $f5.can -text "Close" \
	    -command [list IsoControlCan $t]]
    set sav [button $f5.sav -text "Save Grid" \
	    -command IsoControlSave]
    set sub [button $f5.sub -text "Submit" \
	    -command IsoControl_UpdateColorplane]
    pack $hid $can $sav $sub -side left -pady 5 -expand 1

}


proc IsoControl_SetIsolineColorOK {type t} {
    global isoControl mody_col
    
    if { $type == "OK" } {
	set cID [xcModifyColorGetID]
	set isoControl(isoline_monocolor) "#$mody_col($cID,hxred)$mody_col($cID,hxgreen)$mody_col($cID,hxblue)"
	if [winfo exists .iso] {
	    IsoControl_UpdateColorplane
	} else {
	    UpdatePropertyPlane
	}
    }
    destroy $t
}
proc IsoControl_SetIsolineColor {} {
    global isoControl
    set t [xcToplevel [WidgetName] \
	    "Set Isoline Color" "Isoline Color" . -0 0 1]
    xcModifyColor $t "Set Isoline Color:" $isoControl(isoline_monocolor) \
	    groove left left 100 100 70 5 20
    
    set ok  [DefaultButton [WidgetName $t] -text "OK" \
	    -command [list IsoControl_SetIsolineColorOK OK $t]]
    set can [button [WidgetName $t] -text "Cancel" \
	    -command [list IsoControl_SetIsolineColorOK Cancel $t]]
    pack $ok $can -padx 10 -pady 10 -expand 1
}
proc IsoControl_IsolineColor item {
    global isoControl
    
    if { $isoControl(isoline_color) == "monocolor" } {
	$isoControl(bmc) config -state normal
    } else {
	$isoControl(bmc) config -state disabled
    }
}


proc IsoControl_Animate what {
    global isoControl prop

    if { ![info exists isoControl(stop_playing)] } {
	set isoControl(stop_playing) 0
    }

    #
    # when I will do PLANE123 simultaneously animation, the isoControl(plane)
    # plane must be called as last with IsoControl_SetCurrentSlide routine,
    # because isoControl(current_text_slide) variable is set in this routine
    #
    set do 1
    set item $isoControl(plane)    
    xcDebug "IsoControl_Animate: $item $isoControl(anim_step) $isoControl($item,current_slide) / $isoControl($item,nslide)"

    if $isoControl(anim_apply_to_all) { 
	if { $item == 1 } { 
	    set item {2 3 1} 
	} elseif { $item == 2 } {
	    set item {3 1 2}
	} else {
	    set item {1 2 3}
	}
    }
    switch -exact -- $what {
	stop {
	    set isoControl(stop_playing) [expr {$isoControl(stop_playing) ? 0 : 1}]
	}
	first   {	
	    foreach i $item {
		if { $isoControl($i,current_slide) == 1 && \
			[llength $item] == 1 } {
		    return
		}
		set isoControl($i,current_slide) 1
	    }
	}
	backward {
	    SetWatchCursor
	    while $do {
		set finished 0
		set nitems [llength $item]
		foreach i $item {
		    set do 0
		    if { $isoControl($i,current_slide) == 1 } {
			incr finished
		    }
		    if { $isoControl($i,current_slide) > 1 } {
			IsoControl_SetCurrentSlide $i -$isoControl(anim_step)
			set do 1
		    }		    
		}
		if { $finished == $nitems || $isoControl(stop_playing)} {
		    set do 0
		    set isoControl(stop_playing) 0
		    ResetCursor
		    xcSwapBuffers
		    return
		}  
		IsoControl_UpdateColorplane

		after $isoControl(time_delay)
	    }
	    ResetCursor
	    xcSwapBuffers
	    return
	}
	previous {
	    foreach i $item {
		if { $isoControl($i,current_slide) == 1 && \
			[llength $item] == 1 } {
		    return
		}
		IsoControl_SetCurrentSlide $i -$isoControl(anim_step)
	    }
	}
	next     {
	    foreach i $item {
		if { $isoControl($i,current_slide) == \
			$isoControl($i,nslide) && [llength $item] == 1 } {
		    return
		}
		IsoControl_SetCurrentSlide $i $isoControl(anim_step)
	    }
	}
	forward  {
	    SetWatchCursor
	    while $do {
		set finished 0
		set nitems [llength $item]
		foreach i $item {
		    set do 0
		    if { $isoControl($i,current_slide) == $isoControl($i,nslide) } {
			incr finished
		    }		    
		    if { $isoControl($i,current_slide) < $isoControl($i,nslide) } {
			IsoControl_SetCurrentSlide $i $isoControl(anim_step)
			set do 1
		    }
		}
		if { $finished == $nitems || $isoControl(stop_playing) } {
		    set do 0
		    set isoControl(stop_playing) 0
		    ResetCursor
		    xcSwapBuffers
		    return
		}
		IsoControl_UpdateColorplane	
		
		after $isoControl(time_delay)	
	    }
	    ResetCursor
	    xcSwapBuffers
	    return
	}
	last     {
	    foreach i $item {
		if { $isoControl($i,current_slide) == \
			$isoControl($i,nslide) && [llength $item] == 1 } {
		    return
		}
		set isoControl($i,current_slide) $isoControl($i,nslide)
	    }
	}
    }
    
    set i $isoControl(plane)        
    set isoControl(current_text_slide) "Current slide:  $isoControl($i,current_slide) / $isoControl($i,nslide)"
    
    IsoControl_UpdateColorplane
}


proc IsoControl_SetCurrentSlide {i incr} {
    global isoControl

    xcDebug "IsoControl_SetCurrentSlide: $i, $incr"
    incr isoControl($i,current_slide) $incr
    if { $isoControl($i,current_slide) < 1 } {
	set isoControl($i,current_slide) 1
    } elseif { $isoControl($i,current_slide) > $isoControl($i,nslide) } {
	set isoControl($i,current_slide) $isoControl($i,nslide)
    }
    #
    # update isoControl(current_text_slide) only if isoControl(plane) plane
    # is requested
    #
    if { $i == $isoControl(plane) } {
	set isoControl(current_text_slide) "Current slide:  $isoControl($i,current_slide) / $isoControl($i,nslide)"
	update
    }
    xcDebug "IsoControl_SetCurrentSlide: $i $incr $isoControl($i,current_slide) / $isoControl($i,nslide)"

}



proc IsoControlCan {t {dim 3}} {
    global isoControl isosurf

    #set button [tk_dialog [WidgetName] WARNING "Are You sure to Close PropertyPlane/IsoSurface. All data will be lost.\n\n Really Close???" warning 0 No Yes]
    #if { $button == 0 } {
    #	return
    #}
    #unset isosurf(3Dinterpl_degree_old)
    #unset isoControl(1,nslide)
    #unset isoControl(2,nslide)
    #unset isoControl(3,nslide)
    #xcTraceDelete nxdir
    #xcTraceDelete nydir
    #xcTraceDelete nzdir
    #xcTraceDelete isosurf(3Dinterpl_degree)
    #xc_iso finish
    #.mesa render

    # display-off the isosurface and colorplanes
   
    if { $dim == 3 } {
	set isoControl(close,isosurf) $isoControl(isosurf)
    }
    set isoControl(close,colorplane) $isoControl(colorplane) 
    set isoControl(close,isoline)    $isoControl(isoline)    

    set isoControl(isosurf)    0
    set isoControl(colorplane) 0
    set isoControl(isoline)    0
    if { $dim == 3 } {
	foreach i {1 2 3} {
	    set isoControl(close,$i,colorplane) $isoControl($i,colorplane)
	    set isoControl(close,$i,isoline)    $isoControl($i,isoline)   
	    set isoControl($i,colorplane) 0
	    set isoControl($i,isoline)    0
	}
	UpdateIsosurf	
    } elseif { $dim == 2 } {
	UpdatePropertyPlane
    }
    CancelProc $t
}



proc IsoControlSave {} {
    global xcMisc system
    
    set filehead [file tail [filehead $xcMisc(titlefile)]]
    set filetypes {
	{{XCrySDen Structure File} {.xsf}     }
	{{All Files}               *          }
    }
    cd $system(PWD)
    set sfile [tk_getSaveFile -initialdir [pwd] \
	    -title "Save Calculated Grid of Points" \
	    -defaultextension .xsf \
	    -initialfile $filehead.xsf \
	    -filetypes $filetypes]
    cd $system(SCRDIR)
    if { $sfile == {} } {
	return 0
    }

    set ident1 UNKNOWN
    set ident2 {}
    # Let the user specify some identifier for the datagrid
    OneEntryToplevel [WidgetName] "DataGrid Identifier" Identifier \
	    "Please specify identifier for DataGrid:" 80 ident2 text 30 30
    if { $ident2 == "" } {
	set ident2 DataGrid_generated_by_XCrySDen1.0
    }
    regsub -all { } $ident2 _ ident2
    
    if [winfo exists .iso] {
	# 3D
	set dg DATAGRID3D
    } elseif [winfo exists .iso2D] {
	#2D
	set dg DATAGRID2D
    }
    
    _IsoControlSave $sfile $ident1 $ident2 $dg
}

proc _IsoControlSave {sfile ident1 ident2 dg} {
    global system

    xc_iso save $sfile $ident1

    # take care about the newlines newline
    set content   [ReadFile -nonewline $system(SCRDIR)/xc_struc.$system(PID)]
    set datablock 0
    set out       {}
    foreach line [split $content \n] {
	if [string match BEGIN_BLOCK_* $line] {
	    set datablock 1
	}
	if !$datablock {
	    append out "$line\n"
	}
	if [string match END_BLOCK_* $line] {
	    set datablock 0
	}
    }
    append out "BEGIN_BLOCK_$dg\n"
    append out "$ident2\n"
    append out [ReadFile -nonewline $sfile]
    append out "\nEND_BLOCK_$dg"
    WriteFile $sfile $out w
}



proc ConvertTwoSideVar {{var {}}} {
    global isosurf
    
    if { ! [info exists isosurf(old_twoside_lighting)] } {
	set isosurf(old_twoside_lighting) off
    }

    if { $var == {}} {
	switch -exact -- $isosurf(twoside_lighting) {
	    0   { set isosurf(twoside_lighting) off }
	    1   { set isosurf(twoside_lighting) on }
	}
    } else {
	switch -exact -- $var {
	    on - 1 - true { 
		xc_setGLparam lightmodel -two_side_iso 1 
		
		# (GL_CCW,GL_CW): this is a dirty trick, namely, the
		# two-side lighting is taking for negative-isosurface
		# the back-side as front side and vice versa
		if { $isosurf(old_twoside_lighting) == "off" } {
		    xc_setGLparam isonormal -what isosurf_neg
		}
	    }
	    off - 0 - false { 
		xc_setGLparam lightmodel -two_side_iso 0 
		
		# (GL_CCW,GL_CW): this is a dirty trick, namely, the
		# two-side lighting is taking for negative-isosurface
		# the back-side as front side and vice versa
		if { $isosurf(old_twoside_lighting) == "on" } {
		    xc_setGLparam isonormal -what isosurf_neg
		}
	    }
	}
	# now render the changes
	.mesa render
	
	set isosurf(old_twoside_lighting) $var
    }
}



proc RevertIsoSides what {
    global openGL

    # just in case if openGL(isoside_$what) does not exist
    if ![info exists openGL(isoside_$what)] {
	set openGL(isoside_$what) [xc_getGLparam frontface -what isosurf_$what]
	xcDebug "openGL(isoside_$what):: $openGL(isoside_$what)"
    }

    if { $openGL(isoside_$what) == "CCW" } {
	set openGL(isoside_$what) "CW"
    } else { 
	set openGL(isoside_$what) "CCW"
    }
    
    xc_setGLparam frontface -what isosurf_$what \
	    -frontface $openGL(isoside_$what)
    # now render the changes
    .mesa render
}


proc RevertIsoNormals what {
    global openGL

    xc_setGLparam isonormal -what isosurf_$what
    # now render the changes
    .mesa render
}


proc SurfaceSmoothing {} {
    global isoControl fillEntries

    set isoControl(smooth_nstep)  [xc_iso get smoothsteps]
    set isoControl(smooth_weight) [xc_iso get smoothwieght]

    set t [xcToplevel [WidgetName] "Surface Smoothing" "SurfSmooth" \
	    .iso 20 20 1]

    message $t.m -aspect 800 \
	    -relief groove -bd 2 \
	    -text "Reasonable values for weight are between 0.1 and 1. Lighter weight will require more steps for smoothing, but will perturb the surface less !!!"
    pack $t.m -side top -padx 3m -pady 3m -ipadx 1m -ipady 1m 

    set f [frame $t.f]	        
    
    FillEntries $t {
	"Smoothing steps:" 
	"Smoothing weight:"
    } {isoControl(smooth_nstep) isoControl(smooth_weight)} 17 7
    set foclist $fillEntries
    set varlist {
	{isoControl(smooth_nstep) int} {isoControl(smooth_weight) real}
    }

    button $t.b1 -text "Close"  -command [list CancelProc $t]
    button $t.b2 -text "Update" \
	    -command [list SurfaceSmoothingOK $t $foclist $varlist]
    
    pack $f -side bottom -expand 1 -fill both  -padx 3m -pady 3m
    pack $t.b1 $t.b2 -side left -expand 1 -padx 2m -pady 2m
}
proc SurfaceSmoothingOK {t foclist varlist} {
    global isoControl
    
    if ![check_var $varlist $foclist] {
	return
    }
    xc_iso smoothsteps  $isoControl(smooth_nstep)
    xc_iso smoothweight $isoControl(smooth_weight)
    xc_iso smoothing
    UpdateIsosurf
    return
}



proc IsoControlCommand {} {
    global isoControl prop

    if $prop(pm_isolevel) {
	#
	# render +isolevel and -isolevel isosurfaces
	#
	
	$isoControl(color_button) configure -command \
		[list MultiWidget {} -b_height 2 -testbutton 1 \
		-create_tplw 1 \
		-tplw_args {xcToplevel [WidgetName] "Set OpenGL parameters" \
		"OpenGLPar"} \
		-command { \
		{"Front Side Color\nfor positive values" \
		{SetOpenGLPar _POS_FRONT_COLOR_ ISOSURF}} \
		{"Back Side Color\nfor positive values"  \
		{SetOpenGLPar _POS_BACK_COLOR_  ISOSURF}} \
		{"Front Side Color\nfor negative values" \
		{SetOpenGLPar _NEG_FRONT_COLOR_ ISOSURF}} \
		{"Back Side Color\nfor negative values"  \
		{SetOpenGLPar _NEG_BACK_COLOR_  ISOSURF}} } \
		-bottom_button { \
		{Close CancelProc} {Update UpdateOpenGLPar} } ]

	$isoControl(blend_button) configure -command \
		[list MultiWidget [WidgetName] -b_height 2 -testbutton 2 \
		-create_tplw 1 \
		-tplw_args {xcToplevel [WidgetName] "Set OpenGL parameters" \
		"OpenGLPar"} \
		-command { \
		{"Transparency" \
		{SetOpenGLPar _BLEND_       ISOSURF}} \
		{"Front Side Color\nfor positive values" \
		{SetOpenGLPar _POS_FRONT_COLOR_ ISOSURF}} \
		{"Back Side Color\nfor positive values"  \
		{SetOpenGLPar _POS_BACK_COLOR_  ISOSURF}} \
		{"Front Side Color\nfor negative values" \
		{SetOpenGLPar _NEG_FRONT_COLOR_ ISOSURF}} \
		{"Back Side Color\nfor negative values"  \
		{SetOpenGLPar _NEG_BACK_COLOR_  ISOSURF}} } \
		-bottom_button { \
		{Close CancelProc} {Update UpdateOpenGLPar} } ]
	$isoControl(revert_button1) configure -state normal
	$isoControl(revert_button2) configure -state normal
    } else {
	#
	# render just isolevel isosurface
	#
	$isoControl(color_button) configure -command \
		[list MultiWidget {} -b_height 2 -testbutton 1 \
		-create_tplw 1 \
		-tplw_args {xcToplevel [WidgetName] "Set OpenGL parameters" \
		"OpenGLPar"} \
		-command { \
		{"Front Side Color" \
		{SetOpenGLPar _ONE_FRONT_COLOR_ ISOSURF}} \
		{"Back Side Color"  \
		{SetOpenGLPar _ONE_BACK_COLOR_ ISOSURF}} } \
		-bottom_button { \
		{Close CancelProc} {Update UpdateOpenGLPar} } ]
	
	$isoControl(blend_button) configure -command \
		[list MultiWidget [WidgetName] -b_height 2 -testbutton 2 \
		-create_tplw 1 \
		-tplw_args {xcToplevel [WidgetName] "Set OpenGL parameters" \
		"OpenGLPar"} \
		-command { \
		{"Transparency" \
		{SetOpenGLPar _BLEND_       ISOSURF}} \
		{"Front Side Color" \
		{SetOpenGLPar _ONE_FRONT_COLOR_ ISOSURF}} \
		{"Back Side Color"  \
		{SetOpenGLPar _ONE_BACK_COLOR_  ISOSURF}} } \
		-bottom_button { \
		{Close CancelProc} {Update UpdateOpenGLPar} } ]
	$isoControl(revert_button1) configure -state disabled
	$isoControl(revert_button2) configure -state disabled
    }
}


# procs for updateing the font for thermometer
proc isoControl_thermoFont {} {
    global isoControl

    puts stderr "*** 1. isoControl_thermoFont: $isoControl(cpl_thermoFont)"    

    set font [fontToplevelWidget [WidgetName] \
		  "Sample Font Text" $isoControl(cpl_thermoFont)]

    puts stderr "*** 2. isoControl_thermoFont: $font"

    global isoControl
    if { $font != {} } {
	set isoControl(cpl_thermoFont) $font
    }
}
