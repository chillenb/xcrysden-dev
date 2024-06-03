#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/isoRender.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc xcPrepareIsosurf {command {dim 3D}} {
    global system isostack isosign isodata isofiles vec prop isosurf
    # RENDER ISOSURFACE
    # ----------------------------------------------------------------------
    # to xc_isopoints 3D we submit point1, point2, point3:
    #
    #                      * Z -> point 3
    #                      *
    #                      *  origin -> point 0
    #                    *   *
    #     point 1 <- X *       * Y-> point 2
    # ----------------------------------------------------------------------

    SetWatchCursor

    # from Bohr to Angstroms
    for {set i 0} {$i < 3} {incr i} {
	for {set j 0} {$j < 3} {incr j} {
	    set vec($i,$j) [Bohr2Angs $vec($i,$j)]
	}
	set isosurf(origin,$i) [Bohr2Angs $isosurf(origin,$i)]
    }

    # origin shift::
    set xs $isosurf(origin,0)
    set ys $isosurf(origin,1)
    set zs $isosurf(origin,2)

    ########################################
    # make some nice output
    regsub -all -- ";" $isostack \n isostack_output
    regsub -all -- ";" $isosign  \n isosign_output
    regsub -all -- ";" $isodata  \n isodata_output

    xcDebug "****************************************"
    xcDebug "xc_iso init"
    xcDebug $isostack_output
    xcDebug "xc_iso end isostack\n"
    xcDebug $isosign_output
    xcDebug "xc_iso end isosign\n"
    xcDebug "xc_isofiles $system(SCRDIR)/xc_binVrt.$system(PID) $system(SCRDIR)/xc_bin.$system(PID) $isofiles\n"
    if { $dim == "3D" } {
    xcDebug "xc_isopoints 0 3D \
	    $xs $ys $zs \
	[expr $xs + $vec(0,0)] [expr $ys + $vec(0,1)] [expr $zs + $vec(0,2)] \
	[expr $xs + $vec(1,0)] [expr $ys + $vec(1,1)] [expr $zs + $vec(1,2)] \
	[expr $xs + $vec(2,0)] [expr $ys + $vec(2,1)] [expr $zs + $vec(2,2)]"
    } else {
    xcDebug "xc_isopoints 0 2D \
	    $xs $ys $zs \
	[expr $xs + $vec(0,0)] [expr $ys + $vec(0,1)] [expr $zs + $vec(0,2)] \
	[expr $xs + $vec(1,0)] [expr $ys + $vec(1,1)] [expr $zs + $vec(1,2)]" 
    xcDebug $isodata_output
    xcDebug "xc_isodata make\n"
    xcDebug "****************************************"
    }
    # END_OF_OUTPUT SECTION
    ########################################

    xc_iso init
    eval $isostack
    xc_iso end isostack
    
    # assign singns to all stck-levels 
    eval $isosign
    xc_iso end isosign

    eval xc_isofiles $system(SCRDIR)/xc_binVrt.$system(PID) $system(SCRDIR)/xc_bin.$system(PID) $isofiles
    
    xcDebug "xc_isofiles $system(SCRDIR)/xc_binVrt.$system(PID) $system(SCRDIR)/xc_bin.$system(PID) $isofiles"

    if { $dim == "3D" } {
	xc_isopoints 0 3D \
		$xs $ys $zs \
	 [expr $xs + $vec(0,0)] [expr $ys + $vec(0,1)] [expr $zs + $vec(0,2)] \
	 [expr $xs + $vec(1,0)] [expr $ys + $vec(1,1)] [expr $zs + $vec(1,2)] \
	 [expr $xs + $vec(2,0)] [expr $ys + $vec(2,1)] [expr $zs + $vec(2,2)]
    } elseif { $dim == "2D" } {
	xc_isopoints 0 2D \
	    $xs $ys $zs \
	    [expr $xs + $vec(0,0)] [expr $ys + $vec(0,1)] [expr $zs + $vec(0,2)] \
	    [expr $xs + $vec(1,0)] [expr $ys + $vec(1,1)] [expr $zs + $vec(1,2)]
    }

    eval $isodata
    xc_isodata make
    
    #
    # if isosurf(3Dinterpl_degree) > 1; perform triCubic spline interpolation
    #
    if { $isosurf(3Dinterpl_degree) > 1 } {
	xc_iso interpolate $isosurf(3Dinterpl_degree)
    }

    # get MIN-ISO-VALUE & MAX-ISO-VALUE
    if { $isosurf(spin) == "ALPHA" || \
	    $isosurf(spin) == "BETA" } {
	# in C95 "spin-dependent" properties are ECHD & ECHG, 
	# first record is ALPHA+BETA, second is ALPHA-BETA, so:
	# ALPHA = (record1 + record2) / 2, but in C-code I just use
	# (record1 + record2), so we must correct it here by factor 2.0	    
	set isosurf(minvalue) [expr [xc_iso minvalue] / 2.0]
	set isosurf(maxvalue) [expr [xc_iso maxvalue] / 2.0]
    } else {
	set isosurf(minvalue) [xc_iso minvalue]
	set isosurf(maxvalue) [xc_iso maxvalue]
    }

    set isosurf(rangevalue) [expr $isosurf(maxvalue) - $isosurf(minvalue)]

    ResetCursor
}



proc UpdateIsosurf {{var {}}} {
    global prop isosurf periodic isosign isodata isosurf_struct \
	    isoControl isoControl_struct vec dif_isosurf
    # var.......because this routine is used also by RadioButCmd, and this 
    # routine pass one parameter to command, we must make "var" here a
    # dummy variable
    ###############################################
    # before anything is done, check $prop(isovalue)

    SetWatchCursor

    xcDebug "1,current_slide:: $isoControl(1,current_slide)"
    xcDebug "2,current_slide:: $isoControl(2,current_slide)"
    xcDebug "3,current_slide:: $isoControl(3,current_slide)"
    #
    # ISOSURFACE-SECTION
    #
    xcDebug "\nUpdateIsosurf\n"
    xcDebug " xc_isosurf .mesa 0 \
	    -drawstyle    $isosurf(type_of_isosurf) \
	    -shademodel   $isosurf(shade_model) \
	    -transparency $isosurf(transparency)"

    if { $isoControl(isosurf) == 1 } { 
	if { $var != {} && ![info exists prop(isolevel)] } {
	    ResetCursor
	    return
	}
    }

    if { $var != {} && $prop(isolevel) == {} } {
	ResetCursor
	return
    }

    if { $isoControl(isosurf) == 1 } { 
	if ![check_var {{prop(isolevel) real}} $isosurf(isovalue_entry)] {
	    ResetCursor
	    return
	}
    }

    ##############
    # ISOSURFACE #
    ##############
    set polygonise 0; # do we need to call xc_iso polygonise
    set expand {0 0 0}
    #########################
    # this is needed if spin is changed
    if { $isosurf_struct(isosign)    != $isosign || \
	    $isosurf_struct(isodata) != $isodata || \
	    $isosurf_struct(spin)    != $isosurf(spin) } {	    	    
	SetXC_Iso
	eval $isosign
	xc_iso end isosign
	eval $isodata
	xc_isodata make
	set polygonise 1
    }
    #########################

    if { $isosurf_struct(3Dinterpl_degree) != $isosurf(3Dinterpl_degree) } {
	xc_iso interpolate $isosurf(3Dinterpl_degree)
	#
	# now IsoPlane will have to be updated, because ***plvertex was
	# reallocated with xc_iso interpolate
	#
	set upd {1 2 3}
    }

    if { $isoControl(isosurf) } {
	if { $isosurf_struct(isolevel) != $prop(isolevel) || \
		 $isosurf_struct(pm_isolevel) != $prop(pm_isolevel) || \
		 $isosurf_struct(3Dinterpl_degree) !=  $isosurf(3Dinterpl_degree) || \
		 $isosurf_struct(tessellation_type) != $isosurf(tessellation_type) || \
		 $isosurf_struct(normals_type) != $isosurf(normals_type) || \
		 $isosurf_struct(shade_model) != $isosurf(shade_model) } {
	    if { $isosurf(spin) == "ALPHA" || \
		    $isosurf(spin) == "BETA" } {
		# in C95 "spin-dependent" properties are ECHD & ECHG, 
		# first record is ALPHA+BETA, second is ALPHA-BETA, so:
		# ALPHA = (record1 + record2) / 2, but in C-code I just use
		# (record1 + record2), so we must correct it here by factor 2.0
		if { !$prop(pm_isolevel) || \
			[IsEqual [expr $isosurf(rangevalue) / 1.0e7] \
			$prop(isolevel)] } {
		    #
		    # render just isolevel isosurface
		    #
		    xc_iso isolevel [expr $prop(isolevel) * 2.0]
		} else {
		    #
		    # render +/- isolevel isosurface
		    #
		    xc_iso isolevel \
			    [expr $prop(isolevel) * 2.0] \
			    [expr $prop(isolevel) * -2.0]
		}
	    } else {
		if { !$prop(pm_isolevel) || \
			[IsEqual \
			[expr $isosurf(rangevalue) / 1.0e7] $prop(isolevel)] } {
		    #
		    # render just isolevel isosurface
		    #		
		    xcDebug "ISOSURFACE: +isovalue"
		    xc_iso isolevel $prop(isolevel)
		} else {
		    #
		    # render +/- isolevel isosurface
		    #
		    xcDebug "ISOSURFACE: +/-isovalue"
		    xc_iso isolevel $prop(isolevel) [expr $prop(isolevel) * -1.0]
		    
		}
	    }
	    set polygonise 1
	}
	
	if { $polygonise } { 
	    xc_iso polygonise -algorithm $isosurf(tessellation_type) -shademodel $isosurf(shade_model) -normals $isosurf(normals_type)
	}
	
	xc_isosurf .mesa \
		-drawstyle    $isosurf(type_of_isosurf) \
		-shademodel   $isosurf(shade_model) \
		-transparency $isosurf(transparency)
		
	if { $isosurf(spin) == "ALPHA" || \
		$isosurf(spin) == "BETA" } {
	    # in C95 "spin-dependent" properties are ECHD & ECHG, 
	    # first record is ALPHA+BETA, second is ALPHA-BETA, so:
	    # ALPHA = (record1 + record2) / 2, but in C-code I just use
	    # (record1 + record2), so we must correct it here by factor 2.0	    
	    set isosurf(minvalue) [expr [xc_iso minvalue] / 2.0]
	    set isosurf(maxvalue) [expr [xc_iso maxvalue] / 2.0]
	} else {
	    set isosurf(minvalue) [xc_iso minvalue]
	    set isosurf(maxvalue) [xc_iso maxvalue]
	}
	
	########################################
	# now take care of xc_isoexpand
	# isosurf(expand) ......... none/whole/specify
	#                                         |-> isosurf(expand_[X,Y,Z])
	if { $isosurf(expand) == "none" } {
	    if { $periodic(dim) == 3 } {
		set expand {1 1 1}
	    } elseif { $periodic(dim) == 2 } {
		set expand {1 1 0}
	    } elseif { $periodic(dim) == 1 } {
		set expand {1 0 0}
	    } else {
		set expand {0 0 0}
	    }
	} elseif { $isosurf(expand) == "whole" } {
	    set expand "whole"
	} else {
	    # expand == specify
	    if { $periodic(dim) == 3 } {
		set expand [list $isosurf(expand_X) $isosurf(expand_Y) \
			$isosurf(expand_Z)]
	    } elseif { $periodic(dim) == 2 } {
		set expand [list $isosurf(expand_X) $isosurf(expand_Y) 0]
	    } elseif { $periodic(dim) == 1 } {
		set expand [list $isosurf(expand_X) 0 0]
	    } else {
		set expand {0 0 0}
	    }
	}
	
	xc_isoexpand .mesa 0 -repeattype default \
		-shape parapipedal -expand $expand -render after
    } else {
	# 
	# do not render an isosurface
	#
	xc_isosurf .mesa -isosurf none
    }

    ################################################
    # should PROPERTY-PLANE THERMOMETER be updated #
    ################################################
    set upd_thermo {} 
    foreach i {1 2 3} {
	if { \
		 $isoControl_struct($i,cpl_basis)        != $isoControl($i,cpl_basis) || \
		 $isoControl_struct($i,cpl_function)     != $isoControl($i,cpl_function) || \
		 $isoControl_struct($i,colorplane)       != $isoControl($i,colorplane) || \
		 $isoControl_struct($i,cpl_thermometer)  != $isoControl($i,cpl_thermometer) || \
		 $isoControl_struct($i,cpl_thermoTplw)  != $isoControl($i,cpl_thermoTplw) || \
		 $isoControl_struct($i,cpl_thermoFmt)    != $isoControl($i,cpl_thermoFmt) || \
		 $isoControl_struct($i,cpl_thermoLabel)  != $isoControl($i,cpl_thermoLabel) || \
		 $isoControl_struct($i,cpl_thermoNTics)  != $isoControl($i,cpl_thermoNTics) || \
		 $isoControl_struct($i,cpl_thermoFont)   != $isoControl($i,cpl_thermoFont) || \
		 $isoControl_struct($i,2Dlowvalue)       != $isoControl($i,2Dlowvalue)        || \
		 $isoControl_struct($i,2Dhighvalue)      != $isoControl($i,2Dhighvalue) } {
	    append upd_thermo " $i"
	}
    }
		     
    #############################################
    # should PROPERTY-PLANE be rendered as well #
    #############################################
    #
    # if it is something has changed from the last time, 
    # then all three planes should be updated, else just current plane
    #
    if ![info exists upd] { 
	set upd {} 
	foreach i {1 2 3} {
	    if { \
		     $isoControl_struct($i,cpl_basis)        != $isoControl($i,cpl_basis) || \
		     $isoControl_struct($i,cpl_function)     != $isoControl($i,cpl_function) || \
		     $isoControl_struct($i,colorplane)       != $isoControl($i,colorplane) || \
		     $isoControl_struct($i,isoline)          != $isoControl($i,isoline)	     || \
		     $isoControl_struct($i,colorplane_lighting) != $isoControl($i,colorplane_lighting)	     || \
		     $isoControl_struct($i,cpl_transparency) != $isoControl($i,cpl_transparency)  || \
		     $isoControl_struct($i,cpl_thermometer)  != $isoControl($i,cpl_thermometer) || \
		     $isoControl_struct($i,cpl_thermoTplw)  != $isoControl($i,cpl_thermoTplw) || \
		     $isoControl_struct($i,cpl_thermoFmt)    != $isoControl($i,cpl_thermoFmt) || \
		     $isoControl_struct($i,cpl_thermoLabel)  != $isoControl($i,cpl_thermoLabel) || \
		     $isoControl_struct($i,cpl_thermoNTics)  != $isoControl($i,cpl_thermoNTics) || \
		     $isoControl_struct($i,cpl_thermoFont)   != $isoControl($i,cpl_thermoFont) || \
		     $isosurf_struct($i,2Dexpand)            != $isosurf($i,2Dexpand)  	     || \
		     $isosurf_struct($i,2Dexpand_X)          != $isosurf($i,2Dexpand_X)	     || \
		     $isosurf_struct($i,2Dexpand_Y)          != $isosurf($i,2Dexpand_Y)	     || \
		     $isosurf_struct($i,2Dexpand_Z)          != $isosurf($i,2Dexpand_Z)	     || \
		     $isoControl_struct($i,anim_step)        != $isoControl($i,anim_step)	     || \
		     $isoControl_struct($i,current_slide)    != $isoControl($i,current_slide)     || \
		     $isoControl_struct($i,2Dlowvalue)       != $isoControl($i,2Dlowvalue)        || \
		     $isoControl_struct($i,2Dhighvalue)      != $isoControl($i,2Dhighvalue)	     || \
		     $isoControl_struct($i,2Dnisoline)       != $isoControl($i,2Dnisoline)	     || \
		     $isoControl_struct($i,isoline_color)    != $isoControl($i,isoline_color)     || \
		     $isoControl_struct($i,isoline_width)    != $isoControl($i,isoline_width)     || \
		     $isoControl_struct($i,isoline_monocolor)!= $isoControl($i,isoline_monocolor) || \
		     $isoControl_struct($i,isoline_stipple)  != $isoControl($i,isoline_stipple) } {
		append upd " $i"
	    }
	}
    }
    xcDebug --stderr "UpdateIsosurf -- upd == $upd" 
    foreach i $upd {
	if { $isoControl($i,colorplane) && !$isoControl($i,isoline) } {
	    set what "colorplane"
	} elseif { !$isoControl($i,colorplane) && $isoControl($i,isoline) } {
	    set what "isoline"
	} elseif { $isoControl($i,colorplane) && $isoControl($i,isoline) } {
	    set what "both"
	} else {
	    set what "none"
	}
	if { $isoControl($i,cpl_basis) == {} } {
	    tk_dialog [WidgetName] ERROR "ERROR: You forgot to specify color basis for plane #$i. Please Do it !" error 0 OK
	    ResetCursor
	    return
	}
	if { $isoControl($i,cpl_function) == {} } {
	    tk_dialog [WidgetName] ERROR "ERROR: You forgot to specify scale function for plane #$i. Please Do it !" error 0 OK
	    ResetCursor
	    return
	}
	#
	# this must be the same as in "isosurf.h"
	#
	set COLORBASE_MONO       0
	set COLORBASE_RAINBOW    1
	set COLORBASE_RGB        2
	set COLORBASE_GEOGRAPHIC 3
	set COLORBASE_BLUE_WHITE_RED 4
	set COLORBASE_BLACK_BROWN_WHITE 5
	
	set SCALE_FUNC_LIN      0
	set SCALE_FUNC_LOG      1
	set SCALE_FUNC_LOG10    2
	set SCALE_FUNC_SQRT     3
	set SCALE_FUNC_ROOT3    4
	set SCALE_FUNC_GAUSS    5
	set SCALE_FUNC_SLATER   6

	switch -exact -- $isoControl($i,cpl_basis) {
	    MONOCHROME     { set basis $COLORBASE_MONO }
	    RAINBOW        { set basis $COLORBASE_RAINBOW }
	    RGB	           { set basis $COLORBASE_RGB }
	    GEOGRAPHIC     { set basis $COLORBASE_GEOGRAPHIC }
	    BLUE-WHITE-RED { set basis $COLORBASE_BLUE_WHITE_RED }
	    BLACK-BROWN-WHITE { set basis $COLORBASE_BLACK_BROWN_WHITE }
	}
	switch -exact -- $isoControl($i,cpl_function) {
	    LINEAR     { set func $SCALE_FUNC_LIN    } 
	    LOG        { set func $SCALE_FUNC_LOG    } 
	    LOG10      { set func $SCALE_FUNC_LOG10  } 
	    SQRT       { set func $SCALE_FUNC_SQRT   } 
	    3th-ROOT   { set func $SCALE_FUNC_ROOT3  } 
	    EXP(x)     { set func $SCALE_FUNC_SLATER } 
	    EXP(x^2)   { set func $SCALE_FUNC_GAUSS  } 
	}
	
	if { $isoControl(isoline_color) == "monocolor" } {
	    set linecolor \
		    [concat monocolor [rgb_h2f $isoControl(isoline_monocolor)]]
	} else {
	    set linecolor polycolor
	}
	switch -exact -- $isoControl(isoline_stipple)  {
	    {no stipple}        { set linedash nodash   }	   
	    {stipple negative}  { set linedash negdash  }
	    {full stipple}      { set linedash fulldash }
	}

	xc_iso isoplaneconfig $i \
	    -isoplanemin      $isoControl($i,2Dlowvalue) \
	    -isoplanemax      $isoControl($i,2Dhighvalue) \
	    -isolinecolor     $linecolor \
	    -isolinewidth     $isoControl($i,isoline_width) \
	    -isolinedash      $linedash \
	    -isolinenlevels   $isoControl($i,2Dnisoline) \
	    -isoplanelighting $isoControl($i,colorplane_lighting)
	
	xc_iso isoplane $i $basis $func $what $isoControl($i,current_slide)
	if { $what != "none" } {
	    xc_isoplane .mesa $i -planetype $what \
		    -transparency $isoControl($i,cpl_transparency) \
		    -render after

	    ########################################
	    # now take care of xc_isoexpand
	    # isosurf(expand) ......... none/whole/specify
	    #                                      |-> isosurf(expand_[X,Y,Z])
	    if { $isosurf($i,2Dexpand) == "none" } {
		if { $periodic(dim) == 3 } {
		    set expand2D {1 0 0  0 1 0  0 0 1    0 0 0}
		} elseif { $periodic(dim) == 2 } {
		    set expand2D {1 0 0  0 1 0  0 0 1    0 0 0}
		} elseif { $periodic(dim) == 1 } {
		    set expand2D {1 0 0  0 1 0  0 0 1    0 0 0}
		} else {
		    set expand2D {1 0 0  0 1 0  0 0 1    0 0 0}
		}
	    } elseif { $isosurf($i,2Dexpand) == "whole" } {
		set expand2D "whole"
	    } else {
		# expand == specify
		if { $periodic(dim) == 3 } {
		    set expand2D [list 1 0 0  0 1 0  0 0 1 \
			    $isosurf($i,2Dexpand_X) \
			    $isosurf($i,2Dexpand_Y) \
			    $isosurf($i,2Dexpand_Z)]
		} elseif { $periodic(dim) == 2 } {
		    set expand2D [list 1 0 0  0 1 0  0 0 1 \
			    $isosurf($i,2Dexpand_X) $isosurf($i,2Dexpand_Y) 0]
		} elseif { $periodic(dim) == 1 } {
		    set expand2D [list 1 0 0  0 1 0  0 0 1 \
			    $isosurf($i,2Dexpand_X) 0 0]
		} else {
		    set expand2D {1 0 0  0 1 0  0 0 1    0 0 0}
		}
	    }
	
	    xc_isoexpand .mesa $i -repeattype default \
		    -shape parapipedal -expand $expand2D -render after
	}
    }


    ############################
    # take care of thermometer #
    ############################
    foreach i $upd_thermo {
	if { $isoControl($i,cpl_thermometer) } {
	    # render thermometer !!!	   
	    foreach wid [list .mesa.thermo$i .thermo$i] {
		if { [winfo exists $wid] } {
		    destroy $wid
		}	    
	    }
	    if { $isoControl($i,cpl_thermoTplw) == 1 } {
		set prefix {}
	    } else {
		set prefix .mesa
	    }
	    set tw [thermometerWidget $prefix.thermo$i $i \
			$isoControl($i,cpl_thermoLabel) \
			$isoControl($i,cpl_thermoFont) \
			$isoControl($i,cpl_thermoFmt) \
			$isoControl($i,cpl_basis) \
			$isoControl($i,cpl_function) \
			$isoControl($i,2Dlowvalue) \
			$isoControl($i,2Dhighvalue) \
			$isoControl($i,cpl_thermoNTics) \
			$isoControl($i,cpl_thermoTplw)]
	    if { $isoControl($i,cpl_thermoTplw) == 0 } {
		pack $tw -anchor nw -side top
	    }
	} else {
	    if { [winfo exists .mesa.thermo$i] } {
		destroy .mesa.thermo$i
	    } elseif { [winfo exists .thermo$i] } {
		destroy .thermo$i
	    }
	}
    }
    
    ######################
    # end PROPERTY-PLANE #
    ######################
    
    # update the display
    #xc_display .mesa
    .mesa render

    # now updata isosurf_struct global variable
    Set_UpdateIsosurf_Struct $expand

    ResetCursor
}
    



proc Set_UpdateIsosurf_Struct {{expand {}}} {
    global isosurf_struct isostack isosign isodata prop isosurf \
	    isoControl isoControl_struct

    if [info exists isosurf(3Dinterpl_degree)] {
	set isosurf_struct(3Dinterpl_degree) $isosurf(3Dinterpl_degree)
    } else {
	set isosurf(3Dinterpl_degree) 1
	set isosurf_struct(3Dinterpl_degree) 1
    }
    #############################################
    if [info exists isosurf(spin)] {
	set isosurf_struct(spin) $isosurf(spin)
    } else {
	set isosurf_struct(spin) {}
    }
    #############################################
    if [info exists isostack] {
	set isosurf_struct(isostack) $isostack
    } else {
	set isosurf_struct(isostack) {}
    }
    #############################################
    if [info exists isosign] {	
	set isosurf_struct(isosign) $isosign
    } else {
	set isosurf_struct(isosign) {}
    }
    ##############################################
    if [info exists isodata] {
	set isosurf_struct(isodata) $isodata
    } else {
	set isosurf_struct(isodata) {}
    }
    #############################################
    if [info exists prop(isolevel)] {
	set isosurf_struct(isolevel) $prop(isolevel)
    } else {	
	set isosurf_struct(isolevel) {}
    }
    #############################################
    if [info exists prop(pm_isolevel)] {
	set isosurf_struct(pm_isolevel) $prop(pm_isolevel)
    } else {	
	set isosurf_struct(pm_isolevel) 0
    }

    #############################################
    if [info exists isosurf(shade_model)] {
	set isosurf_struct(shade_model)   $isosurf(shade_model)
    }

    #############################################
    if [info exists isosurf(tessellation_type)] {
	set isosurf_struct(tessellation_type)   $isosurf(tessellation_type)
    }

    #############################################
    if [info exists isosurf(normals_type)] {
	set isosurf_struct(normals_type)   $isosurf(normals_type)
    }

    #############################################
    if [info exists isosurf(expand)] {
	set isosurf_struct(expand)   $isosurf(expand)
    }
    #############################################
    # now here we will do a little different
    if ![info exists isosurf_struct(expand_XYZ)] {
	set isosurf_struct(expand_XYZ) {}
    }	
    if { $expand != {} } {
	set isosurf_struct(expand_XYZ) $expand
    }

    #############################################
    # isoControl array for PROPERTY PLANE       #
    #############################################
    #
    # have we already call IsoControl_InitVar ???
    #
    if ![info exists isoControl(cpl_basis)] {
	IsoControl_InitVar
    }

    foreach i {1 2 3} {
	set isoControl_struct($i,cpl_basis)           $isoControl($i,cpl_basis)
	set isoControl_struct($i,cpl_function)        $isoControl($i,cpl_function)
	set isoControl_struct($i,colorplane)          $isoControl($i,colorplane)
	set isoControl_struct($i,isoline)             $isoControl($i,isoline)
	set isoControl_struct($i,colorplane_lighting) $isoControl($i,colorplane_lighting)
	set isoControl_struct($i,cpl_transparency)    $isoControl($i,cpl_transparency)
	set isoControl_struct($i,cpl_thermometer)     $isoControl($i,cpl_thermometer)
	set isoControl_struct($i,cpl_thermoTplw)     $isoControl($i,cpl_thermoTplw)
	set isoControl_struct($i,cpl_thermoFmt)       $isoControl($i,cpl_thermoFmt)
	set isoControl_struct($i,cpl_thermoLabel)     $isoControl($i,cpl_thermoLabel) 
	set isoControl_struct($i,cpl_thermoNTics)     $isoControl($i,cpl_thermoNTics)
	set isoControl_struct($i,cpl_thermoFont)      $isoControl($i,cpl_thermoFont)
	set isosurf_struct($i,2Dexpand)             $isosurf($i,2Dexpand)  
	set isosurf_struct($i,2Dexpand_X)           $isosurf($i,2Dexpand_X)
	set isosurf_struct($i,2Dexpand_Y)           $isosurf($i,2Dexpand_Y)
	set isosurf_struct($i,2Dexpand_Z)           $isosurf($i,2Dexpand_Z)
	set isoControl_struct($i,anim_step)         $isoControl($i,anim_step)
	set isoControl_struct($i,current_slide)     $isoControl($i,current_slide)
	set isoControl_struct($i,2Dlowvalue)        $isoControl($i,2Dlowvalue)     
	set isoControl_struct($i,2Dhighvalue)       $isoControl($i,2Dhighvalue)
	set isoControl_struct($i,2Dnisoline)        $isoControl($i,2Dnisoline)
	set isoControl_struct($i,isoline_color)     $isoControl($i,isoline_color)
	set isoControl_struct($i,isoline_width)     $isoControl($i,isoline_width)
	set isoControl_struct($i,isoline_monocolor) $isoControl($i,isoline_monocolor)
	set isoControl_struct($i,isoline_stipple)   $isoControl($i,isoline_stipple)
    }

    #
    # ISOOBJ_BASE == 2D::
    #
    set isoControl_struct(cpl_basis)           $isoControl(cpl_basis)
    set isoControl_struct(cpl_function)        $isoControl(cpl_function)
    set isoControl_struct(colorplane)          $isoControl(colorplane)
    set isoControl_struct(isoline)             $isoControl(isoline)
    set isoControl_struct(colorplane_lighting) $isoControl(colorplane_lighting)
    set isoControl_struct(cpl_transparency)    $isoControl(cpl_transparency)
    set isoControl_struct(cpl_thermometer)     $isoControl(cpl_thermometer)
    set isoControl_struct(cpl_thermoTplw)      $isoControl(cpl_thermoTplw)
    set isoControl_struct(cpl_thermoFmt)       $isoControl(cpl_thermoFmt)
    set isoControl_struct(cpl_thermoLabel)     $isoControl(cpl_thermoLabel) 
    set isoControl_struct(cpl_thermoNTics)     $isoControl(cpl_thermoNTics)
    set isoControl_struct(cpl_thermoFont)      $isoControl(cpl_thermoFont)
    set isosurf_struct(2Dexpand)               $isosurf(2Dexpand)  
    set isosurf_struct(2Dexpand_X)             $isosurf(2Dexpand_X)
    set isosurf_struct(2Dexpand_Y)             $isosurf(2Dexpand_Y)
    set isosurf_struct(2Dexpand_Z)             $isosurf(2Dexpand_Z)
    set isoControl_struct(anim_step)           $isoControl(anim_step)
    set isoControl_struct(current_slide)       $isoControl(current_slide)
    set isoControl_struct(2Dlowvalue)          $isoControl(2Dlowvalue)     
    set isoControl_struct(2Dhighvalue)         $isoControl(2Dhighvalue)
    set isoControl_struct(2Dnisoline)          $isoControl(2Dnisoline)
    set isoControl_struct(isoline_color)       $isoControl(isoline_color)
    set isoControl_struct(isoline_width)       $isoControl(isoline_width)
    set isoControl_struct(isoline_monocolor)   $isoControl(isoline_monocolor)
    set isoControl_struct(isoline_stipple)     $isoControl(isoline_stipple)
}

