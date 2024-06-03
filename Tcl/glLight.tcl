#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/glLight.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc glLight {{togl .mesa}} {
    global glLight xcFonts

    if { ! [info exists glLight(nlights)] } {
	set glLight(nlights) 6
    }
    
    set t .gl_light
    if { [winfo exists $t] } {
	return
    }
    xcToplevel $t "Lights Setting" "Lights"

    set nb [NoteBook $t.nb]
    pack $nb -expand 1 -fill both

    for {set i 0} {$i < $glLight(nlights)} {incr i} {
	
	glLight:loadValues $i

	# insert page
	$nb insert $i light$i -text "Light \#.$i"
	set page [$nb getframe light$i]

	set cb [checkbutton $page.cb -text "Enable Light No.$i" -font $xcFonts(big) \
		    -relief raised -bd 1 \
		    -onvalue 1 -offvalue 0 \
		    -variable glLight($i)]
	pack $cb -side top -ipadx 5 -ipady 5 -padx 5 -pady 5 -fill x

	set f1 [frame $page.1 -relief groove -bd 2]
	pack $f1 -side top -expand 1 -fill both -padx 5 -pady 5

	set cb1 [checkbutton $f1.cb -text "Local viewer lightmodel" \
		     -onvalue 1 -offvalue 0 \
		     -relief raised -bd 1 \
		     -variable glLight(lightmodel_local_viewer)]
	pack $cb1 -side top -padx 5 -pady 5 -fill x


	# in frame $f1 first goes the LIGHTMODEL ambient light

	set la [frame $f1.lightmodel_ambient -relief groove -bd 2]
	pack $la -side left -fill both -padx 2 -pady 5 -ipady 3 -expand 1
	setRGBAwidget $la "Lightmodel\nambient light:" \
	    glLight(lightmodel_ambient_R) glLight(lightmodel_ambient_G) \
	    glLight(lightmodel_ambient_B) glLight(lightmodel_ambient_A) \
	    _LIGHTMODEL_AMBIENT_
	
	#
	# in frame $f1 the RGBA widgets for ambient, diffuse, and 
	# specular light will go.
	# 		
	foreach f {ambient diffuse specular} {
	    set frame($f)  [frame $f1.$f -relief groove -bd 2]
	    pack $frame($f) -side left -fill both -padx 2 -pady 5 -ipady 3 -expand 1	    

	    setRGBAwidget $frame($f) "[string totitle $f] light:\n" \
		glLight($i,${f}_R) glLight($i,${f}_G) glLight($i,${f}_B) glLight($i,${f}_A) \
		_LIGHT${i}_[string toupper $f]_
	}

	#
	# if frame $f2 goes the rest ...
	#
	#			     -fract_position {x y z w} 
	#			     -spot_dir       {x,y,z} 	
	#			     -spot_exp       {v} 
	#			     -spot_cutoff    {v} 
	#			     -const_atten    {v} 
	#			     -lin_atten      {v} 
	#			     -quad_atten     {v}


	# LIGHT-POSITION

	set fracpos [frame $page.fracpos -relief groove -bd 2]
	set l1 [label $fracpos.__l -text "Light position (in fractions of molecular size)"]
	pack $fracpos -side top -fill x -padx 5 -pady 5
	pack $l1 -side top -fill x

	set x [Entries $fracpos "X: Y: Z: W:" \
		   [list \
			glLight($i,fract_position_X) \
			glLight($i,fract_position_Y) \
			glLight($i,fract_position_Z) \
			glLight($i,fract_position_W)] \
		   10 1 left -fill x]
	# add here code for entry-checking


	# SPOT-DIRECTION

	set spotdir [frame $page.spotdir -relief groove -bd 2]
	set l2 [label $spotdir.__l -text "Specular spot direction (in fractions of molecular size)"]
	pack $spotdir -side top -fill x -padx 5 -pady 5
	pack $l2 -side top -fill x

	set x [Entries $spotdir "X: Y: Z:" \
		   [list \
			glLight($i,spot_dir_X) \
			glLight($i,spot_dir_Y) \
			glLight($i,spot_dir_Z)] \
		   10 1 left -fill x ]
	# add here code for entry-checking


	# the REST of the coeeficients
	
	set coef [frame $page.coef -relief groove -bd 2 ]
	pack $coef -side top -fill x -padx 5 -pady 5
	FillEntries $coef {
	    "Specular spot exponent:"
	    "Specular spot cuttof:"
	    "Constant light attenuation:"
	    "Linear light attenuation:"
	    "Quadratic light attenuation:"
	} [list \
	       glLight($i,spot_exp) \
	       glLight($i,spot_cutoff) \
	       glLight($i,const_atten) \
	       glLight($i,lin_atten) \
	       glLight($i,quad_atten) \
	      ] \
	    27 10
	# add here code for entry-checking

	#
	# in bottom frame goes the "Close|Update" buttons
	#
	set bottom [frame $page.bottom]
	pack $bottom -side top -fill x -padx 5 -pady 5

	set close  [button $bottom.close  -text "Close"  -command [list CancelProc $t]]
	set save   [button $bottom.save   -text "Save Lights" -command glLight:save]
	set load   [button $bottom.load   -text "Load Lights" -command glLight:load]
	set update [button $bottom.update -text "Update" -command [list glLight:update $togl $i]]
	pack $close $save $load $update -side left -padx 5 -pady 5 -ipadx 3 -ipady 3 -expand 1		    
    }

    # show the page for light #.0

    $nb raise light0
}



proc glLight:update {togl i} {
    global glLight

    # load LIGHTMODEL

    if { $glLight($i) } {
	xcDebug "LIGHTMODEL> xc_setGLparam lightmodel -enable_light $i"
	xc_setGLparam lightmodel -enable_light $i
    } else {
	xcDebug "LIGHTMODEL> xc_setGLparam lightmodel -disable_light $i"
	xc_setGLparam lightmodel -disable_light $i
    }    
    xc_setGLparam lightmodel -ambient \
	[list \
	     $glLight(lightmodel_ambient_R) $glLight(lightmodel_ambient_G) \
	     $glLight(lightmodel_ambient_B) $glLight(lightmodel_ambient_A)]
    xc_setGLparam lightmodel -local_viewer $glLight(lightmodel_local_viewer)

    
    # load LIGHTS

    xc_setGLparam light \
	-light $i \
	-ambient  [list $glLight($i,ambient_R)  $glLight($i,ambient_G)  $glLight($i,ambient_B)  $glLight($i,ambient_A)]  \
	-diffuse  [list $glLight($i,diffuse_R)  $glLight($i,diffuse_G)  $glLight($i,diffuse_B)  $glLight($i,diffuse_A)]  \
	-specular [list $glLight($i,specular_R) $glLight($i,specular_G) $glLight($i,specular_B) $glLight($i,specular_A)] \
	-fract_position [list \
			     $glLight($i,fract_position_X) $glLight($i,fract_position_Y) \
			     $glLight($i,fract_position_Z) $glLight($i,fract_position_W)] \
	-spot_dir    [list $glLight($i,spot_dir_X) $glLight($i,spot_dir_Y) $glLight($i,spot_dir_Z)] \
	-spot_exp    $glLight($i,spot_exp) \
	-spot_cutoff $glLight($i,spot_cutoff) \
	-const_atten $glLight($i,const_atten) \
	-lin_atten   $glLight($i,lin_atten)   \
	-quad_atten  $glLight($i,quad_atten)  

    # now do Togl_PostRedisplay(".mesa") to update the change ...

    # t.k.
    # BEWARE: this is dirty and usefull only for present purpose
    if { $togl != {} } {
	if { $togl == ".mesa" } {
	    $togl render
	} else {
	    # force updating the projection and correspondingly lights position
	    cry_toglzoom $togl 0.0
	}
    }
}



proc glLight:loadValues {i} {
    global glLight
    
    set glLight($i) [xc_getGLparam lightmodel -get [list is_light_enabled $i]]

    GetOpenGLPar _LIGHTMODEL_AMBIENT_ \
	glLight(lightmodel_ambient_R) glLight(lightmodel_ambient_G) \
	glLight(lightmodel_ambient_B) glLight(lightmodel_ambient_A)
    GetOpenGLPar _LIGHTMODEL_LOCAL-VIEWER_ \
	glLight(lightmodel_local_viewer) d1 d2 d3
	
    foreach f {ambient diffuse specular} {
	set F [string toupper $f]
	GetOpenGLPar _LIGHT${i}_${F}_ glLight($i,${f}_R) glLight($i,${f}_G) glLight($i,${f}_B) glLight($i,${f}_A)
    }
    
    GetOpenGLPar _LIGHT${i}_FRACT-POSITION_ \
	glLight($i,fract_position_X) glLight($i,fract_position_Y) \
	glLight($i,fract_position_Z) glLight($i,fract_position_W)
    
    GetOpenGLPar _LIGHT${i}_SPOT-DIR_ \
	glLight($i,spot_dir_X) glLight($i,spot_dir_Y) glLight($i,spot_dir_Z) d1
    
    GetOpenGLPar _LIGHT${i}_SPOT-EXP_    glLight($i,spot_exp)    d1 d2 d3     
    GetOpenGLPar _LIGHT${i}_SPOT-CUTOFF_ glLight($i,spot_cutoff) d1 d2 d3
    GetOpenGLPar _LIGHT${i}_CONST-ATTEN_ glLight($i,const_atten) d1 d2 d3 
    GetOpenGLPar _LIGHT${i}_LIN-ATTEN_   glLight($i,lin_atten)   d1 d2 d3   
    GetOpenGLPar _LIGHT${i}_QUAD-ATTEN_  glLight($i,quad_atten)  d1 d2 d3
}



proc glLight:save {} {
    global glLight system
 
    set file [tk_getSaveFile -defaultextension .lights \
		  -filetypes { 
		      {{XCRYSDEN Lights File}  {.lights}}
		      {{All Files}             {.*}     }
		  } \
		  -initialdir  $system(PWD) \
		  -initialfile my_lights.lights \
		  -title       "Save XCRYSDEN Lights File"]
    
    if { $file == "" } {
	return
    }
    
    puts stderr "Lights File: $file"
    puts stderr "Content:     [array get glLight]"
    WriteFile $file [array get glLight] w
}



proc glLight:load {{file {}}} {
    global glLight system
 
    if { $file != "" } {
	if { ! [file exists $file] } {
	    ErrorDialog "file \"$file\" does not exist !!!"
	    return
	}
    } else {
	set file [tk_getOpenFile -defaultextension .lights \
		      -filetypes { 
			  {{XCRYSDEN Lights File}  {.lights}}
			  {{All Files}             {.*}     }
		      } \
		      -initialdir  $system(PWD) \
		      -title       "Load XCRYSDEN Lights File"]
	
	if { $file == "" } {
	    return
	}
    }

    # load the lights

    array set glLight [ReadFile $file]

    # update the lights

    if { ! [info exists glLight(nlights)] } {
	set glLight(nlights) 6
    }
    
    for {set i 0} {$i < $glLight(nlights)} {incr i} {

	if { $glLight($i) } {
	    glLight:update "" $i
	}
    }
}
