#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/scripting.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################


# ------------------------------------------------------------------------
#****h* XCRYSDEN/Scripting ***
#
# NAME
# scripting -- Scripting facility (v0.1) of XCRYSDEN
#
# COPYRIGHT
# Anton Kokalj (C) 2003
#
# FUNCTION
#
# Often a mechanism for automating a particular job is desirable. For
# example, one might want to produce several plots of molecular
# orbitals of a given molecule. It would be desirable that the display
# parameters are exactly the same for all plots. By using the
# XCRYSDEN's GUI one would need to open each file and repeat all the
# operations for each molecular orbital, and this is largely
# redundant. Once the display parameters for one plot are determined
# one could produce the others using script which automates the task.
# This is the purpose of the "Scripting" class of functions, namely,
# to provide an XCRYSDEN shell-script mechanism for achieving such
# jobs.
#
# AUTHOR
# Anton Kokalj
#
# CREATION DATE
# Sometime in February 2003
#
# NOTES
#
# At this stage the API of Scripting is not yet completely fixed and
# is used for testing the "scripting" possibility, that is, loading
# the tasks from a Tcl script.
#
#****
# ------------------------------------------------------------------------

namespace eval scripting {
    variable scripting

    set scripting(verbosity) 1
}



# ------------------------------------------------------------------------
#****if* Scripting/scripting::source
#
# NAME
# scripting::source
#
# USAGE
# scripting::source scriptfile
#
# PURPOSE

# This proc is used to source the an XCRYSDEN script. The proc is for
# internal use only, namely, the XCRYSDEN uses this proc to sources 
# the scriptFile.

#
# ARGUMENTS
# scriptfile -- name of scripting file
#
# RETURN VALUE
# Undefined.
#****
# ------------------------------------------------------------------------

proc scripting::source {scriptfile} {
    global scriptFile
    
    set scriptFile $scriptfile
        
    # execute script in global level --
    uplevel \#0 {
	cd $system(PWD)
	SetWatchCursor
	source [gunzipFile $scriptFile]
	ResetCursor
    }
}


#------------------------------------------------------------------------
# This function is used to convert the old xcrysden script (xcrysden
# version < 1.6) to new syntax (xcrysden version >= 1.6)
#------------------------------------------------------------------------
proc scripting::oldscript {scriptfile} {
    global system
    
    set toglCmdsRe {(xc_rotate|xc_translate|xc_B1motion|xc_ShiftB1motion|xc_B2motion|xc_Brelease|xc_dump2eps|cry_gl2psPrintTogl|cry_dump2ppm|cry_toglzoom|xc_realtimemovie|xc_setfont|xc_setatomlabel|xc_clearatomlabel|xc_queryfont)}
    
    WriteFile $system(SCRDIR)/script.xcrysden \
        [regsub -all -- ".mesa $toglCmdsRe" [ReadFile $scriptfile] {\1 .mesa}]

    scripting::source $system(SCRDIR)/script.xcrysden
}

	    
# ------------------------------------------------------------------------
#****f* Scripting/scripting::open
#
# NAME
# scripting::open
#
# USAGE
# scripting::open option file ?option file? ...
#
# PURPOSE

# This proc is used in XCRYSDEN scripts to open structure files. The
# usage is very similar to that of xcrysden command, that is, instead
# of "xcrysden option file" one calls "scripting::open option file".

#
# ARGUMENTS
# args -- the command line arguments (those of xcrysden program)
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::open --xsf file.xsf
#****
# ------------------------------------------------------------------------

proc scripting::open {args} {
    parseComLinArg $args
    update
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::save
#
# NAME
# scripting::save
#
# USAGE
# scripting::save ?format? file
#
# PURPOSE
# This proc is used in XCRYSDEN scripts to save the structure
# information of the current structure. So far only the XSF format is
# supported.
#
# ARGUMENTS
# format   -- the format to be used for saved [OPTIONAL]
# filename -- the name of the file to save into
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::save --xsf file.xsf
#****
# ------------------------------------------------------------------------

proc scripting::save {args} {
    
    set len [llength $args] 
    if { $len < 1 ||  $len > 2 } {
	ErrorDialog "wrong number of arguments, $len. Must be ::scripting::save ?format? filename"
	return
    }

    set fmt xsf
    if { $len == 2 } {
	set fmt [lindex $args 0]
	switch $fmt {
	    xsf { DummyProc }
	    default {
		ErrorDialog "wrong format in  ::scripting::save, $fmt. Must be \"xsf\"."
		return
	    }
	}
    }

    set filename [lindex $args end]
    SaveXSFStruct $filename
    return 1
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::exec
#
# NAME
# scripting::exec
#
# USAGE
# scripting::exec option file ?option file? ...
#
# PURPOSE

# This proc is used inside majorScript of scripting::multiScript.
# Instead of calling "xcrysden option file" one calls "scripting::exec
# option file".

#
# ARGUMENTS
# args -- the command line arguments (those of xcrysden program)
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::exec --xsf file.xsf
#****
# ------------------------------------------------------------------------

proc scripting::exec {args} {
    global env system
    variable scripting

    # support the following: cd $dir; scripting::exec ***
    set system(PWD) [pwd]

    set script "namespace eval ::scripting \{\n"
    if { [regexp -- {^--} [lindex $args 0]] } {
	# we have an option, i.e., --option file
	append script "scripting::open $args\n"
    } else {
	# we have a filter-routine, i.e., scripting::filter file
	append script $args\n
    }

    if { [info exists scripting(minorScript)] } {
	append script $scripting(minorScript)
    }
    append script "\n\}\n"

    # write a script file
    set script_file [file join $system(PWD) xc_script.tcl.$system(PID)]
    WriteFile $script_file $script w    

    switch -- [string tolower $scripting(verbosity)] {
	1 - yes - on {
	    puts stderr "From scripting::exec executing a minor-script::"
	    puts stderr "-----------------------------------------------\n$script"
	}
    }

    # load a script file
    #catch {::exec sh $env(XCRYSDEN_TOPDIR)/xcrysden --script $script_file >& $system(SCRDIR)/log}
    xcCatchExec sh $env(XCRYSDEN_TOPDIR)/xcrysden --script $script_file >& $system(SCRDIR)/log

    # check for a Tcl-error 
    set log [ReadFile $system(SCRDIR)/log]
    if { [string match {*Error in startup script:*} $log] } {
	puts stderr $log
	exit 1
    }

    # delete the script file
    if { [file exists $script_file] } {
	file delete $script_file
    }
}
    


# ------------------------------------------------------------------------
#****f* Scripting/scripting::multiScript
#
# NAME
# scripting::multiScript
#
# USAGE
# scripting::multiScript majorScript minorScript
#
# PURPOSE

# This proc is for multiple-task jobs. For example, one might want to
# produce several plots of molecular orbitals of a given molecule. It
# would be desirable that the display parameters are exactly the same
# for all plots. This proc is provided to facilitate such jobs. In
# majorScript one specifies how many times to run XCRYSDEN and what
# files to load, and in minorScript one specify what to do with this
# files. For example the instruction in minorScript would be how to
# render molecular orbitals for each provided file in majorScript (see
# example).

#
# ARGUMENTS
# majorScript -- a major script
# minorScript -- a minor script
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::multiScript {
#    # ------------------------------------------------------------------------
#    # This is the MAJOR script
#    #
#    # Lets open a structures from several files (one at a time), and 
#    # print them.
#    # ------------------------------------------------------------------------
#
#    foreach file {file1.xsf file2.xsf} {
#	scripting::exec --xsf $file
#       file rename -force print.png $file.png
#    } {
# } {
#    # ------------------------------------------------------------------------
#    # This is the MINOR script
#    #
#    # It renders and print a structure
#    # ------------------------------------------------------------------------
#
#    scripting::lighting On
#    scripting::displayMode3D BallStick
#    scripting::zoom +0.5 
#    scripting::rotate x -90 
#    scripting::rotate y -30
#    scripting::printToFile print.png windowdump    
# }
#
#****
# ------------------------------------------------------------------------

proc scripting::multiScript {majorScript minorScript} {
    variable scripting
    global system

    if { [winfo exists .title] } {
	destroy .title
	update
    }
    set scripting(majorScript) $majorScript
    set scripting(minorScript) $minorScript
    
    #puts stderr "DEBUG> multiScript::"
    #puts stderr "$scripting(majorScript)"
    #puts stderr ""

    cd $system(PWD)
    eval $scripting(majorScript)
    exit 0
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::lighting
#
# NAME
# scripting::lighting
#
# USAGE
# scripting::lighting On|Off
#
# PURPOSE

# This proc toggles the lighing Off/On display mode.

#
# ARGUMENTS
# mode -- the lighting mode (must be On or Off)
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::lighting On
#****
# ------------------------------------------------------------------------

proc scripting::lighting {mode} {
    global light
    switch -glob [string totitle $mode] {
	"On"   { set light On;  Lighting On  }
	"Off"  { set light Off; Lighting Off }
	default {
	    ErrorDialog "wrong lighting-mode $mode, should one On or Off"
	    return
	}
    }
}
  

# ------------------------------------------------------------------------
#****f* Scripting/scripting::display
#
# NAME
# scripting::display
#
# USAGE
# scripting::display on|off coordinate-system|atomic-labels|crystal-cells|
#                           unicolor-bonds|wigner-seitz-cell|molecular-surface
# or
# scripting::display on|off perspective ?near-factor? ?far-factor?
# or
# scripting::display as ball|spacefill covalent|vad-der-waals
# or
# scripting::display as crystal-cells rods|lines
# or
# scripting::display as cell-mode primitive|convetional
# or 
# scripting::display as cell-unit cell|asymm
#
# PURPOSE

# This proc is used to set the display property of various objects,
# such as coordinate system, atomic symbols, perspective projection
# etc. The "scripting::display on|off" command will display or
# not-display a requested object, while "scripting::display as" will
# display particular object as requested.

#
# ARGUMENTS
# mode -- the mode of command, must be one of on|off|as
# what -- what object to display
# args -- the rest of arguments (as requested by particular object)
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::display on coordinate-system
# or
# scripting::display as cell-mode primitive
#****
# ------------------------------------------------------------------------

proc scripting::display {mode what args} {
    global check radio
    switch $mode {
	yes - on  - 1 { set mode 1 }
	no  - off - 0 { set mode 0 }
	as { ; }
	default {
	    ErrorDialog "wrong scripting::display mode $mode, should be on or off"
	    return
	}
    }

    if { $mode == "1" || $mode == "0" } {
	# ------------------------------------------------------------------------
	# scripting::display on|off ...
	# ------------------------------------------------------------------------

	switch -glob [string tolower $what] {
	    "coor*" { 
		set check(crds) $mode; 
		CrdSist 
	    }
	    "atom*" - "label*" { 
		set check(labels) $mode
		AtomLabels
	    }
	    "cryst*" - "cell*" {
		set check(frames) $mode
		CrysFrames
	    }
	    "unicol*" - "unibon*" {
		set check(unibond) $mode
		Unibond
	    }
	    "force*" - {
		set check(forces) $mode 
		forceVectors .mesa
	    }
	    "wigner*" - {
		set check(wigner) $mode
		WignerSeitz
	    }
	    "molec*" - "molsurf*" {
		set check(pseudoDens) $mode
		PseudoDensity
	    }
	    "persp*" {
		set check(perspective) $mode
		if { $args != "" } {
		    # args1 = near-factor; args2 = far-factor
		    set i 0
		    foreach var $args {		    
			set var($i) [lindex $args $i]
			if { ! [string is double $var($i)] } {
			    ErrorDialog "expected double, but got $var($i), while executing scripting::display $mode $what $args"
			}
			if { $i == 0 } {
			    xc_newvalue .mesa $mody(L_PERSPECTIVEFOVY) $var($i)
			} else {
			    xc_newvalue .mesa $mody(L_PERSPECTIVEBACK) $var($i)
			}
			incr i
		    }	
		}	    
		Perspective
	    }
	}

    } elseif { $mode == "as" } {
	# ------------------------------------------------------------------------
	# scripting::display as ...
	# ------------------------------------------------------------------------

	switch -glob [string tolower $what] {	    
	    "ball*" {
		switch -glob [string tolower $args] {
		    "cova*" {
			set radio(ball) "Balls based on covalent radii"
			xc_newvalue .mesa $mody(L_BALL_COV)
		    }
		    "van*" {
			set radio(ball) "Balls based on Van der Waals radii"
			xc_newvalue .mesa $mody(L_BALL_VDW)
		    }		    
		}		    
	    }
	    "spacefill*" {
		switch -glob [string tolower $args] {
		    "cova*" {
			set radio(space) "SpaceFill based on covalent radii"
			xc_newvalue .mesa $mody(L_SPACE_COV)
		    }
		    "van*" {
			set radio(space) "SpaceFill based on Van der Waals radii"
			xc_newvalue .mesa $mody(L_SPACE_VDW)
		    }		    
		}		    
	    }
	    "crystal-c" - "crystal-f" - "frame*" {
		switch -glob [string tolower $args] {
		    "rod*" {
			set radio(frames) rods
			DispFramesAs
		    }
		    "line*" {
			set radio(frames) lines
			DispFramesAs
		    }	    
		}
	    }	    		
	    "cell-s*" - "cell-m*" {
		switch -glob [string tolower $args] {
		    "prim*" {
			set radio(cellmode) prim
			CellMode 1
		    }
		    "conv*" {
			set radio(cellmode) conv
			CellMode 1
		    }
		}
	    }
	    "cell-u*" - "unit of rep*" {
		switch -glob [string tolower $args] {
		    "cell*" {
			set radio(unitrep) cell
			CellMode 1
		    }
		    "asym*" {
			set radio(unitrep) asym
			CellMode 1
		    }
		}
	    }
	}
    }
}
	
    

# ------------------------------------------------------------------------
#****f* Scripting/scripting::displayMode3D
#
# NAME
# scripting::displayMode3D
#
# USAGE
# scripting::lighting sticks|pipe&ball|ballstick|spacefill
#
# PURPOSE

# Switch to requested 3D (i.e. Lighting-On) display mode

#
# ARGUMENTS
# mode -- one of allowed 3D display modes
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::displayMode3D ballstick
#****
# ------------------------------------------------------------------------
	    
proc scripting::displayMode3D {mode} {
    switch -glob [string tolower $mode] {
	"stick*"     { DisplayOver3D S }
	"pipe&ball*" { DisplayOver3D PB }
	"ballstick*" { DisplayOver3D BS }
	"spacefill*" { DisplayOver3D SF }
	default {
	    ErrorDialog "wrong 3D display-mode $mode, should one of Stick, Pipe&Ball, BallStick, or SpaceFill"
	    return
	}
    }
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::displayMode2D
#
# NAME
# scripting::displayMode2D
#
# USAGE
# scripting::lighting wireframe|pointline|pipe&ball|ballstick-1|ballstick-2|spacefill
#
# PURPOSE

# Switch to requested 2D (i.e. Lighting-Off) display mode

#
# ARGUMENTS
# mode -- one of allowed 2D display modes
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::displayMode2D ballstick
#****
# ------------------------------------------------------------------------

proc scripting::displayMode2D {mode} {
    switch -glob [string tolower $mode] {
	"wireframe*"  { Display2D WF }
	"pointline*"  { Display2D PL }
	"pipe&ball*"  { Display2D PB } 
	"ballstick-1" { Display2D BS1 } 
	"ballstick*"  { Display2D BS2 }
	"spacefill*"  { Display2D SF }
	default {
	    ErrorDialog "wrong 2D display-mode $mode, should one of WireFrame, PointLine, Pipe&Ball, BallStick-1, BallStick-2, or SpaceFill"
	    return
	}
    }
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::mainWindow
#
# NAME
# scripting::mainWindow
#
# USAGE
# scripting::mainWindow maximize
# or
# scripting::mainWindow resize width height
#
# PURPOSE

# Resize the XCRYSDEN main window.

#
# ARGUMENTS
# action -- mode od resizing (maximize or resize)
# args   -- width and height for resize mode
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::mainWindow resize 500 500
#****
# ------------------------------------------------------------------------

proc scripting::mainWindow {action args} {
    switch -exact -- $action {
	maximize {
	    global maingeom
	    set maingeom Maxi
	    MainGeom
	}
	
	resize {
	    set len [llength $args]
	    if { $len != 2 } {
		ErrorDialog "wrong number of arguments $len, should be \"scripting::mainWindow resize width height"
		return
	    }	    
	    set w [lindex $args 0]
	    set h [lindex $args 0]
	    foreach num $args {
		if { ! [string is integer $num] } {
		    ErrorDialog "expected integer but got $num, while executing w scripting::mainWindow $action $w $h"
		}
	    }
	    PlaceGlobWin 1 $w $h
	}
	
	default {
	    ErrorDialog "wrong subcommand $action, should be maximize or resize, while executing scripting::mainWindow ..."
	}
    }
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::displayWindow
#
# NAME
# scripting::displayWindow
#
# USAGE
# scripting::displayWindow fullscreen
# or
# scripting::displayWindow resize width height
#
# PURPOSE

# Resize the XCRYSDEN display window of XCRYSDEN's main window.

#
# ARGUMENTS
# action -- mode od resizing (fullscren or resize)
# args   -- width and height for resize mode
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::displayWindow resize 500 500
#****
# ------------------------------------------------------------------------

proc scripting::displayWindow {action args} {
    switch -exact -- $action {
	maximize - fullscreen {
	    maximizeDisplay
	}

	resize {
	    set len [llength $args]
	    if { $len != 2 } {
		ErrorDialog "wrong number of arguments $len, should be \"scripting::displayWindow resize width height"
		return
	    }	    
	    set w [lindex $args 0]
	    set h [lindex $args 1]
	    foreach num $args {
		if { ! [string is integer $num] } {
		    ErrorDialog "expected integer but got $num, while executing  scripting::displayWindow $action $w $h"
		}
	    }
	    PlaceGlobWin display_resize $w $h
	}
	
	default {
	    ErrorDialog "wrong subcommand $action, should be maximize or resize, while executing scripting::displayWindow ..."
	}
    }
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::zoom
#
# NAME
# scripting::zoom
#
# USAGE
# scripting::zoom step ?nsteps?
#
# PURPOSE

# This proc zooms the displayed structure. The optional argument
# nsteps determines the number of zooming steps. If it is omitted only
# one zooming is performed. The unit of step is fractional, unless %
# sign is appended to the step number. In this case the step is in
# percentage unit.

#
# ARGUMENTS
# step  -- the zooming step (typical values between 0 and 1)
# nstep -- number of zooming step (default=1)
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::zoom 0.05
# scripting::zoom 5%
#
#****
# ------------------------------------------------------------------------

proc scripting::zoom {step {ntimes 1}} {
    if { ! [string is integer $ntimes] } {
	ErrorDialog "wanted integer, but got $ntimes, while executing scripting::zoom $step $ntimes"
	return
    }
    if { [string range $step end end] == "%" } {
	set step [string trim $step %]
	if { ! [string is double $step] } {
	    ErrorDialog "wanted double, but got $step, while executing scripting::zoom $step $ntimes"
	}
	set step [expr {$step / 100.0}]
    } elseif { ! [string is double $step] } {
	ErrorDialog "wanted double, but got $step, while executing scripting::zoom $step $ntimes"
    }
    for {set i 0} {$i < $ntimes} {incr i} {
	xc_translate .mesa +z $step
	update
    }
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::rotate
#
# NAME
# scripting::rotate
#
# USAGE
# scripting::rotate x|y|z rotationstep ?nsteps?
# or
# scripting::rotate xy|xz|yz rotstep_dir1 rotstep_dir2 ?nsteps?
#
# PURPOSE

# This proc rotates the displayed objects, either around x, y, or z
# axes (x|y|z modes), or around x-and-y, x-and-z, y-and-z axes
# (xy|xz|yz modes). In the latter case the rotste_dir1 is the rotation
# step for the first axis, while rotstep_dir2 is the rotation step for
# second axis. The optional argument nsteps determines the number of
# rotation steps. If it is omitted only one rotation is performed.

#
# ARGUMENTS
# dir  -- direction of rotation (must be one one x, y, z, xy, xz, or yz)
# args -- the rest of arguments as requested by particular mode
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::rotate x 5
# or
# scripting::rotate x 5 10
# or
# scripting::rotate xy 5 5
# or
# scripting::rotate xy 5 5 10
#****
# ------------------------------------------------------------------------

proc scripting::rotate {dir args} {
        
    switch -- $dir {
	xy - xz - yz {
	    set argc [llength $args]
	    if { $argc < 2 && $argc > 4 } {
		ErrorDialog "wrong number of scripting::rotate ++ arguments"
		return
	    }
	    set step1  [lindex $args 0]
	    set step2  [lindex $args 1]
	    set n_step [lindex $args 2]
	    if { $n_step == "" } {
		set n_step 1
	    }
	    if { ! [string is integer $n_step] } {
		ErrorDialog "wanted integer, but got $ntimes, while executing scripting::rotate $dir $args"
		return	    
	    }
	    for {set i 0} {$i < $n_step} {incr i} {
		xc_rotate .mesa ++$dir $step1 $step2
		update
	    }
	}

	x - y - z {
	    set argc [llength $args]
	    if { $argc < 1 && $argc > 2 } {
		ErrorDialog "wrong number of scripting::rotate arguments"
		return
	    }
	    set step   [lindex $args 0]
	    set n_step [lindex $args 1]
	    if { $n_step == "" } {
		set n_step 1
	    }
	    if { ! [string is integer $n_step] } {
		ErrorDialog "wanted integer, but got $ntimes, while executing scripting::rotate $dir $args"
		return	    
	    }
	    for {set i 0} {$i < $n_step} {incr i} {
		xc_rotate .mesa +$dir $step
		update
	    }
	}
    }
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::printToFile
#
# NAME
# scripting::printToFile
#
# USAGE
# scripting::printToFile ?filename? ?windowdump? ?togl?
#
# PURPOSE

# This proc prints to a file the content of the display window. If
# windowdump is specified, the display window is dumped to a file.In
# the windowdump case make sure that no other window obscures the
# display window, as it will appear on the dump. If togl option is
# omitted the XCRYSDEN main display window is assumed.

#
# ARGUMENTS
# filename   -- the name of the file to print to (if not specified or 
#               void, the filename will be queried)
# windowdump -- if non-void, the a display window dump is performed
# togl       -- the pathName of the display window
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::printToFile print.png
#****
# ------------------------------------------------------------------------

proc scripting::printToFile {{filename {}} {windowdump {}} {togl .mesa}} {
    global xcMisc system printSetup light toglEPS
    
    #
    # check if we have converting program
    #
    if { $filename == "" } {
	set filename [printTogl_queryFilename]
	if { $filename == "" } {
	    return
	}
    }

    cd $system(PWD)
    set ext [file extension [file tail $filename]]
    set EXT [string trimleft [string toupper $ext] .]
    
    SetWatchCursor

    global printSetup
    if { ![info exists printSetup(dumpWindow)] } {
        set printSetup(dumpWindow) 0
    }

    # windowdump == 0 --> no window dumping
    
    if { $windowdump != 0 } {
        if { $windowdump != "" || $printSetup(dumpWindow) } {
            dumpWindow $togl $filename
            return
        }
    }

    set fmt GL2PS_EPS    
    set vectorial 0
    switch -exact -- $EXT {
        PS - EPS - PDF - SVG {
            set vectorial 1
            set fmt GL2PS_$EXT
        }
    }
    
    if { $light == "Off" && $vectorial } {
        #
        # ligting-OFF mode : print vectorial
        #
        if { ! [info exists toglEPS(pointsize)] } {
            set toglEPS(pointsize) 2.0
        }
        if { ! [info exists toglEPS(linewidth)] } {
            set toglEPS(linewidth) 2.0
        }
        cry_gl2psPrintTogl $togl $fmt GL2PS_NO_SORT GL2PS_NONE \
            $toglEPS(pointsize) $toglEPS(linewidth) $filename    
    } else {

	printTogl_Antialias begin
	
	if { $EXT == "PPM" } {
	    cry_dump2ppm $togl $filename
	} else {	
	    if { ! [info exists xcMisc(ImageMagick.convert)] } {
		ErrorDialog "cannot print to [file tail $filename], the variable xcMisc(ImageMagick.convert) is not defined in \$HOME/.xcrysden/custom-definition file."
		return
	    }    
	    cry_dump2ppm $togl $system(SCRDIR)/tmp.ppm
	    scripting::_printToFile_imageConvert $system(SCRDIR)/tmp.ppm $filename
	}

	printTogl_Antialias end

    }
    
    ResetCursor    
}

# ------------------------------------------------------------------------
#****if* Scripting/scripting::_printToFile_imageConvert
#
# NAME
# scripting::_printToFile_imageConvert
#
# USAGE
# scripting::printToFile infile outfile
#
# PURPOSE

# This proc performs the PPM to PNG/GIF/JPEG image conversion. It is used 
# internally by printToFile proc.

#
# ARGUMENTS
# infile  -- name of PPM file to convert (input)
# outfile -- name of converted file (output)
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::_printToFile_imageConvert tmp.ppm print.png
#****
# ------------------------------------------------------------------------

proc scripting::_printToFile_imageConvert {infile outfile} {
    global printSetup xcMisc

    if { ! [info exists printSetup(useOptions) ] } {
	# see the printSetup proc
	set printSetup(useOptions) 0
    }
    if { ! [info exists xcMisc(ImageMagick.convertOptions)] } {
	set xcMisc(ImageMagick.convertOptions) ""
    }
    
    if { [info exists xcMisc(ImageMagick.convert)] } {
        if { $printSetup(useOptions) } {
            eval xcCatchExecReturn $xcMisc(ImageMagick.convert) $xcMisc(ImageMagick.convertOptions) $infile $outfile
        } else {
            eval xcCatchExecReturn $xcMisc(ImageMagick.convert) $infile $outfile
        }
    } else {
        ErrorDialog "cannot print to [file tail $outfile], because ImageMagick's \"convert\" program was not found."
        return
    } 
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::load_myParam
#
# NAME
# scripting::load_myParam
#
# USAGE
# scripting::load_myParam
#
# PURPOSE

# This proc loads the myParam values defined in the script.

#
# ARGUMENTS
# None.
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# set myParam(TESSELLATION) 50.0
# set myParam(ATRAD_SCALE)  1.00
# scripting::load_myParam
#
#****
# ------------------------------------------------------------------------

proc scripting::load_myParam {} {
    ::load_myParam
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::buildCrystal
#
# NAME
# scripting::buildCrystal
#
# USAGE
# scripting::buildCrystal nx ?ny? ?nz?
#
# PURPOSE

# This proc builds (i.e. generate) a crystal, that is, nx-cells in 1st
# direction, ny-cells in second direction, and nz-cells in 3rd
# direction.

#
# ARGUMENTS
# nx -- number of cells in 1st crystallographic direction
# ny -- number of cells in 2nd crystallographic direction
# nz -- number of cells in 3rd crystallographic direction
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::buildCrystal 3 3 3
#
#****
# ------------------------------------------------------------------------

proc scripting::buildCrystal {nx {ny 0} {nz 0}} {
    global periodic nxdir nydir nzdir    
    
    foreach num [list $nx $ny $nz] {
	if { ! [string is integer $num] } {
	    ErrorDialog "expected integer but got $num, while executing scripting::buildCrystal $nx $ny $nz"
	    return
	}
	if { $num < 0 } {
	    ErrorDialog "expected positive integer but got $num, while executing scripting::buildCrystal $nx $ny $nz"
	    return
	}
    }
    
    if { $periodic(dim) == 0 } { return }
    if { $periodic(dim) <  2 } { set ny 0 }
    if { $periodic(dim) <  3 } { set nz 0 }
    
    set nxdir $nx
    set nydir $ny
    set nzdir $nz
    
    GenGeomDisplay 1
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::chdir
#
# NAME
# scripting::chdir -- an enhanced version of cd
#
# USAGE
# scripting::chdir newdir
#
# PURPOSE

# This proc is an enhancement of standard "cd" routine. It changes the
# directory to $newdir and updates system(PWD) variable.

#
# ARGUMENTS
# newdir -- new directory to chdir into
#
# RETURN VALUE
# The absolute path of $newdir.
#
# EXAMPLE
# scripting::chdir ../
#
#****
# ------------------------------------------------------------------------

proc scripting::chdir {newdir} {
    global system
    cd $newdir
    set system(PWD) [pwd]
    return $system(PWD)
}
