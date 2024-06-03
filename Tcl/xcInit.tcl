#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/xcInit.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

###############################################################################
# xcrysden TEMPORARY FILES:
# 
# xc_gengeom.$$  --- produced with gengeom with M1_INFO
# xc_struc.$$    --- produced with gengeom with .not. M1_INFO -> used by xcrys
# xc_str2xcr.$$  --- produced by str2xcr program (converted WIEN2k struct file
#                                                 to "xsf" format)
# xc_wnstr.struct--- temporary WIEN2k struct file
# xc_wienstruct.$$-- wien struct file (copy of original WIEN2k struct file)
# xc_inp.$$      --- crystal95's input (also xc_tmp.$$)
# xc_output.$$   --- crystal95's output
# xc_tmp.$$      --- something insignificant  
# xc_bin.$$      --- binary file produced by "xcrys" interpreter; used 
#                    for isosurface evaluation
# xc_binVrt.$$   --- binary file produced by "xcrys" interpreter; used
#                    to stote the current gridvertices
# xc_datagrid.$$ --- binary file where datagrid from DATAGRIDXD is stored
# xc_ndsfp.$$    --- Wigner-Seitz nodes:: nodes for primitive cell mode
# xc_ndsfc.$$    --- Wigner-Seitz nodes:: nodes for conventional cell mode
# xc_klist.$$    --- k-list-file for the kPath program
# xc_rho.$$      --- rho-3Ddatagrid file for WIEN
# here are variables that need to be initialised more than once

proc InitGlobalVar {} {
    global species nxdir nydir nzdir \
	periodic \
	inp n_groupsel groupsel AdvGeom \
	prop \
	XCState \
	dispmode mode2D \
	geng \
	ws \
	colSh \
	sInfo \
	system \
	isoControl \
	openGL \
	undoMenu undoAdvGeom \
	pDen \
	colSh \
	xsfAnim gifAnim \
	light \
	glLight \
	atomLabel \
        toglEPS

    #
    # array xcMisc will be used for various things; it is to prevent to many
    # global variables
    #
    # xcMisc(titlefile) ... name of the file that appears on the main window
    #                       title

    set species {}

    set nxdir 1
    set nydir 1
    set nzdir 1

    # PERIODIC global variable -> data about dimensionality goes here
    # periodic(dim) ... what is the dimensionality of the system
    set periodic(dim)    0 
    set periodic(igroup) 1;    #igroup according to gengeom program

    set dispmode(style)  3D
    set mode2D(WF)  Off
    set mode2D(PL)  Off
    set mode2D(PB)  Off
    set mode2D(BS1) Off
    set mode2D(BS2) Off
    set mode2D(SF)  Off

    set arraylist [list inp AdvGeom XCState prop]

    set varlist [list n_groupsel groupsel]

    foreach array $arraylist {
	if [array exists $array] { unset $array }
    }
    foreach var $varlist {
	if [info exists $var] { unset $var }
    }

    #default #3 gengeom's argument is 1
    if { ! [info exists system(c95_version)] } { set system(c95_version) none }
    set geng(M3_ARGUMENT) [GetGengM3Arg ANGS $system(c95_version)]

    # what to do with IsoXXXX ???
    # close all other toplevels ???

    # WIGNER-SEITZ CELL
    set ws(not_config_yet) 1

    # initialization of mody array
    ModConst

    # this is associated with XCRYSDEN INFO RECORD (core xcrys interpreter)
    if { [info exists sInfo] } { unset sInfo }

    if ![info exists colSh] { 
	set colSh(scheme) atomic
	set colSh(slab_fractional) 1
	set colSh(slab_dir) "-z"
	set colSh(slab_colbas) monochrome
	set colSh(slab_coltyp) combined
	set colSh(slab_alpha) 0.65
	set colSh(dist_x) 0.0
	set colSh(dist_y) 0.0
	set colSh(dist_z) 0.0
	set colSh(dist_r) 1.0
	set colSh(dist_colbas) monochrome
	set colSh(dist_coltyp) combined
	set colSh(dist_alpha) 0.65
    }

    # this should be the same as ISOLINE_MAXLEVEL in isosurf.h file
    set isoControl(max_allowed_2Dnisoline) 100

    ################################
    # initialization of openGL array
    ################################
    set openGL(src_blend_list) { 
	GL_ZERO 
	GL_ONE 
	GL_DST_COLOR 
	GL_ONE_MINUS_DST_COLOR 
	GL_SRC_ALPHA 
	GL_ONE_MINUS_SRC_ALPHA 
	GL_DST_ALPHA 
	GL_ONE_MINUS_DST_ALPHA 
	GL_SRC_ALPHA_SATURATE 
    }
    
    set openGL(dst_blend_list) {
	GL_ZERO
	GL_ONE
	GL_SRC_COLOR
	GL_ONE_MINUS_SRC_COLOR
	GL_SRC_ALPHA
	GL_ONE_MINUS_SRC_ALPHA
	GL_DST_ALPHA
	GL_ONE_MINUS_DST_ALPHA 
    }


    xc_setGLparam frontface -what isosurf_one -frontface CCW
    xc_setGLparam frontface -what isosurf_pos -frontface CCW
    xc_setGLparam frontface -what isosurf_neg -frontface CCW
    
    set pos [xc_getGLparam frontface -what isosurf_pos]
    set neg [xc_getGLparam frontface -what isosurf_neg]
    #xcDebug -stderr "Pos: $pos"
    #xcDebug -stderr "Neg: $neg"
    set openGL(isoside_pos) [lindex $pos 0]
    set openGL(isoside_neg) [lindex $neg 0]

    set openGL(front_ambient_R) 0
    set openGL(front_ambient_G) 0
    set openGL(front_ambient_B) 0
    
    set openGL(front_diffuse_R) 0
    set openGL(front_diffuse_G) 0
    set openGL(front_diffuse_B) 0
    
    set openGL(front_specular_R) 0
    set openGL(front_specular_G) 0
    set openGL(front_specular_B) 0
    
    set openGL(front_emission_R) 0
    set openGL(front_emission_G) 0
    set openGL(front_emission_B) 0
    
    ####
    
    set openGL(back_ambient_R) 0
    set openGL(back_ambient_G) 0
    set openGL(back_ambient_B) 0
    
    set openGL(back_diffuse_R) 0
    set openGL(back_diffuse_G) 0
    set openGL(back_diffuse_B) 0
    
    set openGL(back_specular_R) 0
    set openGL(back_specular_G) 0
    set openGL(back_specular_B) 0
    
    set openGL(back_emission_R) 0
    set openGL(back_emission_G) 0
    set openGL(back_emission_B) 0

    set undoMenu(active_fg)   #ffffff
    set undoMenu(active_bg)   #0000ff
    set undoMenu(default_fg)  #000000
    set undoMenu(default_bg)  #ffffff
    
    set undoAdvGeom(start_index)        0
    set undoAdvGeom(current_index)      0
    set undoAdvGeom(list)               {}
    set undoAdvGeom(redo_start_index)   0
    set undoAdvGeom(redo_current_index) 0
    set undoAdvGeom(redo_list)          {}

    if [info exists pDen(nsurface)] { unset pDen(nsurface) }
    set pDen(type)         gauss
    set pDen(radius)       cov
    set pDen(level)        1.0
    set pDen(cutoff)       1.0
    set pDen(colorscheme)  atomic
    set pDen(drawstyle)    wire
    set pDen(surfacetype)  molsurf
    set pDen(resolution)   0.35
    set pDen(smoothsteps)  0
    set pDen(smoothweight) 0.2
    set pDen(transparent)  0
    set pDen(shademodel)   smooth
    set pDen(monocolor)    {0.8 0.8 0.2}
    set pDen(tessellation) cube
    set pDen(normals)      gradient

    set pDen(t_type)         GAUSSIAN
    set pDen(t_radius)       {Covalent radii}
    set pDen(t_colorscheme)  {Atomic colors}
    set pDen(t_drawstyle)    Wire
    set pDen(t_surfacetype)  {Pseudo density}
    set pDen(t_shademodel)   Smooth
    set pDen(t_tessellation) Cube
    set pDen(t_normals)      Gradient

    set colSh(slabrange_min) 0.00
    set colSh(slabrange_max) 1.00

    set xsfAnim(not_anim) 0

    set gifAnim(create)              0
    set gifAnim(gif_transp)          0
    set gifAnim(gif_minimize)        0
    set gifAnim(gif_global_colormap) 0
    set gifAnim(edit_param)          1
    set gifAnim(movie_format)        mpeg
    set gifAnim(temp_files_dir)      tmp
    set gifAnim(frame_files_format)  PPM
    set gifAnim(ntime_first_frame)   1
    set gifAnim(ntime_last_frame)    1
    set gifAnim(delay)               10
    set gifAnim(loop)                0
    set gifAnim(make_gifAnim)        0


    set atomLabel(fontBrowser)            "Simple Font Browser" 
    set atomLabel(globalFont)             ""
    set atomLabel(globalFont.brightColor) {1.0 1.0 1.0}
    set atomLabel(globalFont.darkColor)   {0.0 0.0 0.0}
    set atomLabel(globalFont.do_display)  1

    set atomLabel(atomFont)             ""
    set atomLabel(atomFont.id)          ""
    set atomLabel(atomFont.label)       ""
    set atomLabel(atomFont.brightColor) {1.0 1.0 1.0}
    set atomLabel(atomFont.darkColor)   {0.0 0.0 0.0}
    set atomLabel(atomFont.do_display)  1

    set light On
    #if { ! [info exists glLight(nlights)] } {
    #	set glLight(nlights) 6
    #}    

    set toglEPS(DRAW_BACKGROUND)    0
    set toglEPS(SIMPLE_LINE_OFFSET) 0
    set toglEPS(SILENT)             0
    set toglEPS(BEST_ROOT)          1
    set toglEPS(OCCLUSION_CULL)     1
    set toglEPS(NO_TEXT)            0
    set toglEPS(LANDSCAPE)          0
    set toglEPS(NO_PS3_SHADING)     0
    set toglEPS(NO_PIXMAP)          0
    set toglEPS(NO_BLENDING)        0
}


proc xcInit {} {
    #global system Const geng

    # take care of the scratch directory
    if { ! [file isdirectory $system(SCRDIR)] } {
	puts stderr "ERROR: SCRATCH directory \"$system(SCRDIR)\" does not exist"
	exit 1
    }
    if { ! [file writable $system(SCRDIR)] } {
        error "SCRATCH directory \"$system(SCRDIR)\" is not writeable"
        exit 1
    }
    
    set system(USER_DIR) $system(TOPDIR)	

    ######################################################################
    # initialize XCRYSDEN'S LIBRARY: variables needed for keeping xcrysden
    # alive will be loaded here
    set xcMisc(status_init_label) "Initializing library ..."
    source $system(TOPDIR)/Tcl/xcInitLib.tcl 
    
    ##########################################################
    # now read the user-custom file; USER MAY CHANGE SOMETHING
    # if ![user-custom-file present] 
    if { [file exists $env(HOME)/.xcrysden/custom-definitions] } {
	source $env(HOME)/.xcrysden/custom-definitions
    } else {
	source $system(TOPDIR)/Tcl/custom-definitions
    }
    # backward compatibility; now in custom-definitions we use 
    # xcMisc(printCommand) for consistency
    if { [info exists xcMisc(printCommand)] } {
	set printCanvas(printCommand) $xcMisc(printCommand)
    }

    #if $xcMisc(debug) {
    #	#debug
    #	#lappend auto_path $xcMisc(dev_dir)
    #}
    
    # ------------------------------------------------------------------------
    # GOTO $system(SCRDIR)    
    cd $system(SCRDIR)
    # ------------------------------------------------------------------------
    # make the 444 core file on $system(SCRDIR)
    catch {exec touch core}
    catch {exec chmod 444 core}
    
    
    ###########################
    # make some subdirectories; so far we need just 1
    if { [file isdirectory $system(SCRDIR)/dir1] } {
	file delete -force $system(SCRDIR)/dir1
    }
    file mkdir $system(SCRDIR)/dir1
    set system(SCRDIR_1) $system(SCRDIR)/dir1
    
    
    set xcMisc(status_init_label) "Checking packages ..."	
    #######################################################################
    # check software packages    
    
    check_package_awk
    check_package_terminal
    check_package_crystal
    find_package_imagemagick
    find_package_gifsicle
    find_package_whirlgif
    find_package_mencoder
    find_package_ppmtompeg
    find_package_babel
    find_package_xwd
    find_package_xsf2_manipulator
    
    determine_movie_encoders

    #######################################################################
    
    ###################
    # FORTRAN UNIT NAME
    set system(ftn_name) [FtnName]
    
    #################################################
    #             IMPORTANT CONSTANTS
    # -----------------------------------------------
    set Const(bohr)      0.529177;  # conversion factor for Ang/Bohr in c95
    
    ##################################################
    #       INITIALIZATION OF GLOBAL VARIABLES
    # ------------------------------------------------
    set xcMisc(status_init_label) "Creating initializing variables ..."
    InitGlobalVar
    
    ##################################################
    # load atom names
    AtomNames
    
    ###################################################
    #        THIS IS FOR GENGEOM PROGRAM
    # usage of "gengeom" program:
    # 
    # gengeom  MODE1  MODE2  MODE3  IGRP  NXDIR  NYDIR  NZDIR  OUTPUT INPUT
    #    0       1      2      3      4     5      6      7      8       9
    #
    set geng(M1_INFO)       0; #INFO mode
    set geng(M1_PRIM)       1; #PRIMITIV CELL; in case of H/R PARAPIPEDAL SHAPE
    set geng(M1_CONV)       2; #CONVENTIONAL CELL; in case of H/R it is PARAP. SHAPE
    set geng(M1_HEXA_SHAPE) 3; #THREEPLE CELL for H/R; HEXAGONAL SHAPE
    set geng(M1_PRIM3)      4; #PRIMITIV cell for H/R; HEXAGONAL SHAPE
    
    set geng(M2_CELL)         1; #CELL is unit of repetition
    set geng(M2_TR_ASYM_UNIT) 2; #translation asymetric unit is unit of repetition
    #default gengeom's #3 argument
    set geng(M3_ARGUMENT)     [GetGengM3Arg ANGS $system(c95_version)]
    
    set geng(IGRP_HEXA)       8; # hexagonal  groups
    set geng(IGRP_TRIG)       9; # trigonal groups
    
    wm iconbitmap . @$system(BMPDIR)/xcrysden.xbm
    #wm iconmask . @$system(BMPDIR)/xcrysden_mask.xbm
    
    set xcMisc(status_init_label) "Building GUI ..."
}





###################################################################
###  MAIN --- MAIN --- MAIN --- MAIN --- MAIN --- MAIN --- MAIN ### 
###################################################################

package provide Tk ; #puts stderr tk_version=$tk_version

# ------------------------------------------------------------------------
# First process the "argc/argv". The order of arguments is
# XCRYSDEN_TOPDIR XCRYSDEN_SCRATCH and the user specified command line
# options
# ------------------------------------------------------------------------

set system(TOPDIR) [lindex $argv 0]
set system(SCRDIR) [lindex $argv 1]

# ------------------------------------------------------------------------
# Load xcrys.dll
# ------------------------------------------------------------------------
if { [file exists $system(TOPDIR)/bin/xcrys.dll] } {    
    load $system(TOPDIR)/bin/xcrys.dll
} elseif { [info exists env(XCRYSDEN_LIB_BINDIR)] && [file exists $env(XCRYSDEN_LIB_BINDIR)/xcrys.dll] } {
    load $env(XCRYSDEN_LIB_BINDIR)/xcrys.dll
}

# ------------------------------------------------------------------------
# some dirty fixes needed to get the program work under CYGWIN
# ------------------------------------------------------------------------
source $system(TOPDIR)/Tcl/cygwin.tcl

# ------------------------------------------------------------------------
# Load Bwidgets package
# ------------------------------------------------------------------------

set bwidget [glob -nocomplain $system(TOPDIR)/external/lib/bwidget-*]
if { $bwidget != "" } {
    set BWidget_dir $bwidget
    lappend auto_path  $BWidget_dir    
}
package require BWidget

# ------------------------------------------------------------------------
# take care of trace/fulltrace utility
# ------------------------------------------------------------------------

if { [info exists env(XCRYSDEN_TRACE)] || [info exists env(XCRYSDEN_FULLTRACE)] } {
    # BWidget needs some special treatment if XCTR is enabled
    foreach file [glob $BWidget_dir/*.tcl] {
	if { $file != "$BWidget_dir/pkgIndex.tcl" && [file exists $file] } {
	    source $file
	}
    }

    # take care of XCTR
    set xctr(recording) 0
    source [file join $system(TOPDIR) Tcl xctr.tcl]
}

#source $system(TOPDIR)/Tcl/parseComLinArg.tcl


# ------------------------------------------------------------------------
# Palette
# ------------------------------------------------------------------------

# this was the 0.3 palette:
#tk_setPalette "#ddd" 
#tk_setPalette "#b5b193"    
#tk_setPalette "#ee9"
#tk_setPalette "#bbb"


# ------------------------------------------------------------------------
# check for stale files in $env(XCRYSDEN_SCRATCH)
# ------------------------------------------------------------------------
proc clean_xcrysden_scratch {} {
    global env
    
    set time [clock seconds]
    set five_days [expr 5*24*3600]
    foreach item [glob -nocomplain -directory $env(XCRYSDEN_SCRATCH) xc_*] {
	set item_time [file atime $item]
	if { [expr {$time - $item_time}] > $five_days } {
	    lappend stale_item $item
	}
    }
    
    if { [info exists stale_item] > 0 } {
	set respond [tk_messageBox -parent . -type yesno -default yes -icon question \
			 -title "Cleaning XCRYSDEN_SCRATCH ?" \
			 -message "Possible stale files older then 10 days exist in XCRYSDEN_SCRATCH directory.\n\nDo you want to delete them?"]

	if { $respond == "yes" } {
	    puts stderr "*** cleaning XCRYSDEN_SCRATCH directory: $stale_item"
	    catch { eval file delete -force $stale_item }
	}
    }
}
clean_xcrysden_scratch


# ------------------------------------------------------------------------
# Welcome images
# ------------------------------------------------------------------------

image create photo kpath -format gif \
    -file $system(TOPDIR)/images/xcrysden_kpath.gif
image create photo welcome -format gif \
    -file $system(TOPDIR)/images/xcrysden-welcome.gif -width 480 -height 320


# ------------------------------------------------------------------------
# make a WELCOME window
# ------------------------------------------------------------------------

proc centerWelcome {thisWin} {
    set w  500
    set h  350
    # root window height/width
    set rh [winfo screenheight $thisWin]     
    set rw [winfo screenwidth $thisWin]
    
    set reqX [expr {($rw-$w)/2}]
    set reqY [expr {($rh-$h)/2}]
    
    wm geometry $thisWin +${reqX}+${reqY}
}

puts stderr "Running on platform : $xcrys(platform)"
puts stderr "   Operating system : $tcl_platform(os)"

if { "[lindex $argv 2]" != "--quiet" } {

    # MACOSX has problems with wm iconify/wm deiconify requests, so
    # don't use them

    if { $tcl_platform(os) != "Darwin" } {
	#catch {wm iconify .}
	catch {wm withdraw .}
	toplevel .title 
	frame .title.f -relief flat -bd 0 -bg #fff
	label .title.f.l -image welcome -anchor center -relief flat -bd 0
	set xcMisc(status_init_label) "Initializing ..."
	label .title.f.l2 -textvariable xcMisc(status_init_label) \
	    -relief flat -bd 0
	pack  .title.f 
	pack .title.f.l .title.f.l2 -side top -fill both -padx 0m -pady 0m
	
	centerWelcome .title
	
	wm overrideredirect .title true 
	update
    }
}
    

# ------------------------------------------------------------------------
# load necessary initialization
# ------------------------------------------------------------------------

eval [info body xcInit]


# ------------------------------------------------------------------------
# start recording the tracing
# ------------------------------------------------------------------------

set xctr(recording) 1

# ------------------------------------------------------------------------
# do renaming such that exit is associated with clean_exit that cleans the TMPDIR
# ------------------------------------------------------------------------
rename exit exit_tcl
rename clean_exit exit

# provide a custom "cd" when in debug mode to trace changing directories

if { [info exists env(XCRYSDEN_DEBUG)] } {
    rename cd cd_tcl

    proc cd_debug {dir} {
	puts stderr "*** "
	puts stderr "*** cd into $dir"
	puts stderr "*** "
	cd_tcl $dir
    }
    
    rename cd_debug cd
}

# create a photo image "container" for printing into PPM files
# (see file xcAppInit.c & togl_ppm.c)

image create photo $xcrys(print.image)

## simple variable debugger
#proc tracer {varname args} {
#    upvar #0 $varname var
#    puts stderr "$varname was updated to be \"$var\"
#          from 
#          [info level 0]
#          from
#          [info level 1]
#          from
#          [info level 2]"
#}
#trace add variable xcMisc(titlefile) write "tracer xcMisc(titlefile)"

# ------------------------------------------------------------------------
# parse command-line options or simple start the Viewer
# ------------------------------------------------------------------------

if { [llength $argv] > 2 } {
    parseComLinArg [lrange $argv 2 end]
} else {
    ViewMol .
}
