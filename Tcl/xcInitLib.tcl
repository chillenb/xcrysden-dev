#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/xcInitLib.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#  Modified by Eric Verfaillie ericverfaillie@yahoo.fr EV                   #
#  may 2004                                                                 #
#  modifications are near EV comments                                       #
#############################################################################

#########################
# CRYSTALXX executables #
#########################
set system(c95_exist)      0
set system(c95_integrals)  crystal
set system(c95_scf)        crystal
set system(c95_scfdir)     crystal
set system(c95_properties) properties

#############
# DEBUGGING #
#############
set xcMisc(debug) 0

###################################
# measurements: distances, angles #
###################################

set select(dist_precision) 4
set select(angl_precision) 3

########################################################################
# some Window Managers doesn't handle proper the "Maxi" command button #
# therefore user may specify necessary offsets in order that Maxi      #
# covers the whole root window starting on upper left corner of the    #
# screen; when I used fvwm2 WM I had to spefify                        # 
# ${rootX}x${rootY}+6+23;                                              #
# Now I am using enlightenment and there is no needed offset. This     #
# definition here is the DEFAULT; each user shoud specify its own      #
# settings in $HOME/.xcrysden/custom-definitions file !!!              #
########################################################################
set xcMisc(wm_rootXshift) 0
set xcMisc(wm_rootYshift) 0

#check_packages_awk

set system(PWD)    [pwd]
set system(PID)    [pid]
set system(TCLDIR) $system(TOPDIR)/Tcl
set system(BMPDIR) $system(TOPDIR)/images
set system(UTIL)   $system(TOPDIR)/util
set system(AWKDIR) $system(TOPDIR)/Awk

if { [info exists env(XCRYSDEN_LIB_BINDIR)] && [file isdirectory $env(XCRYSDEN_LIB_BINDIR)] } {
    # used by system-wide installation made by "make install" or by linux distros
    set system(FORDIR) $env(XCRYSDEN_LIB_BINDIR)
    set system(BINDIR) $env(XCRYSDEN_LIB_BINDIR)
} else {    
    set system(FORDIR) $system(TOPDIR)/bin
    set system(BINDIR) $system(TOPDIR)/bin
}

# make non-writable core
if { ![file exists $system(PWD)/core] } {
    catch {exec touch $system(PWD)/core}
}
catch {exec chmod -w $system(PWD)/core}

#################################################################
# make more save algorithm for assigning the SCRDIR in the future
#
if { [file isdirectory $system(SCRDIR)/xc_$system(PID)] } {
    file delete -force $system(SCRDIR)/xc_$system(PID)
} elseif { [file exists $system(SCRDIR)/xc_$system(PID)] } {
    file delete -force $system(SCRDIR)/xc_$system(PID)
}
file mkdir $system(SCRDIR)/xc_$system(PID)
set system(SCRDIR) $system(SCRDIR)/xc_$system(PID)
#################################################################


#############################################################################
# THIS IS FOR AUTO-LOAD; this is very dengerous if only tclsh has been loaded
lappend auto_path $system(TCLDIR) $system(TCLDIR)/fs

#
# the following lines are executed only if shell is wish
#
if { [string match *Tk* [package names]] } {
    ###########################################################
    # SETTING X RESOURCES - loading database .xcrysden_defaults
    if { [file exists $system(TCLDIR)/Xcrysden_resources] } {
    	option readfile $system(TCLDIR)/Xcrysden_resources 
    }
    if { [file exists $env(HOME)/.xcrysden/Xcrysden_resources] } {
    	option readfile $env(HOME)/.xcrysden/Xcrysden_resources
    }
    # load palette
    set palette [option get . palette Palette]
    if { $palette != "" } {
    	tk_setPalette $palette
    } else {
	tk_setPalette [. cget -background]
    }

    ###########################################################
    # changed: Wed Dec 30 14:33:46 CET 1998
    # now fonts are managed in     ''ProbeResolution''
    ###########################################################
    # FONTS -- FONTS -- FONTS; (AUTO-LOAD must be before FONTS !!!!)
    #button .b -text TONE
    #entry  .e -text TONE
    #	 
    #set xcFonts(normal)       [lindex [.b configure -font] end]
    #set xcFonts(normal_entry) [lindex [.e configure -font] end]
    #set xcFonts(small)        [ModifyFontSize .b 10 \
    #	     {-family helvetica -slant r -weight bold}]
    #set xcFonts(small_entry)  [ModifyFontSize .e 10 \
    #	     {-family helvetica -slant r -weight normal}]
    #destroy .b .e
    	
    ######################################################################
    # COLORS -- COLORS -- COLORS -- COLORS -- COLORS -- COLORS
    set xcColors(disabled_fg) [lindex \
    	    [GetWidgetConfig button -disabledforeground] end]
    set xcColors(enabled_fg)  [lindex \
    	    [GetWidgetConfig scale -foreground] end]
    set xcColors(normal_bg) [lindex \
    	    [GetWidgetConfig button -background] end]
    #xcDebug -stderr "xcColors(normal_bg)=$xcColors(normal_bg)"
    #xcDebug -stderr "xcColors(enabled_fg)=$xcColors(enabled_fg)"
    #xcDebug -stderr "xcColors(disabled_fg)=$xcColors(disabled_fg)"

    ###########################################################
    # CURSORS --- CURSORS --- CURSORS
    set xcCursor(default) [. cget -cursor]
    set xcCursor(watch)   watch
    	
    ####################################
    # set precision to double precission
    #set tcl_precision 17; # this is legacy stuff
    	
    #####################################################
    # CREATE IMAGES --- CREATE IMAGES --- CREATE IMAGES #
    #####################################################
    # rotation images
    set xcMisc(status_init_label) "Creating images ..."

    # K-path selection buttons (no rescaling)

    image create photo printer -format gif \
	-file $system(BMPDIR)/printer.gif
    image create photo rotXmin -format gif \
	-file $system(BMPDIR)/rotXmin.gif
    image create photo rotXplus -format gif \
	-file $system(BMPDIR)/rotXplus.gif
    image create photo rotYmin -format gif \
	-file $system(BMPDIR)/rotYmin.gif
    image create photo rotYplus -format gif \
	-file $system(BMPDIR)/rotYplus.gif
    image create photo rotZmin -format gif \
	-file $system(BMPDIR)/rotZmin.gif
    image create photo rotZplus -format gif \
	-file $system(BMPDIR)/rotZplus.gif
    image create photo zoomUp -format gif \
	-file $system(BMPDIR)/zoomup.gif
    image create photo zoomDown -format gif \
	-file $system(BMPDIR)/zoomdown.gif
    	
    # animation images (no rescaling)

    image create photo first -format gif -file $system(BMPDIR)/first.gif
    image create photo last -format gif -file $system(BMPDIR)/last.gif
    image create photo previous -format gif -file $system(BMPDIR)/previous.gif
    image create photo next -format gif -file $system(BMPDIR)/next.gif
    image create photo backward -format gif -file $system(BMPDIR)/backward.gif
    image create photo forward -format gif -file $system(BMPDIR)/forward.gif    	
    image create photo stop -format gif -file $system(BMPDIR)/stop.gif    	
    image create photo pause -format gif -file $system(BMPDIR)/pause.gif    	
    image create photo snapshot -format gif -file $system(BMPDIR)/camera-photo.gif
	
    # unmapWin images: should be Nx16 dimension (no rescaling)

    image create photo unmap-isosurf -format gif \
	-height 16 -width 80 \
	-file $system(BMPDIR)/unmap-isosurface.gif
    
    image create photo unmap-C95output -format gif \
	-height 16 -width 80 \
	-file $system(BMPDIR)/unmap-C95output.gif
    
    image create photo unmap-empty -format gif \
	-height 16 -width 80 \
	-file $system(BMPDIR)/unmap-empty.gif
    
    image create photo unmap -format gif \
	-height 16 \
	-file $system(BMPDIR)/unmap.gif

    # palette image (no rescaling)

    image create photo colors -format gif -file $system(BMPDIR)/mini-colors.gif

    # nice images (no rescaling)

    #image create photo titleImg -format gif -file $system(BMPDIR)/xcrysden-picture.gif

    # FermiSurface images (no rescaling)

    foreach image {bz cell nocell nocrop solidcell solidwirecell wirecell}  {
	image create photo fs_$image -format gif -file $system(BMPDIR)/$image.gif
    }

    # Viewer Rotate & Orientate buttons

    set b 26; #this is button size
    lappend xcMisc(rescale_image_list) \
	[image create photo up    -format gif \
	     -file $system(BMPDIR)/up1.gif -width $b -height $b]
    lappend xcMisc(rescale_image_list) \
	[image create photo down  -format gif \
	     -file $system(BMPDIR)/down1.gif -width $b -height $b]
    lappend xcMisc(rescale_image_list) \
	[image create photo left  -format gif \
	     -file $system(BMPDIR)/left1.gif -width $b -height $b]
    lappend xcMisc(rescale_image_list) \
	[image create photo right -format gif \
	     -file $system(BMPDIR)/right1.gif -width $b -height $b]
    lappend xcMisc(rescale_image_list) \
	[image create photo center -format gif \
	     -file $system(BMPDIR)/center_nice.gif -width $b -height $b]
    lappend xcMisc(rescale_image_list) \
	[image create photo rotXY    -format gif \
	     -file $system(BMPDIR)/rotXY1.gif -width $b -height $b]
    lappend xcMisc(rescale_image_list) \
	[image create photo rotXZ  -format gif \
	     -file $system(BMPDIR)/rotXZ1.gif -width $b -height $b]
    lappend xcMisc(rescale_image_list) \
	[image create photo rotYZ  -format gif \
	     -file $system(BMPDIR)/rotYZ1.gif -width $b -height $b]
    lappend xcMisc(rescale_image_list) \
	[image create photo rotAB    -format gif \
	     -file $system(BMPDIR)/rotAB.gif -width $b -height $b]
    lappend xcMisc(rescale_image_list) \
	[image create photo rotAC  -format gif \
	     -file $system(BMPDIR)/rotAC.gif -width $b -height $b]
    lappend xcMisc(rescale_image_list) \
	[image create photo rotBC  -format gif \
	     -file $system(BMPDIR)/rotBC.gif -width $b -height $b]

    # displayMode + displayStyle images (rescale image)

    foreach image {
	wireframes_2d
	pointlines_2d
	pipeballs_2d
	ballsticks2_2d
	ballsticks1_2d
	spacefills_2d
	spacefills_3d
	ballsticks_3d
	pipeballs_3d
	sticks_3d
	dm_wire
	dm_solid
	dm_anaglyph    
	dm_stereo
	dm_smooth
	dm_flat
	rep_unit
	rep_asym
    } {
	lappend xcMisc(rescale_image_list) [image create photo $image -format gif -file $system(BMPDIR)/$image.gif]
    }

    ##################################
    # Get Fonts, Images & Other Sizes 
    set xcMisc(status_init_label) "Probing screen dimensions ..."
    ProbeResolution
    SetImageSizes

    # load default options for Togl widgets

    toglOptDefault

    # the below "FoucIn" doesn't work well:
    ## let's select the content of the entry widget upon FocuIn
    #
    #bind Entry <FocusIn> {
    #	%W selection from 0
    #	%W selection to end
    #	%W icursor end
    #}
}
