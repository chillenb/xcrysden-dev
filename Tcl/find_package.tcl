#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/__file__
# ------                                                                    #
# Copyright (c) 2008 by Anton Kokalj                                        #
#############################################################################

proc find_package {package packageVar recognizeCmd} {
    # package -- name or fullpath for package to query/find
    # packageVar -- where to store the found fullpath of package
    # recognizeCmd -- procedure that checks if the found package is the proper one

    upvar #0 $packageVar pkgVar

    # check first if user has explicitly set the package by the packageVar
    if { [info exists pkgVar] } {
	if { [set path [check_package $pkgVar $recognizeCmd]] != "" } {
	    # package $pkgVar is OK
	    set pkgVar $path
	    return
	}
    }

    # try to guess the fullpath of package
    if { [set path [check_package $package $recognizeCmd]] != "" } {
	set pkgVar $path
    } else {
	if { [info exists pkgVar] } {
	    unset pkgVar
	}
    }
}


proc check_package {package recognizeCmd} {

    set path [auto_execok $package]   

    if { $path != "" } {
	if { [$recognizeCmd $path] } {
	    return $path
	}	   
    } else {
	puts stdout "can't find package: $package"
    }
    
    return ""
}



#
# Find: ImageMagick
#
proc find_package_imagemagick {} {
    global xcMisc

    # convert

    find_package convert xcMisc(ImageMagick.convert) query_imagemagick
    if { [info exists xcMisc(ImageMagick.convert)] } { 
	# variable xcMisc(convert) will be used for animated GIF encoding
	set xcMisc(convert) $xcMisc(ImageMagick.convert)
	puts stderr "Package ImageMagick's convert: $xcMisc(ImageMagick.convert)" 
    }

    # import

    find_package import xcMisc(ImageMagick.import) query_imagemagick
    if { [info exists xcMisc(ImageMagick.import)] } { 
	puts stderr "Package ImageMagick's import: $xcMisc(ImageMagick.import)" 
    }
}


proc query_imagemagick {package} {

    catch {set result [exec -- $package -version]}
    
    if { ! [info exists result] } {
	return 0
    }

    if { [string match -nocase *imagemagick* $result] } {	
	return 1
    }

    return 0
}


#
# Find: Gifsicle
#
proc find_package_gifsicle {} {
    global xcMisc
    find_package gifsicle xcMisc(gifsicle) query_gifsicle
    if { [info exists xcMisc(gifsicle)] } { puts stderr "Package Gifsicle: $xcMisc(gifsicle)" }
}
proc query_gifsicle {package} {
    # gifsicle tends to be rather unique name, we don't permorm any checkings
    return 1
}


#
# Find: Whirlgif
#
proc find_package_whirlgif {} {
    global xcMisc
    find_package whirlgif xcMisc(whirlgif) query_whirlgif
    if { [info exists xcMisc(whirlgif)] } { puts stderr "Package Whirlgif: $xcMisc(whirlgif)" }
}
proc query_whirlgif {package} {
    # whirlgif tends to be rather unique name, we don't permorm any checkings
    return 1
}


#
# Find: Mencoder
#
proc find_package_mencoder {} {
    global xcMisc
    find_package mencoder xcMisc(mencoder) query_mencoder
    if { [info exists xcMisc(mencoder)] } { puts stderr "Package Mencoder: $xcMisc(mencoder)" }
}
proc query_mencoder {package} {
    # mencoder seems to be rather unique name, we don't permorm any checkings
    return 1
}


#
# Find: Ppmtompeg
#
proc find_package_ppmtompeg {} {
    global xcMisc
    # backward compatibility ppmtompeg --> mpeg_encode
    if { [info exists xcMisc(mpeg_encode)] } {
	set xcMisc(ppmtompeg) $xcMisc(mpeg_encode) 
    }
    find_package ppmtompeg xcMisc(ppmtompeg) query_ppmtompeg
    if { [info exists xcMisc(ppmtompeg)] } { puts stderr "Package Ppmtompeg: $xcMisc(ppmtompeg)" }
}
proc query_ppmtompeg {package} {
    # ppmtompeg seems to be rather unique name, we don't permorm any checkings
    return 1
}



#
# Find: Babel
#
proc find_package_babel {} {
    global xcMisc
    find_package babel xcMisc(babel) query_babel
    if { [info exists xcMisc(babel)] } { puts stderr "Package Babel: $xcMisc(babel)" }
}


proc query_babel {package} {
    catch {set result [exec $package -H]}
    
    if { ! [info exists result] } {
	return 0
    }

    if { [string match -nocase *babel* $result] } {	
	return 1
    }

    return 0
}


#
# Find: Xwd
#
proc find_package_xwd {} {
    global xcMisc
    find_package xwd xcMisc(xwd) query_xwd
    if { [info exists xcMisc(xwd)] } { puts stderr "Package Xwd: $xcMisc(xwd)" }
}
proc query_xwd {package} {
    # xwd tends to be rather unique name, we don't perform any checkings
    return 1
}


#
# Find: Xsf2_manipulator
#
proc find_package_xsf2_manipulator {} {
    global xcMisc
    find_package xsf2_manipulator xcMisc(xsf2_manipulator) query_xsf2_manipulator
    if { [info exists xcMisc(xsf2_manipulator)] } { puts stderr "Package Xsf2_manipulator: $xcMisc(xsf2_manipulator)" }
}
proc query_xsf2_manipulator {package} {
    # xsf2_manipulator tends to be rather unique name, we don't perform any checkings
    return 1
}
