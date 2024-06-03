#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/scriptingAtomicLabels.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# ------------------------------------------------------------------------
#****c* Scripting/scripting::atomicLabels
#
# NAME
# scripting::atomicLabels
#
# PURPOSE

# This namespace provide the scripting interface for editing
# atomic-labels (fonts, font-colors, label-texts). There are two
# levels of atomic-labels, so-called, global and custom. The global
# labels are kind of default labels with the text of atomic-symbols.
# The global labels are all those atomic-labels, that were not yet
# explicitly edited on per-atom basis. The custom labels are those
# that were edited on per-atom basis. 
#
# Atomic labels have two different font colors, so-called,
# bright-color and dark-color. The bright font color is used for all
# Lighting-On display modes, and for WireFrame-like Lighting Off
# display modes. On the other hand the dark-color is used for
# Lighting-Off BallStick and SpaceFill display modes.

#
# COMMANDS

# -- scripting::atomicLabels::global
# Edits the fonts and font-colors of global atomic labels.

#
# --scripting::atomicLabels::atomID
# Edits the fonts, font-colors of custom atomic label on per-atom basis.

#
# -- scripting::atomicLabels::atomType
# Edits the fonts, font-colors of custom atomic label on per-atom type basis.

#
# -- scripting::atomicLabels::render
# Toggles the display of global or custom atomic labels

#
# -- scripting::atomicLabels::clear
# Clears the custom label according to specified mode, which is either
# all, atomID, or atomType.
#
# SOURCE


namespace eval scripting::atomicLabels {
    variable atomicLabels
    
    # --------------------------------------------------
    # atomicLabels array-variable has these elements
    # --------------------------------------------------
    #
    # atomicLabels(global.font)        
    # atomicLabels(global.brightcolor) 
    # atomicLabels(global.darkcolor)   
    #
    # atomicLabels(custom.label)       
    # atomicLabels(custom.font)        
    # atomicLabels(custom.brightcolor) 
    # atomicLabels(custom.darkcolor)   
    #
    # --------------------------------------------------
    # in addition it also contains the following elements:
    # --------------------------------------------------
    #
    #  atomicLabels(atomID.*)
    #
    # where "*" resperesnt the same elements as for custom.*    
}

#****
# ------------------------------------------------------------------------


proc scripting::atomicLabels::_loadDefaults {} {
    variable atomicLabels
    
    set atomicLabels(global.font)        ""
    set atomicLabels(global.brightcolor) ""
    set atomicLabels(global.darkcolor)   ""

    set atomicLabels(custom.label)       ""
    set atomicLabels(custom.font)        ""
    set atomicLabels(custom.brightcolor) ""
    set atomicLabels(custom.darkcolor)   ""
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::atomicLabels::global
#
# NAME
# scripting::atomicLabels::global
#
# USAGE
# scripting::atomicLabels::global \
#     -xfont       font  \
#     -tkfont      tkfontname \
#     -brightcolor rgb \
#     -darkcolor   rgb 

#
# PURPOSE

# This proc edits the fonts (i.e. font-type and size) and font-colors
# of global atomic labels. One can specify either XLFD font name, i.e.,
# those that look as: -*-itc bookman-demi-r-normal-*-13-*-*-*-*-*-*-*, or
# Tk font name. The Tk font names are those created by the "font create ..."
# Tk command. 

#
# ARGUMENTS
# args -- various configuration "-option value" pairs (see Options)
#
# OPTIONS
# ------------------------------------------------------------------------
#  OPTION        ALLOWED-VALUES + DESCRIPTION
# ------------------------------------------------------------------------
#  -xfont        xfontname  --> X11's XLFD font name
#
#  -tkfont       tkfontname --> Tk font name, previously create by the
#                               "font create ..." Tk command. If both -xfont 
#                               and -tkfont options are specified the value 
#                               of -tkfont option is taken/
#
#  -brightcolor  rgb        --> the red-green-blue specification of bright 
#                               font-color color (see scripting::atomicLabels 
#                               for the description) of bright and dark font 
#                               colors). Color can be specified either in 
#                               hexadecimal forms (i.e. #rrggbb) or as a list
#                               of three floats in range [0,1] (i.e. 
#                               {0.5 0.5 0.5}).
#      
#  -darkcolor    rgb        --> the red-green-blue specification of dark
#                               font-color color (see scripting::atomicLabels 
#                               for the description) of bright and dark font 
#                               colors). Color can be specified either in 
#                               hexadecimal forms (i.e. #rrggbb) or as a list
#                               of three floats in range [0,1] (i.e. 
#                               {0.5 0.5 0.5}).
# ------------------------------------------------------------------------
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::atomicLabels::global \
#     -xfont       -*-itc bookman-demi-r-normal-*-13-*-*-*-*-*-*-* \
#     -brightcolor #ffffff \
#     -darkcolor   #000000
#
#****
# ------------------------------------------------------------------------

proc scripting::atomicLabels::global {args} {
    variable atomicLabels

    _parseOptions global $args

    xc_setfont .mesa \
	$atomicLabels(global.font) \
	$atomicLabels(global.brightcolor) \
	$atomicLabels(global.darkcolor)
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::atomicLabels::atomID
#
# NAME
# scripting::atomicLabels::atomID
#
# USAGE
# scripting::atomicLabels::atomID atomID \
#     -label       string \
#     -xfont       xfont \
#     -tkfont      tkfontname \
#     -brightcolor rgb \
#     -darkcolor   rgb 

#
# PURPOSE

# This proc edits a custom atomic label on per-atom (i.e. atomID)
# basis. It is possible to edit the label-text, font (i.e. font-type
# and size) and font-colors. One can specify either XLFD font name,
# i.e., those that look as: 
# -*-itc bookman-demi-r-normal-*-13-*-*-*-*-*-*-*, or Tk font name. 
# The Tk font names are those created by the "font create ..." Tk 
# command.

#
# ARGUMENTS
# atomID -- the ID (i.e. sequential number) of an atom
# args   -- various configuration "-option value" pairs (see Options)
#
# OPTIONS
# ------------------------------------------------------------------------
#  OPTION        ALLOWED-VALUES + DESCRIPTION
# ------------------------------------------------------------------------
#  -label        string     --> the label-text
#
#  -xfont        xfontname  --> X11's XLFD font name
#
#  -tkfont       tkfontname --> Tk font name, previously create by the
#                               "font create ..." Tk command. If both -xfont 
#                               and -tkfont options are specified the value 
#                               of -tkfont option is taken/
#
#  -brightcolor  rgb        --> the red-green-blue specification of bright 
#                               font-color color (see scripting::atomicLabels 
#                               for the description) of bright and dark font 
#                               colors). Color can be specified either in 
#                               hexadecimal forms (i.e. #rrggbb) or as a list
#                               of three floats in range [0,1] (i.e. 
#                               {0.5 0.5 0.5}).
#      
#  -darkcolor    rgb        --> the red-green-blue specification of dark
#                               font-color color (see scripting::atomicLabels 
#                               for the description) of bright and dark font 
#                               colors). Color can be specified either in 
#                               hexadecimal forms (i.e. #rrggbb) or as a list
#                               of three floats in range [0,1] (i.e. 
#                               {0.5 0.5 0.5}).
# ------------------------------------------------------------------------
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::atomicLabels::atomID 10 \
#     -label       "Fe(spin up)"
#     -xfont       -*-itc bookman-demi-r-normal-*-13-*-*-*-*-*-*-* \
#     -brightcolor #ffffff \
#     -darkcolor   #000000
#
#****
# ------------------------------------------------------------------------

proc scripting::atomicLabels::atomID {atomID args} {
    variable atomicLabels

    _parseOptions    atomID $args
    _checkAtomID     atomID $atomID    
    _mapCustomValues $atomID 

    xc_setatomlabel .mesa $atomID \
	$atomicLabels(atomID.$atomID.label) \
	$atomicLabels(atomID.$atomID.font) \
	$atomicLabels(atomID.$atomID.brightcolor) \
	$atomicLabels(atomID.$atomID.darkcolor)
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::atomicLabels::atomType
#
# NAME
# scripting::atomicLabels::atomType
#
# USAGE
# scripting::atomicLabels::atomType atomType \
#     -label       string \
#     -xfont       xfont  \
#     -tkfont      tkfontname \
#     -brightcolor rgb \
#     -darkcolor   rgb 

#
# PURPOSE

# This proc edits a custom atomic labels on per-atom-type
# (i.e. atomType) basis. It is possible to edit the label-texts, font
# (i.e. font-type and size) and font-colors. One can specify either
# XLFD font name, i.e., those that look as: 
# -*-itc bookman-demi-r-normal-*-13-*-*-*-*-*-*-*, or Tk font name.
# The Tk font names are those created by the "font create ..." Tk
# command.

#
# ARGUMENTS
# atomType -- atomic-number or atomic-symbol
# args     -- various configuration "-option value" pairs (see Options)
#
# OPTIONS
# ------------------------------------------------------------------------
#  OPTION        ALLOWED-VALUES + DESCRIPTION
# ------------------------------------------------------------------------
#  -label        string     --> the label-text
#
#  -xfont        xfontname  --> X11's XLFD font name
#
#  -tkfont       tkfontname --> Tk font name, previously create by the
#                               "font create ..." Tk command. If both -xfont 
#                               and -tkfont options are specified the value 
#                               of -tkfont option is taken/
#
#  -brightcolor  rgb        --> the red-green-blue specification of bright 
#                               font-color color (see scripting::atomicLabels 
#                               for the description) of bright and dark font 
#                               colors). Color can be specified either in 
#                               hexadecimal forms (i.e. #rrggbb) or as a list
#                               of three floats in range [0,1] (i.e. 
#                               {0.5 0.5 0.5}).
#      
#  -darkcolor    rgb        --> the red-green-blue specification of dark
#                               font-color color (see scripting::atomicLabels 
#                               for the description) of bright and dark font 
#                               colors). Color can be specified either in 
#                               hexadecimal forms (i.e. #rrggbb) or as a list
#                               of three floats in range [0,1] (i.e. 
#                               {0.5 0.5 0.5}).
# ------------------------------------------------------------------------
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::atomicLabels::atomType Zn \
#     -label       "Fe(spin up)"
#     -xfont       -*-itc bookman-demi-r-normal-*-13-*-*-*-*-*-*-* \
#     -brightcolor #ffffff \
#     -darkcolor   #000000
#
#****
# ------------------------------------------------------------------------

proc scripting::atomicLabels::atomType {atomType args} {
    variable atomicLabels
    ::global   mody

    _parseOptions atomType $args

    set natoms [xc_getvalue $mody(GET_NATOMS)]
    set nat    [_getAtomicNumber atomType $atomType]

    # search for all atomType atoms and update

    for {set id 1} {$id <= $natoms} {incr id} {
	if { [xc_getvalue $mody(GET_NAT) $id] == $nat } {

	    _mapCustomValues $id

	    xc_setatomlabel .mesa $id \
		$atomicLabels(atomID.$id.label) \
		$atomicLabels(atomID.$id.font) \
		$atomicLabels(atomID.$id.brightcolor) \
		$atomicLabels(atomID.$id.darkcolor)
	}
    }
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::atomicLabels::render
#
# NAME
# scripting::atomicLabels::render
#
# USAGE
# scripting::atomicLabels::render global on|off
# or
# scripting::atomicLabels::render custom on|off
# or
# scripting::atomicLabels::render atomID id on|off
# or
# scripting::atomicLabels::render atomType type on|off
#
# SYNOPSIS
# scripting::atomicLabels::render labelType args
# 
# PURPOSE
# This proc toggles the display of global or custom atomic labels.
#
# scripting::atomicLabels::render global -- toggles the display of all
# global atomic labels
#
# scripting::atomicLabels::render custom -- toggles the display of all
# custom atomic labels. 
#
# scripting::atomicLabels::render atomID id -- toggles the display of
# a particular atomic-label (i.e. that of id-th atom). The argument
# "id" is the an atomic id, i.e., sequential number.
#
# scripting::atomicLabels::render atomType type -- toggles the display
# of all labels of a give atomic-type. The argument "type" is either
# atomic symbol or atomic number.
#
# ARGUMENTS
# labelType -- must be one of globa, custom, atomID or atomType
# args      -- the rest of the arguments (depends on the labelType) 
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::atomicLabels::render atomType Mg off
#
#****
# ------------------------------------------------------------------------

proc scripting::atomicLabels::render {labelType args} {
    variable atomicLabels 
    ::global   mody

    switch -exact -- $labelType {
	"global" {
	    # usage: render global on|off
	    set toggle [_toggleValue $args]
	    xc_newvalue .mesa $mody(SET_GLOBALATOMLABEL_DO_DISPLAY) $toggle
	}

	"custom" {
	    # usage: render custom on|off
	    set antiToggle [_antiToggleValue $args]
	    foreach id [xc_getvalue $mody(GET_ATOMLABEL_ALL_ID)] {
		xc_newvalue .mesa $mody(SET_DO_NOT_DISPLAY_ATOMLABEL) $id $antiToggle
	    }
	}
	
	"atomID" {
	    # usage: render atomID id on|off
	    set id         [lindex $args 0]
	    set antiToggle [_antiToggleValue [lindex $args 1]]
	    set natoms     [xc_getvalue $mody(GET_NATOMS)]
	    _checkAtomID render $id
	    xc_newvalue .mesa $mody(SET_DO_NOT_DISPLAY_ATOMLABEL) $id $antiToggle
	}

	"atomType" {
	    # usage: render atomType type on|off	    
	    set nat        [_getAtomicNumber render [lindex $args 0]]
	    set antiToggle [_antiToggleValue [lindex $args 1]]
	    set natoms     [xc_getvalue $mody(GET_NATOMS)]
	    for {set id 1} {$id <= $natoms} {incr id} {
		if { [xc_getvalue $mody(GET_NAT) $id] == $nat } {
		    xc_newvalue .mesa $mody(SET_DO_NOT_DISPLAY_ATOMLABEL) $id $antiToggle
		}
	    }
	}
	
	default {
	    ErrorIn scripting::atomicLabels::render \
		"wrong type of atomic label $labelType, must be global, custom, atomID, or atomType"
	    return
	}
    }
}
proc scripting::atomicLabels::_toggleValue {toggle} {
    if { $toggle == "on" } {
	return 1
    } elseif { $toggle == "off" } {
	return 0
    } else {
	ErrorIn scripting::atomicLabels::render \
	    "wrong toggle value $toggle, must be \"on\" or \"off\""
	return -code return
    }
}
proc scripting::atomicLabels::_antiToggleValue {toggle} {
    if { $toggle == "on" } {
	return 0
    } elseif { $toggle == "off" } {
	return 1
    } else {
	ErrorIn scripting::atomicLabels::render \
	    "wrong toggle value $toggle, must be \"on\" or \"off\""
	return -code return
    }
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::atomicLabels::clear
#
# NAME
# scripting::atomicLabels::clear
#
# USAGE
# scripting::atomicLabels::clear all
# or
# scripting::atomicLabels::clear atomID id
# or
# scripting::atomicLabels::clear atomType type
#
# SYNOPSIS
# scripting::atomicLabels::clear mode args
# 
# PURPOSE
# This proc clears the custom label according to specified mode, which
# is either all, atomID, or atomType.
#
# scripting::atomicLabels::clear all -- clears all custom labels,
# meaning that all labels become global-labels
#
# scripting::atomicLabels::clear atomID id -- clears the custom label
# of id-th atom. Its label become the global label
#
# scripting::atomicLabels::clear atomType type -- clears all custom
# labels of a given atomic type. These labels become global labels.
#
# ARGUMENTS
# mode -- must be one of all, atomID, or atomType
# args -- either atom-ID or atom-Type for atomID or atomType, respectively
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::atomicLabels::clear atomType Mg
#
#****
# ------------------------------------------------------------------------

proc scripting::atomicLabels::clear {mode args} {
    ::global mody
    
    switch -exact -- $mode {
	"all" {
	    # usage: clear all
	    xc_clearatomlabel .mesa all
	}
	
	"atomID" {
	    # usage: render atomID id on|off
	    _checkAtomID render $args
	    xc_clearatomlabel .mesa $args
	}
	
	"atomType" {
	    # usage: render atomType type on|off	    
	    set nat        [_getAtomicNumber render $args]
	    set natoms     [xc_getvalue $mody(GET_NATOMS)]
	    for {set id 1} {$id <= $natoms} {incr id} {
		if { [xc_getvalue $mody(GET_NAT) $id] == $nat } {
		    xc_clearatomlabel .mesa $id
		}
	    }
	}
	
	default {
	    ErrorDialog "wrong mode \"$mode\" for scripting::atomicLabels::clear command, must be all, atomID, or atomType"
	    return
	}
    }
}



# ------------------------------------------------------------------------
#****if* Scripting/scripting::atomicLabels::_parseOptions
#
# NAME
# scripting::atomicLabels::_parseOptions
#
# USAGE
# scripting::atomicLabels::_parseOptions ?-options value? ...
#
# PURPOSE

# This proc is to parse the options of scripting::atomicLabels::****
# class of commands. For internal use only.

#
# ARGUMENTS
# args -- the "-option value" pairs
#
# RETURN VALUE
# Undefined.
#****
# ------------------------------------------------------------------------

proc scripting::atomicLabels::_parseOptions {cmd args} {
    variable atomicLabels

    #
    # load default values
    #
    _loadDefaults

    set labelType $cmd

    switch -exact -- $labelType {
	"global" {
	    set type global
	    set allowed_options { -xfont -brightcolor -darkcolor }
	}
	"atomID" - "atomType" {
	    set type custom
	    set allowed_options { -label -xfont -brightcolor -darkcolor }
	}
	default {
	    ErrorIn scripting::atomicLabels::_parseOptions \
		"unknown labelType $labelType, must be global, atomID, or atomType"
	    return -code return
	}
    }

    set xfont  ""
    set tkfont ""
    set i 0
    foreach option [lindex $args 0] {
	incr i
	
	if { $i%2 } {
            set tag $option
	} else {	    
            switch -glob -- $tag {
		"-label" { 
		    if { $type == "global" } {
			_unknownOption $cmd $tag $allowed_options
			return -code return
		    }
		    set label $option		
		    set atomicLabels($type.label) $option
		}

		"-xfont" {
		    set xfont $option
		}

		"-tkfont" {
		    set tkfont $option
		}

		"-brightcolor" -    
		"-darkcolor" {
		    set colortype [string trim $tag -]
		    if { [string range $option 0 0] == "#"} {
			# color in hexadecimal #rrggbb form

			set atomicLabels($type.$colortype) [rgb_h2f $option]
		    } else {
			# assuming color in float format, i.e. {1.0 0.5 0.3}

			set atomicLabels($type.$colortype) $option
		    }
		}

		"-render" {
		    switch -exact -- [string tolower $option] {
			yes - on - 1 { set atomicLabels($type.render) 1 }
			no - off - 0 { set atomicLabels($type.render) 0 }
			default {
			    ErrorIn scripting::atomicLabels::$cmd \
				"wrong value $option for -render option, must be on or off"
			    return -code return
			}
		    }
		}
		default {
		    _unknownOption $cmd $tag $allowed_options
		    return -code return
		}
	    }
	}
    }
    if { $i%2 } {
	ErrorIn scripting::atomicLabels::$cmd \
	    "odd number of arguments for \"-option value\" pairs"
	return -code return
    }

    # TkFontname has priority over XFontname

    if { $tkfont != "" } {
	set atomicLabels($type.font) [xcTkFontName2XLFD $tkfont]
    } elseif { $xfont != "" } {
	set atomicLabels($type.font) $xfont
    }
    
    return 1
}


proc scripting::atomicLabels::_mapCustomValues {atomID} {
    variable atomicLabels
    
    foreach elem {label font brightcolor darkcolor} {	
	if { ! [info exists atomicLabels(atomID.$atomID.$elem)] \
		 || $atomicLabels(custom.$elem) != "" } {	   
	    
	    set atomicLabels(atomID.$atomID.$elem) $atomicLabels(custom.$elem)	    
	}
    }
}

proc scripting::atomicLabels::_unknownOption {where tag args} {
    ErrorIn scripting::atomicLabels::$where "unknown option $tag, must be obe of $args"
}

proc scripting::atomicLabels::_checkAtomID {where id} {
    ::global mody
    set natoms [xc_getvalue $mody(GET_NATOMS)]
    
    if { ! [string is integer $id] } {
	ErrorIn scripting::atomicLabels::$where "$id is not atomic ID"
	return -code return
    }
    
    if { $id < 1 || $id > $natoms} {
	ErrorIn scripting::atomicLabels::$where "atomic ID $id out of range, must be within [1,$natoms]"
	return -code return
    }
    
    return 1
}

proc scripting::atomicLabels::_getAtomicNumber {where atomType} {
    set nat [Aname2Nat $atomType]
    if { ! [string is integer $nat] } {
	ErrorIn scripting::atomicLabels::$where "$atomType is not atomic symbol"
	return -code return
    }
    return $nat
}