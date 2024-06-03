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
# Copyright (c) 1996--2019 by Anton Kokalj                                  #
#############################################################################


#
# xcTkFontName2XLFD --
#
# Tries to map TkFontName to XLFD X11 font name, if it does not
# succeed, then returns an empty string.
#

proc xcTkFontName2XLFD {font} {
    global tcl_platform

    if { $tcl_platform(platform) == "windows" } {
    	set fontAttr [font actual $font]
    	set font     [font create]
    	eval {font configure $font} $fontAttr
    	return $font
    }

    puts stderr "*** xcTkFontName2XLFD : font = $font"
    puts stderr "*** xcTkFontName2XLFD : font actual font = [font actual $font]"

    # *** below is for X11 only::

    # --------------------------------------------------
    # construct the font in the following form:
    # --------------------------------------------------
    # -foundry-family-weight-slant-setwidth-addstyle-pixel-point-resx-resy-spacing-width-charset-encoding    
    # ------------------------------------------------------------------------

    # --------------------------------------------------
    # Tk allowed fields
    # --------------------------------------------------
    #          -family name
    #          -size size
    #          -weight weight
    #          -slant slant
    #          -underline boolean
    #          -overstrike boolean

    set family_ helvetica
    foreach opt {family size weight slant} {
	upvar 1 $opt var
	set var  [font actual $font -$opt]
	set $opt $var
	
	# weight: 
	#         normal = normal | regular | medium | book | light
	#         bold   = bold | extrabold | demi | demibold
	#
	# slant: 
	#         italic = i | o

        if { $opt == "family" } {
            set family_ [X11Font_MatchFamily $var]
        }
	if { $opt == "weight" } {
	    if { $var == "normal" } {
		set weightList { medium normal regular book light }
	    } else {
		set weightList { bold extrabold demi demibold }
	    }
	}
	if { $opt == "slant" } {
	    if { $var == "italic" } {
		set slantList { i o }
	    } else {
		set slantList { r }
	    }
	}
    }     

    set family $family_
    
    # a hack for Mac OS X, which doesn't like negative sizes

    global tcl_plaform
    if { $tcl_platform(os) == "Darwin" } {
	if { [string is integer $size] && $size < 0 } {
	    set size [expr $size * (-1)]
	}
    }

    foreach weight $weightList {
	foreach slant $slantList {
	    # example::   "-*-bookman-*      -*     -*-*-64   -*-*-*-*-*-*-*"
	    set XLFD_name "-*-$family-$weight-$slant-*-*-$size-*-*-*-*-*-*-*"
            set result [xc_queryfont .mesa $XLFD_name]

            #puts stderr "*** xc_queryfont .mesa $XLFD_name ; result = $result"
            
            if { $result > 0 } {
		return $XLFD_name
	    }
	}
    }

    # couldn't map tk-font-name --> XLFD name, return a tk-font name
    return $font
}


# dirty fix to map from different Font families to most basic ones supported by X11
proc X11Font_MatchFamily {family} {
    
    if { [regexp -nocase Helvetica|Arial|Verdana|Carlito $family] } {
        set family_ helvetica
    } elseif { [regexp -nocase Serif|Times|Termes|Math $family] } {
        set family_ times
    } elseif { [regexp -nocase Courier|Mono $family] } {
        set family_ courier
    } elseif { [regexp -nocase Fixed $family] } {
        set family_ fixed
    } else {
        set family_ helvetica
    }

    return $family_
}

#togl .mesa
#pack .mesa
#
#foreach font {
#    TkDefaultFont
#    TkTextFont	
#    TkFixedFont	
#    TkMenuFont	
#    TkHeadingFont	
#    TkCaptionFont	
#    TkSmallCaptionFont
#    TkIconFont	
#    TkTooltipFont
#} {
#    puts [xcTkFontName2XLFD TkDefaultFont]
#}
#
