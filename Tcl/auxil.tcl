#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/auxil.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc xcPlace {w1 w2 x y} {
    # this proc place window w2 near w1
    # x,y where to place according to w1

    # query geom of w1
    set geom [wm geometry $w1]
    set x [expr [lindex [split $geom x+-] 2] + $x]
    set y [expr [lindex [split $geom x+-] 3] + $y]
    # now place $w2 to +x +y
    wm geometry $w2 +$x+$y
    
    return
}


    
proc Nat2Aname {nat} {
    global Alist

    # maybe NAT is greater than 100, so
    set nat  [expr $nat % 100]
    set selA [lindex $Alist $nat]
    return $selA
}


proc Aname2Nat {atom} {
    global Alist

    set n 0
    set Atom [string toupper $atom]
    foreach elem $Alist {
	set Elem [string toupper $elem]
	if { $Atom == $Elem } { return $n }
	incr n
    }
    # if we come so far we have an illegal atom name
    return "unknown atom name \"$atom\""
}


proc AnameExt2Nat {atom} {
    # proc names stands for: AtomNameExtended to Nat
    #
    # extended atom name is AtomicSymbolCharacters
    global Alist

    set n    1
    set Nat -1
    foreach elem $Alist {
	if { [string match -nocase $elem* $atom] } { 
	    # first compare if it is two-character match 
	    # example: S vs. Se
	    if { [string equal -nocase -length 2 $elem $atom] } {
		return $n
	    }
	    set Nat $n 
	}
	incr n
    }
    if { $Nat > 0 } {
	return $Nat
    } else {
	# if we come so far we have an illegal atom name
	return 0
    }
}

    
proc AtomNames {} {
    global Alist
    
    set Alist { X \
	    H He Li Be  B  C  N  O  \
	    F Ne Na Mg Al Si  P  S  \
	    Cl Ar  K Ca Sc Ti  V Cr \
	    Mn Fe Co Ni Cu Zn Ga Ge \
	    As Se Br Kr Rb Sr  Y Zr \
	    Nb Mo Tc Ru Rh Pd Ag Cd \
	    In Sn Sb Te  I Xe Cs Ba \
	    La Ce Pr Nd Pm Sm Eu Gd \
	    Tb Dy Ho Er Tm Yb Lu Hf \
	    Ta  W Re Os Ir Pt Au Hg \
	    Tl Pb Bi Po At Rn Fr Ra \
	    Ac Th Pa  U Np Pu Am Cm \
	    Bk Cf Es Fm}
    return $Alist
}


proc exit_pr {{arg {}}} {
    global system

    set button 1
    if { $arg == {} } {
	set button [tk_messageBox -message "Really quit?" \
			-type yesno -icon question]
    } else {
	# exit_pr -silent
	set button yes
    }

    if { $button == "yes" } {
	SetWatchCursor
	if { ![file exists $system(SCRDIR)] } {
	    exit 0
	}

	catch "cd $system(PWD)"
	
	clean_exit
    }
}


proc clean_exit {{returnCode 0}} {
    global system
    if { [file isdirectory $system(SCRDIR)] } {
	xcDebug -stderr "************************************************************************"
	xcDebug -stderr "Deleting directory $system(SCRDIR) ; Please Wait !!!"
	xcDebug -stderr "************************************************************************"
	# we catch deleting, since on some strange NFS systems directory 
	# deletion will fail
	if { ! [catch {file delete -force -- $system(SCRDIR)}] } {
	    xcDebug -stderr "Directory deleted !"
	} else {
	    xcDebug -stderr "Failed to delete the directory !"
	}
    }

    global exit_viewer_win
    if { [info exists exit_viewer_win] } {
	foreach win $exit_viewer_win {
	    if { [winfo exists $win] } {
		bind $win <Destroy> {}
	    }
	}	
    }
    exit_tcl $returnCode
}


#############################################################################
# this is the old routine for manipulating the font; the new one named 
# ModifyFont uses internal Tk "font" command
proc ModifyFontSize {w size {arg {}}} {
    global oldsize
    
    # take care of $arg:: 
    if { $arg != {} } {
	set foundry *
	set family  *
	set weight  *
	set slant   *
	set i 0
	foreach option $arg {
	    if { [regexp {^-} $option] } {
		set tag $option
	    } else {
		switch -- $tag {
		    "-foundry" {set foundry $option}
		    "-family"  {set family  $option}
		    "-weight"  {set weight  $option}
		    "-slant"   {set slant   $option}
		    default    {
			tk_dialog [WidgetName .a] Error \
				"ERROR: Bad option \"$tag\" submited to \n\
				ModifyFontSize procedure"
			error 0 OK
			return fixed 
		    }
		}
	    }
	    incr i
	}
	
	if { $i%2 } {
	    tk_dialog [WidgetName .a] Error \
		    "ERROR: You called ModifyFontSize with an \n\
		    odd number of args !" \
		    error 0 OK
	    return fixed
	}
	set font "-$foundry-$family-$weight-$slant-*-*-$size-*"	
    } else {    
	set font [lindex [$w config -font] 3]
    }
    
    set fontlist [split $font -]    
    #xcDebug "Fontlist:: $fontlist"
    
    if { [llength  $fontlist] > 7 || $arg != {} } {
	# X font name
	# define first four fields; 
	set foundry [lindex [split $font -] 1]
	set family  [lindex [split $font -] 2]
	set weight  [lindex [split $font -] 3]
	set slant   [lindex [split $font -] 4]    
	
	set font    [lindex $font 3]
	
	#puts stdout "TTT:: $font"
	#puts stdout "LLL:: -$foundry-$family-$weight-$slant-*-*-$size-*"
	#first four fields + pixel fields are impotant for us
	set done 1
	if [catch {$w config -font \
		-$foundry-$family-$weight-$slant-*-*-$size-*}] {
	    set done 0
	    set upsize   $size	
	    set downsize $size
	} else {
	    return "-$foundry-$family-$weight-$slant-*-*-$size-*"
	}
	
	for {} {$done != 1} {} {
	    # first try up, than try down
	    set upsize   [expr $upsize + 1]
	    set downsize [expr $downsize - 1]
	    #puts "$upsize $downsize"
	    # maybe we have gone to far
	    if { $downsize == 0 || $upsize > [expr $size + 20]} {
		# use default font
		if [catch {$w config -font $font}] {
		    $w config -font TkFixedFont
		    return fixed
		} else {
		    return $font
		}
	    }
	    if {[catch {$w config -font \
		    -$foundry-$family-$weight-$slant-*-*-$upsize-*}] == 0} {
		return "-$foundry-$family-$weight-$slant-*-*-$upsize-*"
	    } elseif {[catch {$w config -font \
		    -$foundry-$family-$weight-$slant-*-*-$downsize-*}] == 0} {
		return "-$foundry-$family-$weight-$slant-*-*-$downsize-*"
	    }
	}
    } elseif { [llength  $fontlist] == 2 } {
	# maybe alias name is something like: {Helvetica -12 bold}
	set oldsize [lindex [lindex $fontlist 1] 0]
	#puts stdout "oldsize:: $oldsize"
	if [number oldsize posint] {
	    #replace oldsize with size
	    set newfont \
		    [concat [lindex $fontlist 0] -$size [lindex $fontlist 2]]
	    if [catch {$w config -font $newfont}] {		   
		set done 0
		set upsize   $size	
		set downsize $size
	    } else {
		return $newfont
	    }

	    for {} {$done != 1} {} {
		# first try up, than try down
		set upsize   [expr $upsize + 1]
		set downsize [expr $downsize - 1]
		#puts "$upsize $downsize"
		# maybe we have gone to far
		if { $downsize == 0 || $upsize > [expr $size + 20]} {
		    # use default font
		    if [catch {$w config -font $font}] {
			$w config -font TkFixedFont
			return fixed
		    } else {
			return $font
		    }
		}
		set upfont [concat [lindex $fontlist 0] -$upsize \
			[lindex $fontlist 2]]
		set downfont [concat [lindex $fontlist 0] -$downsize \
			[lindex $fontlist 2]]		
		if {[catch {$w config -font $upfont}] == 0 } {
		    return $upfont
		} elseif {[catch {$w config -font $downfont}] == 0 } {
		    return $downfont
		}
	    }	    
	}
    } else {
	# give up; use new Tk font Mechanism
	set font [ModifyFont [$w cget -font] $w -size $size -default 1]
    }
}


#
# this is the new routine and uses Tk font mechanism
#
proc ModifyFont {font window {args {}}} {
    global modifyFont

    # allowed arguments::
    #     -default     ... create default font (0/1) (IGNORED)
    #	  -family      ... name 
    #	  -size        ... size 
    #	  -weight      ... weight 
    #	  -slant       ... slant 
    #	  -underline   ... boolean 
    #	  -overstrike  ... boolean 

    # font actual font ?-displayof window? ?option?
    set default     0
    set family      [font actual $font -displayof $window -family]
    set size        [font actual $font -displayof $window -size]
    set weight      [font actual $font -displayof $window -weight]
    set slant       [font actual $font -displayof $window -slant]
    set underline   [font actual $font -displayof $window -underline]
    set overstrike  [font actual $font -displayof $window -overstrike]    

    # take care of $arg:: 
    if { $args != {} } {
	set i 0
	foreach option $args {
	    if { [regexp {^-} $option] } {
		set tag $option
	    } else {
		switch -- $tag {		    
		    "-default"    {set default    $option}
		    "-family"     {set family     $option}
		    "-size"       {set size       $option}
		    "-weight"     {set weight     $option}
		    "-slant"      {set slant      $option}
		    "-underline"  {set underline  $option}
		    "-overstrike" {set overstrike $option}
		    default    {
			tk_dialog [WidgetName] Error \
				"ERROR: Bad option \"$tag\" submited to \n\
				ModifyFont procedure" error 0 OK
			return fixed 
		    }
		}
	    }
	    incr i
	}
    
	if { $i%2 } {
	    tk_dialog [WidgetName] Error \
		    "ERROR: You called ModifyFont with an \n\
		    odd number of args !" \
		    error 0 OK
	    return fixed
	}
    }

    set new_font [font create]
    
    font configure $new_font \
	-family     $family \
	-size       $size \
	-weight     $weight \
	-slant      $slant \
	-underline  $underline \
	-overstrike $overstrike
    return $new_font
}

proc SetFont {widtype args} {
    set w    [$widtype [WidgetName]]
    set font [$w cget -font]
    xcDebug -debug "SetFont: $font $w $args"
    set Font [eval {ModifyFont $font $w} $args]
    destroy $w
    return $Font
}
    
#
# this routine uses Tk font mechanism
#
proc GetFontAtribute {font window arg} {

    # take care of $arg:: 
    if { [llength $arg] != 1 } {
	tk_dialog [WidgetName] Error \
		"ERROR: You called GetFontAtribute with wrong \
		    number of args !" \
		    error 0 OK
	return fixed
    }
    set tag    [lindex $arg 0]
    set option [lindex $arg 1]
    switch -- $tag {		    
	"-family"     {return [font actual $font -displayof $window -family]}
	"-size"       {return [font actual $font -displayof $window -size]}
	"-weight"     {return [font actual $font -displayof $window -weight]}
	"-slant"      {return [font actual $font -displayof $window -slant]}
	"-underline"  {
	    return [font actual $font -displayof $window -underline]}
	"-overstrike" {
	    return [font actual $font -displayof $window -overstrike]}
	default    {
	    tk_dialog [WidgetName .a] Error \
		    "ERROR: Bad option \"$tag\" submited to \n\
		    GetFontAtribute procedure"
	    error 0 OK
	    return 0
	}
    }
}


# #
# # xcTkFontName2XLFD --
# #
# # Tries to map TkFontName to XLFD X11 font name, if it does not
# # succeed, then returns an empty string.
# #
# 
# proc xcTkFontName2XLFD {font} {
#     global tcl_platform
# 
#     if { $tcl_platform(platform) == "windows" } {
#     	set fontAttr [font actual $font]
#     	set font     [font create]
#     	eval {font configure $font} $fontAttr
#     	return $font
#     }
# 
#     #puts stderr "*** xcTkFontName2XLFD : font = $font"
#     #puts stderr "*** xcTkFontName2XLFD : font actual font = [font actual $font]"
#     
#     # *** below is for X11 only::
# 
#     # --------------------------------------------------
#     # construct the font in the following form:
#     # --------------------------------------------------
#     # -foundry-family-weight-slant-setwidth-addstyle-pixel-point-resx-resy-spacing-width-charset-encoding    
#     # ------------------------------------------------------------------------
# 
#     # --------------------------------------------------
#     # Tk allowed fields
#     # --------------------------------------------------
#     #          -family name
#     #          -size size
#     #          -weight weight
#     #          -slant slant
#     #          -underline boolean
#     #          -overstrike boolean
# 
#     foreach opt {family size weight slant} {
# 	upvar 1 $opt var
# 	set var  [font actual $font -$opt]
# 	set $opt $var
# 	
# 	# weight: 
# 	#         normal = normal | regular | medium | book | light
# 	#         bold   = bold | extrabold | demi | demibold
# 	#
# 	# slant: 
# 	#         italic = i | o
# 	
# 	# if { $opt == "weight" } {
# 	#     if { $var == "normal" } {
# 	# 	set weightList { medium normal regular book light }
# 	#     } else {
# 	# 	set weightList { bold extrabold demi demibold }
# 	#     }
# 	# }
# 	# if { $opt == "slant" } {
# 	#     if { $var == "italic" } {
# 	# 	set slantList { i o }
# 	#     } else {
# 	# 	set slantList { r }
# 	#     }
# 	# }
#     }     
#     set XLFD_name "-*-$family-$weight-$slant-*-*-$size-*-*-*-*-*-*-*"
#     #puts stderr "*** xcTkFontName2XLFD : XLFD_name = $XLFD_name"
#     
#     return $XLFD_name
# 
#     # a hack for Mac OS X, which doesn't like negative sizes
#     
#     #global tcl_plaform
#     #if { $tcl_platform(os) == "Darwin" } {
#     #    if { [string is integer $size] && $size < 0 } {
#     #        set size [expr $size * (-1)]
#     #    }
#     #}
#     #
#     #foreach weight $weightList {
#     #    foreach slant $slantList {
#     #        # example::   "-*-bookman-*      -*     -*-*-64   -*-*-*-*-*-*-*"
#     #        set XLFD_name "-*-$family-$weight-$slant-*-*-$size-*-*-*-*-*-*-*"
#     #        if { [xc_queryfont .mesa $XLFD_name] > 0 } {
#     #    	return $XLFD_name
#     #        }
#     #    }
#     #}
# 
#     # couldn't map tk-font-name --> XLFD name, return tk-font-name
#     #return $font
# }



proc AlwaysOnTopON {lower upperlist} {
    xcDebug -debug "AlwaysOnTopON"
    #there maybe more than one widget to raise
    foreach upper $upperlist {
	xcRaiseRegister $upper $lower
    }
    bind $lower <Button-1> [list xcRaise $lower]
    bind $lower <Button-2> [list xcRaise $lower]
    bind $lower <Button-3> [list xcRaise $lower]    
}


proc xcRaiseRegister {upper lower} {
    global xcRaise
    #
    # parse xcRaise($lower,toplevels)
    #
    if ![info exists xcRaise($lower,toplevels)] {
	set xcRaise($lower,toplevels) {}
    }
    set new_list {}
    foreach win $xcRaise($lower,toplevels) {
	if { [winfo exists $win] && $win != $upper } {
	    append new_list "$win "
	}
    }
    set xcRaise($lower,toplevels) [concat $new_list $upper]
}


proc xcRaise lower {
    global xcRaise
    # xcRaise($lower,toplevels) tells if there are some more toplevels 
    # to raise !!!
    
    foreach widget $xcRaise($lower,toplevels) {
	if { [winfo exists $widget] } {
	    raise $widget $lower
	}
    }
    xcDebug -debug "$lower,xcRaise($lower,toplevels) $xcRaise($lower,toplevels)"
}


proc AlwaysOnTopOFF {{lower {.}}} {
    puts stdout "AlwaysOnTopOFF"
    bind $lower <Button-1> {}
    bind $lower <Button-2> {}
    bind $lower <Button-3> {}
}


proc CancelProc {w {var {}}} {
    upvar $var varn
    if { [winfo exists $w] } {
	AlwaysOnTopOFF
	catch { grab release $w }
	destroy $w
    }
    set varn 0
    #uplevel { return 0 }
    return 0
}

# a simple wrapper to be used with widget "-command"
proc DestroyWid w {
    destroy $w
}

proc winGeom { w } {
    # procedure determines the geometry of $w and return it
    return [wm geometry $w]
}


##############################################################################
# this proc read dimension & group (family) number out of GENGEOM file
proc GetDimGroupXSF {dim group xsfFile} {
    upvar $dim   dm
    upvar $group gr
    set fileID [open $xsfFile r]
    GetDimGroup dm gr $fileID
    close $fileID
}
proc GetDimGroup {dim group fileID} {
    upvar $dim   dm
    upvar $group gr

    # <t.k.>: Thu Jul 13 15:13:57 CEST 2017
    set dm 0
    set gr 1
    # </t.k.>
    set n 0
    set output [split [read $fileID] \n]
    foreach line $output {
	switch -regexp -- $line {
	    {^ *DIM-GROUP} {
		set nn [expr $n + 1]
		set dimgroup [lindex $output $nn]
		set dm       [lindex $dimgroup 0]
		set gr       [lindex $dimgroup 1]
		xcDebug "GET-DIM-GROUP:: [lindex $dimgroup 0] [lindex $dimgroup 1]"
		return
	    }
	    {^ *POLYMER} {
		set dm 1
		set gr 1
		return
	    }
	    {^ *SLAB} {
		set dm 2
		set gr 1
		return
	    }
	    {^ *CRYSTAL} {
		set dm 3
		set gr 1
		return
	    }
	}
	incr n
    }
}


# this proc is synonym for CellMode
proc GenGeomDisplay {{update 0}} {
    xcDebug "In GenGeomDisplay"
    CellMode $update
}


##############################################################################
# conversion between angstrom & bohr
proc Bohr2Angs var {
    global Const
    return [expr $var * $Const(bohr)]
}


proc Angs2Bohr var {
    global Const
    return [expr $var / $Const(bohr)]
}


###########################################################
#this proc generate a widget name that do not already exist
proc WidgetName {{w {}}} {
    set i 0
    for {} {1} {} {
	if [winfo exist $w.a$i] {
	    incr i
	} else {
	    return $w.a$i
	}
    }
}


##############################################################################
proc GetWidgetConfig {widget_com option} {
    
    for {set i 1} {1} {incr i} {
	if ![winfo exists .gwc$i] {
	    set w .gwc$i
	    break
	}
    }

    $widget_com $w
    set res [$w cget $option]
    if {    $option == "-background" || \
	    $option == "-bg" || \
	    $option == "-foreground" || \
	    $option == "-fg" || \
	    $option == "-activebackground" || \
	    $option == "-activeforeground" || \
	    $option == "-highlightbackground" || \
	    $option == "-hightlightcolor" || \
	    $option == "-disabledforeground" || \
	    $option == "-insertbackground" || \
	    $option == "-selectbackground" || \
	    $option == "-selectcolor" || \
	    $option == "-selectforeground" || \
	    $option == "-troughcolor" } {
	if { [string range $res 0 0] != "#" } { 
	    set norm [lindex [winfo rgb . white] 0]
	    set rgb  [winfo rgb . $res]
	    set res  [format "#%02x%02x%02x" \
		    [expr 256 * [lindex $rgb 0] / $norm] \
		    [expr 256 * [lindex $rgb 1] / $norm] \
		    [expr 256 * [lindex $rgb 2] / $norm]]
	}
    }
    destroy $w
    return $res
}


##############################################################################
# return the filehead out of filename (file.poss --> file)
proc FileHead file {
    return [file rootname $file]
}


##############################################################################
# return the possix out of filename (file.poss --> poss)
proc FilePossix file {
    set filename [split $file .]
    set nfield   [llength $filename]
    set possix [lrange $filename [expr $nfield - 1] [expr $nfield - 1]]

    return $possix
}


#####################################
# used by, for example, DefaultButton
proc DummyProc {{args {}}} {
    return
}


###############################################################
# if numbers specified in $args differ .leq. $limit -> return 1
#                                                 else return 0
proc IsEqual {limit args} {

    set oldnum [lindex $args 1]
    foreach num $args {
	if { [expr abs($oldnum - $num)] >= $limit } {
	    return 0
	}
	set oldnum $num
    }
    
    return 1
}

###############################################################################
# is a lower_or_equal to b (within $limit)
proc IsLEQ {limit a b} {

    if { $a <= [expr $b + $limit] } {
	return 1
    } else {
	return 0

    }    
}


proc xcPause sec {
    set iter [lindex [time { for {set i 1} {$i <= 10} {incr i} {update} }] 0]
    set count [expr int(1e7 * $sec / $iter)]
    xcDebug "xcPause:: $count"
    for {set i 1} {$i < $count} {incr i} {update}
}


proc rgb_h2d rgb {

    set len [string length $rgb]
    # len can be 4,7,10,13 

    set i [expr $len / 3]
    set norm 1
    for {set n 1} {$n <= $i} {incr n} {
	set norm [expr $norm * 16]
    }

    set r [string range $rgb 1 $i]
    set g [string range $rgb [expr 1 + $i] [expr 2 * $i]] 
    set b [string range $rgb [expr 1 + 2 * $i] end]
    
    set r [h2df $r]
    set g [h2df $g]
    set b [h2df $b]
    
    return [list $r $g $b]
}



proc rgb_h2f rgb {
    set len [string length $rgb]
    # len can be 4,7,10,13 

    set i [expr $len / 3]
    set norm 1
    for {set n 1} {$n <= $i} {incr n} {
	set norm [expr $norm * 16]
    }

    set r [string range $rgb 1 $i]
    set g [string range $rgb [expr 1 + $i] [expr 2 * $i]] 
    set b [string range $rgb [expr 1 + 2 * $i] end]

    set r [h2df $r $norm]
    set g [h2df $g $norm]
    set b [h2df $b $norm]
    
    return [list $r $g $b]
}


proc rgb_f2h {rgba} {
    set r [d2h [expr round([lindex $rgba 0] * 255)]]
    set g [d2h [expr round([lindex $rgba 1] * 255)]]
    set b [d2h [expr round([lindex $rgba 2] * 255)]]

    return #${r}${g}${b}
}


proc rgb_f2d {rgba} {
    # f is clamped float in range [0--1]
    # returns decimal-list {255 255 255}
    set r [expr round([lindex $rgba 0] * 255)]
    set g [expr round([lindex $rgba 1] * 255)]
    set b [expr round([lindex $rgba 2] * 255)]
    return [list $r $g $b]
}


proc rgb_ac_f2h rgba {
    # same as rgb_f2h, just to get color a little briter
    
    set r [expr round([lindex $rgba 0] * 280)]
    set g [expr round([lindex $rgba 1] * 280)]
    set b [expr round([lindex $rgba 2] * 280)]
    if { $r > 255 } {set r 255}
    if { $g > 255 } {set g 255}
    if { $b > 255 } {set b 255}
    set r [d2h $r]
    set g [d2h $g]
    set b [d2h $b]

    return #${r}${g}${b}
}


proc rgb_d2f {rgba} {
    # f is clamped float in range [0--1]
    # BEWARE: assuming rgba (INPUT) as {255 255 255}
    set r [expr double([lindex $rgba 0]) / 255.0]
    set g [expr double([lindex $rgba 1]) / 255.0]
    set b [expr double([lindex $rgba 2]) / 255.0]

    return [list $r $g $b]
}


proc h2df {h {norm 1}} {
    # usage: h2df #rrggbb     --> returns decimal-list, i.e., {255 255 255}
    # usage: h2df #rrggbb 255 --> returns float-list, i.e. {1.0 1.0 1.0}
    set d   0
    set len [expr [string length $h] - 1]
    for {set i $len} {$i >= 0} {incr i -1} {
	set j [expr $len - $i]
	switch -regexp -- [set a [string range $h $j $j]] {
	    [fF] {set a 15}
	    [eE] {set a 14}
	    [dD] {set a 13}
	    [cC] {set a 12}
	    [bB] {set a 11}
	    [aA] {set a 10}
	}
	set d [expr $d + $a * [xcOnPower 16 $i]]
    }
    if { $norm > 1.0 } {
	return [expr double($d) / double($norm-1)]
    } else {
	return [expr int($d)]
    }
}


proc d2h {num} {
    set n1 [expr int( $num / 16 )]
    set n2 [expr int($num) - $n1 * 16]

    switch -exact -- $n1 {
	15 {set n1 f}
	14 {set n1 e}
	13 {set n1 d}
	12 {set n1 c}
	11 {set n1 b}
	10 {set n1 a}
    }

    switch -exact -- $n2 {
	15 {set n2 f}
	14 {set n2 e}
	13 {set n2 d}
	12 {set n2 c}
	11 {set n2 b}
	10 {set n2 a}
    }
    return [format "%s%s" $n1 $n2]
}

proc d2f {d} {
    # BEWARE: assuming d (INPUT) in range [0,255]
    # f is clamped float in range [0--1]
    return [expr double([lindex $d 0]) / 255.0]
}

proc xcOnPower {a n} {
    set res 1
    for {set i 1} {$i <= $n} {incr i} {
	set res [expr $res * $a]
    }
    return $res
}


#####################
# set cursor to watch
proc SetWatchCursor {} {
    global xcCursor
    foreach t [winfo children .] {
	catch { puts stderr "SetWatchCursor: $t" }
	if {"[info commands $t]" != {} } {
	    $t config -cursor $xcCursor(watch)
	}
    }
    . config -cursor $xcCursor(watch)
    CursorUpdate
}

proc CursorUpdate {} {
    global xcCursor
    if { [info exists xcCursor(dont_update)] } {
	if { ! $xcCursor(dont_update) } {
	    update
	}
    } else {
	update
    }
}

#######################
# set cursor to default
proc ResetCursor {} {
    global xcCursor
    if { [info exists xcCursor(dont_update)] } {
	if { $xcCursor(dont_update) } {
	    return
	}
    }
    foreach t [winfo children .] {
	if {"[info commands $t]" != {} } {
	    $t config -cursor $xcCursor(default)
	}
    }
    CursorUpdate
    . config -cursor $xcCursor(default)
    #CursorUpdate
}


proc xcSwapBuffers {} {
    if { [winfo exists .mesa] } {
	update
	xc_swapbuffer .mesa
    }
}


##############################################################################
# 
# Purpose: find out what is the name of fortran units (without number)    
# Return:  the name of fortran UNIT
proc FtnName {} {
    global system

    # create an empty $system(SCRDIR)/fort_unit/ directory
    set pwd [pwd]
    cd $system(SCRDIR)
    if { [file exists fort_unit] } {
	file delete -force fort_unit
    }
    file mkdir fort_unit

    # cd to dirt_unit and run a simple fortran test

    cd fort_unit
    xcCatchExecReturn $system(FORDIR)/ftnunit
    update
    set file [glob -nocomplain *]
    regsub {\.99} $file {} ftn_name

    # delete the fort_unit directory
    cd ..
    file delete -force fort_unit

    cd $pwd    
    return $ftn_name

    #
    # this was the old routine
    #
    #set cwd [pwd]
    #cd $system(SCRDIR)
    #exec $system(FORDIR)/ftnunit
    #update
    #set file [file tail [lindex [glob -nocomplain $system(SCRDIR)/*99] 0]]
    #regsub 99 $file {} file
    #exec rm -f ${file}99
    #cd $cwd
    #return $file
}


#
# capitalizes the word
#
proc capitalize word {
    set w1 [string toupper [string range $word 0 0]]
    set w2 [string range $word 1 end]
    return [format %s%s $w1 $w2]
}


#
# return the filehead (filename without extension)
#
proc filehead {filename} {
    return [file rootname $filename]
}


proc WriteFile {filename content {flag w}} {
    global tcl_platform
    set fID [open $filename $flag]
    if { $tcl_platform(platform) == "windows" } {
	fconfigure $fID -translation {auto lf}
    }
    puts $fID $content
    flush $fID
    close $fID
}

proc ReadFile {filename {arg {}}} {
    # Usage: ReadFile filename   OR   ReadFile -nonewline filename
    if { $arg != {} } {
	set filename $arg
    }
    set fID [open $filename r] 
    if { $arg != {} } {
	set output [read -nonewline $fID]
    } else {
	set output [read $fID]
    }
    close $fID
    return $output
}

proc GetAbsoluteFileName file {
    global system

    # try this:
    return [file normalize [file join $system(PWD) $file]]
    
    # if filename starts with / or ~ the absolute file name is assumed,
    # otherwise absolute filename should be: $system(PWD)/$file

    #if { $file == "." } {
    #    set file $system(PWD)
    #}
    #set file [string trimright $file /]
    #set c0 [string index $file 0]
    #if { $c0 == "/" || $c0 == "~" } {
    #    return $file
    #} else {
    #    return [file normalize [file join $system(PWD) $file]]
    #}
}

#-----------------------------------------
# convert angstrom unit to fractional unit
proc GetFracCoor {coor} {
#-----------------------------------------
    global system

    set x [lindex $coor 0]
    set y [lindex $coor 1]
    set z [lindex $coor 2]

    xcDebug -debug "exec $system(BINDIR)/fracCoor \
	    $system(SCRDIR)/xc_struc.$system(PID) $x $y $z"

    if { [catch {set coor [exec $system(BINDIR)/fracCoor $system(SCRDIR)/xc_struc.$system(PID) $x $y $z]} errmsg] } {
	ErrorDialog "error occured while executing \"fracCoor\" program.\n\nError Message:\n$errmsg"
	xcDebug -debug "GetFracCoor: $coor"
	return {0.0 0.0 0.0}
    }

    xcDebug -debug "GetFracCoor: $coor"	
    return $coor
}

# -----------------------------------------------
# convert coordinates from Angstrom to $unit unit
proc coorToUnit {unit x y z} {
    # unit must be one of: angs bohr prim conv alat
    global Const 
    
    switch -- $unit {
	bohr {
	    return [list [expr $x / $Const(bohr)] [expr $y / $Const(bohr)] [expr $z / $Const(bohr)]]
	}
	prim - conv {
	    return [xc_fractcoor -ctype $unit -coor [list $x $y $z]]
	}
	alat {
	    global mody
	    set alat [xc_getvalue $mody(GET_ALAT)]
	    return [list [expr $x / $alat] [expr $y / $alat] [expr $z / $alat]]
	}
	angs - default {
	    return [list $x $y $z]
	}
    }
}


##############################################################################
# DEBUGING
proc xcDebug {line {args {}}} {
    global xcMisc

    set channel stdout
    if { $line == "-stderr" } {
	set channel stderr
	set line [string trim $args \{\}]
    } elseif { $line == "-debug" && $xcMisc(debug) == 1 } {
	set channel stderr
	set line [string trim $args \{\}]
    }
    if ![catch {puts $channel $line}] {
	flush $channel
    }
}


proc xcEditFile {file {foreground 0}} {
    global env system
    
    if { [info exists env(EDITOR)] && [info exists system(term)] } {
	if { $foreground != 0 } {
	    exec $system(term) -e $env(EDITOR) $file
	} else {
	    exec $system(term) -e $env(EDITOR) $file &
	}
    } else {	
	if { $foreground != 0 } {
	    tkwait window [defaultEditor $file]
	} else {
	    defaultEditor $file
	}
    }
}



proc xcDeleteAllChildren {wlist} {
    
    foreach w $wlist {
	if ![winfo exists $w] continue
	set children [winfo children $w]
	if { $children != "" } {
	    foreach child $children {
		xcDeleteAllChildren $child
		catch [destroy $child]
	    }
	}
    }	
}


proc gunzipFile {file} {
    global system

    xcDebug -debug "gunzipFile: $file"

    ####################
    set gunzipName $file
    ####################

    set name [file tail $file]
    
    if { [file extension $name] == ".gz" } {
        
        set here [pwd]
        cd $system(SCRDIR)

	# maybe file is already located in $system(SCRDIR); if not copy it there
	if { [file dirname $file] != $system(SCRDIR) && $file != $name } {
	    file copy -force $file $name
	}

	catch {exec -- gzip -df $name}
	set gunzipName [file rootname $name]
	if { ![file exists $gunzipName] } {
	    tk_dialog [WidgetName] "ERROR" \
		"ERROR: error when gunzip-ing file $file" warning 0 OK
	    uplevel 1 { return }
	}
	set gunzipName $system(SCRDIR)/$gunzipName

        cd $here
    }
    
    return $gunzipName
}

# Purpose: clean a welcome window
proc destroyWelcome {} {
    if { [winfo exists .title] } {
	# destroy WELCOME window
	destroy .title
    }
}


proc ErrorDialogInfo {text {errMsg {}}} {
    destroyWelcome
    #error $text
    set id [tk_dialog [WidgetName] ERROR "ERROR: $text." error 0 OK ErrorInfo]
    if { $id == 1 } {
	tkwait window [xcDisplayVarText $errMsg "Error Info"]
    }
}


# Purpose: do exec and report an error upon failure
# Return:  0 on success, 1 on failure
proc xcCatchExec {args} {
    destroyWelcome
    xcDebug -stderr "Executing: $args"
    if { [catch {eval exec $args} errMsg] } {
	ErrorDialogInfo "while executing\nexec $args" $errMsg
	return 1
    } 
    return 0
}

# same as xcCatchExec but with redirection of stdout/stderr !!!
proc xcCatchExecRedirectStdErr {args} {
    destroyWelcome
    xcDebug -stderr "Executing: $args"
    if { [catch {eval exec $args 2> /dev/null} errMsg] } {
	ErrorDialogInfo "while executing\nexec $args" $errMsg
	return 1
    } 
    return 0
}

proc xcCatchExecReturn {args} {
    destroyWelcome
    xcDebug -stderr "Executing: $args"

    SetWatchCursor
    if { [catch {eval exec $args} errMsg] } {
	ErrorDialogInfo "while executing\nexec $args" $errMsg
        ResetCursor
	uplevel 1 {
	    return 1
	}
    }
    ResetCursor

    return 0
}


# same as xcCatchExecReturn but with redirection of stdout/stderr !!!
proc xcCatchExecReturnRedirectStdErr {args} {
    destroyWelcome
    xcDebug -stderr "Executing: $args"
    if { [catch {eval exec $args 2> /dev/null} errMsg] } {
	ErrorDialogInfo "while executing\nexec $args" $errMsg
	uplevel 1 {
	    return 1
	}
    }
    return 0
}

proc ErrorDialog {text {errMsg {}}} {
    destroyWelcome
    set text "ERROR: $text."
    if { $errMsg != "" } {
	append text "\n\nError Mesage:\n$errMsg"
    }
    tk_messageBox -title ERROR -message $text -type ok -icon error
}

proc WarningDialog {text {warnMsg {}}} {
    if { [winfo exists .title] } {
	# destroy WELCOME window
	destroy .title
    }
    set text "WARNING: $text"
    if { $warnMsg != "" } {
	append text "\n\nWarning Mesage:\n$warnMsg"
    }
    tk_messageBox -title WARNING -message $text -type ok -icon warning
}

proc ErrorIn {where text} {
    tk_messageBox -title ERROR -message "ERROR: $text\n\nThis error was triggered from $where procedure" -type ok -icon error
}


#
# xcSkipEmptyLines --
#
# Purpose: skip empty lines from the variable content
proc xcSkipEmptyLines {text} {
    foreach line [split $text \n] {
	if { [regexp -- {\w} $line] } {
	    append out [format "%s\n" $line]
	}
    }
    return $out
}


# ------------------------------------------------------------------------
# evaluate the Tcl commands within the catch command and if error occurs
# prints the errorMsg. If errorMsg is void, than prints the error message
# returned by the Catch command.
# ------------------------------------------------------------------------
proc xcCatchEval {cmd {errorMsg {}}} {

    if { [catch {eval $cmd} _errorMsg] } {
	if { $errorMsg == "" } {
	    set errorMsg $_errorMsg
	} else {
            append errorMsg "\n$_errorMsg"
        }
	ErrorDialog "An ERROR occured while executing:\n$cmd\n\nERROR MESSAGE: $errorMsg"
    }
}


# ------------------------------------------------------------------------
#****f* Scripting/repeat
#
# NAME
# repeat
#
# USAGE
# repeat ntimes script
#
# PURPOSE

# This proc is for repetitive execution of a script supllied by
# "script" argument. For example:
# 
# repeat 10 { puts "Hello !!!" }
#
# will print "Hello !!!" 10-times. The repeat is nothing else then
# simplified "for" loop. Above example could be also achieved by:
#
# for {set i 0} {$i < 10} {incr i} {
#    puts "Hello !!!"
# }
#

#
# SIDE EFFECTS
# Inside repeat scripts, the "repeat" variable have the value of the current
# repeat-iteration. For example:
#
# repeat 4 { puts "This is the $repeat. iteration !!!" }
#
# will print:
#
# This is the 1. iteration
# This is the 2. iteration
# This is the 3. iteration
# This is the 4. iteration

#
# ARGUMENTS
# ntimes -- how many times to execute a script
# script -- script to execute
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# repeat 10 { 
#    scripting::rotate x 5
#    scripting::makeMovie::makeFrame
# }
#
# SOURCE

proc repeat {ntimes script} {
    global repeat_script repeat

    set repeat_script $script
    for {set repeat 1} {$repeat <= $ntimes} {incr repeat} {
	uplevel 1 {eval $repeat_script}
    }
}
#****
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
#****f* Scripting/wait
#
# NAME
# wait
#
# USAGE
# wait ms
#
# PURPOSE
# This proc is similar to the Tcl command "after ms". However before
# waiting for period of ms milliseconds, it updates all the events
# (the after command does not make the update before waiting !!!)
#
# ARGUMENTS
# ms -- waiting time in milliseconds
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# wait 500
#
# SOURCE

proc wait {ms} {
    if { ! [string is integer $ms] } {
	ErrorIn wait "expected integer, but got $ms"
	return
    }
    update
    after $ms
}
#****
# ------------------------------------------------------------------------

proc positiveInteger {string} {
    if { ![string is integer $string] } {
	return 0
    } elseif { $string <= 0 } {
	return 0
    } else {
	return 1
    }
}

proc nonnegativeInteger {string} {
    if { ![string is integer $string] } {
	return 0
    } elseif { $string < 0 } {
	return 0
    } else {
	return 1
    }
}


# return 1 if string is the OpenGL RGBA color spec, 0 otherwise
proc rgba {string} {
    if { [llength $string] != 4 } {
	return 0
    }
    for {set i 0} {$i < 4} {incr i} {
	set v [lindex $string $i]
	if { ! [string is double $v] } {
	    return 0
	} elseif { $v > 1.0 || $v < 0.0 } {
	    return 0
	}
    }
    return 1
}


proc allowedValue {value allowedValues} {
    # returns 1 if $value is among item in $allowedValues list

    foreach item $allowedValues {
	if { $value == $item } {
	    return 1
	}
    }
    return 0
}

proc destroyWelcomeWindow {} {
    if { [winfo exists .title] } {
	destroy .title
    }
}


proc xcTempFile {name} {
    global system
    return $system(SCRDIR)/$name.$system(PID)
}


#
# evalInScratch -- evaluate the script in SCRATCH, i.e. $system(SCRDIR), directory
#
proc evalInScratch {script} {
    global system

    set here [pwd]
    cd $system(SCRDIR)
    uplevel 1 [list eval $script]
    cd $here
}

#
# evalInDir -- evaluate the script in $dir directory
#
proc evalInDir {dir script} {
    set here [pwd]
    cd $dir
    uplevel 1 eval $script
    cd $here
}


# evalInPWD -- This is a workaround routine: the code does many times
# "cd $system(SCRDIR)", hence the real pwd is lost. There is a global
# $system(PWD), but for example user might change in scripting-scripts
# the cd then his [pwd] is lost as system(PWD) was not updated. This
# routine execute the code either in [pwd], but if [pwd] ==
# $system(SCRDIR), then it executes the code in $system(PWD)
#
proc evalInPWD {script} {
    global system
    set here [pwd]
    if { $here != $system(SCRDIR) } {
	cd $here
    } else {
	cd $system(PWD)
    }
    uplevel 1 eval $script
    cd $here
}


#------------------------------------------------------------------------
#****f* auxil/putsFlush
#  NAME
#    putsFlush -- Tcl "puts" + "flush"
#  USAGE
#    putsFlush ?-nonewline? ?channelId? string
#
#  DESCRIPTION
#    Identical to Tcl's puts, but invoke the flush immediately after.
#    See puts man-page of Tcl.
#********
#------------------------------------------------------------------------
			
proc putsFlush {args} {
    update; update idletask
    # puts ?-nonewline? ?channelId? string
    set ind 0
    set flags "" 
    if { [lindex $args $ind] == "-nonewline" } {
	set flags "-nonewline"
	incr ind
    }
    if { [llength [lrange $args $ind end]] == 1 } {
	set channel stdout
    } else {
	set channel [lindex $args $ind]
	incr ind
    }

    eval puts $flags $channel [lrange $args $ind end]
    flush $channel
}

#
# Tcl's file copy will copy the link instead of the file. If link has
# a relative filename value, that's will be a mass: correct for this.
#
proc fileCopy {src dst} {

    catch {set file [file readlink $src]}
    
    if { [info exists file] } {
	global system
	return [file copy -force [file join $system(PWD) $file] $dst]
    } else {
	return [file copy -force $src $dst]
    }
}

proc lineRead {var file script} {
    # PURPOSE
    #   Read entire file line-by-line and at each line execute a
    #   script at one level up.
    # ARGUMENTS
    # * var    -- name of variable where the content of line will be stored
    # * file   -- name of file to read
    # * script -- script to execute when line is read 
    #
    # CREDITS
    #   Based on fileutils::foreachLine from tcllib (almost verbatim).
    # SOURCE
    upvar $var line

    set fid    [open $file r]
    set code   0
    set result {}

    while { ! [eof $fid] } {
        gets $fid line
        set code [catch {uplevel 1 $script} result]
        if {($code != 0) && ($code != 4)} { 
            break 
        }
    }
    close $fid

    if { ($code == 0) || ($code == 3) || ($code == 4) } {
        return $result
    }
    if { $code == 1 } {
        global errorCode errorInfo
        return \
            -code      $code      \
            -errorcode $errorCode \
            -errorinfo $errorInfo \
            $result
    }
    return -code $code $result
}




# Purpose: returns all the descendents of the given window (including
# itself)
proc getAllDescendantWid {w} {
    global getAllDescendantWid_list
    
    if { [info exists getAllDescendantWid_list] } {
	set getAllDescendantWid_list ""
    }
    
    return [getAllDescendantWid_ $w]
}
proc getAllDescendantWid_ {wlist} {
    global getAllDescendantWid_list
    
    foreach w $wlist {	
	if { ![winfo exists $w] } continue
	
	lappend getAllDescendantWid_list $w
	
	set children [winfo children $w]
	
	if { $children != "" } {
	    foreach child $children {
		getAllDescendantWid_ $child		
	    }
	}
    }
    return $getAllDescendantWid_list
}


# set a variable only if it does not exist
proc ifset {varName value} {
    upvar 1 $varName var

    if { ! [info exists var] } {
        uplevel 1 $script
    } 
}
