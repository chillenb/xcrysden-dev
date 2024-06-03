#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/colormenu.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc ColorMenu {can mcolor} {
    global mody radio mesa_bg xcColors cm

    set cm(can) $can

    ###########################################################################
    # COLOR MENU
    set mesa_bg(black)      {0.0 0.0 0.0 1.0}
    set mesa_bg(white)      {1.0 1.0 1.0 1.0}
    set mesa_bg(red)        {0.6 0.0 0.0 1.0}
    set mesa_bg(green)      {0.0 0.6 0.0 1.0}
    set mesa_bg(blue)       {0.0 0.0 0.6 1.0}
    set mesa_bg(darkcyan)   {0.0 0.362 0.410 1.0}
    set mesa_bg(lightbrown) {0.699 0.684 0.258 1.0}
    set mesa_bg(default) [concat [rgb_h2f $xcColors(normal_bg)] 1.0]

    set mesa_bg(hx_black)      [rgb_f2h $mesa_bg(black)]        
    set mesa_bg(hx_white)      [rgb_f2h $mesa_bg(white)]   
    set mesa_bg(hx_red)        [rgb_f2h $mesa_bg(red)]     
    set mesa_bg(hx_green)      [rgb_f2h $mesa_bg(green)]   
    set mesa_bg(hx_blue)       [rgb_f2h $mesa_bg(blue)]    
    set mesa_bg(hx_darkcyan)   [rgb_f2h $mesa_bg(darkcyan)]
    set mesa_bg(hx_lightbrown) [rgb_f2h $mesa_bg(lightbrown)]
    set mesa_bg(hx_default)    [rgb_f2h $mesa_bg(default)] 

    set mesa_bg(hx_ac_black)      [rgb_ac_f2h $mesa_bg(black)]        
    set mesa_bg(hx_ac_white)      [rgb_ac_f2h $mesa_bg(white)]   
    set mesa_bg(hx_ac_red)        [rgb_ac_f2h $mesa_bg(red)]     
    set mesa_bg(hx_ac_green)      [rgb_ac_f2h $mesa_bg(green)]   
    set mesa_bg(hx_ac_blue)       [rgb_ac_f2h $mesa_bg(blue)]    
    set mesa_bg(hx_ac_darkcyan)   [rgb_ac_f2h $mesa_bg(darkcyan)]
    set mesa_bg(hx_ac_lightbrown) [rgb_ac_f2h $mesa_bg(lightbrown)] 
    set mesa_bg(hx_ac_default)    [rgb_ac_f2h $mesa_bg(default)] 

    if ![info exists radio($can,bg)] {
	set radio($can,bg) $mesa_bg(hx_ac_black)
    }

    $mcolor add radiobutton -label "black" \
	    -command [list eval {xc_newvalue $cm(can) \
	    $mody(L_BACKGROUND) } $mesa_bg(black)]\
	    -variable radio($can,bg) \
            -value $mesa_bg(hx_ac_black) \
	    -background $mesa_bg(hx_black) -foreground #ffffff \
	    -activebackground $mesa_bg(hx_ac_black) \
	    -activeforeground #ffffff
    # Date: Mon Mar 22 19:39:01 CET 1999
    # t.k: I stoped here; update all values to $mesa_bg(hx_ac_XXX)
    ##############################################################
    $mcolor add radiobutton -label "white" \
	    -command [list eval {xc_newvalue $cm(can) \
	    $mody(L_BACKGROUND)} $mesa_bg(white)] \
	    -variable radio($can,bg) \
            -value $mesa_bg(hx_ac_white) \
	    -background $mesa_bg(hx_white) -foreground #000000 \
	    -activebackground $mesa_bg(hx_ac_white)
    $mcolor add radiobutton -label "red" \
	    -command [list eval {xc_newvalue $cm(can) \
	    $mody(L_BACKGROUND)} $mesa_bg(red)]\
	    -variable radio($can,bg) \
            -value $mesa_bg(hx_ac_red) \
	    -background $mesa_bg(hx_red) -foreground #ffffff \
	    -activebackground $mesa_bg(hx_ac_red)
    $mcolor add radiobutton -label "green" \
	    -command [list eval {xc_newvalue $cm(can) \
	    $mody(L_BACKGROUND)} $mesa_bg(green)]\
	    -variable radio($can,bg) \
            -value $mesa_bg(hx_ac_green) \
 	    -background $mesa_bg(hx_green) -foreground #ffffff \
	    -activebackground $mesa_bg(hx_ac_green)
    $mcolor add radiobutton -label "blue" \
	    -command [list eval {xc_newvalue $cm(can) \
	    $mody(L_BACKGROUND)} $mesa_bg(blue)]\
	    -variable radio($can,bg) \
            -value $mesa_bg(hx_ac_blue) \
	    -background $mesa_bg(hx_blue) -foreground #ffffff \
	    -activebackground $mesa_bg(hx_ac_blue)
    $mcolor add radiobutton -label "dark cyan" \
	    -command [list eval {xc_newvalue $cm(can) \
	    $mody(L_BACKGROUND)} $mesa_bg(darkcyan)]\
	    -variable radio($can,bg) \
            -value $mesa_bg(hx_ac_darkcyan) \
	    -background $mesa_bg(hx_darkcyan) -foreground #ffffff \
	    -activebackground $mesa_bg(hx_ac_darkcyan)
    $mcolor add radiobutton -label "light brown" \
	    -command [list eval {xc_newvalue $cm(can) \
	    $mody(L_BACKGROUND)} $mesa_bg(lightbrown)]\
	    -variable radio($can,bg) \
            -value $mesa_bg(hx_ac_lightbrown) \
	    -background $mesa_bg(hx_lightbrown) -foreground #ffffff \
	    -activebackground $mesa_bg(hx_ac_lightbrown)
    $mcolor add radiobutton -label "default" \
	    -command [list eval {xc_newvalue $cm(can) \
	    $mody(L_BACKGROUND)} $mesa_bg(default)] \
	    -variable radio($can,bg) \
            -value $mesa_bg(hx_ac_default) \
	    -background $mesa_bg(hx_default) \
	    -activebackground $mesa_bg(hx_ac_default)
    $mcolor add separator 
    $mcolor add command -label "Custom ..." \
	    -command [list xcMesaChangeBg $can]
    trace variable radio($can,bg) w xcTrace
}