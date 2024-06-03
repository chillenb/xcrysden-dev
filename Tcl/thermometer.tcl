#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/thermometer.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

set tc(GEOGRAPHIC) {
    {   0./255.    0./255.  204./255.  1.0 } 
    { 150./255.  153./255.  255./255.  1.0 } 
    { 192./255.  216./255.  255./255.  1.0 } 
    {   0./255.  192./255.    3./255.  1.0 } 
    { 103./255.  255./255.   51./255.  1.0 } 
    { 206./255.  255./255.   51./255.  1.0 } 
    { 255./255.  255./255.  124./255.  1.0 } 
    { 223./255.  191./255.  143./255.  1.0 } 
    { 157./255.   48./255.  152./255.  1.0 } 
    { 255./255.    0./255.  255./255.  1.0 } 
    { 255./255.  255./255.  255./255.  1.0 }
}

set tc(RAINBOW) {
    { 1.0  0.0  0.0  1.0 } 
    { 1.0  1.0  0.0  1.0 } 
    { 0.0  1.0  0.0  1.0 } 
    { 0.0  1.0  1.0  1.0 } 
    { 0.0  0.0  1.0  1.0 } 
    { 0.97 0.0  1.0  1.0 }
}

set tc(RGB) {
    { 1.0  0.0  0.0  1.0 } 
    { 0.0  1.0  0.0  1.0 } 
    { 0.0  0.0  1.0  1.0 }
}

set tc(MONOCHROME) {
    { 0.1  0.1  0.1  1.0 } 
    { 1.0  1.0  1.0  1.0 }
};

set tc(BLUE-WHITE-RED) {
  {0.0  0.0  1.0  1.0} 
  {0.6  0.6  1.0  1.0} 
  {1.0  1.0  1.0  1.0} 
  {1.0  0.6  0.6  1.0} 
  {1.0  0.0  0.0  1.0}
}

set tc(BLACK-BROWN-WHITE) {
  {0.00  0.00  0.00  0.75} 
  {0.48  0.14  0.02  0.75} 
  {0.78  0.63  0.07  0.75} 
  {1.00  1.00  1.00  0.75}
}

proc funcLinear  {x} { return $x }
proc funcLog     {x} { return [expr log($x)] }
proc funcExp     {x} { return [expr exp($x)] }
proc funcSqrt    {x} { return [expr sqrt($x)] }
proc funcPow2    {x} { return [expr $x*$x] }
proc func3thRoot {x} { return [expr exp(log($x)/3.0)] }
proc funcPow3    {x} { return [expr $x*$x*$x] }
proc funcExp2    {x} { return [expr exp($x*$x)] }    
proc funcSqrtLog {x} { return [expr sqrt(log(exp($x)))] }

proc thermometerSetIntv {termoVar scaleFunc colorBasis from to nstep} {    
    global tc    
    upvar #0 $termoVar termo
    
    if { $scaleFunc == "LINEAR"   } { 
	set function     funcLinear  
	set invFunction  funcLinear
    }
    if { $scaleFunc == "LOG10"    } { 
	set function     funcLog     
	set invFunction  funcExp
    }
    if { $scaleFunc == "SQRT"     } { 
	set function     funcSqrt    
	set invFunction  funcPow2
    }
    if { $scaleFunc == "3th-ROOT" } { 
	set function     func3thRoot 
	set invFunction  funcPow3
    }
    if { $scaleFunc == "EXP(x)"   } { 
	set function     funcExp     
	set invFunction  funcLog
    }
    if { $scaleFunc == "EXP(x^2)" } { 
	set function     funcExp2    
	set invFunction  funcSqrtLog
    }

    set shift 0.0
    if { ( \
	       $scaleFunc == "LOG10" || \
	       $scaleFunc == "SQRT" || \
	       $scaleFunc == "3th-ROOT" || \
	       $scaleFunc == "EXP(x^2)" ) && $from < 1.0e-7 } {
	set shift [expr -$from + 1.0e-7]
	set from  1.0e-7
	set to    [expr $shift + $to]
    }

    set min   [$function $from]
    set max   [$function $to]
    set range [expr $max - $min]

    #
    # first the basic intervals
    #
    set i    0
    set N    [llength $tc($colorBasis)]
    set N1   [expr $N - 1]
    set STEP [expr double($range) / double($N1)]
    set termo(STEP) $STEP
    foreach color $tc($colorBasis) {
	for {set j 0} {$j < 3} {incr j} {
	    set termo(basic_col,$i,$j) [lindex $color $j]
	}
	set termo(basic_intv,$i) [expr double($i) * $STEP]
	incr i
    }
    
    #
    # now the actual intervals
    #    
    set nstep1 [expr $nstep - 1]
    set step   [expr double($range) / double($nstep1)]
    for {set i 0} {$i < $nstep} {incr i} {
	set termo(intv,$i) [expr double($i) * $step]
	set termo(tic,$i)  [$invFunction \
				[expr $min + $termo(intv,$i) - $shift]]
    }
}


proc thermometerWidget {w iplane labeltext font fmt colorbasis scalefunc from to nstep {toplevel 0}} {
    global tc
    ## IN/OUT: w            name of the widget
    ## IN:     iplane       planes index (i.e. 1/2/3)
    ## IN:     labeltext    text of label
    ## IN:     font         display font-name
    ## IN:     fmt          format for the numeric-labels
    ## IN:     colorbasis   which colorbasis
    ## IN:     from         starting value
    ## IN:     to           ending   value
    ## IN:     nstep        number of steps

    global $w
    upvar #0 $w termo
    
    if { $toplevel == 0 } {
	frame $w -relief solid -bd 2 -bg \#ffffff
    } else {
	xcToplevel $w "Thermometer: $iplane" "Thermometer: $iplane"
    }

    set frame_atrib "-bg \#ffffff -relief flat -bd 0"    
    set con [eval frame $w.f $frame_atrib]
    pack $con -fill both -expand 1
    
    set canvas   [eval frame $con.f -bg \#ffffff -relief flat -bd 2]    
    
    if { $font == {} } {
	set font fixed
    }
    label $con.l \
	-text $labeltext \
	-font $font \
	-bg   \#ffffff \
	-relief solid \
	-bd 2
    
    pack $con.l -side top -fill y -padx 5 -pady 5 -expand 1 -ipadx 5    
    set fontsize [font actual $font -size]
    if { $fontsize < 0 } {
	set scale    [tk scaling]
	set fontsize [expr round(abs($fontsize*$scale))]
    }
    
    # calculate the intervals
    thermometerSetIntv $w $scalefunc $colorbasis $from $to $nstep
    
    set N    [llength $tc($colorbasis)]
    set N1   [expr $N - 1]

    #
    # now plot the Thermometer
    #
    set nstep1  [expr $nstep - 1]
    for {set i 0} {$i < $nstep} {incr i} {
    	set value $termo(intv,$i)
    
    	for {set j 0} {$j < $N1} {incr j} {
    	    set j1 [expr $j + 1]
	    
    	    if { $value >= $termo(basic_intv,$j) && 
    		 $value <= $termo(basic_intv,$j1) } {

		set f [expr ($value - $termo(basic_intv,$j)) / $termo(STEP)]
		
    		for {set k 0} {$k < 3} {incr k} {		    
    		    set delta [expr \
    				   $termo(basic_col,$j1,$k) - \
    				   $termo(basic_col,$j,$k)]
    		    set col($i,$k) [expr \
    					$termo(basic_col,$j,$k) + $f * $delta]
    		}		
    		set rgb($i) [format "\#%02x%02x%02x" \
    				 [expr int($col($i,0)*255.0)] \
    				 [expr int($col($i,1)*255.0)] \
    				 [expr int($col($i,2)*255.0)] ]		
    		break
    	    }
    	}
    	if { ![info exists rgb($i)] } {
    	    # if we come here, the assign the "last" color
    	    set rgb($i) [format "\#%02x%02x%02x" \
    			     [expr int($termo(basic_col,$N1,0)*255.0)] \
    			     [expr int($termo(basic_col,$N1,1)*255.0)] \
    			     [expr int($termo(basic_col,$N1,2)*255.0)] ]
    	}
    
    	set container [eval frame $canvas.f$i $frame_atrib]
    	set cof [frame $container.c \
    		     -width      $fontsize \
    		     -height     $fontsize \
    		     -background $rgb($i) \
    		     -relief     solid \
    		     -bd         2]
    	set text [format $fmt $termo(tic,$i)]
    	set lab [eval label $container.l \
    		     -text $text \
    		     -justify left \
    		     -font $font \
    		     -width 10 \
    		     -bg \#ffffff]
    	
    	pack $cof -padx 3 -pady 1 -side left -anchor e
    	pack $lab -padx 3 -pady 1 -side left -fill x -anchor e
	pack $container -side top -fill x \
	    -padx 0 -pady 0 -ipadx 0 -ipady 0 -anchor e
    }
    pack $canvas -padx 5 -pady 5 -ipadx 2 -ipady 2 -anchor s -fill both    

    if { $toplevel == 1 } {
	button $w.b -text "Print Thermo" \
	    -font $font -command [list dumpWindow $con]
	pack $w.b -side top -fill x -expand 1
    }
    return $w
}
