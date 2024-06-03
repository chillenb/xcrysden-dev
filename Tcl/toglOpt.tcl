#------------------------------------------------------------------------
# Set default options for Togl widget
#------------------------------------------------------------------------
proc toglOptDefault {} {
    global toglOpt tcl_platform

    set toglOpt(rgba)           true  
    set toglOpt(redsize)        1     
    set toglOpt(greensize)      1     
    set toglOpt(bluesize)       1     
    set toglOpt(double)         true  
    set toglOpt(depth)          true  
    set toglOpt(depthsize)      1     
    set toglOpt(accum)          true
    set toglOpt(accumredsize)   1     
    set toglOpt(accumgreensize) 1     
    set toglOpt(accumbluesize)  1     
    set toglOpt(accumalphasize) 1     
    set toglOpt(alpha)          false 
    set toglOpt(alphasize)      1     
    set toglOpt(stencil)        false 
    set toglOpt(stencilsize)    1     
    set toglOpt(auxbuffers)     0     
    set toglOpt(overlay)        false 
    set toglOpt(stereo)         none 
    set toglOpt(time)           100

    # on Mac OS X/X11, Cygwin, and Windows-WSL xcrysden crashes if
    # toglOpt(accum) = true ...
    
    if { [string match -nocase *Darwin* $tcl_platform(os)] \
             || [string match -nocase *CYGWIN* $tcl_platform(os)] \
             || [string match -nocase *Microsoft* $tcl_platform(osVersion)] \
         } {
	# disable accumulation buffer
        
        set toglOpt(accum) false
    }

    # Togl2.0 fails to print-to-file for Cygwin and Windows-WSL; use
    # Imagemagick for printing instead
    
    if { [string match -nocase *CYGWIN* $tcl_platform(os)] \
             || [string match -nocase *Microsoft* $tcl_platform(osVersion)] \
         } {
        global printSetup
        set printSetup(dumpWindow) 1
    }   
}
