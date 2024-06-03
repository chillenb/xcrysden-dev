proc check_stereo {} {
    global stereo_visual

    set togl [togl [WidgetName .] -stereo native \
		  -rgba           true  \
		  -redsize        1     \
		  -greensize      1     \
		  -bluesize       1     \
		  -double         true  \
		  -depth          true  \
		  -depthsize      1     \
		  -accum          true  \
		  -accumredsize   1     \
		  -accumgreensize 1     \
		  -accumbluesize  1     \
		  -accumalphasize 1     \
		  -alpha          false \
		  -alphasize      1     \
		  -stencil        false \
		  -stencilsize    1     \
		  -auxbuffers     0     \
		  -overlay        false]
    pack $togl
    update

    set stereo [lindex [$togl configure -stereo] end]
    #set stereo [xc_stereo]	
    
    if { $stereo != "none" } {
    	set stereo_visual true
    } else {
    	set stereo_visual false
    }
 
    destroy $togl
}
