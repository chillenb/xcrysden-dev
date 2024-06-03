proc FS_displayBandWidths {bandWidthFile {spin ""}} {
    global wn

    set fID [open $bandWidthFile r]
    
    set nbands [gets $fID]
    set wn($spin,nbands) $nbands
    set    text "Number of bands: $nbands\n"
    append text "---------------\n\n"
    #append text "# Band  No.           Min-ENE         Max-ENE\n"
    while { ! [eof $fID] } {
	set line [gets $fID]
	if { [llength $line] == 3 } {
	    set i [lindex $line 0]
	    set e [lindex $line 1]
	    set E [lindex $line 2]
	    append text [format "Band No.: %3d;     Min-ENE: %12.6f     Max-ENE: %12.6f\n" $i $e $E]
	}
    }
    close $fID

    #
    # display the bandwidths in a Text widget
    #
    return [xcDisplayVarText $text "Band widths"]
}
