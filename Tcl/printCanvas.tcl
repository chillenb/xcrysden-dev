#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/printCanvas.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc xcPrintCanvas {can {args {}}} {
    global printCanvas system
    # first query canvas dimension
    # WARNING:: if canvas is scrolled we will get just viewable dimesions;
    #           in that case user can specified scroll-dimensions!!!
    if [info exists printCanvas(printCommand)] {
	set com $printCanvas(printCommand)
    }
    if [array exists printCanvas] {
	unset printCanvas
    }
    if [info exists com] {
	set printCanvas(printCommand) $com
    }

    set printCanvas(gnuplot_print) 0
    if { $can == "none" } {
	set printCanvas(c_w) 300
	set printCanvas(c_h) 300
    } else {
	set printCanvas(c_w) [$can cget -width]
	set printCanvas(c_h) [$can cget -height]
    }
    set i 0
    foreach option $args {
	incr i
	# odd cycles are tags, even options
        if { $i%2 } {
            set tag $option
        } else {
	    #
	    # practically non of the option is supported yet
	    #
            switch -exact -- $tag {
		-canvas_scrollwidth  { set printCanvas(c_w)        $option }
		-canvas_scrollheight { set printCanvas(c_h)        $option }
		-colormap            { set printCanvas(colormap)   $option }
		-colormode           { set printCanvas(colormode)  $option }
		-file                { set printCanvas(file)       $option } 
		-fontmap             { set printCanvas(fontmap)    $option }
		-height              { set printCanvas(height)     $option }
		-pageanchor          { set printCanvas(pageanchor) $option }
		-pageheight          { set printCanvas(pageheight) $option }
		-pagewidth           { set printCanvas(pagewidth)  $option }
		-pagex		     { set printCanvas(pagex)      $option }
		-pagey		     { set printCanvas(pagey)      $option }
		-rotate		     { set printCanvas(rotate)     $option }
		-width		     { set printCanvas(width)      $option }
		-x		     { set printCanvas(x)	   $option }
		-y                   { set printCanvas(y)          $option }
		-gnuplot_print       { set printCanvas(gnuplot_print) $option }
		default { tk_dialog [WidgetName] Error \
			"ERROR: Bad xcPrintCanvas configure option $tag" \
			error 0 OK }
	    }
	}
    }

    if { $i%2 } {
	tk_dialog [WidgetName] Error "ERROR: You called xcPrintCanvas with an odd number of args !" \
		error 0 OK
	return 0
    }

    if { $can == "none" && $printCanvas(gnuplot_print) != 1 } {
	tk_dialog [WidgetName] Error "ERROR: xcPrintCanvas was misused; code: none && !gnuplot_print" error 0 OK
	return 0
    }

    #
    # how much is one screen-pixel
    #
    #set sp2cm [winfo fpixels $can 1c]
    #set sp2in [winfo fpixels $can 1i]
    #puts "$sp2cm $sp2in"

    #
    # initialize variables
    #
    set printCanvas(margin) 2c
    if ![info exists printCanvas(paperSize)] {
	set printCanvas(paperFormat) A4
	set printCanvas(paperSize) {209.9m 297m}
    }
    if ![info exists printCanvas(paperOrient)] {
	set printCanvas(paperOrient) Portrait
    }
    if ![info exists printCanvas(printTo)] {
	set printCanvas(printTo) Printer
    }
    if ![info exists printCanvas(printCommand)] {
	set printCanvas(printCommand) lpr
    }
    
    set oldgrab [grab current]
    set t [xcToplevel [WidgetName] Print Print . 100 100 1]
    catch { grab $t }
    set printCanvas(toplevel) $t
    set f1 [frame $t.f1]
    set f2 [frame $t.f2 -relief raised -bd 2]
    set f3 [frame $t.f3 -relief raised -bd 2]
    
    pack $f1 $f2 $f3 -side top -expand 1 -fill both -ipady 5
    
    #
    # upper frame - f1: left -> canvas; right -> options
    #
    set f1l [frame $f1.l -relief raised -bd 2]
    set f1r [frame $f1.r]
    pack $f1l $f1r -side left -expand 1 -fill both -ipady 5
    # f1-left
    set c [canvas $f1l.c -width 200 -height 200 -relief sunken -bd 2]
    set printCanvas(can) $c
    if { $can == "none" } {
	set can $c
    }
    # create background rectangle
    $c create rectangle -2 -2 210 210 \
	-fill    "#999" \
	-outline "#999" \
	-stipple  gray25 \
	-tags "bg"
    
    set f1l1 [frame $f1l.1]
    set min [button $f1l1.min \
		 -bitmap "@$system(BMPDIR)/xcPrintCanvas_left.xbm" \
		 -command [list xcPrintDrawPic $c lt]]
    set maj [button $f1l1.maj \
		 -bitmap "@$system(BMPDIR)/xcPrintCanvas_right.xbm" \
		 -command [list xcPrintDrawPic $c gt]]
    set max [button $f1l1.max -text "Max" \
		 -command [list xcPrintDrawPic $c max]]
    pack $c $f1l1 -side top -expand 1  -padx 10
    pack $min $maj $max -side left -expand 1 -padx 10

    # f1-right: up
    set f1ru [frame $f1r.up -relief raised -bd 2]
    pack $f1ru -side top  -fill both -ipady 5 -expand 1
    RadioButCmd $f1ru {Paper Size:} printCanvas(paperFormat) \
	xcPrintPaperFormat top left 1 0 3 \
	A4 {8.5"x11"}
    RadioButCmd $f1ru {Paper Orientation:} printCanvas(paperOrient) \
	xcPrintOrient top left 1 0 3 \
	Portrait Landscape
    set f1ruu [frame $f1ru.fu -relief groove -bd 2]
    pack $f1ruu -side top -padx 3 -pady 3 -expand 1 -fill both
    set f1ru1 [frame $f1ru.f1]
    set f1ru2 [frame $f1ru.f2]
    set f1ru3 [frame $f1ru.f3]
    pack $f1ru1 $f1ru2 $f1ru3 -in $f1ruu -side top -fill both -expand 1
    set lt [label $f1ru1.lt -text "Picture Size:"]
    pack $lt -side top
    set ls [label $f1ru1.ls -text "Current Size:"]
    set es [label $f1ru1.es \
		-textvariable printCanvas(currentSize)]
    pack $ls $es -side left -pady 5

    set sw [label $f1ru2.ss -text "Set Width:" -width 11]
    set ew [entry $f1ru2.ew -textvariable printCanvas(cw_m) -width 3]
    set lw [label $f1ru2.lw -text "mm"]
    bind $ew <Return> [list xcPrintSetSize $c $ew]
    
    if { $printCanvas(gnuplot_print) != 1 } {	
	set upd [button $f1ru2.upd -text "Update" \
		     -width 4 \
		     -command [list xcPrintSetSize $c $ew]] 
	pack $sw $ew $lw $upd -side left -pady 5
    } else {
	set sh [label $f1ru3.ss -text "Set Height:" -width 11]
	set eh [entry $f1ru3.ew -textvariable printCanvas(ch_m) -width 3]
	set lh [label $f1ru3.lw -text "mm"]
	set upd [button $f1ru3.upd -text "Update" \
		     -width 4 \
		     -command [list xcPrintSetSize $c $ew $eh]]
	bind $eh <Return> [list xcPrintSetSize $c $ew $eh]
	pack $sw $ew $lw -side left -pady 5
	pack $sh $eh $lh $upd -side left -pady 5
    }
    
    # f1-right: down
    set f1rd [frame $f1r.dn -relief raised -bd 2]
    pack $f1rd -side top -fill x -ipady 3    
    RadioButCmd $f1rd {Print to:} printCanvas(printTo) \
	xcPrintPrintTo top left 1 1 3 \
	Printer File

    #
    # frame f2: Print Command
    #
    set f21 [frame $f2.1]
    set f22 [frame $f2.2]
    set f22l [frame $f22.l]
    set f22r [frame $f22.r]
    pack $f21 $f22 -side top -expand 1 -fill both
    pack $f22l -side left -expand 1 -fill both
    pack $f22r -side left -fill both
    
    set printCanvas(printWid) $f21
    set printCanvas(fileWid)  $f22
    set printCanvas(printEntry) \
	[FillEntries $f21 {"Print Command:"} printCanvas(printCommand) \
	     15 {} top left]
    set printCanvas(fileEntry) \
	[FillEntries $f22l {"Filename:"} printCanvas(file) \
	     15 {} top left]
    set types {
	{{All Files}       *            }
	{{PS  Files}       {.ps}        }
	{{EPS Files}       {.eps}       }
    }
    set fn [button $f22r.b1 -text "Filename" \
		-command [list xcPrintFilename $types]]
    pack $fn -side top -padx 5

    #
    # frame 3: Cancel & Print
    #
    set cancel [button $f3.can -text Cancel -command [list CancelProc $t]]
    set print  [DefaultButton $f3.print -text Print \
	    -done_var printCanvas(done)]
    pack $cancel $print -side left -padx 5 -expand 1

    #
    # initial setup
    #
    xcPrintPrintTo $printCanvas(printTo)
    xcPrintDrawPaper $c
    xcPrintDrawPic $c max

    tkwait variable printCanvas(done)
    #
    # now print postscript
    #
    if { $printCanvas(paperOrient) == "Portrait" } {
	set printCanvas(rotate) 0
    } else {
	set printCanvas(rotate) 1
    }

    #
    # maybe we want to print using gnuplot
    #
    if $printCanvas(gnuplot_print) {
	GnuplotPrint
    } else {
	xcDebug "PictureSize:: [xcPrintPix2MM $printCanvas(pic_w)]m"
	if { $printCanvas(printTo) == "Printer" } {
	    set ps [$can postscript \
			-pagewidth [xcPrintPix2MM $printCanvas(pic_w)]m \
			-width  $printCanvas(c_w) \
			-height $printCanvas(c_h) \
			-rotate $printCanvas(rotate)]
	    if { [catch {eval exec $printCanvas(printCommand) {<<\n$ps}} error] } {
		tk_dialog [WidgetName] ERROR \
		    "ERROR: An Error has occurred when trying to print: $error" error 0 OK
		destroy $t
		return 0
	    }
	} else {
	    # parse $printCanvas(file)
	    set file [GetAbsoluteFileName $printCanvas(file)]
	    $can postscript \
		    -pagewidth [xcPrintPix2MM $printCanvas(pic_w)]m \
		    -width  $printCanvas(c_w) \
		    -height $printCanvas(c_h) \
		    -file   $file \
		    -rotate $printCanvas(rotate)
	}
    }

    catch { grab release $t }
    destroy $t

    if { $oldgrab != "" } {
	catch { grab $oldgrab }
    }
}


proc xcPrintDrawPaper {can} {
    global printCanvas
    
    set ration [expr 190 / [winfo fpixels $can 297m]]

    if { $printCanvas(paperOrient) == {Portrait} } {
	set pw [expr $ration * \
		[winfo fpixels $can [lindex $printCanvas(paperSize) 0]]]
	set ph [expr $ration * \
		[winfo fpixels $can [lindex $printCanvas(paperSize) 1]]]
    } else {
	set pw [expr $ration * \
		[winfo fpixels $can [lindex $printCanvas(paperSize) 1]]]
	set ph [expr $ration * \
		[winfo fpixels $can [lindex $printCanvas(paperSize) 0]]]
    }

    set bbox [xcPrintBBox $pw $ph]
    set x1 [lindex $bbox 0]
    set y1 [lindex $bbox 1]
    set x2 [lindex $bbox 2]
    set y2 [lindex $bbox 3]
    $can delete paper
    $can create rectangle \
	    [expr $x1 + 3] [expr $y1 + 3] \
	    [expr $x2 + 3] [expr $y2 + 3] \
	    -fill "#000" \
	    -outline "#000" \
	    -tags "paper back"
    $can create rectangle $x1 $y1 $x2 $y2 \
	    -fill "#ffffff" \
	    -outline "#000" \
	    -tags "paper front"
}

proc xcPrintDrawPic {can {mode {}}} {
    global printCanvas
    

    set ration [expr 190 / [winfo fpixels $can 297m]]

    switch -exact -- $mode { 
	gt {
	    # magnify by 5%
	    set printCanvas(pic_w) [expr $printCanvas(pic_w) * 1.05]
	    set printCanvas(pic_h) [expr $printCanvas(pic_h) * 1.05]
	}
	lt {
	    # shrink
	    set printCanvas(pic_w) [expr $printCanvas(pic_w) / 1.05]
	    set printCanvas(pic_h) [expr $printCanvas(pic_h) / 1.05]
	}
	max {
	    set pw [expr $ration * \
		    [winfo fpixels $can [lindex $printCanvas(paperSize) 0]]]
	    set ph [expr $ration * \
		    [winfo fpixels $can [lindex $printCanvas(paperSize) 1]]]
	    
	    set maxw [expr $pw - $ration * [winfo fpixels $can \
		    $printCanvas(margin)]]
	    set maxh [expr $ph - $ration * [winfo fpixels $can \
		    $printCanvas(margin)]]
	    set pic_w [expr $ration * \
		    [winfo fpixels $can $printCanvas(c_w)]]
	    set pic_h [expr $ration * \
		    [winfo fpixels $can $printCanvas(c_h)]]
	    
	    set wr [expr $pic_w / $pw]
	    set hr [expr $pic_h / $ph]
	    # is it higher or wider
	    if { $wr > $hr } {
		# it is wider
		set picr [expr $maxw / $pic_w]
	    } else {
		set picr [expr $maxh / $pic_h]
	    }
	    set printCanvas(pic_w) [expr $picr * $pic_w] 
	    set printCanvas(pic_h) [expr $picr * $pic_h] 
	}
    }

    # update printCanvas(currentSize)
    set w [xcPrintPix2MM $printCanvas(pic_w)]
    set h [xcPrintPix2MM $printCanvas(pic_h)]
    xcDebug "$printCanvas(pic_w) $printCanvas(pic_h)"
    set printCanvas(currentSize) [format "%3.0f%2s x %3.0f%2s" $w mm $h mm]
    set bbox [xcPrintBBox $printCanvas(pic_w) $printCanvas(pic_h)]
    $can delete picture
    xcDebug "Pic: $bbox $mode"
    eval {$can create rectangle} $bbox {\
	    -fill "#0f0" \
	    -outline "#000" \
	    -stipple  gray25 \
	    -tags "picture"}
}
    

proc xcPrintPaperFormat item {
    global printCanvas
    xcDebug "xcPrintPaperFormat: $item"
    if { $item == {A4} } {
	set printCanvas(paperSize) {209.9m 297m}
    } elseif { $item == {8.5"x11"} } { 
	set printCanvas(paperSize) {8.5i 11i}
    }
    xcPrintDrawPaper $printCanvas(can)
    xcPrintDrawPic   $printCanvas(can)
}


proc xcPrintOrient item {
    global printCanvas
    xcDebug "xcPrintOrient: $item"
    set printCanvas(paperOrient) $item
    xcPrintDrawPaper $printCanvas(can)
    xcPrintDrawPic   $printCanvas(can)
}


proc xcPrintPrintTo item {
    global printCanvas

    set disabled_color [lindex \
	    [GetWidgetConfig button -disabledforeground] end]
    set enabled_color  [lindex \
	    [GetWidgetConfig button -foreground] end]

    set fileLabel  "[string trimright $printCanvas(fileEntry) entry1]lab1"
    set printLabel "[string trimright $printCanvas(printEntry) entry1]lab1"

    if { $item == {Printer} } {
	xcDisableAll $printCanvas(fileWid)
	xcEnableAll  $printCanvas(printWid)
	$printCanvas(fileEntry)  configure -relief flat
	$printCanvas(printEntry) configure -relief sunken
	$fileLabel  configure -foreground $disabled_color
	$printLabel configure -foreground $enabled_color
    } elseif { $item == {File} } {
	xcDisableAll $printCanvas(printWid)
	xcEnableAll  $printCanvas(fileWid)
	$printCanvas(fileEntry)  configure -relief sunken
	$printCanvas(printEntry) configure -relief flat
	$fileLabel  configure -foreground $enabled_color
	$printLabel configure -foreground $disabled_color
    }
}

proc xcPrintSetSize {can args} {
    global printCanvas

    if { $printCanvas(gnuplot_print) != 1 } {
	set varlist {{printCanvas(cw_m) posreal}}
    } else {
	set varlist {{printCanvas(cw_m) posreal} {printCanvas(ch_m) posreal}}
    }
    
    if [check_var $varlist $args] {
	set pic_w [xcPrintMM2Pix $printCanvas(cw_m)]
	set r [expr $pic_w / double($printCanvas(pic_w))]
	set printCanvas(pic_w) [expr $r * $printCanvas(pic_w)]
	if { $printCanvas(gnuplot_print) != 1 } {
	    set printCanvas(pic_h) [expr $r * $printCanvas(pic_h)]
	} else { 
	    set pic_h [xcPrintMM2Pix $printCanvas(ch_m)]
	    set r [expr $pic_h / double($printCanvas(pic_h))]
	    set printCanvas(pic_h) [expr $r * $printCanvas(pic_h)]
	}
	xcPrintDrawPic   $can
    }
}


proc xcPrintPix2MM x {
    return [expr 297. * double($x) / 190.]
}


proc xcPrintMM2Pix x {
    return [expr 190. * double($x) / 297.]
}


proc xcPrintBBox {pw ph} {
    set x1 [expr (200 - $pw) / 2]
    set y1 [expr (200 - $ph) / 2]
    set x2 [expr $x1 + $pw]
    set y2 [expr $y1 + $ph]

    return [list $x1 $y1 $x2 $y2]
}


proc xcPrintFilename types {
    global printCanvas system

    set pwd [pwd]
    cd $system(PWD)
    set printCanvas(file) [tk_getSaveFile \
	    -title {Specify Filename} \
	    -filetypes $types \
	    -parent $printCanvas(toplevel)]    
    cd $pwd
}

#lappend auto_path "/home/tone/src/xcrysden0.0/src"
#lappend auto_path [pwd]
#
#set system(TOPDIR) [pwd]
#set system(SCRDIR) /tmp
#set system(PWD)    [pwd]
#set system(PID)    [pid]
#set system(FORDIR) $system(TOPDIR)/F
#set system(TCLDIR) $system(TOPDIR)
#set system(BINDIR) $system(TOPDIR)
#set system(BMPDIR) $system(TOPDIR)/bitmap
#button .b -text TONE
#entry  .e -text TONE
#set xcFonts(normal)       [lindex [.b configure -font] end]
#set xcFonts(normal_entry) [lindex [.e configure -font] end]
#set xcFonts(small)        [ModifyFontSize .b 10 \
#	 {-family helvetica -slant r -weight bold}]
#set xcFonts(small_entry)  [ModifyFontSize .e 10 \
#	 {-family helvetica -slant r -weight normal}]
#destroy .b .e
#
#set ft [frame .ft -relief sunken -bd 2]
#pack $ft -side top -expand true -fill y
#
#set c [canvas $ft.canv -yscrollcommand [list $ft.yscroll set]]
#set scb [scrollbar $ft.yscroll -orient vertical -command [list $c yview]]
#pack $scb -side right -fill y
#pack $c -side left -fill both -expand true
#	 
#$c config -width 300 -height 100 \
#	 -scrollregion "0 0 300 300"
#xcPrintCanvas $c 
