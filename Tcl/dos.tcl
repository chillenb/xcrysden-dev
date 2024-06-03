#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/dos.tcl
# ------                                                                    #
# Copyright (c) 1996--2004 by Anton Kokalj                                  #
#############################################################################


proc DOSSGraph { NPRO } {
    global prop ftn25 grafdata grafsize graf system
    
    set fileID [open $prop(dir)/$prop(file)25 r]
         
    #
    # query GrapherID
    #
    set gID [NextGrapherID]
    xcDebug "GrapherID:: -- $gID --"
    set Ymax 0.0
    if { $prop(type_of_run) == "UHF" } {
	set N_segm 2
    } else {
	set N_segm 1
    }

    # first PROJECTED DOS, if any
    set grafdata($gID,N_text) 0
    for {set j 1} {$j <= $N_segm} {incr j} {
	set i 1; # don't change
	for {} {$i <= $NPRO} {incr i} {
	    xcDebug "DOSSGraph:: i=$i, j=$j"
	    ReadFTN25 DOSS $fileID $gID $j $i	
	    if { $ftn25(Ymax) > $Ymax } { set Ymax $ftn25(Ymax) }
	    if { $prop(type_of_run) == "UHF" && abs($ftn25(Ymin)) > $Ymax } {
		set Ymax abs($ftn25(Ymin))
	    }
	}
	# now TOTAL DOS
	ReadFTN25 DOSS $fileID $gID $j $i
    }

    for {set i 1} {$i <= $NPRO} {incr i} {
	# someday makes this better
	set atoms ""
	set AOs ""
	if { $prop(N,$i) < 0 } {
	    append atoms "$prop(NDM,$i) "	    
	    set grafdata($gID,text$i)  \
		    "DOS projected to atoms: $atoms"
	    incr grafdata($gID,N_text)
	} else {
	    append AOs "$prop(NDM,$i) "	    
	    set grafdata($gID,text$i)  \
		    "DOS projected to AOs: $AOs"
	    incr grafdata($gID,N_text)
	}
    }


    if { $ftn25(Ymax) > $Ymax } { set Ymax $ftn25(Ymax) }
    if { $prop(type_of_run) == "UHF" && abs($ftn25(Ymin)) > $Ymax } {
	set Ymax abs($ftn25(Ymin))
    }
    set grafdata($gID,text$i) "Total DOS"
    incr grafdata($gID,N_text)
    close $fileID
    

    # Yoffset is 0.1 --> 1.1 
    set Yoffset                  1.1
    set dX                       $ftn25(X4)
    set grafdata($gID,N_graf)    [expr $NPRO + 1]; # N of graphs
    set grafdata($gID,type)      DOSS
    set grafdata($gID,X_title,[expr $NPRO + 1])  "E \[a.u.\]"
    set grafsize($gID,canW)      600
    if { $NPRO == 0 } { set grafsize($gID,canH) 400 }
    if { $NPRO == 1 } { set grafsize($gID,canH) 500 }
    if { $NPRO == 2 } { set grafsize($gID,canH) 600 }

    if { $NPRO <= 2 } {
	set grafsize($gID,Yscroll) 0
    } else {
	set grafsize($gID,Yscroll) [expr ($NPRO + 1) * 200]
    }	

    for {set g 1} {$g < $grafdata($gID,N_graf)} {incr g} {
	set grafdata($gID,Xtick_text_draw,$g)   0
    }

    set ftn25(X2) [FTN25_MustBeNumber $ftn25(X2) 0.0 \
	    "An Error has occured while reading unit 25; Energy corresponding to the first point can not be read. A value of 0.0 will be taken for first energy point instead !"]
    set ftn25(X4) [FTN25_MustBeNumber $ftn25(X4) 1.0 \
	    "An Error has occured while reading unit 25; Energy increment can not be read. A value of 1.0 will be taken for energy increment instead !"]
    set ftn25(X5) [FTN25_MustBeNumber $ftn25(X5) 0.0 \
	    "An Error has occured while reading unit 25; Fermi energy can not be read. A value of 0.0 will be taken for Fermi energy instead !"]
    for {set g 1} {$g <= $grafdata($gID,N_graf)} {incr g} {
	set grafdata($gID,Ytick_text_draw,$g) 0
	set grafdata($gID,Y_title,$g)   "   DOS\n(arb. un.)"
	set grafdata($gID,Ymax,$g)      [expr $Yoffset * $Ymax]
	set grafdata($gID,Yoffset,$g)   0.0
	set grafdata($gID,Xoffset,$g)   0.05	 
	set grafdata($gID,Xmin,$g)      $ftn25(X2); # min Energy
	set grafdata($gID,Xmax,$g)      \
		[expr $ftn25(X2) + $ftn25(NCOL) * $ftn25(X4)]
	set grafdata($gID,Xline1,$g)    $ftn25(X5);       # Fermi Energy
	set grafdata($gID,Xline2,$g)    0.5
	set grafdata($gID,N_Xline,$g)   1;                
	set grafdata($gID,N_point,$g)   $ftn25(NCOL);    # N of point per graph
	set grafdata($gID,Xline1_text,$g)        "E_Fermi"
	set grafsize($gID,Xline1_textstipple,$g) @$system(BMPDIR)/dotH.bmp
	set grafsize($gID,Xline1_fill,$g)        "#ff0000"
	set grafsize($gID,Xline1_width,$g)       1
	set grafsize($gID,Xline1_stipple,$g)     @$system(BMPDIR)/dotV.bmp
	set grafsize($gID,margin_X1)             $grafsize(margin_X1)
	set grafsize($gID,margin_X2)             $grafsize(margin_X2)
	set grafsize($gID,margin_Y1)             $grafsize(margin_Y1)
	set grafsize($gID,margin_Y2)             $grafsize(margin_Y2)
	set grafsize($gID,Yaxe_title_offset)     15
	if { $prop(type_of_run) == "UHF" } {
	    set grafdata($gID,Ymin,$g)           -$grafdata($gID,Ymax,$g) 
	    set grafdata($gID,N_Yline,$g)        1
	    set grafdata($gID,Yline1,$g)         0.0
	    set grafsize($gID,Yline1_fill,$g)    "#000000"
	    set grafsize($gID,Yline1_width,$g)   1
	    set grafsize($gID,Yline1_stipple,$g) {}
	    set grafdata($gID,N_segment,$g)      2
	} else {
	    # RHF
	    set grafdata($gID,Ymin,$g)        0.0
	    set grafdata($gID,N_Yline,$g)     0
	    set grafdata($gID,N_segment,$g)   1
	}
		
	set graf($gID,RorigX,$g) 0.0
	set graf($gID,RorigY,$g) \
		[expr double($NPRO + 1 - $g) / double($NPRO + 1)]
	set graf($gID,RsizeX,$g) 1.0
	set graf($gID,RsizeY,$g) \
		[expr 1.0 / double($NPRO + 1) - [XGraph2nMRel $gID 5]]
	xcDebug "origin:: $graf($gID,RorigX,$g) $graf($gID,RorigY,$g)"
	xcDebug "size::   $graf($gID,RsizeX,$g) $graf($gID,RsizeY,$g)"
	
	set grafdata($gID,text${g}_X) 0.02
	set grafdata($gID,text${g}_Y) \
		[expr double($g - 0.9) / double($NPRO + 1)]

	# X vales of points
	for {set i 1} {$i <= $grafdata($gID,N_point,$g)} {incr i} {
	    set grafdata($gID,X$i,$g) \
		    [expr $grafdata($gID,Xmin,$g) + $dX * ($i - 1)]
	}
	# take care of X ticks
	set grafdata($gID,N_mXtick,$g) 2
	set grafdata($gID,N_MXtick,$g) 5
	set dx [expr abs(($grafdata($gID,Xmax,$g) - \
		$grafdata($gID,Xmin,$g)) / \
		($grafdata($gID,N_MXtick,$g) - 1))]
	for {set i 1} {$i <= $grafdata($gID,N_MXtick,$g)} {incr i} {
	    set grafdata($gID,Xtick$i,$g) \
		    [expr $grafdata($gID,Xmin,$g) + ($i - 1) * $dx]
	
	    if { abs($grafdata($gID,Xtick$i,$g)) < 0.01 || \
		    abs($grafdata($gID,Xtick$i,$g)) >= 100 } {
		set grafdata($gID,Xtick${i}_text,$g) \
			[format "%6.3E" $grafdata($gID,Xtick$i,$g)]
	    } else {
		set grafdata($gID,Xtick${i}_text,$g) \
			[format "%6.3f" $grafdata($gID,Xtick$i,$g)]
	    }
	}
	# take care of X ticks
	set grafdata($gID,N_mYtick,$g) 2
	set grafdata($gID,N_MYtick,$g) 5
	set dy [expr abs(($grafdata($gID,Ymax,$g) - \
		$grafdata($gID,Ymin,$g)) / \
		($grafdata($gID,N_MYtick,$g) - 1))]
	for {set i 1} {$i <= $grafdata($gID,N_MYtick,$g)} {incr i} {
	    set grafdata($gID,Ytick$i,$g) \
		    [expr $grafdata($gID,Ymin,$g) + ($i - 1) * $dy]
	
	    if { abs($grafdata($gID,Ytick$i,$g)) < 0.01 || \
		    abs($grafdata($gID,Ytick$i,$g)) >= 100 } {
		set grafdata($gID,Ytick${i}_text,$g) \
			[format "%6.3E" $grafdata($gID,Ytick$i,$g)]
	    } else {
		set grafdata($gID,Ytick${i}_text,$g) \
			[format "%6.3f" $grafdata($gID,Ytick$i,$g)]
	    }
	}
    }
    
    xcDebug "Ymax, Ymin:: $grafdata($gID,Ymax,1), $grafdata($gID,Ymin,1)"
}
