#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/band.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc BANDGraph { NLINE {N_proj 1}} {
    global properties prop ftn25 grafdata grafsize graf system

    # ftn25(X1) ... Emin
    # ftn25(X2) ... Emax
    # ftn25(X3) ... not used
    # ftn25(X4) ... dK
    # ftn25(X5) ... Efermi
    set fileID [open $prop(dir)/$prop(file)25 r]
    
    #
    # query GrapherID
    #
    set gID [NextGrapherID]
    xcDebug "GrapherID:: -- $gID --"

    set grafsize($gID,canW) 600
    if { $prop(type_of_run) == "UHF" } {
	set N_graf 2
	set grafsize($gID,canH) 600
	set grafdata($gID,N_text) 2
	set grafdata($gID,Xtick_text_draw,1) 0
	set grafdata($gID,text1) "Alpha Band Structure"
	set grafdata($gID,text2) "Beta  Band Structure"
    } else {
	set N_graf 1
	set grafsize($gID,canH)   400
	set grafdata($gID,N_text) 0
    }
    
    set grafdata(Xtick_draw)        0
    set grafdata(Xtick_text_draw)   1
    set grafdata(Ytick_draw)        1
    set grafdata(Ytick_text_draw)   1
    set grafdata($gID,N_graf)       $N_graf
    set grafsize($gID,margin_X1)    $grafsize(margin_X1)
    set grafsize($gID,margin_X2)    $grafsize(margin_X2)
    set grafsize($gID,margin_Y1)    $grafsize(margin_Y1)
    set grafsize($gID,margin_Y2)    $grafsize(margin_Y2)

    for {set i 1} {$i <= $N_graf} {incr i} {
	set grafdata($gID,N_MXtick,$i)     [expr $NLINE + 1]
	set grafdata($gID,N_mXtick,$i)     0
	set grafdata($gID,Xtick1,$i)       0.0
	set grafdata($gID,Xtick1_text,$i)  $properties(TICK1)
	set grafdata($gID,Xmin,$i)         0.0
	set grafdata($gID,N_Xline,$i)      [expr $NLINE - 1]; #lines at K-nodes
    }

    #
    # this will be for projected band structure
    #
    for {set ig 1} {$ig <= $N_graf} {incr ig} {
	for {set ip 1} {$ip <= $N_proj} {incr ip} {
	    xcDebug "what1:: # of proj.:: $ip"
	    set grafdata($gID,N_point,$ig)      0
	    for {set i 1} {$i <= $NLINE} {incr i} {
		ReadFTN25 BAND $fileID $gID $ip $ig
		set grafdata($gID,X0,$ig) -$ftn25(X4); # used just here    
		for {set ni 1} {$ni <= $ftn25(NCOL)} {incr ni} {
		    set nxi  [expr $grafdata($gID,N_point,$ig) + $ni]
		    set nxi1 [expr $grafdata($gID,N_point,$ig) + $ni - 1]
		    set grafdata($gID,X$nxi,$ig) \
			    [expr $grafdata($gID,X$nxi1,$ig) + $ftn25(X4)]
		}
		set grafdata($gID,N_point,$ig) \
			[expr $grafdata($gID,N_point,$ig) + $ftn25(NCOL)]
		if { $i < $NLINE } {
		    set grafdata($gID,Xline$i,$ig) $grafdata($gID,X$nxi,$ig)
		    set grafsize($gID,Xline${i}_width,$ig)   2
		    set grafsize($gID,Xline${i}_fill,$ig)    "#000"
		    set grafsize($gID,Xline${i}_stipple,$ig) {}
		}
		set ii [expr $i + 1]	
		set grafdata($gID,Xtick$ii,$ig) $grafdata($gID,X$nxi,$ig)
		set grafdata($gID,Xtick${ii}_text,$ig) $properties(TICK$ii)
		
		if { $i == 1 && $ig == 1} {
		    set Ymin $ftn25(X1)
		    set Ymax $ftn25(X2)
		} elseif { $ftn25(X1) < $Ymin } {
		    set Ymin $ftn25(X1) 
		} elseif { $ftn25(X2) > $Ymax } {
		    set Ymax $ftn25(X2) 
		}
	    }
	}
    }

    for {set ig 1} {$ig <= $N_graf} {incr ig} {
	set yrange                       [expr $Ymax - $Ymin]
	set grafdata($gID,Yoffset,$ig)   0.05
	set grafdata($gID,Xoffset,$ig)   0.0

	set grafdata($gID,Ymin,$ig) \
		[expr $Ymin - $yrange * $grafdata($gID,Yoffset,$ig)]
	set grafdata($gID,Ymax,$ig) \
		[expr $Ymax + $yrange * $grafdata($gID,Yoffset,$ig)]

	set grafdata($gID,text${ig}_X) 0.02
	set grafdata($gID,text${ig}_Y) \
		[expr double($ig - 0.95) / double($N_graf)]
	
	set grafdata($gID,N_segment,$ig) [expr $ftn25(NROW) * $N_proj]
	set grafdata($gID,X_title,$ig)   ""
	set grafdata($gID,Y_title,$ig)   "E / a.u."
	set grafdata($gID,Xmax,$ig)      $grafdata($gID,X$nxi,$ig)
	set grafdata($gID,Yline1,$ig)    $ftn25(X5); # Fermi Energy
	set grafdata($gID,N_Yline,$ig)   1;                
	set grafdata($gID,Yline1_text,$ig)    "E_Fermi"
	set grafsize($gID,Yline1_stipple,$ig) @$system(BMPDIR)/dotH.bmp
	set grafsize($gID,Yline1_fill,$ig)    "#ff0000"
	set grafsize($gID,Yline1_width,$ig)   1
	set grafdata($gID,N_MYtick,$ig) 6
	set grafdata($gID,N_mYtick,$ig) 2
	set graf($gID,RorigX,$ig) 0.0
	set graf($gID,RorigY,$ig) \
		[expr double($N_graf - $ig) / double($N_graf)]
	set graf($gID,RsizeX,$ig) 1.0
	set graf($gID,RsizeY,$ig) \
		[expr 1.0 / double($N_graf) - [XGraph2nMRel $gID 5]]
	# Y tics
	set dy [expr abs(($grafdata($gID,Ymax,$ig) - \
		$grafdata($gID,Ymin,$ig)) / \
		($grafdata($gID,N_MYtick,$ig) - 1))]
	for {set i 1} {$i <= $grafdata($gID,N_MYtick,$ig)} {incr i} {
	    set grafdata($gID,Ytick$i,$ig) \
		    [expr $grafdata($gID,Ymin,$ig) + ($i - 1) * $dy]
	    set grafdata($gID,Ytick${i}_text,$ig) \
		    [TickFormat $grafdata($gID,Ytick$i,$ig) $dy]
	    #puts stdout "TICKS: $grafdata($gID,Ytick$i,$ig)"
	}
    }

    close $fileID
}


