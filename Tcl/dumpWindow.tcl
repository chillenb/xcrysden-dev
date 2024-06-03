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
# Copyright (c) 2008 by Anton Kokalj                                        #
#############################################################################


proc dumpWindow {window {outfile ""}} {
    raise $window
    dumpWindowByID [winfo id $window] $outfile
}

proc dumpWindowByID {id {outfile ""}} {
    global xcMisc system printSetup
    
    SetWatchCursor
    
    update
    update idletask
    
    if { $outfile == "" } {
	set outfile [dumpWindow_queryFilename]
    }

    if { [info exists xcMisc(ImageMagick.import)] } {

        if { [info exists printSetup(useOptions) ] } {
            set tmpFile [file join $system(SCRDIR) [file tail $outfile]]
        } else {
            set tmpFile $outfile
        }

	xcCatchExecReturn $xcMisc(ImageMagick.import) -window $id $tmpFile

        if { [info exists printSetup(useOptions) ] } {
            scripting::_printToFile_imageConvert $tmpFile $outfile
        }
    } elseif { [info exists xcMisc(xwd)] } {
	
	xcCatchExecReturn $xcMisc(xwd) -id $id -out $system(SCRDIR)/tmp.xwd
	
	if { [info exists xcMisc(ImageMagick.convert)] } {
	    # convert by ImageMagick
	    scripting::_printToFile_imageConvert $system(SCRDIR)/tmp.xwd $outfile
	    
	} else {	    
	    # try to convert by netbmp
	    
	    set head [file rootname $outfile]
	    file copy -force $system(SCRDIR)/test.xwd $head.xwd
	    
	    set ext [string trim [file extension $outfile] .]
	    switch $ext {
		jpg { set ext jpeg }
		tif { set ext tiff }
		eps { set ext ps }
	    }
	    
	    set xwdtopnm [auto_execok xwdtopnm]
	    set pnmtoext [auto_execok pnmto$ext]
	    
	    if { $xwdtopnm != "" && $pnmtoext != "" } {
		catch {exec $xwdtopnm $head.xwd > $head.pnm}
		catch {exec $pnmtoext $head.pnm > $head.$ext}
	    }
	}
    } else {
	ErrorDialog "cannot make windowDump: neither \"import\" nor \"xwd\" programs are available"
    }
    
    ResetCursor
}


proc dumpWindow_queryFilename {} {
    global xcMisc system
    
    set head [file rootname [file tail $xcMisc(titlefile)]]
    
    set deffile $head.png
    set defext  .png
    set filetypes {
	{{PNG}        {.png} }
	{{JPEG}       {.jpg .jpeg} }
	{{GIF}        {.gif} }
	{{TIFF}       {.tif .tiff}}
	{{EPS}        {.eps .ps} }
	{{All Files}  *      }
    }
    
    set sfile [tk_getSaveFile -initialdir $system(PWD) \
		   -title             "Print to File" \
		   -defaultextension  $defext \
		   -initialfile       $deffile \
		   -filetypes         $filetypes]
    
    if { $sfile == {} } {
	return -code return
    }
    return $sfile
}
