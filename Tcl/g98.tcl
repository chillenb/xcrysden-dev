#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/g98.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc g98Cube {{file {}} {viewmol_exists {}}} {
    global g98 system

    # delete welcome window
    set title_off 0
    if {[winfo exists .title]} { 
	destroy .title 
	set title_off 1
    }

    if { $file == {} } {
	# get-file
	set file [tk_getOpenFile -defaultextension .cube \
		      -filetypes { 
			  {{All Files}                  {.*}  }
		      } -initialdir $system(PWD) \
		      -title "Open Gaussian CUBE File"]
	if { $file == "" } {
	    return
	}
    } else {
	if { ![file exists $file] } {	    
	    ErrorDialog "WARNING: File \"$file\" does not exist !!!"
	    return
	}
    }

    # now create a Viewer ...

    if { $viewmol_exists == "" } {
	if { $title_off } {		
	    wm deiconify .
	}
	ViewMol .
    }

    # parse the cube file ...
    
    set files [lsort -increasing [g98Cube:cube2xsf $file]]

    # ------------------------------------------------------------------------
    # Now we have either one or several XSF files. If there is only one,
    # simply load it, otherwise ask for which one to load !!!
    # ------------------------------------------------------------------------

    set nfiles [llength $files]
    
    # ... and load the structure 

    if { $nfiles == 0 } {
	ErrorDialog "An ERROR occured: cube2xsf has failed !!!"
	return
    } elseif { $nfiles == 1 } {
	set file $files
    } else {
	# make a selection-list, where user will select a directory !!!
	set file [g98Cube:selectFile $files]
    }
	
    xsfOpen $g98(cube_dir)/$file .mesa
}


proc g98Cube:cube2xsf {file} {
    global g98 system

    SetWatchCursor    

    # uncompress the file is necessary

    set file [gunzipFile $file]


    # make a new directory where converting files will be placed

    set dir $system(SCRDIR)/g98cube
    file mkdir $dir
    cd $dir
    set g98(cube_dir) [pwd]

    # execute cube2xsf
    set cw [DisplayUpdateWidget "Calculating" "Converting GAUSSIAN CUBE file to XSF"]
    update
    xcCatchExec $system(BINDIR)/cube2xsf 1 $file
    destroy $cw

    ResetCursor

    return [glob -nocomplain *.xsf]
}


proc g98Cube:selectFile {files} {
    global g98cube

    set result ""

    set t [xcToplevel .g98cube "Select File" "Select File"]

    set f1 [frame $t.1 -class StressText]
    set f2 [frame $t.2]
    set f3 [frame $t.3]    
    pack $f1 $f2 $f3 -side top -padx 10 -pady 10 -ipadx 10 -ipady 10


    set m [message $f1.m -text "Your CUBE file contain several molecular orbitals, therefore several XSF files were created.\n\nPlease select a particular XSF file that you would like to analyze." -width 400 -justify left]
    pack $m -side top -fill both -expand 1 -padx 3 -pady 3

    set ind 1
    set g98cube(xsf_file) [lindex $files 0]
    
    set irow 1
    set icol 1
    foreach file $files {
	set r [radiobutton $f2.r$ind \
		   -variable g98cube(xsf_file) \
		   -text     "File: $file" \
		   -value    $file \
		   -relief   ridge \
		   -bd       2 \
		   -anchor   w]
	grid $r -row $irow -column $icol -padx 3 -pady 3 -sticky we
	incr icol
	if { $icol > 3 } {
	    set icol 1
	    incr irow
	}
	incr ind
    }
    
    set b [button $f3.b -text "Continue" -command [list CancelProc $t]]
    pack $b -side top -ipadx 5 -ipady 5

    tkwait window $t
    return $g98cube(xsf_file)
}
