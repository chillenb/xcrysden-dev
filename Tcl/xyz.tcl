#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/xyz.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc xyzOpen {{file {}} {can .mesa} {update 0}} {
    global system

    if { $file == {} } {
	# get-file
	set file [tk_getOpenFile -defaultextension .xyz \
		      -filetypes { 
			  {{XYZ Files}                  {.xyz} }
			  {{All Files}                  {.*}   }
		      } -initialdir $system(PWD) \
		      -title "Open XYZ File"]
	if { $file == "" } {
	    return
	}
    } else {
	if { ![file exists $file] } {	    
	    ErrorDialog "File \"$file\" does not exist !!!"
	    return
	}
    }

    reloadRegister $file XYZ
    UpdateWMTitle $file
    
    # uncompress the file is necessary
    
    set file    [gunzipFile $file]
    set head    [file rootname [file tail $file]]
    set xsfFile $system(SCRDIR)/$head.xsf

    # open the XYZ(read) and XSF(write) file

    set xyzID [open $file r]
    set xsfID [open $xsfFile w]

    
    # --------------------------------------------------
    # first-pass read (estimate number of animsteps)
    # --------------------------------------------------

    gets $xyzID line
    set natoms [lindex $line 0]
    if { ![string is integer $natoms] } {
	ErrorDialog "error parsing XYZ file: \"$file\" !!!"
	return
    }
    set nline 0
    seek $xyzID 0 start 
    while {[gets $xyzID line] > -1} {	
	incr nline
    }

    # calculate number of animsteps

    set animsteps [expr {int($nline / ($natoms + 2))}]
    if { $animsteps > 1 } {
	puts $xsfID "ANIMSTEPS  $animsteps"
    }


    # --------------------------------------------------
    # second-pass read: write XSF file
    # --------------------------------------------------

    seek $xyzID 0 start 
    set  iframe 0
    while {[gets $xyzID line] > -1} {	
	incr iframe
	if { [llength $line] == 0 } {
	    # skip empty line
	    continue
	}
	set natoms [lindex $line 0]
	if { ![string is integer $natoms] } {
	    ErrorDialog "error parsing XYZ file: \"$file\" !!!"
	    return
	}

	# read comment line

	gets $xyzID lin

	# read atoms

	if { $animsteps > 1} {
	    puts $xsfID "ATOMS  $iframe"
	} else {
	    puts $xsfID "ATOMS"
	}
	for {set ia 0} {$ia < $natoms} {incr ia} {
	    if {[gets $xyzID line] < 0} {
		ErrorDialog "error parsing XYZ file: \"$file\" !!!"
		return
	    }	
	    puts $xsfID $line
	}
    }
    
    # close files
    close $xyzID
    close $xsfID

    # load the structure ...

    xsfOpenParsed $xsfFile $can $update
}
