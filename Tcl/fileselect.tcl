#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/tonefile.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc fileselectResources {} {
    # path is used to enter the file name
    option add *Fileselect*path.relief		sunken	startup
    option add *Fileselect*path.background	white	startup
    option add *Fileselect*path.foreground	black	startup

    # Text for the OK and Cancel buttons
    option add *Fileselect*ok*text		OK	startup
    option add *Fileselect*ok*underline		0	startup
    option add *Fileselect*cancel.text		Cancel	startup
    option add *Fileselect*cancel.underline 	0	startup
    # Size of the listbox
    option add *Fileselect*list.width		20	startup
    option add *Fileselect*list.height		15	startup
    # Size of the directory
    option add *Fileselect*menub.width		30	startup
}

# fileselect returns the selected pathname, or {}
proc fileselect {{why "File Selection"} {default {}} {mustExist 1} \
	{seldir {0}}} {
    global fileselect system

    ########################################
    # CD to $system(PWD)
    xcDebug "pwd = [pwd]"
    cd $system(PWD)
    ########################################

    # if seldir -> select a directory
    set fileselect(seldir) $seldir
    puts stdout "fileselect(seldir):: $fileselect(seldir)"
    
    set t [toplevel .fileselect -bd 4 -class Fileselect]
    xcPlace . .fileselect 50 100
    wm title $t "File Browser"
    wm iconname $t "Open Structure"
    wm transient $t
    focus $t

    # if toplevel is destroyed we must cd back to $system(SCRDIR)
    bind $t <Destroy> { cd $system(SCRDIR) }
    fileselectResources
    
    message $t.msg \
	    -aspect 1000 \
	    -text $why \
	    -relief groove	   
    $t.msg config -font [ModifyFontSize $t.msg 18]
    pack $t.msg -side top -expand 1 -padx 10 -pady 10
    
    # a custom optionMenu where UpDirs will be
    set tfrm [frame $t.tfrm]
    pack $tfrm -expand 1 -fill x -pady 10
    set dlab [label $tfrm.lbl -text "Directory:" -padx 0]
    set fileselect(optMenu) [menubutton $tfrm.menub -relief raised \
	    -text "Listing ..." -menu $tfrm.menub.m -indicatoron 1 -padx 0]
    pack $dlab -side left
    pack $fileselect(optMenu) -side left -fill x -expand 1 
    set fileselect(menu) [menu $tfrm.menub.m -tearoff 0]

    # Create an entry for the pathname
    # The value is kept in fileselect(path)
    frame $t.top
    if $seldir {
	label $t.top.l -text "Directory:" -padx 0
    } else {
	label $t.top.l -text "File:" -padx 0
    }
    set e [entry $t.top.path \
	    -textvariable fileselect(path)]
    pack $t.top -side top -fill x
    pack $t.top.l -side left
    pack $t.top.path -side right -fill x -expand true
    
    # Create two listboxes to hold the directory contents
    # one for files; one for subdirs
    set tm [frame $t.mid -relief groove -bd 2]
    pack $tm -side top -fill both -expand 1 -pady 15 -ipadx 3 -ipady 3
    set lb1 [ScrolledListbox2 $tm.lb1]
    if !$seldir { set lb2 [ScrolledListbox2 $tm.lb2] }
    # Create the OK and Cancel buttons
    # The OK button has a rim to indicate it is the default
    frame $t.buttons 
    frame $t.buttons.ok -bd 2 -relief sunken
    set ok [button $t.buttons.ok.b \
	    -command [list fileselectOK $seldir]]
    set can [button $t.buttons.cancel \
	    -command fileselectCancel]
    
    # Pack the list, scrollbar, and button box
    # in a horizontal stack below the upper widgets
    pack $t.buttons -side bottom -expand 1 -pady 5 
    pack $t.buttons.ok $t.buttons.cancel \
	    -side left -padx 10 -pady 5 -expand 1
    pack $t.buttons.ok.b -padx 4 -pady 4 -expand 1
    puts stdout "1."
    if $seldir {
	fileselectBindings $t $e $lb1 {} $ok $can 
    } else {
	fileselectBindings $t $e $lb1 $lb2 $ok $can 
    }
    puts stdout "2."
    # Initialize variables and list the directory
    if {[string length $default] == 0} {
	set fileselect(path) {}
	if $fileselect(seldir) { set fileselect(path) $system(PWD) }
	set diry $system(PWD)
    } else {
	set fileselect(path) [file tail $default]
	puts stdout $fileselect(path)
	set diry [file dirname $default]
    }
    set fileselect(dir) {}
    set fileselect(done) 0
    set fileselect(mustExist) $mustExist
    
    # Wait for the listbox to be visible so
    # we can provide feedback during the listing 
    if $seldir {
	tkwait visibility .fileselect.mid.lb1.list
    } else {
	tkwait visibility .fileselect.mid.lb2.list
    }
    fileselectList $diry

    catch { grab $t }
    
    tkwait variable fileselect(done)
    destroy $t
    return $fileselect(path)
}


proc fileselectBindings { t e lb1 lb2 ok can } {
    global fileselect
    # t - toplevel
    # e - name entry
    # lb - listbox
    # ok - OK button
    # can - Cancel button
    
    # Elimate the all binding tag because we
    # do our own focus management
    
    if $fileselect(seldir) {
	set binlist [list $e $lb1 $ok $can]
    } else {
	set binlist [list $e $lb1 $lb2 $ok $can]
    }
    foreach w $binlist {
	    bindtags $w [list $t [winfo class $w] $w]
	}
	# Dialog-global cancel binding
	bind $t <Control-c> fileselectCancel

	# Entry bindings
	bind $e <Return> fileselectOK
	bind $e <space> fileselectComplete

	# A single click, or <space>, puts the name in the entry
	# A double-click, or <Return>, selects the name
	bind $lb1 <space> "fileselectTake %W ; focus $e"
	bind $lb1 <Button-1> \
		"fileselectClick %W %y ; focus $e"
	bind $lb1 <Return> "fileselectTake %W ; fileselectOK"
	bind $lb1 <Double-Button-1> \
		"fileselectClick %W %y ; fileselectOK"
	if !$fileselect(seldir) {
	    bind $lb2 <space> "fileselectTake %W ; focus $e"
	    bind $lb2 <Button-1> \
		    "fileselectClick %W %y ; focus $e"
	    bind $lb2 <Return> "fileselectTake %W ; fileselectOK"
	    bind $lb2 <Double-Button-1> \
		"fileselectClick %W %y ; fileselectOK"
	}

	# Focus management.  	# <Return> or <space> selects the name.
	bind $e <Tab> "focus $lb1 ; $lb1 select set 0"
	if !$fileselect(seldir) { 
	    bind $lb1 <Tab> "focus $lb2; $lb2 select set 0"
	    bind $lb2 <Tab> "focus $e" 
	} else {
	    bind $lb1 <Tab> "focus $e"
	}
	

	# Button focus.  Extract the underlined letter
	# from the button label to use as the focus key.
	foreach but [list $ok $can] {
		set char [string tolower [string index  \
			[$but cget -text] [$but cget -underline]]]
		bind $t <Alt-$char> "focus $but ; break"
	}
	bind $ok <Tab> "focus $can"
	bind $can <Tab> "focus $ok"

	# Set up for type in
	focus $e
}


proc PreOptMenu {name} {
    global fileselect dir 

    set dir $name
    set fileselect(dir) $name    
    set path $fileselect(path)
    set fileselect(path) ""
    fileselectOK
    set fileselect(path) $path 
    if $fileselect(seldir) {
	set fileselect(path) $name
    }
    puts stdout "PreOptMenu:: $path"
    OptMenu $name
    puts stdout "OptMenu:> $name"
}


proc OptMenu {dirnames} {
    global dir fileselect

    set dirnames [split $dirnames /]
    set dirnames [lrange $dirnames 1 end]
    set dirnames [Print $dirnames]
    set sel [lindex $dirnames [expr [llength $dirnames] - 1]]
    destroy $fileselect(optMenu)
    menubutton $fileselect(optMenu) -text $sel -menu $fileselect(menu) \
	    -indicatoron 1 -relief raised
    pack $fileselect(optMenu) -side left -fill x -expand true
    menu $fileselect(menu) -tearoff 0

    foreach name $dirnames {
	$fileselect(menu) add command -label $name -command \
		[list PreOptMenu $name]
    }
}


proc Print {names} {
    set prev {}
    set lis {}
    # $names could be ""
    puts stdout "NAMES> $names"
    set lis "/ "
    foreach el $names {
	if { $el != "" && $el != "."} {
	    append lis "$prev/$el "
	    set prev "$prev/$el"
	}
    } 
    return [format "%s" $lis]
}


proc fileselectList { diry {files {}} } {
    global fileselect dir

    set dir $diry
    # Update the directory display
    destroy $fileselect(optMenu)
    puts stdout "DIR> $dir"
    set dirnames [lrange [split $dir /] 1 end]
    set dirnames [Print $dirnames]
    set sel [lindex $dirnames [expr [llength $dirnames] - 1]]
    menubutton $fileselect(optMenu) -text $sel -menu $fileselect(menu) \
	    -indicatoron 1 -relief raised
    pack $fileselect(optMenu) -side right -fill x -expand 1
    menu $fileselect(menu) -tearoff 0

    foreach name $dirnames {
	$fileselect(menu) add command -label $name -command \
		[list PreOptMenu $name]
    }

    .fileselect.mid.lb1.list delete 0 end
    if !$fileselect(seldir) { .fileselect.mid.lb2.list delete 0 end }
    
    set fileselect(dir) $dir
    puts stdout "fileselectList.1 $fileselect(dir)"
    if ![file isdirectory $dir] {
	.fileselect.mid.lb1.list insert 0 "Bad Directory"
	return
    }
    if !$fileselect(seldir) {
	.fileselect.mid.lb2.list insert 0 Listing...
    } else {
	.fileselect.mid.lb1.list insert 0 Listing...
    }
    update idletasks
    if !$fileselect(seldir) {
	.fileselect.mid.lb2.list delete 0
    } else {
	.fileselect.mid.lb1.list delete 0
    }
    if {[string length $files] == 0} {
	# List the directory and add an
	# entry for the parent directory	    
	if {[catch {glob -nocomplain $fileselect(dir)/*}]} {
	    tk_messageBox -message "Can not change to directory \"\
		    $fileselect(dir)\" \
		    Permission denied." -type ok -icon warning
	} else {
	    set files [glob -nocomplain $fileselect(dir)/*]
	    puts stdout "@ $fileselect(dir)"
	    if { $fileselect(dir) != "//." && $fileselect(dir) != "/." && \
		$fileselect(dir) != "/"} {
	    .fileselect.mid.lb1.list insert end ../
	    }
	}
    }
    # Sort the directories in lb1 & files in lb2
    set dirs {}
    set others {}
    foreach f [lsort $files] {
	if [file isdirectory $f] {
	    lappend dirs [file tail $f]/
	} else {
	    lappend others [file tail $f]
	}
    }
    foreach f $dirs {
	.fileselect.mid.lb1.list insert end $f
    }
    if !$fileselect(seldir) {
	foreach f $others { 
	    .fileselect.mid.lb2.list insert end $f
	}
    }
}


proc fileselectOK {{okbut {0}}} {
	global fileselect 

    # $fileselect(seldir) == 1 when opening directory

    if { $fileselect(path) == "/" } {
	if $fileselect(seldir) {
	    set fileselect(path) /
	} else {
	    set fileselect(path) {}
	}
	set fileselect(dir) $fileselect(path)
    }
    # Handle the parent directory specially & current dir (./)
    if {[regsub {^\.\./?} $fileselect(path) {} newpath] != 0} {	
	set fileselect(path) $newpath	
	set fileselect(dir) [file dirname $fileselect(dir)]
	puts stdout "fileselect(dir).1> $fileselect(dir)"
	fileselectOK
	return
    }
    if { $fileselect(seldir) && \
	    [regsub {/\.\.} $fileselect(path) {} newpath] != 0} {
	
	set fileselect(path) [file dirname $newpath]
	set fileselect(dir)  [file dirname $fileselect(dir)]
    }
    # this is to prevent the //./
    if {[regsub {^//\./?} $fileselect(dir) / newdir] != 0} {
	set fileselect(dir) $newdir
	puts stdout "fileselect(dir).2> $fileselect(dir)"
	fileselectOK
	return
    }
   
    # if fileselect(path) != '/' then
    if { $fileselect(path) != "/" } {
	set path [string trimright $fileselect(dir)/$fileselect(path) /]
    } else {
	set path ""
    }
    
    if { [file isdirectory $path] && $fileselect(seldir) == 0 } {
	set fileselect(path) {}
	puts stdout "back to fileselectList"
	fileselectList $path
	return
    }
    if { [file isdirectory $fileselect(path)] && $fileselect(seldir) == 1 && \
	    $okbut == 0 } {
	fileselectList $fileselect(path)
	return
    }
    if { $fileselect(seldir) == 1 && [file isdirectory $fileselect(path)] && \
	    $okbut == 1 } {
	set fileselect(path) $fileselect(path)
	set fileselect(done) 1
	return
    } elseif { [file exists $path] && $fileselect(seldir) == 0 } {
	set fileselect(path) $path
	set fileselect(done) 1
	return
    }

    # Neither a file or a directory.
    # See if glob will find something
    if [catch {glob $path} files] {
	# No, perhaps the user typed a new
	# absolute pathname
	if [catch {glob $fileselect(path)} path] {
	    puts stdout "#"
	    tk_messageBox -message "File/Directory \"$fileselect(path)\" \
		    does not exist !" -icon error -type ok
	    raise .fileselect .
	    set fileselect(dir) [file dirname $fileselect(dir)]
	    set fileselect(path) {}
	    fileselectOK
	    return
	} 
	#else
	# OK - try again
	# set fileselect(dir) [file dirname $fileselect(dir)]
	# set fileselect(path) [file tail $fileselect(path)]
	# puts stdout "here"
	# fileselectOK
	#return
    } else {
	# Ok - current directory is ok,
	# either select the file or list them.
	if {[llength [split $files]] == 1} {
	    if { $fileselect(seldir) && $fileselect(dir) == "/" } {
		set fileselect(path) "/"
	    } else {
		set fileselect(path) $files
	    }
	    puts stdout "go back to fileselectOK"
	    fileselectOK
	} else {
	    set fileselect(dir) [file dirname [lindex $files 0]]
	    puts stdout "OK.1> $files"
	    puts stdout "OK.2> $fileselect(dir)"
	    #append the directorys to files; also ../
	    if { $fileselect(dir) != "//." && $fileselect(dir) != "/"} {
		append files " ../"
	    }
	    set allfiles [glob -nocomplain $fileselect(dir)/*]
	    foreach f $allfiles {
		#set f [file tail $f]
		if {[file isdirectory $f]} {
		    # if dir is not yet on the list-->append
		    if { [lsearch -exact $files $f] == -1} {
			append files " $f"
		    }
		}
	    }
	    fileselectList $fileselect(dir) $files
	}
    }
}


proc fileselectCancel {} {
	global fileselect
	set fileselect(done) 1
	set fileselect(path) {}
}


proc fileselectClick { lb y } {
    # Take the item the user clicked on
    global fileselect
    set fileselect(path) [$lb get [$lb nearest $y]]
    if $fileselect(seldir) {
	set fileselect(path) [string trimright $fileselect(path) /]
	set fileselect(path) $fileselect(dir)/$fileselect(path)
    }
    
}


proc fileselectTake { lb } {
	# Take the currently selected list item
	global fileselect
	set fileselect(path) [$lb get [$lb curselection]]
}


proc fileselectComplete {} {
    global fileselect

    # Do file name completion
    # Nuke the space that triggered this call
    set fileselect(path) [string trim $fileselect(path) \t\ ]
    
    # Figure out what directory we are looking at
    # dir is the directory
    # tail is the partial name
    if {[string match /* $fileselect(path)]} {
	set dir [file dirname $fileselect(path)]
	set tail [file tail $fileselect(path)]
    } elseif [string match ~* $fileselect(path)] {
	if [catch {file dirname $fileselect(path)} dir] {
	    return	;# Bad user
	}
	set tail [file tail $fileselect(path)]
    } else {
	set path $fileselect(dir)/$fileselect(path)
	set dir [file dirname $path]
	set tail [file tail $path]
    }
    # See what files are there
    set files [glob -nocomplain $dir/$tail*]	
    if {[llength [split $files]] == 1} {
	# Matched a single file
	set fileselect(dir) $dir
	set fileselect(path) [file tail $files]
    } else {
	if {[llength [split $files]] > 1} {
	    # Find the longest common prefix
	    set l [expr [string length $tail]-1]
	    set miss 0
	    # Remember that files has absolute paths
	    set file1 [file tail [lindex $files 0]]
	    while {!$miss} {
		incr l
		if {$l == [string length $file1]} {
		    # file1 is a prefix of all others
		    break
		}
		set new [string range $file1 0 $l]
		foreach f $files {
		    if ![string match $new* [file tail $f]] {
			set miss 1
			incr l -1
			break
		    }
		}
	    }
	    set fileselect(path) [string range $file1 0 $l]
	}
	#append the directorys to files; also ../
	puts stdout ">>> $fileselect(dir)"
	if { $fileselect(dir) != "//." && $fileselect(dir) != "/" } {
	    append files " ../"
	}
	set allfiles [glob -nocomplain $dir/*]
	foreach f $allfiles {
	    #set f [file tail $f]
	    if {[file isdirectory $f]} {
		# if dir is not yet on the list-->append
		if { [lsearch -exact $files $f] == -1} {
		    append files " $f"
		}
	    }
	}
	puts "FILES= $files"
	fileselectList $dir $files
    }
}
#lappend auto_path /net/surf/users/tone/prog/XCrys/TCL
#set pwd [pwd]
#lappend auto_path /home/tone/prog/XCrys/Mesa
#puts stdout "END:: [fileselect {Open Directory} {} 1 1]" 





