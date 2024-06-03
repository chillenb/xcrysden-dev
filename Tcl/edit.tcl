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

proc defaultEditor {{file ""}} {
    global edit


    set t [xcToplevel [WidgetName] "File: $file" "Editor"]
    
    set edit(saveFile,$t) $file

    # scrolled text
    
    set edit(sw,$t) [ScrolledWindow $t.sw -auto both -scrollbar both]
    set edit(tw,$t) [text $edit(sw,$t).text -width 80 -heigh 25]
    #pack $edit(sw,$t) $edit(tw,$t) -side top -expand 1 -fill both
    pack $edit(sw,$t) -side top -expand 1 -fill both
    $edit(sw,$t) setwidget $edit(tw,$t)

    if { $file != "" } {
	if { [file readable $file] } {
	    set content [ReadFile $file]
	    $edit(tw,$t) insert 1.0 $content
	} else {
	    ErrorDialog "file \"$file\" does not exist or is not readable"	    
	}
    }

    set bb [frame $t.f]
    pack $bb -side top -fill x

    # buttons

    foreach c {
	cancel save saveAs close
    } text {
	Cancel Save "Save As" "Save & Close"
    } {
	set b [string tolower $c]
	button $bb.$b -text $text -command [list defaultEditor_$c $t]
	pack $bb.$b -side left -expand 1
    }

    tkwait visibility $edit(tw,$t)
    focus $edit(tw,$t)
    return $t
}

proc defaultEditor_cancel {t} {
    set button [tk_messageBox \
		    -message "Changes will be lost. Really close this window?" \
                    -type yesno -icon question]
    if { $button == "yes" } {
	defaultEditor_done $t
    }
}

proc defaultEditor_close {t} {
    defaultEditor_save $t
    defaultEditor_done $t

}

proc defaultEditor_done {t} {
    global edit

    array unset edit *,$t
    destroy $t
}

proc defaultEditor_save {t} {
    global edit
    if { $edit(saveFile,$t) != "" } {
	set content [$edit(tw,$t) get 1.0 end]
	WriteFile $edit(saveFile,$t) $content
    } else {
	defaultEditor_saveAs $t
    }
}

proc defaultEditor_saveAs {t} {
    global edit system

    set filetypes {
        {{Text Files}     {.txt .text}}
        {{Shell Scripts}  {.sh}}
	{{Input Files}    {.in .inp}}
        {{All Files}      * }
    }

    set saveFile [tk_getSaveFile \
		      -initialdir $system(PWD) \
		      -title      "Save File As" \
		      -defaultextension "" \
		      -filetypes $filetypes]

    # maybe Cancel button was pressed    
    if { $saveFile == "" } {     
        return
    }

    # now save the file
    set edit(saveFile,$t) $saveFile
    defaultEditor_save $t
    wm title $t "File: [file tail $saveFile]"
}

#set system(TOPDIR) [pwd]
#set system(PWD) [pwd]
#lappend auto_path .
#set BWidget_dir $system(TOPDIR)/external/lib/bwidget1.8.0
#lappend auto_path  $BWidget_dir
#package require BWidget
#defaultEditor
