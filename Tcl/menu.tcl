#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/menu.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2014 by Anton Kokalj                                   #
#############################################################################


proc ViewMolMenu {w can} {
    global radio check

    # load constant values for xc_newvalue, xc_resetvar, xc_getvalue 
    # and xc_getdefault

    # $w --- widget to place MENU in
    # $can - widget where STRUCTURE is displayed

    # colors image
    #image create photo colors -format gif -file $system(BMPDIR)/jpaint.gif

    menubutton $w.vmcolor -image colors       -menu $w.vmcolor.menu 
    menubutton $w.vmfile  -text  "File"       -menu $w.vmfile.menu -underline 0
    menubutton $w.vmdis   -text  "Display"    -menu $w.vmdis.menu -underline 0
    menubutton $w.vmmod   -text  "Modify"     -menu $w.vmmod.menu -underline 0
    menubutton $w.vmadvg  -text  "AdvGeom"    -menu $w.vmadvg.menu -underline 0 
    menubutton $w.vmpro   -text  "Properties" -menu $w.vmpro.menu -underline 0 
    menubutton $w.vmdat   -text  "Tools"      -menu $w.vmdat.menu -underline 0
    menubutton $w.vmhelp  -text  "Help"       -menu $w.vmhelp.menu -underline 0

    pack $w.vmcolor -side left -padx 3
    pack $w.vmfile  $w.vmdis $w.vmmod $w.vmadvg $w.vmpro $w.vmdat $w.vmhelp \
	    -side left -padx 10 -pady 3 
    #################################################################    

    set check(crds)       0
    set check(labels)     0
    set check(frames)     1
    set check(unibond)    0
    set check(perpective) 0

    set radio(cellmode) "conv"
    set radio(ball)     "Balls based on covalent radii"
    set radio(space)    "SpaceFill based on covalent radii"
    set radio(frames)   rods
    set radio(hexamode) "parapipedal"
    set radio(unitrep)  "cell"


    ########################################################################
    set mcolor [menu $w.vmcolor.menu]
    set mfile  [menu $w.vmfile.menu]
    set mmod   [menu $w.vmmod.menu]
    set mdis   [menu $w.vmdis.menu]
    set madvg  [menu $w.vmadvg.menu]
    set mpro   [menu $w.vmpro.menu]
    set mdat   [menu $w.vmdat.menu]
    set mhelp  [menu $w.vmhelp.menu]

    mainMenu $can $mcolor $mfile $mmod $mdis $madvg $mpro $mdat $mhelp

    ### accelerators

    bind . p {ToggleMenuCheckbutton perspective Perspective}
    bind . d {ToggleMenuCheckbutton depthcuing  DepthCuing}
    bind . a {ToggleMenuCheckbutton antialias   AntiAlias}
    bind . y {ToggleMenuCheckbutton crds        CrdSist}
    bind . s {ToggleMenuCheckbutton labels      AtomLabels}
    bind . c {ToggleMenuCheckbutton frames      CrysFrames}
    bind . u {ToggleMenuCheckbutton unibond     Unibond}
    bind . h {ToggleMenuCheckbutton Hbonds      Hbonds .mesa}
    bind . f {ToggleMenuCheckbutton forces      forceVectors .mesa}
    bind . w {ToggleMenuCheckbutton wigner      WignerSeitz}
    bind . m {ToggleMenuCheckbutton pseudoDens  PseudoDensity}

    bind . <F1> {DisplayMode3D}
    bind . <F2> {DisplayMode2D}
    bind . <F3> {ToggleMenuRadiobutton cellmode prim CellMode 1}
    bind . <F4> {ToggleMenuRadiobutton cellmode conv CellMode 1}
}


proc mainMenu {can mcolor mfile mmod mdis madvg mpro mdat mhelp} {
    global mody ball radio check species mesa_bg xcColors system \
	    undoMenu undoAdvGeom light

    #------------------------------------------------------------------------
    # COLOR MENU (i.e. Palette)
    #------------------------------------------------------------------------

    ColorMenu $can $mcolor


    #------------------------------------------------------------------------
    # FILE MENU
    #------------------------------------------------------------------------

    $mfile add command -label "New CRYSTAL Input" -command new_file 
    $mfile add separator
    #$mfile add command -label "Open Structure" \
    #	     -command [list OpenXSFStruct $can]
    $mfile add cascade -label "Open Structure ..."        -menu $mfile.opstr
    $mfile add cascade -label "Open PWscf ..."            -menu $mfile.pwscf 
    $mfile add cascade -label "Open CRYSTAL ..." -menu $mfile.crystal
    $mfile add cascade -label "Open WIEN2k ..."           -menu $mfile.wien
    #$mfile add command -label "Band Path Selection - Test" \
    #	      -command Bz_MakeToplevel
    $mfile add separator
    $mfile add command -label "Close" -command CloseCase -accelerator "Crlt-w"
    bind . <Control-w> CloseCase

    $mfile add separator

    $mfile add command -label "Save XSF Structure" -command SaveXSFStruct -accelerator "Ctrl-s"
    $mfile add command -label "Save Current State and Structure" -command saveState -accelerator "Ctrl-t"
    $mfile add command -label "Save Current State" -command {saveState state} -accelerator "Ctrl-u"
    bind . <Control-s> SaveXSFStruct
    bind . <Control-t> saveState
    bind . <Control-Alt-s> saveState
    bind . <Control-u> { saveState state }

    $mfile add command -label "Save CRYSTAL Input" \
	    -command SaveCrystalInput
    $mfile add command -label "Save WIEN2k Struct File" -command wnSaveSFile
    $mfile add separator
    $mfile add command -label "Print to File" -command {printTogl .mesa} -accelerator "Ctlr-p"
    bind . <Control-p> {printTogl .mesa}
    #$mfile add command -label "Print " -command PrintStructure
    $mfile add command -label "Print Setup" -command printSetup -accelerator "Ctlr-Alt-p"
    bind . <Control-Alt-p> printSetup
    

    $mfile add separator
    $mfile add cascade -label "Utilities ..." -menu $mfile.util
    $mfile add cascade -label "XCrySDen Examples ..." -menu $mfile.examples

    $mfile add separator
    $mfile add command -label "Exit" -command exit_pr -underline 1 -accelerator "Ctrl-q"
    trace variable species w xcTrace

    # Open Structure ...

    set opstr [menu $mfile.opstr -tearoff 0]
    $opstr add command -label "Open XSF (XCrySDen Structure File)" \
	-command [list xsfOpenMenu $can]
    $opstr add command -label "Open AXSF (Animation XCrySDen Structure File)" \
	-command [list xsfAnimOpenMenu $can]
    $opstr add command -label "Open BXSF (i.e. Fermi Surface Files)" \
	-command [list bxsfOpenMenu $can]
    $opstr add separator
    $opstr add command -label "Open XCrySDen Scripting File" \
	-command [list scriptOpenMenu $can]    
    $opstr add separator
    $opstr add command -label "Open XYZ" -command xyzOpen
    $opstr add command -label "Open PDB" \
	-command [list OpenXYZPDB $can pdb]
    $opstr add separator
    $opstr add command -label "Open Gaussian Z-Matrix File" \
	-command [list gzmat_menu {Gaussian Z-Matrix File}]
    $opstr add command -label "Open Gaussian Output File" \
	-command [list addOption:hardcoded \
		      [list sh $system(TOPDIR)/scripts/g98toxsf.sh] {} \
		      "Open Gaussian Output File"]
    $opstr add command -label "Open Gaussian Cube File" \
	-command [list g98Cube {} viewmol_exists]

    $opstr add separator

    $opstr add command -label "Open Orca Output File" \
	-command [list addOption:hardcoded \
		      [list sh $system(TOPDIR)/scripts/orca2xsf.sh] {} \
		      "Open Orca Output File"]

    $opstr add separator
    
    $opstr add command -label "Open PWscf Input File" \
	-command [list openExtStruct 3 crystal external \
		      [list sh $system(TOPDIR)/scripts/pwi2xsf.sh] \
		      pwi2xsf.xsf_out \
		      {PWSCF Input} ANGS \
		      -preset pwInputPreset]
    $opstr add command -label "Open PWscf Output File" \
	-command [list openExtStruct 3 crystal external \
		      [list sh $system(TOPDIR)/scripts/dummy.sh] \
		      pwo2xsf.xsf \
		      {PWSCF Output} ANGS \
		      -preset pwOutputPreset]

    $opstr add separator

    $opstr add command -label "Open FHI98MD \"inp.ini\" File" \
	    -command [list openExtStruct 3 crystal external34 \
	    $system(BINDIR)/fhi_inpini2ftn34 \
	    $system(ftn_name).34 {FHI98MD "ini.inp"} BOHR \
	    -preset [list fhiPreset inpini]]
    $opstr add command -label "Open FHI98MD \"coord.out\" File" \
	    -command [list openExtStruct 3 crystal external \
	    $system(BINDIR)/fhi_coord2xcr \
	    fhi_coord.xcr {FHI98MD "coord.out"} BOHR \
	    -preset [list fhiPreset coord]]

    # Open PWscf ...
    set pwscf [menu $mfile.pwscf -tearoff 0]
    $pwscf add command -label "Open PWscf Input File" \
	-command [list openExtStruct 3 crystal external \
		      [list sh $system(TOPDIR)/scripts/pwi2xsf.sh] \
		      pwi2xsf.xsf_out \
		      {PWSCF Input} ANGS \
		      -preset pwInputPreset]
    $pwscf add command -label "Open PWscf Output File" \
	-command [list openExtStruct 3 crystal external \
		      [list sh $system(TOPDIR)/scripts/dummy.sh] \
		      pwo2xsf.xsf \
		      {PWSCF Output} ANGS \
		      -preset pwOutputPreset]

    # Open CRYSTAL-95/98/03/06 ...
    set crystal [menu $mfile.crystal -tearoff 0]
    $crystal add command -label "Open CRYSTAL Input"      -command OpenCyrstalInput
    $crystal add command -label "Open CRYSTAL Properties (e.g. fort.9)" -command PropC95

    # Open WIEN2k

    set wien [menu $mfile.wien -tearoff 0]
    $wien add command -label "Open WIEN2k Struct File" \
	    -command wnOpenSFile
    $wien add command -label "Render pre-Calculated Density" \
	    -command wnOpenRenderDensity
    $wien add command -label "Calculate & Render Density" \
	    -command wnOpenCalcAndRenderDensity
    $wien add command -label "Select k-path" \
	    -command wnOpenKPath
    $wien add command -label "Fermi Surface" \
	    -command wnOpenFS

    set util [menu $mfile.util -tearoff 0]
    $util add command -label "Periodic Table of Elements" \
	    -command [list ptable .]

    #
    # XCRYSDEN examples !!!
    #
    set exam [menu $mfile.examples -tearoff 0]
    $exam add command -label "XSF Files" \
	-command [list xsfOpenMenu $can $system(TOPDIR)/examples/XSF_Files]
    $exam add command -label "BXSF (i.e. Fermi Surface) Files" \
	-command {
	    set file [tk_getOpenFile -defaultextension .bxsf \
			  -filetypes { 
			      {{All Files}       {.*}   }
			      {{Band XSF Files}  {.bxsf}}
			  } -initialdir $system(TOPDIR)/examples/FermiSurface \
			  -title "Open Fermi Surface (BXSF File)"]
	    if { $file != "" } {
		genFSInit $file
	    }
	}   
    #--------------------
    $exam add separator
    #--------------------

    # SCRIPTING Files ...
    $exam add command -label "Scripting Files" \
	-command {
	    set file [tk_getOpenFile -defaultextension .tcl \
			  -filetypes {
			      {{XCrySDen Scripting Files}  {.tcl} }
			      {{XCrySDen Scripting Files}  {.xcrysden} }
			      {{All Files}  {.*}   }
			  } -initialdir $system(TOPDIR)/examples/Scripting \
			  -title "Open Tcl Script File"]
	    if { $file != "" } {
		scripting::source $file
	    }
	}

    #--------------------
    $exam add separator
    #--------------------
    
    $exam add command -label "CRYSTAL Input Files" \
	-command {
	    OpenCyrstalInput $system(TOPDIR)/examples/CRYSTALxx_input_files 
	}
    $exam add command -label "WIEN2k Struct Files" \
	-command {
	    set file [tk_getOpenFile -defaultextension .struct \
			  -filetypes {
			      {{All Files}            {.*}   }
			      {{WIEN2k Struct Files}  {.stuct}}
			  } -initialdir $system(TOPDIR)/examples/WIEN_struct_files \
			  -title "Open WIEN2k Struct File"]
	    if { $file != "" } {
		reloadRegister $file WIEN2k
                wnOpenSFile $file
	    }
	}
    $exam add cascade -label "PWSCF Files ..." -menu $exam.pwscf
    $exam add cascade -label "FHI98MD Files ..." -menu $exam.fhi
    $exam add cascade -label "GAUSSIAN Files ..." -menu $exam.gaus

    # ORCA examples 
    $exam add command -label "Orca Output Files ..." \
	-command {
	    set file [tk_getOpenFile -defaultextension .out \
			  -filetypes { 
			      {{All Files}       {.*}  }
			      {{Output Files}    {.out}}
                          } \
                          -initialdir $system(TOPDIR)/examples/Orca_files \
			  -title "Open ORCA Output File"]
	    if { $file != "" } {
                reloadRegister $file addOption:reload
		addOption:hardcoded [list sh $system(TOPDIR)/scripts/orca2xsf.sh] \
		    $file "Open ORCA Output File" viewmol_exists
	    }
	}
    
    #--------------------
    $exam add separator
    #--------------------

    $exam add command -label "PDB Files" \
	-command [list OpenXYZPDB $can pdb $system(TOPDIR)/examples/PDB]
    $exam add command -label "XYZ Files" \
	-command {
	    set file [tk_getOpenFile -defaultextension .xyz \
			  -filetypes {
			      {{All Files}  {.*}   }
			      {{XYZ Files}  {.xyz} }
			  } -initialdir $system(TOPDIR)/examples/XYZ \
			  -title "Open XYZ File"]
	    if { $file != "" } {
		xyzOpen $file
	    }
	}
    
    # CASCADES
    #---------------------
    # FHI98MD examples ...
    set fhi [menu $exam.fhi -tearoff 0]
    $fhi add command -label "Open FHI98MD \"inp.ini\" File" \
	-command {
	    set file [tk_getOpenFile -defaultextension .ini \
			  -filetypes { 
			      {{All Files}             {.*}  }
			      {{FHI98MD inp.ini Files} {.ini}} } \
			  -initialdir $system(TOPDIR)/examples/FHI98MD_files \
			  -title "Open FHI98MD \"inp.ini\" File"]
	    if { $file != "" } {
                reloadRegister $file openExtStruct:reload                

		openExtStruct 3 crystal external34 \
		    $system(BINDIR)/fhi_inpini2ftn34 \
		    $system(ftn_name).34 {FHI98MD "ini.inp"} BOHR \
		    -preset [list fhiPreset inpini] -file $file
	    }
	}
    
    $fhi add command -label "Open FHI98MD \"coord.out\" File" \
	-command {
	    set file [tk_getOpenFile -defaultextension .out \
			  -filetypes { 
			      {{All Files}               {.*}  }
			      {{FHI98MD coord.out Files} {.out}} } \
			  -initialdir $system(TOPDIR)/examples/FHI98MD_files \
			  -title "Open FHI98MD \"coord.out\" File"]
	    if { $file != "" } {
                reloadRegister $file openExtStruct:reload
                
		openExtStruct 3 crystal external \
		    $system(BINDIR)/fhi_coord2xcr \
		    fhi_coord.xcr {FHI98MD "coord.out"} BOHR \
		    -preset [list fhiPreset coord] -file $file
	    }
	}

    # PWSCF examples ...
    set pwscf [menu $exam.pwscf -tearoff 0]
    $pwscf add command -label "Open PWSCF Input File" \
	-command {
	    set file [tk_getOpenFile -defaultextension .inp \
			  -filetypes { 
			      {{All Files}          {.*}  }
			      {{PWSCF Input Files}  {.inp}} } \
			  -initialdir $system(TOPDIR)/examples/PWSCF_files \
			  -title "Open PWSCF Input File"]
	    if { $file != "" } {
                reloadRegister $file openExtStruct:reload
                
		openExtStruct 3 crystal external \
		    [list sh $system(TOPDIR)/scripts/pwi2xsf.sh] \
		    pwi2xsf.xsf_out \
		    {PWSCF Input} ANGS \
		    -preset pwInputPreset -file $file
	    }
	}

    $pwscf add command -label "Open PWSCF Output File" \
	-command {
	    set file [tk_getOpenFile -defaultextension .out \
			  -filetypes { 
			      {{All Files}          {.*}  }
			      {{PWSCF Output Files} {.out}} } \
			  -initialdir $system(TOPDIR)/examples/PWSCF_files \
			  -title "Open PWSCF Output File"]
	    if { $file != "" } {
                reloadRegister $file openExtStruct:reload
                
		openExtStruct 3 crystal external \
		    [list sh $system(TOPDIR)/scripts/dummy.sh] \
		    pwo2xsf.xsf \
		    {PWSCF Output} ANGS \
		    -preset pwOutputPreset -file $file
	    }
	}
    
    # GAUSSIAN examples ...
    set gaus [menu $exam.gaus -tearoff 0]
    $gaus add command -label "Gaussian Z-Matrix File" \
	-command {
	    set file [tk_getOpenFile -defaultextension .com \
			  -filetypes { 
			      {{All Files}               {.*}  }
			      {{GAUSSIAN Z-Matrix Files} {.gzmat}}
			      {{GAUSSIAN Input Files}    {.com}}
			      {{GAUSSIAN Input Files}    {.inp}} 
			      {{GAUSSIAN Input Files}    {.g98}} } \
			  -initialdir $system(TOPDIR)/examples/GAUSSIAN_files \
			  -title "Open Gaussian Z-Matrix File"]
	    if { $file != "" } {
                reloadRegister $file GZMAT
		gzmat $file
	    }
	}

    $gaus add command -label "Gaussian Output File" \
	-command {
	    set file [tk_getOpenFile -defaultextension .log \
			  -filetypes { 
			      {{All Files}                {.*}  }
			      {{GAUSSIAN Output Files}    {.log}}
			      {{GAUSSIAN Output Files}    {.out}} 
			      {{GAUSSIAN Output Files}    {.g98_out}} } \
			  -initialdir $system(TOPDIR)/examples/GAUSSIAN_files \
			  -title "Open Gaussian Output File"]
	    if { $file != "" } {
                reloadRegister $file addOption:reload
		addOption:hardcoded [list sh $system(TOPDIR)/scripts/g98toxsf.sh] \
		    $file "Open Gaussian Output File" viewmol_exists
	    }
	}

    $gaus add command -label "Gaussian Cube File" \
	-command {
	    set file [tk_getOpenFile -defaultextension .cube \
			  -filetypes { 
			      {{All Files}                {.*}  }
			      {{GAUSSIAN Cube   Files}    {.cube}} } \
			  -initialdir $system(TOPDIR)/examples/GAUSSIAN_files \
			  -title "Open Gaussian Cube File"]
	    if { $file != "" } {
		g98Cube $file viewmol_exists
	    }
	}



    #------------------------------------------------------------------------
    # DISPLAY MENU
    #------------------------------------------------------------------------

    #$mdis add command -label "WireFrame" -command WireFrame
    #$mdis add command -label "Points&Lines" -command PointLine
    $mdis add radiobutton -label "Lighting On" -command {Lighting On} \
	-variable light -value On -accelerator "F1"
    $mdis add radiobutton -label "Lighting Off" -command {Lighting Off} \
	-variable light -value Off -accelerator "F2"
    #----------
    $mdis add separator
    #----------
    $mdis add checkbutton -label "Coordinate System" -variable check(crds) \
	-command CrdSist    -accelerator "y"
    $mdis add checkbutton -label "Atomic Symbols" -variable check(labels) \
	-command AtomLabels -accelerator "s"
    $mdis add checkbutton -label "Crystal Cells" -variable check(frames) \
	-command CrysFrames -accelerator "c"
    $mdis add checkbutton -label "Unicolor Bonds" -variable check(unibond) \
	-command Unibond    -accelerator "u"
    $mdis add checkbutton -label "H-Bonds" -variable check(Hbonds) \
	-command [list Hbonds $can] -accelerator "h"
    $mdis add checkbutton -label "Forces" -variable check(forces) \
	-command [list forceVectors $can] -accelerator "f"
    $mdis add checkbutton -label "Wigner-Seitz Cells" -variable check(wigner) \
	-command WignerSeitz  -accelerator "w"
    $mdis add checkbutton -label "Molecular Surface" \
	-variable check(pseudoDens) \
	-command PseudoDensity -accelerator "m"
    #----------
    $mdis add separator
    #----------
    $mdis add checkbutton -label "Perspective Projection" \
	-variable check(perspective) \
	-command Perspective \
	-accelerator "p"
    $mdis add checkbutton -label "Depth Cuing" \
	-variable check(depthcuing) \
	-command DepthCuing \
	-accelerator "d"
    $mdis add checkbutton -label "Anti-Aliasing" \
	-variable check(antialias) \
	-command AntiAlias \
	-accelerator "a"
    #----------
    $mdis add separator
    #----------
    $mdis add cascade -label "Balls As ..." -menu $mdis.ball
    $mdis add cascade -label "SpaceFill As ..." -menu $mdis.space
    $mdis add cascade -label "Crystal Cells As ..." -menu $mdis.frames
    #----------
    $mdis add separator
    #----------
    $mdis add radiobutton -label "Primitive Cell Mode" \
	-command [list CellMode 1]\
	-variable radio(cellmode) \
	-value "prim" -accelerator "F3"
    $mdis add radiobutton -label "Conventional Cell Mode" \
	-command [list CellMode 1] \
	-variable radio(cellmode) \
	-value "conv" -accelerator "F4"
    $mdis add cascade -label "Hexagonal/Rhombohedral ..." -menu $mdis.1
    #----------
    $mdis add separator
    #----------
    $mdis add cascade -label "Unit of Repetition ..." -menu $mdis.unitr

    # BALLS MENU
    set mdisball [menu $mdis.ball -tearoff 0]
    $mdisball add radiobutton -label "Balls based on covalent radii" \
	    -command [list xc_newvalue .mesa $mody(L_BALL_COV)] \
	    -variable radio(ball)
    $mdisball add radiobutton -label "Balls based on Van der Waals radii" \
	    -command [list xc_newvalue .mesa $mody(L_BALL_VDW)] \
	    -variable radio(ball)
    # SPACEFILL MENU
    set mdisspace [menu $mdis.space -tearoff 0]
    $mdisspace add radiobutton -label "SpaceFill based on covalent radii" \
	    -command [list xc_newvalue .mesa $mody(L_SPACE_COV)] \
	    -variable radio(space)
    $mdisspace add radiobutton -label \
	    "SpaceFill based on Van der Waals radii" \
	    -command [list xc_newvalue .mesa $mody(L_SPACE_VDW)] \
	    -variable radio(space)
    # FRAMES/CELLS MENU
    set mdisframes [menu $mdis.frames -tearoff 0]
    $mdisframes add radiobutton \
	-label "Display crystal cells as lines when lighting" \
	-value lines \
	-variable radio(frames) \
	-command DispFramesAs 
    $mdisframes add radiobutton \
	-label "Display crystal cells as rods when lighting" \
	-value rods \
	-variable radio(frames) \
	-command DispFramesAs 
    #HEXAGONAL/RHOMBOHEDRAL MENU
    set mdis1 [menu $mdis.1 -tearoff 0]  
    $mdis1 add radiobutton \
	    -label "Parapipedal-shaped cell" \
	    -command [list CellMode 1]\
	    -variable radio(hexamode) \
	    -value "parapipedal"
    $mdis1 add radiobutton \
	    -label "Hexagonal-shaped cell" \
	    -command [list CellMode 1]\
	    -variable radio(hexamode) \
	    -value "hexagonal"
    #UNIT OF REPETITION
    set mdis3 [menu $mdis.unitr -tearoff 0]
    $mdis3 add radiobutton -label "Unit cell" \
	    -command [list CellMode 1]\
	    -variable radio(unitrep) \
            -value "cell"
    $mdis3 add radiobutton -label "Translational asymmetric unit" \
	    -command [list CellMode 1]\
	    -variable radio(unitrep) \
	    -value "asym"

    ###########################################################################
    # TEMPORARY
    #$mdis add command -label {Print Dimenisonality} \
    #	    -command {
    #	global periodic; 
    #	xcDebug -stderr "\nPERIODICITY:: $periodic(dim)"
    #}
    

    #------------------------------------------------------------------------
    # MODIFY MENU
    #------------------------------------------------------------------------

    #--
    # disable "Atomic Symbols/Fonts" for now, because they don't work
    $mmod add command -label "Atomic Symbols/Fonts" -command ModAtomLabels -accelerator "Shift-s"
    bind . <S> ModAtomLabels
    #$mmod entryconfigure "Atomic Symbols/Fonts" -state disabled
    #--
    
    $mmod add separator
    $mmod add command -label "Atomic Color" -command ModAtomCol -accelerator "Shift-a"
    bind . <A> ModAtomCol
    $mmod add command -label "Unibond Color" -command ModUnibondCol -accelerator "Shift-u"
    bind . <U> ModUnibondCol
    $mmod add command -label "Crystal Cell Color" -command ModCellCol -accelerator "Shift-c"
    bind . <C> ModCellCol
    $mmod add command -label "Coordinate System Color" -command ModXYZCol -accelerator "Shift-z"
    bind . <Z> ModXYZCol

    $mmod add separator

    $mmod add command -label "Atomic Radius" -command ModAtomRad -accelerator "Shift-r"
    bind . <R> ModAtomRad

    $mmod add command -label "Ball Factor" \
	    -command [list ModFactor Ball L_BALLF D_BALLF "Ball Factor"] -accelerator "Shift-b"
    bind . <B> [list ModFactor Ball L_BALLF D_BALLF "Ball Factor"]

    $mmod add command -label "Ball/Stick Ratio" \
	    -command [list ModFactor Stick L_RODF D_RODF "Ball/Stick Ratio"] -accelerator "Shift-q"
    bind . <Q> [list ModFactor Stick L_RODF D_RODF "Ball/Stick Ratio"]

    $mmod add cascade -label "LineWidth ..." -menu $mmod.1
    $mmod add command -label "Point Radius" \
	    -command [list ModFactor Point L_PLRADIUS D_PLRADIUS \
	    "Point Radius"] -accelerator "Shift-d"
    bind . <D> [list ModFactor Point L_PLRADIUS D_PLRADIUS "Point Radius"]

    $mmod add command -label "Crystal Cell's Rod Factor" \
	    -command [list ModFactor FrameRod L_FRAMERODF D_FRAMERODF \
	    "Crystal Cell's Rod Factor"] -accelerator "Shift-o"
    bind . <O> [list ModFactor FrameRod L_FRAMERODF D_FRAMERODF "Crystal Cell's Rod Factor"]


    $mmod add separator

    $mmod add command -label "Tessellation Factor" -command \
	[list ModFactor Tessellation L_TESSELLATION D_TESSELLATION "Tessellation Factor"] -accelerator "Shift-t"
    bind . <T> [list ModFactor Tessellation L_TESSELLATION D_TESSELLATION "Tessellation Factor"]

    $mmod add command -label "Lighting Parameters" -command glLight -accelerator "Shift-l"
    bind . <L> glLight

    $mmod add command -label "Perspective Parameters" -command ModPerspective -accelerator "Shift-p"
    bind . <P> ModPerspective

    #$mmod add command -label "Advanced OpenGL Parameters" -command glModParam -accelerator "Shift-x"
    $mmod add command -label "Material/Fog/Antialias Parameters" -command glModParam -accelerator "Shift-x"
    bind . <X> glModParam
    $mmod add separator

    $mmod add command -label "Force Settings" \
	    -command [list forceVectorsSet $can] -accelerator "Shift-f" 
    bind . <F> [list forceVectorsSet $can]

    $mmod add command -label "H-bonds Settings" \
	    -command [list HbondsSetting $can] -accelerator "Shift-h" 
    bind . <H> [list HbondsSetting $can]

    $mmod add command -label "Wigner-Seitz Cell Settings" \
	    -command SetWignerSeitzInit -accelerator "Shift-w"
    bind . <W> SetWignerSeitzInit

    $mmod add command -label "Molecular Surface Settings" \
	    -command PseudoDenConfigWid -accelerator "Shift-m"
    bind . <M> PseudoDenConfigWid

    $mmod add command -label "Animation Controls" \
	    -command [list xsfAnimWid $can] -accelerator "Shift-z"
    bind . <Z> [list xsfAnimWid $can]

    # screen factor is no more used !!!
    #$mmod add separator
    #$mmod add command -label "Screen Factor" \
    #	     -command [list ModFactor Screen L_SCRF D_SCRF]
    $mmod add separator
    $mmod add command -label "Number of Units Drawn" -command NumCellDrawn -accelerator "Shift-n"
    bind . <N> NumCellDrawn


    set mmod1 [menu $mmod.1 -tearoff 0]
    $mmod1 add command -label "WireFrame's line width" \
	    -command [list ModFactor Wire L_WFLINEWIDTH D_WFLINEWIDTH \
	    "WireFrame's line width"]
    $mmod1 add command -label "PointLines's line width" \
	    -command [list ModFactor Pointline L_PLLINEWIDTH D_PLLINEWIDTH \
	    "PointLine's line width"]
    $mmod1 add command -label "Crystal Cell's line width" \
	    -command [list ModFactor Frame L_FRAMELINEWIDTH D_FRAMELINEWIDTH \
	    "Crystal cell's line width"]
    $mmod1 add separator
    $mmod1 add command -label "Lighting-Off outline width" \
	    -command [list ModFactor Outline L_OUTLINEWIDTH D_OUTLINEWIDTH \
	    "Lighting-Off outline width"] 
    $mmod1 add command -label "Lighting-On Wire line width" \
	    -command [list ModFactor Wire3D L_WF3DLINEWIDTH D_WF3DLINEWIDTH \
	    "Lighting-On wire line width"] 

    # TEMPORARY
    #global tcl_platform
    #if { $tcl_platform(platform) == "windows" } {
    #	# disable the modification of font on windows as it does not
    #	# work ...
    #	
    #	$madvg entryconfig "Atomic Symbols/Fonts" -state disabled
    #}
    #/

    #------------------------------------------------------------------------
    # ADVGEOM MENU
    #------------------------------------------------------------------------

    $madvg add cascade -label "Undo ..." -menu $madvg.undo
    $madvg add cascade -label "Redo ..." -menu $madvg.redo
    $madvg add separator    
    $madvg add command -label "Rotate Cartesian Frame" -command RotFrame
    $madvg add cascade -label "Atoms & Cell Manipulation ..." -menu $madvg.1
    $madvg add command -label "Cut a SLAB" -command CutSlab
    $madvg add cascade -label "Cut a non-Periodical Structure ..." \
	    -menu $madvg.2
    $madvg add separator
    $madvg add command -label "Add an Option Manually" \
	-command cxxAdvGeom.manualOption
    $madvg add command -label "Edit Manually" \
	-command cxxAdvGeom.manualEdit
    $madvg add command -label "View Input Script" \
	-command cxxAdvGeom.viewScript
    $madvg add separator
    $madvg add cascade -label "Multi-Slab ..." -menu $madvg.3

    # UNDO/REDO
    set undoAdvGeom(menu) $madvg
    $madvg entryconfig "Undo ..." -state disabled
    $madvg entryconfig "Redo ..." -state disabled
    set undoMenu(undo)  [menu $madvg.undo -tearoff 0]
    set undoMenu(redo)  [menu $madvg.redo -tearoff 0]
    bind $undoMenu(undo) <Motion> [list UndoMenuMotion $undoMenu(undo) @%y]
    bind $undoMenu(redo) <Motion> [list UndoMenuMotion $undoMenu(redo) @%y]

    set madvg1 [menu $madvg.1 -tearoff 0]
    $madvg1 add command -label "Substitute Atom" -command AtomSubs
    $madvg1 add command -label "Remove Atom" -command AtomRemo
    $madvg1 add command -label "Insert Atom" -command AtomInse
    $madvg1 add command -label "Displace Atom" -command AtomDisp
    $madvg1 add separator
    $madvg1 add command -label "SuperCell" -command SuperCell
    $madvg1 add command -label "Elastic Deformation" -command Elastic

    set madvg2 [menu $madvg.2 -tearoff 0]
    $madvg2 add command -label "Cut a Cluster" -command CutCluster
    $madvg2 add command -label "Cut a Molecule" -command CutMol
    $madvg2 entryconfig "Cut a Molecule" -state disabled

    set madvg3 [menu $madvg.3 -tearoff 0]
    $madvg3 add command -label "Create a Multi-Slab" -command wnMultiSlab
    $madvg3 add command -label "Change Multi-Slab Vacuum Thickness" \
	    -command [list wnMultiSlab change]


    #------------------------------------------------------------------------
    # PROPERTIES MENU
    #------------------------------------------------------------------------

    $mpro add command -label "Get INFO" \
	    -command [list PropC95Cmd INFO]
    $mpro add command -label "Display Bandwidths" \
	    -command [list PropC95Cmd BWID]
    $mpro add command -label "Density of States" \
	    -command [list PropC95Cmd DOSS]
    $mpro add command -label "Band Structure" \
	    -command [list PropC95Cmd BAND]
    $mpro add separator
    $mpro add cascade -label "Isosurfaces ..."            -menu $mpro.1
    $mpro add cascade -label "Properties on Planes ..."   -menu $mpro.2
    #$mpro add cascade -label "Properties on Surfaces ..." -menu $mpro.3
    
    ###########################
    #### THIS IS TEMPORAL  ####
    ###########################
    #$mpro add separator
    #$mpro add command -label "Charge Density Maps" \
    #	     -command [list Iso3DInit ECHD "Charge Density Maps" ECHD \
    #	     "Charge Density Maps" 5]
    #$mpro add command -label "Charge Density Gradient" \
    #	     -command [list Iso3DInit ECHG "Charge Density Gradient" ECHG \
    #	     "Charge Density Gradient" 5]
    #$mpro add command -label "Charge Density Plane" \
    #	      -command [list Iso3DInit ECHP "Charge Density Plane" ECHP \
    #	      "Charge Density Plane" 5]    
    #$mpro add command -label "Electrostatic potential map" \
    #	     -command [list Iso3DInit POTM "Electrostatic potential map" POTM \
    #	     "Electrostatic potential map" 5]        
    ############################
    
    ##############################
    # ISOSURFACES
    set mpro1 [menu $mpro.1 -tearoff 0]
    $mpro1 add command -label "Charge Density" \
	    -command [list SetIsoGridSpace \
	    "Charge Density 3D Map - Grid Specification" "Grid" {ECHG\n0}]
    $mpro1 add command -label "Electrostatic Potential" \
	    -command [list SetIsoGridSpace \
	    "Electrostatic Potential 3D Map - Grid Specification" "Grid" POTM]
    $mpro1 add separator
    $mpro1 add command -label "Difference Maps" -command DiffIsoSurf_Widget
    ##############################
    # PROPERTIES ON PLANES
    set mpro2 [menu $mpro.2 -tearoff 0]
    $mpro2 add command -label "Charge Density" \
	    -command [list SetIsoGridSpace \
	    "Charge Density 2D Map - Grid Specification" "Grid" {ECHG\n0} 2D]
    $mpro2 add command -label "Electrostatic Potential" \
	    -command [list SetIsoGridSpace \
	    "Electrostatic Potential 2D Map - Grid Specification" "Grid" \
	    POTM 2D]
    $mpro2 add separator
    $mpro2 add command -label "Difference Maps" \
	    -command [list DiffIsoSurf_Widget 2D]

    ##############################
    # PROPERTIES ON SURFACES
    set mpro3 [menu $mpro.3 -tearoff 0]


    #------------------------------------------------------------------------
    # TOOLS MENU
    #------------------------------------------------------------------------

    $mdat add command -label "Color Scheme" -command ColorScheme
    $mdat add command -label "Data Grid" -command DataGrid
    $mdat add command -label "k-path Selection" -command kPath
    $mdat add command -label "Movie Maker" -command [list MovieMaker $can]
    $mdat add separator
    $mdat add command -label "Periodic Table of Elements" -command [list ptable .]

    #------------------------------------------------------------------------
    # HELP MENU
    #------------------------------------------------------------------------

    $mhelp add command -label "About" -command xcAbout

}


proc OpenXYZPDB {can {format xyz} {file {}}} {
    global fileselect mode2D mode3D light xcMisc species system

    if { $file == {} } {
	set dir $system(PWD)
    } elseif { [file isdirectory $file] } {
	set dir $file
	set file {}
    } else {
	set dir $system(PWD)
    }

    if { $species == {} } { set species structure }

    if { $file == {} } {
	set title "Open [string toupper $format] File"
	set filetypes {
	    { {XSF File (XCrySDen Structure File)}  {.xsf} }
	    { {XYZ File}                 {.xyz} }	    
	    { {PDB File}                 {.pdb} }
	    { {All Files}                *      }
	}
	if { $format == "xyz" } {
	    set filetypes {
		{ {XYZ File}                 {.xyz} }
		{ {PDB File}                 {.pdb} }
		{ {XSF File (XCrySDen Structure File)}  {.xsf} }
		{ {All Files}                *      }
	    }
	} elseif { $format == "pdb" } {
	    set filetypes {
		{ {PDB File}                 {.pdb} }
		{ {XYZ File}                 {.xyz} }
		{ {XSF File (XCrySDen Structure File)}  {.xsf} }
		{ {All Files}                *      }
	    }
	}
	set file [tk_getOpenFile -defaultextension .$format \
		-filetypes $filetypes \
		-initialdir $dir \
		-title $title]
	if { $file == "" } {
	    return
	}    
    }
    if { ![file exists $file] } {
	tk_dialog .update \
		"WARNING !!!" "WARNING: File \"$file\" does not exist !!!" \
		warning 0 OK
	return
    }

    reloadRegister $file PDB
    set file [gunzipXSF $file]

    ResetDispModes    
    if { [catch {xc_openstr $format $file $can PL}] } {
	tk_dialog [WidgetName] ERROR \
		"ERROR: An Error occured, while reading file $file" error 0 OK
	CloseCase
	return
    }
    Get_sInfoArray
    DisplayDefaultMode
    # append 'render' to XCState(state) if "render" is not defines yet
    xcAppendState render
    xcUpdateState
    
    #
    # now update the title of the program 
    #
    set xcMisc(titlefile) $file
    wm title . "XCrySDen: [file tail $file]"

    # the following copy is used if someone wants to save the structure into the file
    file copy -force $file [file join $system(SCRDIR) xc_struc.$system(PID)]
    return
}



proc OpenXSFStruct {can {file {}}} {
    global fileselect mode2D mode3D light xcMisc species system

    if { $species == {} } { set species structure }

    set update_title 0
    if { $file == "" } { 
	set update_title 1
	#fileselect "Open Structure"
        set fileselect(path) [tk_getOpenFile -defaultextension .xsf \
                                  -filetypes { 
                                      {{All Files}                  {.*}  }
                                      {{XSF Files}                  {.xsf}}
                                      {{GZipped XSF Files}          {.xsf.gz}}
                                      {{Animated XSF Files}         {.axsf}}
                                      {{GZipped Animated XSF Files} {.axsf.gz}}
                                  } -initialdir $system(PWD) -title "Open XSF Structure"]
	if { $fileselect(path) != "" } {
	    set file $fileselect(path)
	} else {
	    puts stderr "WARNING:: \fileselect(path) = \"\""
	    flush stderr
	    return
	}
    }
    
    if ![file exists $file] {
	tk_dialog .update \
		"WARNING !!!" "WARNING: File \"$file\" does not exist !!!" \
		warning 0 OK
	return
    }

    set file [gunzipXSF $file]

    ResetDispModes
    if [catch {xc_openstr xcr $file $can PL}] {
	tk_dialog [WidgetName] ERROR \
		"ERROR: An Error occured, while reading file $file" error 0 OK
	CloseCase
	return
    }
    Get_sInfoArray
    DisplayDefaultMode

    # append 'render' to XCState(state) if "render" is not defines yet
    xcAppendState render
    xcUpdateState

    CrysFrames
    #
    # now update the title of the program 
    # (only if it does not yet contain the filename in the title)
    #
    #if { [string first : [wm title .]] == -1 } {
    #	wm title . "XCrySDen: [file tail $file]"
    #}

    # update the title of the program if $update_title
    if { $update_title } {
	set xcMisc(titlefile) $file
	wm title . "XCrySDen: [file tail $file]"
	# the following copy is used if someone wants to save file as 
	# WIEN2k STRUCT FILE
	file copy -force $file $system(SCRDIR)/xc_struc.$system(PID)
    }
    return
}



proc UpdateStruct {can file {fmt xsf}} {
    global colSh select system

    #
    # first check if we are in selection mode
    #
    if { [info exists select(selection_mode)] } {
	if { $select(selection_mode) } {
	    tk_dialog [WidgetName] "WARNING !!!" "WARNING: structure can not be updated while in selection mode. First exit from selection mode, then repeat operation" warning 0 OK
	    return
	}
    }
    
    if { ![file exists $file] } {
	tk_dialog .update \
	    "WARNING !!!" "WARNING: File \"$file\" does not exist !!!" \
	    warning 0 OK
	return
    }

    if { [catch {xc_updatestr $file $can $fmt} errMsg] } {
	ErrorDialog "an error occured while trying to update a structure; Error reading file $file\nError Message: $errMsg"
        #InitGlobalVar    
	#ResetDispModes
	#DisplayDefaultMode
	#xcUpdateState
	#wm title . "XCrySDen"
	#set xcMisc(titlefile) {}
	#array set check {crds 0  labels 0  frames 0  wigner 0}
	#array set radio {cellmode prim}
	return 0
    }

    # save $file to xc_struc.$system(PID) so that SaveXSFStruct will work
    file copy -force $file $system(SCRDIR)/xc_struc.$system(PID)

    # there is a bug associated with crystal-frames in xcrys interpreter;
    # this is to recover from the bug
    CrysFrames
    
    if { [info exists colSh] } {
	if { $colSh(scheme) != "atomic" } {
	    #
	    # update colorscheme
	    #
	    ColorSchemeUpdate .mesa
	}
    }
    return 1
}
  


proc SaveXSFStruct {{sfile {}}} {
    global system xcMisc

    # this option (SaveXSFStruct is possible only if structure is opened, that 
    # means that render is active
    if { ![xcIsActive render] } { return }

    if { $sfile == "" } {
	set filetypes {
	    {{XCrySDen Structure File} {.xsf} }
	    {{All Files}          *           }
	}
	
	set sfile [tk_getSaveFile -initialdir $system(PWD) \
		       -title "Save Structure in XSF Format" \
		       -defaultextension .xsf \
		       -filetypes $filetypes]
	
	# maybe Cancel button was pressed
	if { $sfile == {} } { return }
    }

    file copy -force $system(SCRDIR)/xc_struc.$system(PID) $sfile

    # update the title of the program
    wm title . "XCrySDen: [file tail $sfile]"
    set xcMisc(titlefile) $sfile
}



proc SaveCrystalInput {} {
    global system xcMisc

    # this option (SaveCrystalInput is posible only if 
    # XCState == c95_openinput || XCState == c95_newinput
    if ![xcIsActive c95] {
	return
    } else {
	if { ![xcIsActive newinput] && ![xcIsActive openinput] } {
	    return
	}
    }

    set sfile [tk_getSaveFile -initialdir $system(PWD) \
	    -title "Save Crystal Input"]
    # maybe Cancel button was pressed
    if { $sfile == {} } { return }
    set geoInput [MakeInput]
    set fileID [open $sfile w]
    puts $fileID $geoInput
    flush $fileID
    close $fileID

    # update the title of the program
    wm title . "XCrySDen: [file tail $sfile]"
    set xcMisc(titlefile) $sfile
    return
}

	

# this proc communicate with gengeom program
proc CellMode {{update {0}}} {
    global nxdir nydir nzdir system radio geng periodic geng

    if { $periodic(dim) == 0 } { 
	GenGeom $geng(M1_PRIM) 1 $geng(M3_ARGUMENT) 1  1 1 1 $system(SCRDIR)/xc_struc.$system(PID)
	if { $update } { 
	    puts stderr "UpdateStruct .mesa $system(SCRDIR)/xc_struc.$system(PID)"
	    UpdateStruct .mesa $system(SCRDIR)/xc_struc.$system(PID)
	}
	return 
    }

    # usage of "gengeom" program:
    # 
    # gengeom  MODE1  MODE2  MODE3  IGRP  NXDIR  NYDIR  NZDIR  OUTPUT  INPUT
    #    0       1      2      3      4     5      6      7      8       9
    #

    # MODE1 .... radio(cellmode), radio(hexamode)
    # MODE2 .... radio(unitrep)
    # MODE3 .... still free

    ########################################
    # CD TO $system(SCRDIR)
    cd $system(SCRDIR)
    ########################################

    puts stdout "cellmode, hexamode, unitrep:: \
	    $radio(cellmode)   $radio(hexamode)   $radio(unitrep)"
    puts stdout "igroup:: $periodic(igroup)"
    flush stdout
    # if igroup < 7 disable Hexagonal/Rhomohedral entry, else enable
    # if igroup == 7 and "prim" -> omly parapipedal shape is possible
    xcHexaRhomboEntryState

    if { $radio(unitrep) == "cell"} { 
	set MODE2 $geng(M2_CELL) 
    } elseif  {  $radio(unitrep) == "asym"} {
	# radio(hexamode) must be set to parapipedal
	set radio(hexamode) "parapipedal"
	# we should disable all hexagonal-shaped menus !!!!!
	puts stdout "TTT_L"
	flush stdout
	.menu.vmdis.menu entryconfig "Hexagonal/Rhombohedral ..." \
		-state disabled
	set MODE2 $geng(M2_TR_ASYM_UNIT) 
    }

    set MODE3 $geng(M3_ARGUMENT)

    puts stdout "$geng(M1_PRIM) $MODE2 $MODE3"
    puts stdout "igroup>> $periodic(igroup)" 

    if { $radio(cellmode) == "prim" } {
	if { $periodic(igroup) < 8 } {
	    GenGeom $geng(M1_PRIM) $MODE2 $MODE3 $periodic(igroup) \
		$nxdir $nydir $nzdir $system(SCRDIR)/xc_struc.$system(PID)
	}
	if { $periodic(dim) == 3 && $periodic(igroup) >= 8 } {
	    if { $radio(hexamode) == "parapipedal" } {
		GenGeom $geng(M1_PRIM) $MODE2 $MODE3 $periodic(igroup) \
		    $nxdir $nydir $nzdir \
		    $system(SCRDIR)/xc_struc.$system(PID)
	    } elseif { $radio(hexamode) == "hexagonal" && \
			   $radio(unitrep) == "cell" } {
		GenGeom $geng(M1_PRIM3) $MODE2 $MODE3 $periodic(igroup) \
		    $nxdir $nydir $nzdir \
		    $system(SCRDIR)/xc_struc.$system(PID)
	    } else {
		puts stdout "exec $system(BINDIR)/gengeom \
			$geng(M1_PRIM) $MODE2 $MODE3 \
			$periodic(igroup) \
			$nxdir $nydir $nzdir \
			$system(SCRDIR)/xc_struc.$system(PID)"	
		# this is hexagonal & tr_asym_unit -> that's not allowed
		tk_dialog .diag "Error" \
		    "Error in program; Please report to author!!!\n\
			Code: CellMode hexagonal & tr_asym_unit" \
		    error 0 OK
		return
	    }
	}
    } elseif { $radio(cellmode) == "conv" } {
	# dim must be 3
	if { $periodic(dim) < 3 } {
	    #
	    # conventional cell mode is supported only for CRYSTALs.
	    # Revert to PRIMITIVE cell mode !!!
	    set radio(cellmode) prim
	    CellMode 1
	    return
	}
	if { $periodic(igroup) < 7 || $radio(hexamode) == "parapipedal" } {
	    # just parapipedal shaped cells
	    GenGeom $geng(M1_CONV) $MODE2 $MODE3 $periodic(igroup) \
		$nxdir $nydir $nzdir \
		$system(SCRDIR)/xc_struc.$system(PID)
	} elseif { $radio(hexamode) == "hexagonal" } {
	    GenGeom $geng(M1_HEXA_SHAPE) $MODE2 $MODE3 $periodic(igroup) \
		$nxdir $nydir $nzdir \
		$system(SCRDIR)/xc_struc.$system(PID)
	} else {
	    tk_dialog .diag "Error" \
		"Error in program; Please report to author!!!\n\
		    Code: CellMode conv & shape" \
		error 0 OK
	    return 
	}	
    }	
    if { $update } { 
	UpdateStruct .mesa $system(SCRDIR)/xc_struc.$system(PID)
    }
}



# proc check the group and enable/disable Hexagonal/Rhombohedral entry in
# Display menu 
proc xcHexaRhomboEntryState {} {
    global periodic radio
    
    xcDebug "xcHexaRhomboEntryState: igroup == $periodic(igroup)"
    if { $periodic(igroup) < 7 } {
	.menu.vmdis.menu entryconfig "Hexagonal/Rhombohedral ..." \
		-state disabled
    } else {
	.menu.vmdis.menu entryconfig "Hexagonal/Rhombohedral ..." \
		-state normal
	#
	# t.k.: this is temporal
	global xcMisc
	if ![info exists xcMisc(0.2.x_hexa/rhombo_debug)] {
	    .menu.vmdis.menu entryconfig "Hexagonal/Rhombohedral ..." \
		    -state disabled
	}
    }	
    if { $periodic(igroup) == 7 } {
	if { $radio(cellmode) == "prim" } {
	    .menu.vmdis.menu.1 entryconfig 2 -state disabled
	} elseif { $radio(cellmode) == "conv" } {
	    .menu.vmdis.menu.1 entryconfig 2 -state normal
	    #
	    # t.k.: this is temporal
	    global xcMisc
	    if ![info exists xcMisc(0.2.x_hexa/rhombo_debug)] {
		.menu.vmdis.menu entryconfig "Hexagonal/Rhombohedral ..." \
			-state disabled
	    }
	}   
    }
}



# constants for xc_getdefault, xc_newvalue, xc_resetvar
proc ModConst {} {
    global mody
    
    #WARNING: this definitions here must be identical to definitions in
    #         "struct.h"
    
    #this definitions here are used for xc_resetvar <var>

    set mody(R_ATRAD)                0
    set mody(R_RCOV)                 10000
    set mody(R_RBALL)                1
    set mody(R_RROD)                 2
    set mody(R_ATCOL)                3
    set mody(R_WFLINEWIDTH)          4
    set mody(R_WF3DLINEWIDTH)        10004
    set mody(R_OUTLINEWIDTH)         10005
    set mody(R_PLLINEWIDTH)          5
    set mody(R_PLRADIUS)             6
    set mody(R_SCRF)                 7
    set mody(R_ALL)                  8
    set mody(R_FRAMECOL)             9
    set mody(R_FRAMELINEWIDTH)       10
    set mody(R_FRAMERODF)            11
    set mody(R_BACKGROUND)           23
    set mody(R_TESSELLATION)         24
    set mody(R_UNIBOND)              25
    set mody(R_UNIBONDCOLOR)         26
    set mody(R_PERSPECTIVE)          27
    set mody(R_PERSPECTIVEBACK)      28
    set mody(R_PERSPECTIVEFOVY)      29
    set mody(R_PERSPECTIVEFRONT)     10029
    set mody(R_FOG)                  30
    set mody(R_ANTIALIAS)            31
    set mody(R_AMBIENT_BY_DIFFUSE)   32
    set mody(R_CURRENTFILEFORMAT)    33
    set mody(R_HBOND)                34
    set mody(R_FORCE_RODTHICKF)      35
    set mody(R_FORCE_ARRTHICKF)      36
    set mody(R_FORCE_ARRLENF)        37
    set mody(R_FORCE_COLOR)          38
    set mody(R_XYZ_AXIS_COLOR)       39
    set mody(R_XYZ_XYPLANE_COLOR)    40

    #this definitions here are used for xc_newvalue <var> ?values ....?

    set mody(L_SPACE_COV)            0 
    set mody(L_SPACE_VDW)            1 
    set mody(L_BALL_COV)             2
    set mody(L_BALL_VDW)             3
    set mody(L_RCOV_ONE)             10004 
    set mody(L_ATRAD_ONE)            4 
    set mody(L_ATRAD_SCALE)          5 
    set mody(L_BALLF)                6 
    set mody(L_RODF)                 7 
    set mody(L_ATCOL_ONE)            8 
    set mody(L_WFLINEWIDTH)          9 
    set mody(L_WF3DLINEWIDTH)        10009
    set mody(L_OUTLINEWIDTH)         10010
    set mody(L_PLLINEWIDTH)          10 
    set mody(L_PLRADIUS)             11
    set mody(L_SCRF)                 12
    set mody(L_COV_SCALE)            13
    set mody(L_XYZ_ON)               14
    set mody(L_LABELS_ON)            15
    set mody(L_FRAME_ON)             16
    set mody(L_FRAMECOL)             17
    set mody(L_FRAMELINEWIDTH)       18
    set mody(L_FRAMERODF)            19
    set mody(L_LINEFRAME)            20
    set mody(L_TR_XTRANSL)           21
    set mody(L_TR_YTRANSL)           22
    set mody(L_BACKGROUND)           23
    set mody(L_TESSELLATION)         24
    set mody(L_UNIBOND)              25
    set mody(L_UNIBONDCOLOR)         26
    set mody(L_PERSPECTIVE)          27
    set mody(L_PERSPECTIVEBACK)      28
    set mody(L_PERSPECTIVEFOVY)      29
    set mody(L_PERSPECTIVEFRONT)     10029
    set mody(L_FOG)                  30
    set mody(L_ANTIALIAS)            31
    set mody(L_AMBIENT_BY_DIFFUSE)   32
    set mody(L_CURRENTFILEFORMAT)    33
    set mody(L_HBOND)                34
    set mody(L_FORCE_RODTHICKF)      35
    set mody(L_FORCE_ARRTHICKF)      36
    set mody(L_FORCE_ARRLENF)        37
    set mody(L_FORCE_COLOR)          38
    set mody(L_XYZ_AXIS_COLOR)       39
    set mody(L_XYZ_XYPLANE_COLOR)    40

    # this definitions here are used for xc_getdefault <var>
 
    set mody(D_SCRF)                 0
    set mody(D_BALLF)                1
    set mody(D_RODF)                 2
    set mody(D_COVF)                 3
    set mody(D_ALL)                  4
    set mody(D_WFLINEWIDTH)          5
    set mody(D_WF3DLINEWIDTH)        10005
    set mody(D_OUTLINEWIDTH)         10006
    set mody(D_PLLINEWIDTH)          6
    set mody(D_PLRADIUS)             7
    set mody(D_ATCOL_ONE)            8
    set mody(D_ATRAD_SCALE)          9
    set mody(D_ATRAD_ONE)            10
    set mody(D_RCOV_ONE)             10004 
    set mody(D_FRAMECOL)             11
    set mody(D_FRAMELINEWIDTH)       12
    set mody(D_FRAMERODF)            13
    set mody(D_MAXSTRUCTSIZE)        14
    set mody(D_BACKGROUND)           23
    set mody(D_TESSELLATION)         24
    set mody(D_UNIBOND)              25
    set mody(D_UNIBONDCOLOR)         26
    set mody(D_PERSPECTIVE)          27
    set mody(D_PERSPECTIVEBACK)      28
    set mody(D_PERSPECTIVEFOVY)      29
    set mody(D_PERSPECTIVEFRONT)     10029
    set mody(D_FOG)                  30
    set mody(D_ANTIALIAS)            31
    set mody(D_AMBIENT_BY_DIFFUSE)   32
    set mody(D_CURRENTFILEFORMAT)    33
    set mody(D_HBOND)                34
    set mody(D_FORCE_RODTHICKF)      35
    set mody(D_FORCE_ARRTHICKF)      36
    set mody(D_FORCE_ARRLENF)        37
    set mody(D_FORCE_COLOR)          38
    set mody(D_XYZ_AXIS_COLOR)       39
    set mody(D_XYZ_XYPLANE_COLOR)    40

    # this definitions are used fot xc_getvalue
    set mody(GET_NATOMS)       100
    set mody(GET_NAT)          101
    set mody(GET_SS_MINX)      115
    set mody(GET_SS_MINY)      116
    set mody(GET_SS_MINZ)      117
    set mody(GET_SS_MAXX)      118
    set mody(GET_SS_MAXY)      119
    set mody(GET_SS_MAXZ)      120
    set mody(GET_AT_MINX)      10115
    set mody(GET_AT_MINY)      10116
    set mody(GET_AT_MINZ)      10117
    set mody(GET_AT_MAXX)      10118
    set mody(GET_AT_MAXY)      10119
    set mody(GET_AT_MAXZ)      10120
    
    set mody(GET_ATOMLABEL_LABEL)       121
    set mody(GET_ATOMLABEL_BRIGHTCOLOR) 122
    set mody(GET_ATOMLABEL_DARKCOLOR)   123
    set mody(GET_ATOMLABEL_DO_DISPLAY)  124
    set mody(GET_ATOMLABEL_ALL_ID)      125
    set mody(GET_GLOBALATOMLABEL_BRIGHTCOLOR) 126
    set mody(GET_GLOBALATOMLABEL_DARKCOLOR)   127
    set mody(GET_GLOBALATOMLABEL_DO_DISPLAY)  128

    set mody(GET_FOG_COLORMODE)   130
    set mody(GET_FOG_COLOR)       131
    set mody(GET_FOG_MODE)        132
    set mody(GET_FOG_DENSITY)     133
    set mody(GET_FOG_ORT_START_F) 134
    set mody(GET_FOG_ORT_END_F)   135
    set mody(GET_FOG_PERSP_F1)    136
    set mody(GET_FOG_PERSP_F2)    137

    set mody(GET_ATOMLABEL_LABEL)       121
    set mody(GET_ATOMLABEL_BRIGHTCOLOR) 122
    set mody(GET_ATOMLABEL_DARKCOLOR)   123

    set mody(SET_ATOMLABEL_DO_DISPLAY)       200
    set mody(SET_GLOBALATOMLABEL_DO_DISPLAY) 201
    set mody(SET_DO_NOT_DISPLAY_ATOMLABEL)   202

    set mody(GET_ANTIALIAS_DEGREE) 140
    set mody(GET_ANTIALIAS_OFFSET) 141

    set mody(GET_ALAT)             150

    set mody(SET_FOG_COLORMODE)   210
    set mody(SET_FOG_COLOR)       211
    set mody(SET_FOG_MODE)        212
    set mody(SET_FOG_DENSITY)     213
    set mody(SET_FOG_ORT_START_F) 214
    set mody(SET_FOG_ORT_END_F)   215
    set mody(SET_FOG_PERSP_F1)    216
    set mody(SET_FOG_PERSP_F2)    217
    
    set mody(SET_ANTIALIAS_DEGREE) 220
    set mody(SET_ANTIALIAS_OFFSET) 221
}
    


proc ModFactor { fac l_const d_const {text {}} } {    
    upvar #0 $fac fact
    global fact mody

    if { $text == "" } { set text "$fac Factor" }

    set top .factor
    # just in case in window already exists
    if { [winfo exists $top] } { return } 
    xcToplevel $top $text $text . 0 0 1
    #grab $top
    
    # place $top according to "."
    xcPlace . $top 100 50

    set f [frame $top.f -relief raised -bd 2]
    set f2 [frame $top.f2 -relief raised -bd 2]
    pack $f $f2 -side top -expand 1 -fill both -ipadx 10 -ipady 10

    # TOP FRAME
    puts stdout "GET VALUE:: [xc_getvalue $mody($d_const)]"
    set fact(f) [xc_getvalue $mody($d_const)]
    Entries $f [list "$text:"] fact(f) 8 0
    bind $f.frame.entry1 <Return> [list ModFacFOK $top fact $l_const]
    set b1 [button $f.frame.b1 -text "Default" \
	    -command [list ModFacDef fact $d_const]]
    pack $b1 -side left -side left -padx 5 -pady 5 -expand 1
    focus $f.frame.entry1

    # BOTTOM FRAME
    set clear  [button $f2.cl  -text "Clear"  -command [list ModFacClear fact $f.frame.entry1] ]
    set update [button $f2.upd -text "Update" -command [list ModFacFOK $top fact $l_const update] ]    
    set ok     [button $f2.ok  -text "Close"  -command [list ModFacFOK $top fact $l_const] ]    

    bind $ok <Return> [list ModFacFOK $top fact $l_const]
    pack $clear $update $ok -expand 1 -side left -padx 5
}



proc ModFacFOK {top fac l_const {update ""}} {
    upvar #0 $fac fact
    global mody

    xc_newvalue .mesa $mody($l_const) $fact(f)
    
    if { $update == "" } {
	if { [winfo exists $top] } { 
	    #grab release $top
	    destroy $top
	}    
    }
    return
}



proc ModFacDef {fac d_const} {
    upvar #0 $fac fact
    global mody

    puts stdout "d_const:: $d_const"
    flush stdout
    set fact(f) [xc_getdefault $mody($d_const)]
    puts stdout "set fact(f) [xc_getdefault $mody($d_const)]"
    return
}



proc ModFacClear {fac w} {
    upvar #0 $fac fact
    global fact

    set fact(f) ""
    focus $w
    return
}



proc CrdSist {} {
    global mody check

    xc_newvalue .mesa $mody(L_XYZ_ON) $check(crds)
}



proc AtomLabels {} {
    global mody check

    xc_newvalue .mesa $mody(L_LABELS_ON) $check(labels)
}



proc CrysFrames {} {
    global mody check periodic

    if { $periodic(dim) == 0 } { 
	set check(frames) 0
	return 
    }

    DispFramesAs
    xc_newvalue .mesa $mody(L_FRAME_ON) $check(frames)
}



proc Unibond {} {
    global mody check

    xc_newvalue .mesa $mody(L_UNIBOND) $check(unibond)
}



proc DispFramesAs {} {
    global mody check radio

    # if crystal frames are not taken On, return
    if { $mody(L_FRAME_ON) == 0 } { return }
    # also if we are in 2D mode return !!!!!!!!!!
    
    # OK, proceed
    switch -glob -- $radio(frames) {
	line* { xc_newvalue .mesa $mody(L_LINEFRAME) 1 }
	rod*  { xc_newvalue .mesa $mody(L_LINEFRAME) 0 }
    }
}



# this proc is for modifying crystal cell's color
proc ModCellCol {} {
    global mody cellcol

    set top .cellcol
    # just in case in window already exists
    if { [winfo exists $top] } { return } 
    toplevel $top
    wm title $top "Crystal Cell's  Color"
    wm iconname $top "Crystal Cell's  Color"
    #grab $top

    # place $top according to "."
    xcPlace . $top 100 50

    # there will be three frames left, right and bottom
    set t [frame $top.t -relief raised -bd 2]
    set b [frame $top.b -relief raised -bd 2]

    pack $b -side bottom -padx 0 -pady 0 -ipadx 0 -ipady 0 \
	    -fill both -expand 1
    pack $t -side top -padx 0 -pady 0 -ipadx 0 -ipady 0 \
	    -fill both -expand 1

    # TOP-FRAME --- TOP-FRAME 
    # before we select any atom from above Listbox, a default will be H
    set color [xc_getvalue $mody(D_FRAMECOL)]
    set cellcol(red)   [lindex $color 0]
    set cellcol(green) [lindex $color 1]
    set cellcol(blue)  [lindex $color 2]     

    set cellcol(hxred) [ d2h [expr int($cellcol(red) * 255)] ] 
    set cellcol(hxgreen) [ d2h [expr int($cellcol(green) * 255)] ] 
    set cellcol(hxblue) [ d2h [expr int($cellcol(blue) * 255)] ] 
    
    set fr  [frame $t.1 -relief sunken -bd 2]
    set col [frame $fr.col  -bd 0 \
	    -bg "#$cellcol(hxred)$cellcol(hxgreen)$cellcol(hxblue)" \
	    -width 100 -height 100]
    
    scale $t.red -from 0 -to 1 -length 100 -variable cellcol(red) \
	    -orient horizontal -label "Red:" -tickinterval 1.0 \
	    -digits 3 -resolution 0.01 -showvalue true -width 10 \
	    -command CellColor
    scale $t.green -from 0 -to 1 -length 100 \
	    -variable cellcol(green) \
	    -orient horizontal -label "Green:" -tickinterval 0.5 \
	    -digits 3 -resolution 0.01 -showvalue true -width 10 \
	    -command CellColor
    scale $t.blue -from 0 -to 1 -length 100 \
	    -variable cellcol(blue) \
	    -orient horizontal -label "Blue:" -tickinterval 0.5 \
	    -digits 3 -resolution 0.01 -showvalue true -width 10 \
	    -command CellColor

    pack $fr -side top -fill both -expand 1 -padx 30 -pady 15 -ipadx 0 -ipady 0
    pack $col -side top -fill both -expand 1 -padx 0 -pady 0 
    pack $t.red $t.green $t.blue -side top -fill both -expand 1 \
	    -ipadx 0 -ipady 0  -pady 0

    # BOTTOM FRAME --- BOTTOM FRAME
    # OK, Default & Cancel Button
    set col [button $b.upd -text "Update Color" -width 9 \
	    -command CellColUpD]
    set def [button $b.def -text "Default Color" -width 10 \
	    -command [list CellColLoad default]]
    set clo [button $b.clo -text "Close" -command [list CellColOK $top] -width 3]
    
    pack $col $def $clo -side left -expand 1 -pady 7 -padx 2 
    
    return
}



proc CellColor { {trash {}} } {
    global cellcol

    set cellcol(hxred) [ d2h [expr int($cellcol(red) * 255)] ]
    set cellcol(hxgreen) [ d2h [expr int($cellcol(green) * 255)] ]
    set cellcol(hxblue) [ d2h [expr int($cellcol(blue) * 255)] ] 
    
    .cellcol.t.1.col configure \
	    -bg "#$cellcol(hxred)$cellcol(hxgreen)$cellcol(hxblue)"
}



proc CellColLoad { {type {}} } {
    global mody cellcol
    
    if { $type == "default"} {
	set color [xc_getdefault $mody(D_FRAMECOL)]	
    } else {
	puts stdout "GETVALUE::"
	set color [xc_getvalue $mody(D_FRAMECOL)]
    }

    set cellcol(red) [lindex $color 0]
    set cellcol(green) [lindex $color 1]
    set cellcol(blue) [lindex $color 2] 

    CellColor
}



proc CellColOK {top} {
    global cellcol
    
    if { [winfo exists $top] } { 
	#grab release $top
	destroy $top
    }
    return
}



proc CellColUpD {} {
    global mody cellcol

    # puts stdout "xc_newvalue .mesa $mody(L_ATCOL_ONE) $atcol(nat) \
    #                                   $atcol(red) $atcol(green) $atcol(blue)"
    # flush stdout
    xc_newvalue .mesa $mody(L_FRAMECOL) \
	    $cellcol(red) $cellcol(green) $cellcol(blue)

    return
}



proc ModUnibondCol {} {
    global mody
    
    set top .unibondcol
    if [winfo exists $top] { return } 
    set con [xcUpdateWindow \
		 -name  $top \
		 -title "Unibond Color" \
		 -cancelcom  [list ModUnibondCol_Cancel $top] \
		 -updatecom  ModUnibondCol_Update \
		 -closecom   [list ModUnibondCol_Close $top]]
    
    set hxcol [rgb_f2h [xc_getvalue $mody(D_UNIBONDCOLOR)]]
    xcModifyColor $con "Unibond Color" \
	$hxcol groove left left 100 100 70 5 20
}
proc ModUnibondCol_Update {} {
    global mody mody_col
    set cID [xcModifyColorGetID]
    xc_newvalue .mesa $mody(L_UNIBONDCOLOR) \
	$mody_col($cID,red) $mody_col($cID,green) $mody_col($cID,blue) 1.0
}
proc ModUnibondCol_Close {top} {
    ModUnibondCol_Update
    ModUnibondCol_Cancel $top
}
proc ModUnibondCol_Cancel {top} {
    destroy $top
}



proc ModXYZCol {} {
    global mody xyz
    
    set top .xyz
    if { [winfo exists $top] } { return } 

    set con [xcUpdateWindow \
		 -name  $top \
		 -title "Coordinate System Color" \
		 -cancelcom  [list ModXYZCol_Cancel $top] \
		 -updatecom  ModXYZCol_Update \
		 -closecom   [list ModXYZCol_Close $top]]
    
    set f1 [frame $con.1]
    set f2 [frame $con.2]
    pack $f1 $f2 -side left -padx 3m -pady 3m
    
    set axis_hxcol    [rgb_f2h [xc_getvalue $mody(D_XYZ_AXIS_COLOR)]]
    set xyplane_hxcol [rgb_f2h [xc_getvalue $mody(D_XYZ_XYPLANE_COLOR)]]
    
    xcModifyColor $f1 "Axis Color" \
	$axis_hxcol groove left left 100 100 70 5 20
    set xyz(axisID) [xcModifyColorGetID]
    
    xcModifyColor $f2 "XY-Plane Color" \
	$xyplane_hxcol groove left left 100 100 70 5 20
    set xyz(xyplaneID) [xcModifyColorGetID]
}
proc ModXYZCol_Update {} {
    global mody mody_col xyz

    set cID $xyz(axisID)
    xc_newvalue .mesa $mody(L_XYZ_AXIS_COLOR) \
	$mody_col($cID,red) $mody_col($cID,green) $mody_col($cID,blue) 1.0

    set cID $xyz(xyplaneID)
    xc_newvalue .mesa $mody(L_XYZ_XYPLANE_COLOR) \
	$mody_col($cID,red) $mody_col($cID,green) $mody_col($cID,blue) 1.0
}
proc ModXYZCol_Close {top} {
    ModXYZCol_Update
    ModXYZCol_Cancel $top
}
proc ModXYZCol_Cancel {top} {
    destroy $top
}

		    

proc CloseCase {} {
    global xcMisc radio check system

    #
    # it is much saver to start new XCRYSDEN application, then
    # to do a complete clean-up, since the is a great danger that something is
    # still left
    #
    cd $system(PWD)
    global env tcl_platform
    exec sh $system(TOPDIR)/xcrysden --quiet &
    update
    after idle {
	exit_pr 1
    }
}


proc xcMesaChangeBg mesa {
    global radio radio mesa_bg

    if ![info exists mesa_bg($mesa,last)] {
	set mesa_bg($mesa,last) $mesa_bg(default)
    }

    switch -exact -- $radio($mesa,bg) {
	black    {set col $mesa_bg(black)}
	white    {set col $mesa_bg(white)}
	red      {set col $mesa_bg(red)}
	green    {set col $mesa_bg(green)}
	blue     {set col $mesa_bg(blue)}
	darkcyan {set col $mesa_bg(darkcyan)}
	Default  {set col $mesa_bg(default)}
	default  {set col $mesa_bg($mesa,last)}
    }
    
    set c1 [d2h [expr [lindex $col 0] * 255]]
    set c2 [d2h [expr [lindex $col 1] * 255]]
    set c3 [d2h [expr [lindex $col 2] * 255]]
    xcDebug #${c1}${c2}${c3}
    set t [xcToplevel [WidgetName] "Set Canvas Background" "SetBg" . 0 0 1]
    xcModifyColor $t "Set Canvas Background" #${c1}${c2}${c3} groove \
	    left left 100 100 70 5 20

    proc xcMesaChangeBgOK {type t mesa} {
	global mody_col mody radio mesa_bg
	
	if { $type == "OK" } {
	    set cID [xcModifyColorGetID]
	    set mesa_bg($mesa,last) [list $mody_col($cID,red) \
		    $mody_col($cID,green) $mody_col($cID,blue) 1.0]
	    xc_newvalue $mesa $mody(L_BACKGROUND) \
		    $mody_col($cID,red) $mody_col($cID,green) \
		    $mody_col($cID,blue) 1.0
	    set radio($mesa,bg) [rgb_f2h [list $mody_col($cID,red) \
		    $mody_col($cID,green) $mody_col($cID,blue)]]
	}
	destroy $t
    }

    set ok  [DefaultButton [WidgetName $t] -text "OK" \
	    -command [list xcMesaChangeBgOK OK $t $mesa]]
    set can [button [WidgetName $t] -text "Cancel" \
	    -command [list xcMesaChangeBgOK Cancel $t $mesa]]
    pack $ok $can -padx 10 -pady 10 -expand 1
}


proc Perspective {} {
    global mody check select

    if { [info exists select(selection_mode)] } {
	if { $select(selection_mode) } {
	    return
	}
    }

    SetWatchCursor
    # we zoom for perspective for 33%. This makes the picture to
    # appear almost the same size for ortho and perspective
    # projections
    if { $check(perspective) } {
	xc_translate .mesa +z 0.33
    } else {
	xc_translate .mesa -z 0.33
    }
    xc_newvalue .mesa $mody(L_PERSPECTIVE) $check(perspective)
    ResetCursor
}
#proc TogglePerspective {} {
#    global check
#    if { $check(perspective) } {
#	set check(perspective) 0
#	Perspective
#    } else {
#    	set check(perspective) 1
#	Perspective
#    }
#}

proc DepthCuing {} {
    global mody check
    SetWatchCursor
    xc_newvalue .mesa $mody(L_FOG) $check(depthcuing)
    ResetCursor
}
proc AntiAlias {} {
    global mody check
    SetWatchCursor
    xc_newvalue .mesa $mody(L_ANTIALIAS) $check(antialias)
    ResetCursor
}


proc ToggleMenuCheckbutton {which args} {    
    global check

    SetWatchCursor
    if { $check($which) } {
	set check($which) 0
	eval $args
    } else {
    	set check($which) 1
	eval $args
    }
    ResetCursor
}

proc ToggleMenuRadiobutton {which value args} {
    global radio periodic

    if { $periodic(dim) > 0 } { 
	set radio($which) $value
	SetWatchCursor
	eval $args
	ResetCursor
    }
}


proc ModPerspective {} {
    global mody persp

    set top .perspective
    if { [winfo exists $top] } { return } 
    set con [xcUpdateWindow \
		 -name  $top \
		 -title "Perspective Settings" \
		 -cancelcom  [list ModPerspective_Cancel $top] \
		 -updatecom  ModPerspective_Update \
		 -closecom   [list ModPerspective_Close $top]]
    
    set f1 [frame $con.f1]
    set f2 [frame $con.f2]
    set f3 [frame $con.f3]
    pack $f1 $f2 $f3 -side top -fill both -expand 1
    set persp(fovy)  [xc_getvalue $mody(D_PERSPECTIVEFOVY)]
    set persp(front) [xc_getvalue $mody(D_PERSPECTIVEFRONT)]
    set persp(back)  [xc_getvalue $mody(D_PERSPECTIVEBACK)]

    set near [button $f1.b -text Deafult -command {
	global persp mody
	set persp(fovy) [xc_getdefault $mody(D_PERSPECTIVEFOVY)]
    }]	      
    set size [button $f2.b -text Default -command {
	global persp mody
	set persp(front) [xc_getdefault $mody(D_PERSPECTIVEFRONT)]
    }]
    set far [button $f3.b -text Default -command {
	global persp mody
	set persp(back) [xc_getdefault $mody(D_PERSPECTIVEBACK)]
    }]
    pack $near $far $size -side right -padx 3 -expand 1

    FillEntries $f1 {"Perspective fovy factor:"}  persp(fovy) 24 8
    FillEntries $f2 {"Perspective front factor:"} persp(front) 24 8
    FillEntries $f3 {"Perspective back factor:"}  persp(back) 24 8
}
proc ModPerspective_Cancel {top} {
    destroy $top
}
proc ModPerspective_Update {} {
    global mody persp
    SetWatchCursor
    xc_newvalue .mesa $mody(L_PERSPECTIVEFOVY) $persp(fovy)
    xc_newvalue .mesa $mody(L_PERSPECTIVEBACK)  $persp(back)
    xc_newvalue .mesa $mody(L_PERSPECTIVEFRONT) $persp(front)
    ResetCursor
}
proc ModPerspective_Close {top} {
    ModPerspective_Update
    ModPerspective_Cancel $top
}


proc scriptOpenMenu {{togl .mesa}} {
    global system

    set file [tk_getOpenFile -defaultextension .xcrysden \
		  -filetypes { 
		      {{All Files}                          {.*}  }
		      {{XCrySDen Scripting File}            {.xcrysden}}
		      {{XCrySDen Scripting File}            {.tcl}}
		      {{GZipped XCrySDen Scripting File}    {.xcrysden.gz}}
		  } \
		  -initialdir $system(PWD) \
		  -title "Open XCrySDen Scripting File"]
    if { $file == "" } {
	return
    }
    scripting::source $file
}
