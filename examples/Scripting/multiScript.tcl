# ------------------------------------------------------------------------
#****** ScriptingExamples/multiScript.tcl ***
#
# NAME
# multiScript.tcl -- a simple example of multi-job script (i.e. multiScript)
#
# USAGE
# xcrysden --script multiScript.tcl
#
# COPYRIGHT
# Anton Kokalj (C) 2003
#
# PURPOSE

# This is a scripting example that shows how to produce a multi-jobs.
# It uses the scripting::multiScript facility. Many times, one might
# want to produce several plots of molecular orbitals of a given
# molecule. It would be desirable that the display parameters are
# exactly the same for all plots. This is one such example, where CO
# HOMO and LUMO molecular-orbitals are printed. 
#
# The multiScript uses major and minor script. In the first script the
# loading of the files is specified, while in the latter script the
# operations of what to do with each file are defined. The usage of
# multiScript is:
#
# scripting::multiScript scriptMajor scriptMinor
#
# scripting::multiScript scriptMajor scriptMinor


#
# WARNINGS
# Inside the major-script (i.e. scriptMajor) the scripting::exec
# should be used instead of scripting::open.

#
# AUTHOR
# Anton Kokalj
#
# CREATION DATE
# Sometime in February 2003
# 
# SOURCE

scripting::multiScript {

    # ------------------------------------------------------------------------
    # This is the MAJOR script 
    #
    # It opens two files, one at a time, and then for each opened
    # file proceed according to minor script instructions
    # ------------------------------------------------------------------------
    global env
    set dir [file join $env(XCRYSDEN_TOPDIR) examples XSF_Files]
    
    set file1 [file join $dir CO_homo.xsf.gz]
    set file2 [file join $dir CO_lumo.xsf.gz]
    
    foreach file [list $file1 $file2] {
	# load the $file XSF
	scripting::exec --xsf $file
	
	# rename print.png to $file.png
	if { [file exists print.png] } {
	    set filehead [file tail [file rootname $file]]
	    file rename -force print.png $filehead.png
	}
    }

} {

    # ------------------------------------------------------------------------
    # This is the MINOR script
    #
    # It produces and prints the isosurface + colorplane
    # ------------------------------------------------------------------------


    # ------------------------------------------------------------------------
    # display the structure in appropriate display-mode
    # ------------------------------------------------------------------------
    
    scripting::lighting On
    scripting::displayMode3D Pipe&Ball

    # ------------------------------------------------------------------------
    # zoom and rotate the structure 
    # ------------------------------------------------------------------------

    scripting::zoom -0.5
    scripting::rotate x -90 

    # ------------------------------------------------------------------------
    # load the 3D scalar field
    # ------------------------------------------------------------------------
    
    scripting::scalarField3D::load

    # ------------------------------------------------------------------------
    # configure, i.e., specify how to render the scalar field
    #
    # for the usage see, for example, isosurface+colorplane+print.tcl file
    # ------------------------------------------------------------------------

    scripting::scalarField3D::configure \
	-isosurface           1 \
	-interpolation_degree 2 \
	-isolevel             0.1 \
	-plusminus            1 \
	-basalplane           2 \
	-colorbasis           BLUE-WHITE-RED \
	-scalefunction        LINEAR \
	-expand2D             specify \
	-expand2D_X 	      1 \
	-expand2D_Y           1 \
	-expand2D_Z           1 \
	-colorplane           1 \
	-isoline              1 \
	-colorplane_lighting  0 \
	-cpl_transparency     0 \
	-cpl_thermometer      1 \
	-2Dlowvalue           -0.1 \
	-2Dhighvalue          +0.1 \
	-2Dnisoline           11 \
	-anim_step            1 \
	-current_slide        25 \
	-isoline_color        monocolor \
	-isoline_width        3 \
	-isoline_monocolor    \#000000

    # ------------------------------------------------------------------------
    # hide the isosourface control window
    # ------------------------------------------------------------------------
    
    wm withdraw .iso

    # ------------------------------------------------------------------------
    # render the 3D scalar field as requested by 
    # scripting::scalarField3D::configure
    # ------------------------------------------------------------------------
    
    scripting::scalarField3D::render
    
    
    # # ------------------------------------------------------------------------
    # # revert the isosurface normals (this should be done after rendering of 
    # # isosurface)
    # # ------------------------------------------------------------------------
    # 
    # scripting::scalarField3D::configure -revertnormal {pos neg}

    # ------------------------------------------------------------------------
    # now lets print to file what we have on the display window
    # ------------------------------------------------------------------------

    scripting::printToFile print.png 

    # we've done all; exit
    exit 0
}

#******
