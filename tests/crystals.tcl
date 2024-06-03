scripting::multiScript {
    global env

    foreach job {
	{scripting::filter::pwscfInput  $env(XCRYSDEN_TOPDIR)/examples/PWSCF_files/EthAl001-2x2.inp 2 {1 13   2 6   3 1}}
	
	{scripting::filter::pwscfOutput -oc $env(XCRYSDEN_TOPDIR)/examples/PWSCF_files/EthAl001-2x2.out 2 {1 13   2 6   3 1}}
	
	{scripting::filter::crystalInput $env(XCRYSDEN_TOPDIR)/examples/CRYSTALxx_input_files/argonite.r1}
	{scripting::filter::crystalInput $env(XCRYSDEN_TOPDIR)/examples/CRYSTALxx_input_files/Pt322.r1}
	{scripting::filter::crystalInput $env(XCRYSDEN_TOPDIR)/examples/CRYSTALxx_input_files/polymer.r1}

	{scripting::filter::fhiInpini $env(XCRYSDEN_TOPDIR)/examples/FHI98MD_files/GaAs_inp.ini 3 {Gallium 31   Arsenic 33}}
	{scripting::filter::fhiCoord  $env(XCRYSDEN_TOPDIR)/examples/FHI98MD_files/GaAsSurface_coord.out 2 {
	    Gallium 31	    Arsenic 33	    hy_1.25 1	    hy_0.75 1}
	}

	{--wien_struct $env(XCRYSDEN_TOPDIR)/examples/WIEN_struct_files/ferrocen.struct}

	{--xsf $env(XCRYSDEN_TOPDIR)/examples/XSF_Files/fcc-410-1x1.xsf}	
    } {
	eval scripting::exec $job
    }
} {

    #if { [xcIsActive c95] } {
    #	global crystalInput
    #	CalStru
    #	foreach t $crystalInput(two_toplevels) {
    #	    destroy $t
    #	}	
    #}

    scripting::zoom 0.3

    scripting::lighting On
    scripting::buildCrystal 2 2 2

    scripting::display on atomic-labels

    scripting::displayMode3D Stick
    scripting::rotate xy +5 +1 5
    scripting::displayMode3D Pipe&Ball
    scripting::rotate xy +5 +1 5

    scripting::display on unicolor-bonds

    scripting::displayMode3D BallStick
    scripting::rotate xy +5 +1 5
    scripting::displayMode3D SpaceFill
    scripting::rotate xy +5 +1 5
    
    scripting::lighting Off
    scripting::displayMode2D WireFrame  
    scripting::rotate xy +5 +1 5
    scripting::displayMode2D PointLine  
    scripting::rotate xy +5 +1 5
    scripting::displayMode2D Pipe&Ball  

    scripting::display off unicolor-bonds

    scripting::rotate xy +5 +1 5
    scripting::displayMode2D BallStick-1 
    scripting::rotate xy +5 +1 5
    scripting::displayMode2D BallStick-2  
    scripting::rotate xy +5 +1 5
    scripting::displayMode2D SpaceFill  
    scripting::rotate xy +5 +1 5
    exit
}