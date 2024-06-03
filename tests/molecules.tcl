scripting::multiScript {
    global env

    foreach job {
	{--xyz      $env(XCRYSDEN_TOPDIR)/examples/XYZ/mol2.xyz}
	{--pdb      $env(XCRYSDEN_TOPDIR)/examples/GAUSSIAN_files/benzene-6CH3-OCH3.g98}
	{--xsf      $env(XCRYSDEN_TOPDIR)/examples/XSF_Files/1symb.xsf}
	{--gzmat    $env(XCRYSDEN_TOPDIR)/examples/GAUSSIAN_files/benzene-6CH3-OCH3.g98}
	{--g98_out  $env(XCRYSDEN_TOPDIR)/examples/GAUSSIAN_files/benzene.g98_out}
    } {
	eval scripting::exec $job
    }
} {

    scripting::lighting On
    scripting::display on atomic-labels

    scripting::displayMode3D Stick
    scripting::rotate xy +5 +1 5
    scripting::displayMode3D Pipe&Ball

    scripting::display on unicolor-bonds

    scripting::rotate xy +5 +1 5
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
    scripting::rotate xy +5 +1 5

    scripting::display off unicolor-bonds

    scripting::displayMode2D BallStick-1 
    scripting::rotate xy +5 +1 5
    scripting::displayMode2D BallStick-2  
    scripting::rotate xy +5 +1 5
    scripting::displayMode2D SpaceFill  
    scripting::rotate xy +5 +1 5

    exit
}