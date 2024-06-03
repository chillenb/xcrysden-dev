# check which type of coorinates exists in a pw.x output file

BEGIN {
    initial      = 0;
    dynamics     = 0;
    intermediate = 0;
    optimized    = 0;
}

/Cartesian axes/ || /Carthesian axes/  { initial = 1; }
/   Entering Dynamics/                 { dynamics = 1; } 
/ATOMIC_POSITIONS/ || /Search of equilibrium positions/ || /Entering Dynamics/ { intermediate = 1; }
/Begin final coordinates/ || /Final estimate of positions/ || /   Final energy/ { optimized = 1; }

END {
    if (initial)      { printf "%s ", "initial"; }
    if (dynamics)     { printf "%s ", "dynamics"; }
    if (intermediate) { printf "%s ", "intermediate"; }
    if (optimized)    { printf "%s ", "optimized"; }
}
