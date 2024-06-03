c     parameters for order of arguments
c     IARG_IGRP:
c               1 ..... TRIGONAL  systems
c               2 ..... HEXAGONAL systems
      PARAMETER (
     *     IARG_MODE1 =1,
     *     IARG_MODE2 =2,
     *     IARG_MODE3 =3,
     *     IARG_IGRP  =4,   !IGROUP for TRIGONAL/HEXAGONAL systems
     *     IARG_NXDIR =5,
     *     IARG_NYDIR =6,
     *     IARG_NZDIR =7,
     *     IARG_OUTPUT=8,
     *     IARG_INPUT =9 )
     

c     M1 ....... mode1
c     M2 ....... mode2
      PARAMETER (
     *     M1_INFO=0,  !INFO mode
     *     M1_PRIM=1,  !PRIMITIV CELL; in case of H/R PARAPIPEDAL SHAPE
     *     M1_CONV=2,  !CONVENTIONAL CELL; in case of H/R it is PARAPIPEDAL SHAPE
     *     M1_HEXA_SHAPE=3,!CONV CELL for H/R;  HEXAGONAL SHAPE
     *     M1_PRIM3=4, !PRIMITIV cell for H/R; HEXAGONAL SHAPE
     *     M2_CELL=1,  !CELL is unit of repetition
     *     M2_TR_ASYM_UNIT=2) !BASIS (traslational asymmetric part) is unit of repetition

C     example: IMODE3 = 13 --> M3 = 1, M4 = 3 
      PARAMETER (
     *     M3_ANGS=1,   !input is given in angs.
     *     M3_BOHR=2,   !input is given in bohrs
     *     M4_CR95=1,   !CRYSTAL version is CRYSTAL95
     *     M4_CR98=2)   !CRYSTAL version is CRYSTAL98

