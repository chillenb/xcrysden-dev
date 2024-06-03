/*
 
 *****************************************************************************
 * Author:                                                                   *
 * ------                                                                    *
 *  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  *
 *  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    *
 *  Jozef Stefan Institute                          Fax: x 386 1 477 3811    *
 *  Jamova 39, SI-1000 Ljubljana                                             *
 *  SLOVENIA                                                                 *
 *                                                                           *
 * Source: $XCRYSDEN_TOPDIR/C/cells.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#define XC_CPP_NO_STDIO
#include "struct.h"

#define CON13 0.333333333333333
#define CON23 0.666666666666667

float pc[4][3] = {
  {0.0, 0.0, 0.0},

  {0.0, 0.0, 0.0},
  {0.0, 0.0, 0.0},
  {0.0, 0.0, 0.0}
};
     
float ac[4][3] = {
  {0.0, 0.0, 0.0},
  {0.0, 0.5, 0.5},

  {0.0, 0.0, 0.0},
  {0.0, 0.0, 0.0}
};

float bc[4][3] = {
  {0.0, 0.0, 0.0},
  {0.5, 0.0, 0.5},

  {0.0, 0.0, 0.0},
  {0.0, 0.0, 0.0}
};

float cc[4][3] = {
  {0.0, 0.0, 0.0},
  {0.5, 0.5, 0.0},

  {0.0, 0.0, 0.0},
  {0.0, 0.0, 0.0}
};

float fc[4][3] = {
  {0.0, 0.0, 0.0},
  {0.5, 0.5, 0.0},
  {0.0, 0.5, 0.5},
  {0.5, 0.0, 0.5}
};

float ic[4][3] = {
  {0.0, 0.0, 0.0},
  {0.5, 0.5, 0.5},

  {0.0, 0.0, 0.0},
  {0.0, 0.0, 0.0}
};
 
float rc[4][3] = {
  {0.0, 0.0, 0.0},
  {CON23, CON13, CON13},
  {CON13, CON23, CON23},

  {0.0, 0.0, 0.0}
};

float hc[4][3] = {
  {0.0, 0.0, 0.0},
  {CON23, CON13, 0.0},
  {CON13, CON23, 0.0},

  {0.0, 0.0, 0.0}
};

float tnrc[4][3] = {
  {0.0, 0.0, 0.0},
  {CON23, CON13, 0.0},
  {CON13, CON23, 0.0},

  {0.0, 0.0, 0.0}
};

/* 
   load the positions within all cell types; for more detailed descriptions
   look at gengeom.f
*/
void 
CellTypes(void)
{
  int i, j;
  for (i=0; i<4; i++)
    for (j=0; j<3; j++) {
      xcr.cellpos[XCR_CELL_PC][i][j]   = pc[i][j];
      xcr.cellpos[XCR_CELL_AC][i][j]   = ac[i][j];
      xcr.cellpos[XCR_CELL_BC][i][j]   = bc[i][j];
      xcr.cellpos[XCR_CELL_CC][i][j]   = cc[i][j];
      xcr.cellpos[XCR_CELL_FC][i][j]   = fc[i][j];
      xcr.cellpos[XCR_CELL_IC][i][j]   = ic[i][j];
      xcr.cellpos[XCR_CELL_RC][i][j]   = rc[i][j];
      xcr.cellpos[XCR_CELL_HC][i][j]   = hc[i][j];
      xcr.cellpos[XCR_CELL_TNRC][i][j] = tnrc[i][j];
    }

  xcr.npos[XCR_CELL_PC]   = XCR_CELL_IPC;   
  xcr.npos[XCR_CELL_AC]   = XCR_CELL_IAC;  
  xcr.npos[XCR_CELL_BC]   = XCR_CELL_IBC;  
  xcr.npos[XCR_CELL_CC]   = XCR_CELL_ICC;  
  xcr.npos[XCR_CELL_FC]   = XCR_CELL_IFC;  
  xcr.npos[XCR_CELL_IC]   = XCR_CELL_IIC;  
  xcr.npos[XCR_CELL_RC]   = XCR_CELL_IRC;  
  xcr.npos[XCR_CELL_HC]   = XCR_CELL_IHC;  
  xcr.npos[XCR_CELL_TNRC] = XCR_CELL_ITNRC;
}

