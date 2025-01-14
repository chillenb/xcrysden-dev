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
 * Source: $XCRYSDEN_TOPDIR/C/atoms.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <GL/gl.h>


char *element[] = {
  "X",
  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O", 
  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S", 
  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", 
  "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", 
  "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", 
  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", 
  "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", 
  "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", 
  "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", 
  "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", 
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", 
  "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", 
  "Bk", "Cf", "Es", "Fm"};


double rcovdef[MAXNAT + 1] = { 
  0.38000, 
  0.38000,0.38000,1.23000,0.89000,0.91000,
  0.77000,0.75000,0.73000,0.71000,0.71000,
  1.60000,1.40000,1.25000,1.11000,1.00000,
  1.04000,0.99000,0.98000,2.13000,1.74000,
  1.60000,1.40000,1.35000,1.40000,1.40000,
  1.40000,1.35000,1.35000,1.35000,1.35000,
  1.30000,1.25000,1.15000,1.15000,1.14000,
  1.12000,2.20000,2.00000,1.85000,1.55000,
  1.45000,1.45000,1.35000,1.30000,1.35000,
  1.40000,1.60000,1.55000,1.55000,1.41000,
  1.45000,1.40000,1.40000,1.31000,2.60000,
  2.00000,1.75000,1.55000,1.55000,1.55000,
  1.55000,1.55000,1.55000,1.55000,1.55000,
  1.55000,1.55000,1.55000,1.55000,1.55000,
  1.55000,1.55000,1.45000,1.35000,1.35000,
  1.30000,1.35000,1.35000,1.35000,1.50000,
  1.90000,1.80000,1.60000,1.55000,1.55000,
  1.55000,2.80000,1.44000,1.95000,1.55000,
  1.55000,1.55000,1.55000,1.55000,1.55000,
  1.55000,1.55000,1.55000,1.55000,1.55000};

double rvdw[MAXNAT + 1] = { 
  0.38000, 	
  1.17000,1.40000,1.80000,2.10000,2.10000,
  1.70000,1.58000,1.52000,1.47000,1.60000,
  2.30000,1.70000,2.10000,1.70000,1.80000,
  1.80000,1.78000,1.90000,2.80000,2.10000,
  2.10000,2.10000,2.10000,2.10000,2.10000,
  2.10000,2.10000,1.60000,1.40000,1.40000,
  1.90000,2.10000,1.85000,1.90000,1.85000,
  2.00000,2.10000,2.10000,2.10000,2.10000,
  2.10000,2.10000,2.10000,2.10000,2.10000,
  1.60000,1.70000,1.60000,1.90000,1.00000,
  2.10000,2.06000,1.96000,2.16000,2.10000,
  2.10000,2.10000,2.10000,2.10000,2.10000,
  2.10000,2.10000,2.10000,2.10000,2.10000,
  2.10000,2.10000,2.10000,2.10000,2.10000,
  2.10000,2.10000,2.10000,2.10000,2.10000,
  2.10000,2.10000,1.75000,1.70000,1.50000,
  2.00000,2.10000,2.10000,2.10000,2.10000,
  2.10000,3.50000,3.50000,3.50000,2.10000,
  2.10000,2.10000,2.10000,2.10000,2.10000,
  2.10000,2.10000,2.10000,2.10000,2.10000};

GLfloat DefAtCol[MAXNAT + 1][3] = {
  {0.7700, 1.0000, 0.4200},  /* 0 */
  {1.0000, 1.0000, 1.0000},  /* 1 */
  {0.0000, 1.0000, 1.0000},  /* 2 */
  {0.9500, 0.9500, 0.9500},  /* 3 */
  {0.7079, 0.7079, 0.7079},  /* 4 */
  {0.8779, 0.7079, 0.5579},  /* 5 */
  {0.3000, 0.3000, 0.3000},  /* 6 */
  {0.4500, 0.7200, 1.0000},  /* 7 */
  {0.7000, 0.0000, 0.0000},  /* 8 */
  {0.6100, 0.9000, 0.0000},  /* 9 */
  {0.0000, 0.9500, 0.9500},  /* 10 */
  {0.0000, 0.9500, 0.9500},  /* 11 */
  {0.7800, 0.8900, 0.7800},  /* 12 */
  {0.8500, 0.8500, 0.8500},  /* 13 */
  {1.0000, 0.8300, 0.6600},  /* 14 */
  {0.9500, 0.6500, 0.0000},  /* 15 */
  {0.9500, 0.9500, 0.4500},  /* 16 */
  {0.7100, 1.0000, 0.0000},  /* 17 */
  {0.0000, 1.0000, 1.0000},  /* 18 */
  {0.0000, 0.9500, 0.9500},  /* 19 */
  {0.0000, 0.9500, 0.9500},  /* 20 */
  {0.0000, 0.9500, 0.9500},  /* 21 */
  {0.6000, 0.6000, 0.6000},  /* 22 */
  {0.6445, 0.8047, 0.8569},  /* 23 */
  {0.6445, 0.8047, 0.8569},  /* 24 */
  {0.6445, 0.8047, 0.8569},  /* 25 */
  {0.7500, 0.1500, 0.0000},  /* 26 */
  {0.6445, 0.8047, 0.8569},  /* 27 */
  {0.6445, 0.8047, 0.8569},  /* 28 */
  {0.8200, 0.4500, 0.1400},  /* 29 */
  {0.7200, 0.7500, 0.7500},  /* 30 */
  {0.9500, 0.0000, 0.9500},  /* 31 */
  {0.9500, 0.0000, 0.9500},  /* 32 */
  {0.9500, 0.9500, 0.0000},  /* 33 */
  {0.9500, 0.9500, 0.0000},  /* 34 */
  {0.6500, 0.3000, 0.0000},  /* 35 */
  {0.0000, 1.0000, 1.0000},  /* 36 */
  {0.0000, 0.9500, 0.9500},  /* 37 */
  {0.0000, 0.9500, 0.9500},  /* 38 */
  {0.7500, 0.7500, 0.7500},  /* 39 */
  {0.7500, 0.7500, 0.7500},  /* 40 */
  {0.7500, 0.7500, 0.7500},  /* 41 */
  {0.7500, 0.7500, 0.7500},  /* 42 */
  {0.7500, 0.7500, 0.7500},  /* 43 */
  {0.7500, 0.7500, 0.7500},  /* 44 */
  {0.7500, 0.7500, 0.7500},  /* 45 */
  {1.0000, 1.0000, 1.0000},  /* 46 */
  {1.0000, 1.0000, 1.0000},  /* 47 */
  {0.7500, 0.7500, 0.7500},  /* 48 */
  {0.7500, 0.7500, 0.7500},  /* 49 */
  {0.7500, 0.7500, 0.7500},  /* 50 */
  {0.7500, 0.7500, 0.7500},  /* 51 */
  {0.7500, 0.7500, 0.7500},  /* 52 */
  {0.5000, 0.0000, 0.5000},  /* 53 */
  {0.0000, 1.0000, 1.0000},  /* 54 */
  {0.0000, 0.9500, 0.9500},  /* 55 */
  {0.0000, 0.9500, 0.9500},  /* 56 */
  {0.7500, 0.7500, 0.7500},  /* 57 */
  {0.7500, 0.7500, 0.7500},  /* 58 */
  {0.7500, 0.7500, 0.7500},  /* 59 */
  {0.7500, 0.7500, 0.7500},  /* 60 */
  {0.7500, 0.7500, 0.7500},  /* 61 */
  {0.7500, 0.7500, 0.7500},  /* 62 */
  {0.7500, 0.7500, 0.7500},  /* 63 */
  {0.7500, 0.7500, 0.7500},  /* 64 */
  {0.7500, 0.7500, 0.7500},  /* 65 */
  {0.7500, 0.7500, 0.7500},  /* 66 */
  {0.7500, 0.7500, 0.7500},  /* 67 */
  {0.7500, 0.7500, 0.7500},  /* 68 */
  {0.7500, 0.7500, 0.7500},  /* 69 */
  {0.7500, 0.7500, 0.7500},  /* 70 */
  {0.7500, 0.7500, 0.7500},  /* 71 */
  {0.7500, 0.7500, 0.7500},  /* 72 */
  {0.7500, 0.7500, 0.7500},  /* 73 */
  {0.7500, 0.7500, 0.7500},  /* 74 */
  {0.7500, 0.7500, 0.7500},  /* 75 */
  {0.7500, 0.7500, 0.7500},  /* 76 */
  {0.7500, 0.7500, 0.7500},  /* 77 */
  {0.7500, 0.7500, 0.7500},  /* 78 */
  {1.0000, 0.8500, 0.0000},  /* 79 */
  {0.7500, 0.7500, 0.7500},  /* 80 */
  {0.7500, 0.7500, 0.7500},  /* 81 */
  {0.7500, 0.7500, 0.7500},  /* 82 */
  {0.7500, 0.7500, 0.7500},  /* 83 */
  {0.7500, 0.7500, 0.7500},  /* 84 */
  {0.7500, 0.7500, 0.7500},  /* 85 */
  {0.0000, 1.0000, 1.0000},  /* 86 */
  {0.7500, 0.7500, 0.7500},  /* 87 */
  {0.7500, 0.7500, 0.7500},  /* 88 */
  {0.7500, 0.7500, 0.7500},  /* 89 */
  {0.7500, 0.7500, 0.7500},  /* 90 */
  {0.7500, 0.7500, 0.7500},  /* 91 */
  {0.7500, 0.7500, 0.7500},  /* 92 */
  {0.7500, 0.7500, 0.7500},  /* 93 */
  {0.7500, 0.7500, 0.7500},  /* 94 */
  {0.7500, 0.7500, 0.7500},  /* 95 */
  {0.7500, 0.7500, 0.7500},  /* 96 */
  {0.7500, 0.7500, 0.7500},  /* 97 */
  {0.7500, 0.7500, 0.7500},  /* 98 */
  {0.9314, 0.8755, 0.8010},  /* 99 */
  {0.9314, 0.8755, 0.8010}  /* 100 */
};

GLfloat DefAtColOld[MAXNAT + 1][3] = {
  {0.7700, 1.0000, 0.4200},  /* this is color of selected atom */
  {0.0000, 0.9500, 0.9500},
  {0.9500, 0.9500, 0.9500},
  {0.9500, 0.9500, 0.9500},
  {0.7079, 0.7079, 0.7079},
  {0.7079, 0.7079, 0.7079},
  {0.9500, 0.9500, 0.0000},
  {0.6445, 0.8047, 0.8569},
  {0.7000, 0.0000, 0.0000}, /* changed */
  {0.8345, 0.9500, 0.9500},
  {0.9500, 0.9500, 0.9500},
  {0.0000, 0.9500, 0.9500},
  {0.0000, 0.9500, 0.9500},
  {0.9500, 0.0000, 0.9500},
  {0.0000, 0.9500, 0.9500},
  {0.9500, 0.9500, 0.0000},
  {0.9500, 0.9500, 0.4500},
  {0.7100, 1.0000, 0.0000},
  {0.9500, 0.9500, 0.9500},
  {0.0000, 0.9500, 0.9500},
  {0.0000, 0.9500, 0.9500},
  {0.0000, 0.9500, 0.9500},
  {0.6000, 0.6000, 0.6000},
  {0.6445, 0.8047, 0.8569},
  {0.6445, 0.8047, 0.8569},
  {0.6445, 0.8047, 0.8569},
  {0.9500, 0.0000, 0.0000},
  {0.6445, 0.8047, 0.8569},
  {0.6445, 0.8047, 0.8569},
  {0.8200, 0.4500, 0.1400},
  {0.9500, 0.0000, 0.9500},
  {0.9500, 0.0000, 0.9500},
  {0.9500, 0.0000, 0.9500},
  {0.9500, 0.9500, 0.0000},
  {0.9500, 0.9500, 0.0000},
  {0.9500, 0.0000, 0.0000},
  {0.7500, 0.7500, 0.7500},
  {0.0000, 0.9500, 0.9500},
  {0.0000, 0.9500, 0.9500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {1.0000, 1.0000, 1.0000}, /* changed */
  /*  {0.7500, 0.7500, 0.7500},  silver */
  {1.0000, 1.0000, 1.0000}, /* silver-white */
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.0000, 0.9500, 0.9500},
  {0.0000, 0.9500, 0.9500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {1.0000, 0.8500, 0.0000},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.7500, 0.7500, 0.7500},
  {0.9314, 0.8755, 0.8010},
  {0.9314, 0.8755, 0.8010}
};

GLclampf DefFrameCol[3]   = { 0.88, 1.00, 0.67 };
GLclampf DefBg[4]         = { 0.00, 0.00, 0.00, 1.00 };
GLclampf DefUnibondCol[4] = { 1.00, 1.00, 1.00, 1.00 };

GLfloat DefAtCol_Ambient_by_Diffuse = 1.0;
