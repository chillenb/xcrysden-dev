/*
 
 This code was taken with permission from: 
   http://astronomy.swin.edu.au/~pbourke/opengl/redblue/.
   This web page is due to:
        P a u l   B o u r k e
        Astrophysics and Supercomputing
        e  pbourke@swin.edu.au
        w  astronomy.swin.edu.au/~pbourke
        p  +61 -3 9214 8624
        f  +61 -3 9214 8797
        l  Mail number 31
            PO Box 218 Hawthorn
            Swinburne University of Technology
            Victoria 3122, Australia

 The code was modifies by:
 
 *****************************************************************************
   Modified by Daniel van Vugt (daniel@computing.edu.au)
   I have removed all accumulation buffer code, and without needing blending.
   This works because when you're writing to different colour planes, it's
   effectively the same as having two separate framebuffers. In fact, blending
   would just slow it all down and add some ugly bugs.
   My changes are marked with a "DVV".

   -----

   Anaglyphs implemented in OpenGL without stereo buffer support
   Adds the two images in the accumulation buffer
   Note
   1. All objects must be drawn as greyscale!
   2. This is written for illustrative purposes...not efficiency!

 * ------------------------------------------------------------------------- *

 * For the purpose of XCRYSDEN modified by:                                  *
 * Eric Verfaillie ericverfaillie@yahoo.fr EV                                *
 * May 2004                                                                  *

 *****************************************************************************

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FALSE 0
#define TRUE (!FALSE)


typedef struct {
   unsigned char r,g,b,a;
} PIXELA;
typedef struct {
   double r,g,b;
} COLOUR;


#define SIGN(x) (x < 0 ? (-1) : 1)
#define MODULUS(p) (sqrt(p.x*p.x + p.y*p.y + p.z*p.z))
#define CROSSPROD(p1,p2,p3) \
   p3.x = p1.y*p2.z - p1.z*p2.y; \
   p3.y = p1.z*p2.x - p1.x*p2.z; \
   p3.z = p1.x*p2.y - p1.y*p2.x

typedef struct {
   XYZ vp;              /* View position           */
   XYZ vd;              /* View direction vector   */
   XYZ vu;              /* View up direction       */
   XYZ pr;              /* Point to rotate about   */
   double focallength;  /* Focal Length along vd   */
   double aperture;     /* Camera aperture         */
   double eyesep;       /* Eye separation          */
   int screenheight,screenwidth;
} CAMERA;

/*void HandleDisplay(struct Togl *togl);
void CameraHome();
void Normalise(XYZ *);
XYZ  CalcNormal(XYZ,XYZ,XYZ);
*/

#define DTOR            0.0174532925
#define RTOD            57.2957795
#define TWOPI           6.283185307179586476925287

#define PID2            1.570796326794896619231322
#define ESC 27
