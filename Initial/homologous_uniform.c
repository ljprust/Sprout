
#include "../defs.h"

static double x_zero = 0.0;
static double y_zero = 0.0;
static double z_zero = 0.0;
static double gam    = 0.0;
static double D      = 0.0;

void setICParams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
   x_zero = theDomain->theParList.Lx/2.*0.;
   y_zero = theDomain->theParList.Ly/2.*0.;
   z_zero = theDomain->theParList.Lz/2.*0.;
   if( theDomain->theParList.Num_x!=1 ) D += 1.;
   if( theDomain->theParList.Num_y!=1 ) D += 1.;
   if( theDomain->theParList.Num_z!=1 ) D += 1.;
}

void initial( double * prim , double * xi , double t ){
   
   double x = xi[0]-x_zero;
   double y = 0.;
   double z = 0.;
   if(D>1.) y = xi[1]-y_zero;
   if(D>2.) z = xi[2]-z_zero;

   prim[RHO] = 1.0;
   prim[PPP] = 1.0e-4/gam;
   prim[UU1] = x/t;
   prim[UU2] = y/t;
   prim[UU3] = z/t;
   prim[XXX] = 0.0;

}
