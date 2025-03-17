
#include "../defs.h"
#include <stdbool.h>
#include <stdio.h>

static double vmax   = 0.0;
static double Eej    = 0.0;
static double Mej    = 0.0;
static double vwind  = 0.0;
static double Mdot   = 0.0;
static double t0     = 0.0;
static double Msun   = 0.0;
static double yr     = 0.0;
static double day    = 0.0;
static double Lz     = 0.0;
static int    quadrant = 0;
static double nPower = 0.0;
static double deltaPower = 0.0;
static double ramPressureFactor = 0.0;

void setICParams( struct domain * theDomain ){
   // constants
   Msun   = 2.0e33;
   yr     = 365.25*24.0*3600.0; // sec
   day    = 24.0*3600.0;
   Lz     = theDomain->theParList.Lz;

   // ejecta parameters
   Eej    = theDomain->theParList.E_ejecta;
   Mej    = Msun * theDomain->theParList.M_ejecta;
   t0     = theDomain->theParList.t_min;
   vmax   = theDomain->theParList.v_max;
   ramPressureFactor = theDomain->theParList.Ram_Pressure_Factor;

   // power laws
   deltaPower = theDomain->theParList.delta_power;
   nPower     = theDomain->theParList.n_power;

   // CSM parameters
   vwind  = 1.0e5 * theDomain->theParList.v_wind;
   Mdot   = Msun/yr * theDomain->theParList.Mdot_wind;

   // model a quadrant of the cube or just an octant
   quadrant = theDomain->theParList.useQuadrant;
}

void initial( double * prim , double * xi , double t , bool debug ){

   double x, y, z, r, r0;
   double K, vt, rt, rhoprefactor, rhoOut, rhoIn;

   x = xi[0];
   y = xi[1];
   z = xi[2];

   if( quadrant>0 ) {
      z = z - Lz/2.0;
   }

   r = sqrt( x*x + y*y + z*z );
   r0 = vmax*t0;

   K = (nPower-3.0)*(3.0-deltaPower)/4.0/3.14159/(nPower-deltaPower);
   vt = sqrt((nPower-5.0)*(5.0-deltaPower)/(nPower-3.0)/(3.0-deltaPower)*2.0*Eej/Mej);
   rt = vt*t0;
   rhoprefactor = K*Mej/rt/rt/rt;
   rhoOut = rhoprefactor*pow(r/rt,-nPower);
   rhoIn  = rhoprefactor*pow(r/rt,-deltaPower);

   if(r <= rt) { // delta power law
      prim[RHO] = rhoIn;
      prim[UU1] = x/r0 * vmax;
      prim[UU2] = y/r0 * vmax;
      prim[UU3] = z/r0 * vmax;
      prim[XXX] = 1.0;
   } else if( r <= r0 ) { // n power law
      prim[RHO] = rhoOut;
      prim[UU1] = x/r0 * vmax;
      prim[UU2] = y/r0 * vmax;
      prim[UU3] = z/r0 * vmax;
      prim[XXX] = 1.0;
   } else { // ISM
      prim[RHO] = Mdot/4.0/3.14159/r/r/vwind;
      prim[UU1] = 0.0;
      prim[UU2] = 0.0;
      prim[UU3] = 0.0;
      prim[XXX] = 0.0;
   }

   prim[PPP] = ramPressureFactor*vmax*vmax*prim[RHO];
}
