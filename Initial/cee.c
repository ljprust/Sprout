
#include "../defs.h"
#include <stdbool.h>
#include <stdio.h>

static double Eej    = 0.0;
static double Mej    = 0.0;
static double Msun   = 0.0;
static double vmax   = 0.0;
static double rhoISM = 0.0;
static double t0     = 0.0;
static double yr     = 0.0;
static double day    = 0.0;
static double Lz     = 0.0;
static double tinput = 0.0;
static double tdelay = 0.0;
static double Rdisk  = 0.0;
static double Rsun   = 0.0;
static double ramPressureFactor = 0.0;
static double Rgas   = 0.0;
static bool   quadrant = false;

void setICParams( struct domain * theDomain ){
   // constants
   yr     = 365.25*24.0*3600.0; // sec
   day    = 24.0*3600.0; // sec
   Rsun   = 7.0e10; // cm
   Msun   = 2.0e33; // g
   Rgas   = 8.314e7; // cgs

   // domain size
   Lz     = theDomain->theParList.Lz;

   // ejecta parameters
   t0     = theDomain->theParList.t_min;
   vmax   = theDomain->theParList.v_max;
   ramPressureFactor = theDomain->theParList.Ram_Pressure_Factor;
   Eej    = 1.0e51;
   Mej    = 1.0*Msun;

   // time of CEE data input
   tinput = 500.0*day;
   tdelay = 1000.0*yr;
   Rdisk  = 500.0*Rsun;

   // CSM parameters
   rhoISM = theDomain->theParList.rho_ISM;

   // model a quadrant of the cube or just an octant
   quadrant = false;
}

void initial( double * prim , double * xi , double t , bool debug ){

   double x, y, z, r;
   double vx, vy, vz;
   bool isEjecta = false;
   double scaleFactor, rhoCEE, diskHeight, rhoSunny, v0sq, temp;

   // scale CEE rho to current time
   scaleFactor = 1.0; // pow((tdelay-tinput)/tinput,3.0)

   x = xi[0];
   y = xi[1];
   z = xi[2];
   r = sqrt(x*x+y*y+z*z);

   // shift z values if we're doing a quadrant
   if( quadrant ) z = z - Lz/2.0;

   // determine fluid velocities
   vx = x/t0;
   vy = y/t0;
   vz = z/t0;

   diskHeight = Rsun*(95.0*log10(r/Rsun)-125.0);
   rhoCEE = 0.01*pow(r/10.0/Rsun,-4.0)*pow(1.0+pow(125.0*Rsun/r,3.5),-1.05)
          * exp(-x*x/2.0/diskHeight/diskHeight);

   v0sq = 4.0/3.0*Eej/Mej;
   rhoSunny = pow(3.0/4.0/3.14159, 1.5) * pow(Mej, 2.5)/pow(Eej, 1.5) /t0/t0/t0 * exp(-r*r/t0/t0/v0sq);

   temp = 4.5e4/(r/100.0/Rsun);

   // define primitives
   if( r < vmax*t0 ) { // ejecta
      prim[RHO] = rhoSunny;
      prim[UU1] = vx;
      prim[UU2] = vy;
      prim[UU3] = vz;
      prim[PPP] = 0.7e14*pow(rhoSunny,1.666666667);
      prim[XXX] = 1.0; // tracks ejecta fraction
   } else if ( r < Rdisk ) { // disk
      prim[RHO] = rhoCEE;
      prim[UU1] = 0.0;
      prim[UU2] = 0.0;
      prim[UU3] = 0.0;
      prim[PPP] = rhoCEE*Rgas*temp;
      prim[XXX] = 0.0;
   } else if (rhoCEE*scaleFactor>rhoISM) { // CEE outflow
      prim[RHO] = rhoCEE/scaleFactor;
      prim[UU1] = 0.0;
      prim[UU2] = 0.0;
      prim[UU3] = 0.0;
      prim[PPP] = rhoCEE*Rgas*temp;
      prim[XXX] = 0.0;
   } else { // ISM
      prim[RHO] = rhoISM;
      prim[UU1] = 0.0;
      prim[UU2] = 0.0;
      prim[UU3] = 0.0;
      prim[PPP] = rhoISM*Rgas*100.0;
      prim[XXX] = 0.0;
   }

   // set pressure to small fraction of ram pressure
   //prim[PPP] = ramPressureFactor*vmax*vmax*prim[RHO];
}

