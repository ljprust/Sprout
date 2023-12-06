
#include "../defs.h"

static double vmax   = 0.0;
static double rhoISM = 0.0;
static double Eej    = 0.0;
static double Mej    = 0.0;
static double vwind  = 0.0;
static double Mdot   = 0.0;
static double t0     = 0.0;
static double Msun   = 0.0;
static double yr     = 0.0;
static double day    = 0.0;

void setICParams( struct domain * theDomain ){
   Msun   = 2.0e33;
   yr     = 365.25*24.0*3600.0; // sec
   day    = 24.0*3600.0;

   Eej    = 1.31e51; // 1.0e51;
   Mej    = 2.25*Msun;
   t0     = 50.0*day;
   vmax   = 2.5926e9; // 1.72e9;
   vwind  = 10.0e5;
   Mdot   = 4.0e-5*Msun/yr;
   rhoISM = 5.0e-25; // 1.6e-24;
}

void initial( double * prim , double * xi , double t ){

   double x, y, z, r;
   double v0, r0, vr, rhoSunny, pres;

   x = xi[0];
   y = xi[1];
   z = xi[2]*0.0;
   r = sqrt( x*x + y*y + z*z );

   r0 = vmax*t0;
   v0 = sqrt(4.0/3.0*Eej/Mej);
   vr = vmax*r/r0;
   rhoSunny = pow(3.0/4.0/3.14159, 1.5) * pow(Mej, 2.5)/pow(Eej, 1.5) /t0/t0/t0 * exp(-vr*vr/v0/v0);
   
   if( r <= r0 ) {
      prim[RHO] = rhoSunny;
      prim[UU1] = x/r0 * vmax;
      prim[UU2] = y/r0 * vmax;
      prim[UU3] = z/r0 * vmax;
      prim[XXX] = 1.0;
   } else {
      prim[RHO] = rhoISM;
      prim[UU1] = 0.0;
      prim[UU2] = 0.0;
      prim[UU3] = 0.0;
      prim[XXX] = 0.0;
   }
   prim[PPP] = 1e-5*vmax*vmax*prim[RHO];
}
