
#include "../defs.h"
#include <stdbool.h>

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
static bool   kasen  = false;
static double fh     = 0.0;
static double mpower = 0.0;
static double thetah = 0.0;
static double thetap = 0.0;
static double kasenA = 0.0;

void setICParams( struct domain * theDomain ){
   Msun   = 2.0e33;
   yr     = 365.25*24.0*3600.0; // sec
   day    = 24.0*3600.0;

   Eej    = 1.0e51;
   Mej    = 1.0*Msun;
   t0     = 100.0*day; // r0 = 1.728e16 cm
   vmax   = 2.0e9;
   vwind  = 10.0e5;
   Mdot   = 4.0e-5*Msun/yr;
   rhoISM = 1.0e-20; // 1.7e-24;

   // Kasen fit params
   kasen  = true;
   fh     = 0.1;
   mpower = 8.0;
   thetah = 30.0;
   thetap = 15.0;
   kasenA = 1.8;
}

void initial( double * prim , double * xi , double t ){

   double x, y, z, r;
   double v0, r0, vr, rhoSunny;
   double kasenFactor, theta;

   x = xi[0];
   y = xi[1];
   z = xi[2];
   r = sqrt( x*x + y*y + z*z );

   r0 = vmax*t0;
   v0 = sqrt(4.0/3.0*Eej/Mej);
   vr = vmax*r/r0;
   rhoSunny = pow(3.0/4.0/3.14159, 1.5) * pow(Mej, 2.5)/pow(Eej, 1.5) /t0/t0/t0 * exp(-vr*vr/v0/v0);
   
   theta = acos(x/r)*180.0/3.14159;
   kasenFactor = fh+(1.0-fh)*pow(theta/thetah,mpower)/(1.0+pow(theta/thetah,mpower)) * (1.0+kasenA*exp(-pow(theta/thetah-1.0,2.0)/pow(thetap/thetah,2.0)));

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
   
   if( kasen ) {
      prim[RHO] *= kasenFactor;
   }

   prim[PPP] = 1.0e-6*prim[RHO]; // 1.0e-3*vmax*vmax*prim[RHO];
}
