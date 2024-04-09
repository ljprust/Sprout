
#include "../defs.h"
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

static double vmax     = 0.0;
static double rhoISM   = 0.0;
static double Eej      = 0.0;
static double Mej      = 0.0;
static double vwind    = 0.0;
static double Mdot     = 0.0;
static double t0       = 0.0;
static double Msun     = 0.0;
static double yr       = 0.0;
static double day      = 0.0;
static double Lz       = 0.0;
static bool   quadrant = false;
static double thetamin = 0.0;

double outerfit( double th, double thpeak, double thbow, double rad );
double innerfit( double th, double threcomppeak, double rad );
double transfit( double th, double thetaL, double thetaR, double rhoL, double rhoR );
double rmax( double theta, double rad0 );
double afit( double radius );
double dfit( double radius );
double c1fit( double radius );
double c2fit( double radius );
double c3fit( double radius );
double thetabowfit( double radius );
double thetapeakfit( double radius );
double thetarecompfit( double radius );
double thetarecomppeakfit( double radius );

void setICParams( struct domain * theDomain ){
   // constants
   Msun   = 2.0e33;
   yr     = 365.25*24.0*3600.0; // sec
   day    = 24.0*3600.0;
   Lz     = theDomain->theParList.Lz;

   // ejecta parameters
   Eej    = 0.97e51; // 1.0e51;
   Mej    = 1.789623e33; // 1.0*Msun;
   t0     = 10.0*yr; // 100.0*day; // r0 = 1.728e16 cm
   vmax   = 2.53e9; // 2.0e9;
   thetamin = 0.01;

   // CSM parameters
   vwind  = 10.0e5;
   Mdot   = 4.0e-5*Msun/yr;
   rhoISM = 6.31e-25; // 1.7e-24; (constant)

   // model a quadrant of the cube or just an octant
   quadrant = true;

   double afit( double radius );
}

double afit( double radius ) {
  double AU = 1.496e13;
  double rad = radius/AU;
  double a = 1.03125262-3.59159214*rad+5.7664357*rad*rad-4.43428767*rad*rad*rad+1.41344816*rad*rad*rad*rad;
  return a;
}

double dfit( double radius ) {
  double d = exp(-0.25425228*log(radius)+8.36636415);
  return d;
}

double c1fit( double radius ) {
  double AU = 1.496e13;
  double rad = radius/AU;
  double c1 = -0.46218021-3.2527453*rad+2.22623454*rad*rad+0.31*rad*rad*rad-0.0600881*rad*rad*rad*rad;
  return c1;
}

double c2fit( double radius ) {
  double c2 = -0.50302244*log(radius)+15.95891016;
  return c2;
}

double c3fit( double radius ) {
  double c3 = -0.06285498*log(radius)+1.93430793;
  return c3;
}

double thetabowfit( double radius ) {
  double th = exp(-0.29336977*log(radius)+8.56898219);
  return th;
}

double thetapeakfit( double radius ) {
  double th = exp(-0.30071382*log(radius)+8.75640844);
  return th;
}

double thetarecompfit( double radius ) {
  double th = exp(-0.19975411*log(radius)+4.84977536);
  return th;
}

double thetarecomppeakfit( double radius ) {
  double th = exp(-0.24091959*log(radius)+6.01677967);
  return th;
}

double outerfit( double th, double thpeak, double thbow, double rad ) {
   double a, d, result;
   a = afit(rad);
   d = dfit(rad);
   result = a + 20.0 * sin( th/thbow*1.015*3.14159 ) * pow( th/thpeak, 16.0*d )
          + d * pow( th/thpeak, 2.75*d );
   if (th > thpeak) {
      result = MAX(result, 1.0);
   }
   return result;
}

double innerfit( double th, double threcomppeak, double rad ) {
   double c1, c2, c3, y;
   c1 = c1fit(rad);
   c2 = c2fit(rad);
   c3 = c3fit(rad);
   y = log(th/threcomppeak);
   return exp( c1 + c2*y + c3*y*y );
}

double transfit( double th, double thetaL, double thetaR, double rhoL, double rhoR ) {
   double rhomax, rhomin, tanharg;
   rhomax = MAX(rhoL, rhoR);
   rhomin = MIN(rhoL, rhoR);
   tanharg = -3.0*((th-thetaL)/(thetaR-thetaL)*2.0-1.0);
   return (tanh(tanharg) + 1.0) * 0.5 * (rhomax-rhomin) + rhomin;
}

double rmax( double theta, double rad0 ) {
   double theta1, theta2, theta3, rad1, rad2, rad3, maxr;
   theta1 = 0.205;
   theta2 = 0.28;
   theta3 = 0.6;
   rad1 = 1.44516;
   rad2 = 1.5484;
   rad3 = 1.2258;

   if ( theta<theta1 ) {
      maxr = rad1 + 2.4563*theta*theta;
   } else if ( theta<theta2 ) {
      maxr = -4.30108*theta + 2.430108;
   } else if ( theta<theta3 ) {
      maxr = -2.20514*(theta-theta2)*(theta-theta2) + rad3;
   } else {
      maxr = 1.0;
   }
   maxr *= rad0;
   return maxr;
}

void initial( double * prim , double * xi , double t , bool debug ){

   double x, y, z, r;
   double v0, vr, rhoSunny, rhoNorm;
   double thetaDeg, theta;
   double r0, r_max;
   double thetabow, thetapeak, thetarecomp, thetarecomppeak;

   r0 = 0.0;

   x = xi[0];
   y = xi[1];
   z = xi[2];
   if( quadrant ) {
      z = z - Lz/2.0;
   }

   r = sqrt( x*x + y*y + z*z );
   theta = acos(z/r);
   thetaDeg = theta*180.0/3.14159;
   r_max = rmax(theta, r0);
   theta = MAX(theta, thetamin);

   thetabow        = thetabowfit(r);
   thetapeak       = thetapeakfit(r);
   thetarecomp     = thetarecompfit(r);
   thetarecomppeak = thetarecomppeakfit(r);

   v0 = sqrt(4.0/3.0*Eej/Mej);
   vr = r/t0;
   rhoSunny = pow(3.0/4.0/3.14159, 1.5) * pow(Mej, 2.5)/pow(Eej, 1.5) /t0/t0/t0 * exp(-vr*vr/v0/v0);

   if( r > r_max ) {
      prim[UU1] = 0.0;
      prim[UU2] = 0.0;
      prim[UU3] = 0.0;
      prim[XXX] = 0.0;
   } else {
      prim[UU1] = x/t0;
      prim[UU2] = y/t0;
      prim[UU3] = z/t0;
      prim[XXX] = 1.0;
   }

   if( r > r_max ) { // ambient medium
      prim[RHO] = rhoISM;
   } else if ( theta > 1.5707963267948966 ) { // unshocked
      prim[RHO] = rhoSunny;
   } else if ( theta >= thetabow ) { // also unshocked
      prim[RHO] = rhoSunny;
   } else if ( theta >= thetarecomp ) { // in rarefaction wave
      rhoNorm = outerfit(theta, thetapeak, thetabow, r);
      prim[RHO] = rhoNorm*rhoSunny;
   } else if ( theta <= thetarecomppeak ) { // within recompression
      rhoNorm = innerfit(theta, thetarecomppeak, r);
      prim[RHO] = rhoNorm*rhoSunny;
   } else { // transition region
      double rhoL, rhoR;
      rhoL = innerfit(thetarecomppeak, thetarecomppeak, r);
      rhoR = outerfit(thetarecomp, thetapeak, thetabow, r);
      rhoNorm = transfit(theta, thetarecomppeak, thetarecomp, rhoL, rhoR);
      prim[RHO] = rhoNorm*rhoSunny;
   }

   if( thetaDeg < 30.0 && vr < 2.0e8 ) {
      prim[YYY] = 1.0;
   } else {
      prim[YYY] = 0.0;
   }

   prim[PPP] = 1.0e-5*vmax*vmax*prim[RHO];
}

