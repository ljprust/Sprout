
#include "../defs.h"
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

static double vmax     = 0.0;
static double rhoISM   = 0.0;
static double rhoNormFloor = 0.0;
static double Eej      = 0.0;
static double Mej      = 0.0;
static double t0       = 0.0;
static double tfit     = 0.0;
static double Msun     = 0.0;
static double yr       = 0.0;
static double day      = 0.0;
static double AU       = 0.0;
static double Lz       = 0.0;
static bool   quadrant = false;
//static double thetamin = 0.0;
static double c1param1 = 0.0;
static double c1param2 = 0.0;
static double c1param3 = 0.0;
static double c1param4 = 0.0;
static double c1param5 = 0.0;
static double c2param1 = 0.0;
static double c2param2 = 0.0;
static double c3param1 = 0.0;
static double c3param2 = 0.0;
static double c4param1 = 0.0;
static double c4param2 = 0.0;
static double c4param3 = 0.0;
static double c4param4 = 0.0;
static double c4param5 = 0.0;
static double c5param1 = 0.0;
static double c5param2 = 0.0;
static double bowparam1 = 0.0;
static double bowparam2 = 0.0;
static double peakparam1 = 0.0;
static double peakparam2 = 0.0;
static double recompparam1 = 0.0;
static double recompparam2 = 0.0;
static double recomppeakparam1 = 0.0;
static double recomppeakparam2 = 0.0;
static double theta1 = 0.0;
static double theta2 = 0.0;
static double theta3 = 0.0;
static double theta4 = 0.0;
static double r1 = 0.0;
static double r2 = 0.0;
static double scalebowshock = 0.0;
static double Pratio = 0.0;
static bool   removeTransitionRegion = false;

double outerfit( double th, double thpeak, double thbow, double rad );
double innerfit( double th, double threcomppeak, double rad );
double transfit( double th, double thetaL, double thetaR, double rhoL, double rhoR );
double rmax( double theta, double rad0 );
double c1fit( double radius );
double c2fit( double radius );
double c3fit( double radius );
double c4fit( double radius );
double c5fit( double radius );
double thetabowfit( double radius );
double thetapeakfit( double radius );
double thetarecompfit( double radius );
double thetarecomppeakfit( double radius );

void setICParams( struct domain * theDomain ){
   // constants
   Msun   = 2.0e33;
   yr     = 365.25*24.0*3600.0; // sec
   day    = 24.0*3600.0;
   AU     = 1.496e13;
   Lz     = theDomain->theParList.Lz;

   // ejecta parameters
   Eej    = 0.97e51;
   Mej    = 1.789623e33;
   t0     = 10.0*yr;
   tfit   = 10000.0;
   vmax   = 2.318e9; // CHANGE THIS!!!
   //thetamin = 0.01;
   rhoNormFloor = 0.01;

   // ISM parameters
   rhoISM = 6.31e-25;
   Pratio = 1.0e-5;

   // model a quadrant of the cube or just an octant
   quadrant = true;

   // fit parameters

   removeTransitionRegion = false;

   // same as gabe
   c1param1 = -0.46218021;
   c1param2 = -3.2527453;
   c1param3 = 2.22623454;
   c1param4 = 0.31;
   c1param5 = -0.0600881;
/*
   // sqrt3floor
   c1param1 = 0.85005667;
   c1param2 = -7.16547361;
   c1param3 = 8.2780451;
   c1param4 = -4.03173747;
   c1param5 = 0.9874284;

   // times3floor
   c1param1 = -0.22036875;
   c1param2 = -2.05216039;
   c1param3 = 1.5731635;
   c1param4 = 0.0;
   c1param5 = 0.0;
*/
   c2param1 = -0.50302244;
   c2param2 = 15.95891016;
/*
   // times3floor
   c2param1 = -0.31235143;
   c2param2 = 10.49238271;
*/
   c3param1 = -0.06285498;
   c3param2 = 1.93430793;
/*
   // times3floor
   c3param1 = -0.02075052;
   c3param2 = 0.74863516;
*/
   c4param1 = 1.03125262;
   c4param2 = -3.59159214;
   c4param3 = 5.7664357;
   c4param4 = -4.43428767;
   c4param5 = 1.41344816;

   c5param1 = -0.25425228;
   c5param2 = 8.36636415;

   scalebowshock = 1.0; // same as gabe
   //scalebowshock = 0.73; // sqrt3floor
   //scalebowshock = 0.55; // times3floor

   bowparam1 = -0.29336977;
   bowparam2 = 8.56898219;

   peakparam1 = -0.30071382;
   peakparam2 = 8.75640844;

   recompparam1 = -0.19975411;
   recompparam2 = 4.84977536;

   recomppeakparam1 = -0.24091959;
   recomppeakparam2 = 6.01677967;

   // rmax fit

   // same as gabe
   theta1 = 0.20;
   theta2 = 0.35;
   theta3 = 0.43;
   theta4 = 0.73;
   r1     = 1.72;
   r2     = 1.25;
/*
   // sqrt3floor
   theta1 = 0.11;
   theta2 = 0.17;
   theta3 = 0.32;
   theta4 = 0.45;
   r1     = 1.58;
   r2     = 1.30;

   // times3floor
   theta1 = 0.10;
   theta2 = 0.16;
   theta3 = 0.29;
   theta4 = 0.36;
   r1     = 1.29;
   r2     = 1.19;
*/
}

double c1fit( double radius ) {
  double rad = radius/AU;
  double c1 = c1param1 + c1param2*rad + c1param3*rad*rad
            + c1param4*rad*rad*rad + c1param5*rad*rad*rad*rad;
  return c1;
}

double c2fit( double radius ) {
  double c2 = c2param1*log(radius)+c2param2;
  return c2;
}

double c3fit( double radius ) {
  double c3 = c3param1*log(radius)+c3param2;
  return c3;
}

double c4fit( double radius ) {
  double rad = radius/AU;
  double c4 = c4param1 + c4param2*rad + c4param3*rad*rad
            + c4param4*rad*rad*rad + c4param5*rad*rad*rad*rad;
  return c4;
}

double c5fit( double radius ) {
  double c5 = exp(c5param1*log(radius)+c5param2);
  return c5;
}

double thetabowfit( double radius ) {
  double th = scalebowshock * exp(bowparam1*log(radius)+bowparam2);
  return th;
}

double thetapeakfit( double radius ) {
  double th = scalebowshock * exp(peakparam1*log(radius)+peakparam2);
  return th;
}

double thetarecompfit( double radius ) {
  double th = exp(recompparam1*log(radius)+recompparam2);
  return th;
}

double thetarecomppeakfit( double radius ) {
  double th = exp(recomppeakparam1*log(radius)+recomppeakparam2);
  return th;
}

double outerfit( double th, double thpeak, double thbow, double rad ) {
   double c4, c5, result;
   c4 = c4fit(rad);
   c5 = c5fit(rad);
   result = c4 + 20.0 * sin( th/thbow*1.015*3.14159 ) * pow( th/thpeak, 16.0*c5 )
          + c5 * pow( th/thpeak, 2.75*c5 );
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
   double maxr, slope, intercept;

   if ( theta<theta1 ) {
      maxr = r1;
   } else if ( theta<theta2 ) {
      slope     = (r2-r1)/(theta2-theta1);
      intercept = r1 - slope*theta1;
      maxr      = slope*theta + intercept;
   } else if ( theta<theta3 ) {
      maxr = r2;
   } else if ( theta<theta4 ) {
      slope     = (1.0-r2)/(theta4-theta3);
      intercept = r2 - slope*theta3;
      maxr      = slope*theta + intercept;
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
   double r0, r_max, rfit;
   double thetabow, thetapeak, thetarecomp, thetarecomppeak;

   x = xi[0];
   y = xi[1];
   z = xi[2];
   if( quadrant ) z = z - Lz/2.0;

   r = sqrt( x*x + y*y + z*z );
   r0 = vmax*t0;
   rfit = r*tfit/t0;
   theta = acos(z/r);
   thetaDeg = theta*180.0/3.14159;
   r_max = rmax(theta, r0);
   //theta = MAX(theta, thetamin);

   thetabow        = thetabowfit(rfit);
   thetapeak       = thetapeakfit(rfit);
   thetarecomp     = thetarecompfit(rfit);
   thetarecomppeak = thetarecomppeakfit(rfit);

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
      rhoNorm = outerfit(theta, thetapeak, thetabow, rfit);
      rhoNorm = MAX(rhoNorm,rhoNormFloor);
      prim[RHO] = rhoNorm*rhoSunny;
   } else if ( removeTransitionRegion && theta >= thetarecomppeak ) {
      // remove transition region, fill instead with outer fit
      rhoNorm = outerfit(theta, thetapeak, thetabow, rfit);
      rhoNorm = MAX(rhoNorm,rhoNormFloor);
      prim[RHO] = rhoNorm*rhoSunny;
   } else if ( theta <= thetarecomppeak ) { // within recompression
      rhoNorm = innerfit(theta, thetarecomppeak, rfit);
      rhoNorm = MAX(rhoNorm,rhoNormFloor);
      prim[RHO] = rhoNorm*rhoSunny;
   } else { // transition region
      double rhoL, rhoR;
      rhoL = innerfit(thetarecomppeak, thetarecomppeak, rfit);
      rhoR = outerfit(thetarecomp, thetapeak, thetabow, rfit);
      rhoNorm = transfit(theta, thetarecomppeak, thetarecomp, rhoL, rhoR);
      rhoNorm = MAX(rhoNorm,rhoNormFloor);
      prim[RHO] = rhoNorm*rhoSunny;
   }

   if( thetaDeg < 30.0 && vr < 2.0e8 ) {
      prim[YYY] = 1.0;
   } else {
      prim[YYY] = 0.0;
   }

   prim[PPP] = Pratio * vmax*vmax*prim[RHO];
}
