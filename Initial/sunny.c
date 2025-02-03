
#include "../defs.h"
#include <stdbool.h>
#include <stdio.h>

#define NTHETAINPUT 224
#define NRADINPUT 1000

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
static double Lz     = 0.0;
static bool   readIC = false;
static bool   quadrant = false;
static double bowshockscale = 0.0;
static double xInput[NTHETAINPUT][NRADINPUT];
static double yInput[NTHETAINPUT][NRADINPUT];
static double rhoInput[NTHETAINPUT][NRADINPUT];
static double vrInput[NTHETAINPUT][NRADINPUT];
static double tracerInput[NTHETAINPUT][NRADINPUT];

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

   // CSM parameters
   vwind  = 10.0e5;
   Mdot   = 4.0e-5*Msun/yr;
   rhoISM = 6.31e-25; // 1.7e-24; (constant)

   // Kasen fit parameters
   kasen  = true;
   bowshockscale = 0.55;
   fh     = 0.1;
   mpower = 8.0;
   thetah = 30.0*bowshockscale;
   thetap = 15.0;
   kasenA = 1.8;

   // read ICs from file
   readIC = false;

   // model a quadrant of the cube or just an octant
   quadrant = true;

   if (readIC) {
      FILE *xInputFile;
      FILE *yInputFile;
      FILE *rhoInputFile;
      FILE *vrInputFile;
      FILE *tracerInputFile;

      char filename_x[256];
      char filename_y[256];
      char filename_rho[256];
      char filename_vr[256];
      char filename_tracer[256];

      sprintf(filename_x,      "sproutinput_x.txt");
      sprintf(filename_y,      "sproutinput_y.txt");
      sprintf(filename_rho,    "sproutinput_rho.txt");
      sprintf(filename_vr,     "sproutinput_vr.txt");
      sprintf(filename_tracer, "sproutinput_tracer.txt");

      xInputFile = fopen(filename_x,"r");
      printf("opened %s\n", filename_x);
      yInputFile = fopen(filename_y,"r");
      printf("opened %s\n", filename_y);
      rhoInputFile = fopen(filename_rho,"r");
      printf("opened %s\n", filename_rho);
      vrInputFile = fopen(filename_vr,"r");
      printf("opened %s\n", filename_vr);
      tracerInputFile = fopen(filename_tracer,"r");
      printf("opened %s\n", filename_tracer);

      int i, j;
      for( i=0 ; i<NTHETAINPUT ; ++i ){
         for( j=0 ; j<NRADINPUT ; ++j ){
            fscanf( xInputFile,      "%lf", &xInput[i][j]      );
            fscanf( yInputFile,      "%lf", &yInput[i][j]      );
            fscanf( rhoInputFile,    "%lf", &rhoInput[i][j]    );
            fscanf( vrInputFile,     "%lf", &vrInput[i][j]     );
            fscanf( tracerInputFile, "%lf", &tracerInput[i][j] );
         }
      }

      fclose(xInputFile);
      fclose(yInputFile);
      fclose(rhoInputFile);
      fclose(vrInputFile);
      fclose(tracerInputFile);
      printf("done reading input data, closed files\n");
   }
}

void initial( double * prim , double * xi , double t , bool debug ){

   double x, y, z, r;
   double v0, r0, vr, rhoSunny;
   double kasenFactor, thetaDeg, theta;
   double xForReading, yForReading;
   int i, j, minIndex_i, minIndex_j;
   double dist2;
   double minDist2 = 1.0e100;
   double rhoRead, vrRead, tracerRead;

   x = xi[0];
   y = xi[1];
   z = xi[2];

   if( quadrant ) {
      z = z - Lz/2.0;
   }

   r = sqrt( x*x + y*y + z*z );
   theta = acos(z/r);

   if (readIC) {
      xForReading = z;
      yForReading = sqrt( x*x + y*y );

      for( i=0 ; i<NTHETAINPUT ; ++i ){
         for( j=0 ; j<NRADINPUT ; ++j ){
            dist2 = ( xForReading - xInput[i][j] ) * ( xForReading - xInput[i][j] )
                  + ( yForReading - yInput[i][j] ) * ( yForReading - yInput[i][j] );
            if(dist2 < minDist2) {
               minDist2 = dist2;
               minIndex_i = i;
               minIndex_j = j;
            }
         }
      }

      //if (debug) printf("found neighbor with indices %d %d\n",minIndex_i,minIndex_j);

      rhoRead    = rhoInput[minIndex_i][minIndex_j];
      vrRead     = vrInput[minIndex_i][minIndex_j];
      tracerRead = tracerInput[minIndex_i][minIndex_j];
      //if (debug) printf("read rho %5.3e vr %5.3e tracer %5.3e\n",rhoRead,vrRead,tracerRead);
      //if (debug) printf("x y z minDist2 %5.3e %5.3e %5.3e %5.3e\n",x,y,z,minDist2);
      //if (debug) printf("xForReading yForReading xInput yInput %5.3e %5.3e %5.3e %5.3e\n", xForReading, yForReading, xInput[minIndex_i][minIndex_j], yInput[minIndex_i][minIndex_j]);
   }

   r0 = vmax*t0;
   v0 = sqrt(4.0/3.0*Eej/Mej);
   vr = vmax*r/r0;
   rhoSunny = pow(3.0/4.0/3.14159, 1.5) * pow(Mej, 2.5)/pow(Eej, 1.5) /t0/t0/t0 * exp(-vr*vr/v0/v0);
   
   thetaDeg = theta*180.0/3.14159;
   kasenFactor = fh+(1.0-fh)*pow(thetaDeg/thetah,mpower)/(1.0+pow(thetaDeg/thetah,mpower)) * (1.0+kasenA*exp(-pow(thetaDeg/thetah-1.0,2.0)/pow(thetap/thetah,2.0)));

   if( readIC ) {
      prim[RHO] = rhoRead;
      prim[UU1] = vrRead * x/r;
      prim[UU2] = vrRead * y/r;
      prim[UU3] = vrRead * z/r;
      prim[XXX] = tracerRead;
   } else if( r <= r0 ) {
      prim[RHO] = rhoSunny;
      if( kasen ) {
         prim[RHO] *= kasenFactor;
      }
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

   if( thetaDeg < 30.0 && vr < 2.0e8 ) {
      prim[YYY] = 1.0;
   } else {
      prim[YYY] = 0.0;
   }

   prim[PPP] = 1.0e-5*vmax*vmax*prim[RHO];
}
