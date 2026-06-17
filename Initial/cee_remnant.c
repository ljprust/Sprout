
#include "../defs.h"
#include <stdbool.h>
#include <stdio.h>

#define NINPUT 819200

// This initial conditions file is designed to import the ejecta from Boos et al. 2024 into Sprout.

static double vmax   = 0.0;
static double rhoISM = 0.0;
static double t0     = 0.0;
static double tinput = 0.0;
static double yr     = 0.0;
static double Lx     = 0.0;
static double ramPressureFactor = 0.0;
static bool   quadrant = false;
static double vrInput[NINPUT];
static double rhoInput[NINPUT];
static double tracerInput[NINPUT];
static double xInput[NINPUT];
static double zInput[NINPUT];

void setICParams( struct domain * theDomain ){
   // constants
   yr     = 365.25*24.0*3600.0; // sec

   // domain size
   Lx     = theDomain->theParList.Lx;

   // ejecta parameters
   t0     = theDomain->theParList.t_min;
   vmax   = theDomain->theParList.v_max;
   ramPressureFactor = theDomain->theParList.Ram_Pressure_Factor;

   // CSM parameters
   rhoISM = theDomain->theParList.rho_ISM;

   // epoch of intial data
   tinput = theDomain->theParList.t_input;

   // model a quadrant of the cube or just an octant
   quadrant = true;

   ////// READ DATA FROM TXT FILES //////

   FILE *vrInputFile;
   FILE *rhoInputFile;
   FILE *tracerInputFile;
   FILE *xInputFile;
   FILE *zInputFile;

   char filename_vr[256];
   char filename_rho[256];
   char filename_tracer[256];
   char filename_x[256];
   char filename_z[256];

   sprintf(filename_vr,  "sproutinput_vr.txt");
   sprintf(filename_rho, "sproutinput_rho.txt");
   sprintf(filename_tracer, "sproutinput_tracer.txt");
   sprintf(filename_x, "sproutinput_x.txt");
   sprintf(filename_z, "sproutinput_z.txt");

   vrInputFile  = fopen(filename_vr, "r");
   rhoInputFile = fopen(filename_rho,"r");
   tracerInputFile = fopen(filename_tracer,"r");
   xInputFile = fopen(filename_x, "r");
   zInputFile = fopen(filename_z, "r");

   int i;
   for( i=0 ; i<NINPUT ; ++i ){
      fscanf( vrInputFile,  "%lf", &vrInput[i]  );
      fscanf( rhoInputFile, "%lf", &rhoInput[i] );
      fscanf( tracerInputFile, "%lf", &tracerInput[i] );
      fscanf( xInputFile, "%lf", &xInput[i] );
      fscanf( zInputFile, "%lf", &zInput[i] );
   }

   fclose(vrInputFile);
   fclose(rhoInputFile);
   fclose(tracerInputFile);
   fclose(xInputFile);
   fclose(zInputFile);

   printf("Done reading input data, closed files\n");
}

void initial( double * prim , double * xi , double t , bool debug ){

   double x, y, z, rcyl;
   double vx, vy, vz;
   int i;
   int minIndex = 0;
   double dist2;
   double minDist2 = 1.0e100;
   double rhoRead, vrRead, tracerRead;
   bool isEjecta = false;
   double scaleFactor;

   // scale rho to current time
   scaleFactor = tinput/t0;

   x = xi[0];
   y = xi[1];
   z = xi[2];
   rcyl = sqrt(y*y + z*z);

   // shift z values if we're doing a quadrant
   if( quadrant ) {
      x = x - Lx/2.0;
   }

   // determine fluid velocities
   vx = x/t0;
   vy = y/t0;
   vz = z/t0;

   for( i=0 ; i<NINPUT ; ++i ){
      dist2 = ( x    - xInput[i]*scaleFactor ) * ( x    - xInput[i]*scaleFactor )
            + ( rcyl - zInput[i]*scaleFactor ) * ( rcyl - zInput[i]*scaleFactor );
      if(dist2 < minDist2) {
         minDist2 = dist2;
         minIndex = i;
      }
   }

   // set density to value of nearest neighbor
   rhoRead = rhoInput[minIndex]*scaleFactor*scaleFactor*scaleFactor;

   // set abundance to that of nearest neighbor
   tracerRead = tracerInput[minIndex];

   // check if nearest neighbor is ejecta cell (vr > 1 km/s)
   vrRead = vrInput[minIndex];
   if( vrRead > 1.0e5 ) isEjecta = true;

   // various debug messages
   //if ( debug || false ) {
   //   printf("found neighbor with index %d\n",minIndex);
   //   printf("read rho %5.3e\n",rhoRead);
   //   printf("x y z minDist2 %5.3e %5.3e %5.3e %5.3e\n",x,y,z,minDist2);
   //}

   // make sure ejecta doesn't touch +z boundary
   //if( quadrant && vz > 4.5e9 ) isEjecta = false;

   // define primitives
   if( isEjecta ) {
      prim[RHO] = rhoRead;
      prim[UU1] = vx;
      prim[UU2] = vy;
      prim[UU3] = vz;
      prim[XXX] = tracerRead; // tracks ejecta fraction
   } else {
      prim[RHO] = rhoISM;
      prim[UU1] = 0.0;
      prim[UU2] = 0.0;
      prim[UU3] = 0.0;
      prim[XXX] = 0.0;
   }

   // set pressure to small fraction of ram pressure
   prim[PPP] = ramPressureFactor*vmax*vmax*prim[RHO];
}

