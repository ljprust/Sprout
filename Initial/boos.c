
#include "../defs.h"
#include <stdbool.h>
#include <stdio.h>

#define XINPUT 90
#define ZINPUT 180

// This initial conditions file is designed to import the ejecta from Boos et al. 2024 into Sprout.

static double vmax   = 0.0;
static double rhoISM = 0.0;
static double Eej    = 0.0;
static double Mej    = 0.0;
static double vwind  = 0.0;
static double Mdot   = 0.0;
static double t0     = 0.0;
static double tinput = 0.0;
static double Msun   = 0.0;
static double yr     = 0.0;
static double day    = 0.0;
static double Lz     = 0.0;
static double ramPressureFactor = 0.0;
static bool   quadrant = false;
static double vxInput[XINPUT][ZINPUT];
static double vzInput[XINPUT][ZINPUT];
static double rhoInput[XINPUT][ZINPUT];

void setICParams( struct domain * theDomain ){
   // constants
   Msun   = 2.0e33;
   yr     = 365.25*24.0*3600.0; // sec
   day    = 24.0*3600.0;
   Lz     = theDomain->theParList.Lz;

   // ejecta parameters
   Eej    = 1.0e51;
   Mej    = 1.0*Msun;
   t0     = 10.0*yr;
   vmax   = 3.0e9; // 2.53e9;
   ramPressureFactor = 1.0e-5;

   // CSM parameters
   vwind  = 10.0e5;
   Mdot   = 0.0; // 4.0e-5*Msun/yr;
   rhoISM = 6.31e-25; // 1.7e-24; (constant)

   // epoch of intial data
   tinput = 53.059; // sec

   // model a quadrant of the cube or just an octant
   quadrant = true;

   ////// READ DATA FROM TXT FILES //////

   FILE *vxInputFile;
   FILE *vzInputFile;
   FILE *rhoInputFile;

   char filename_vx[256];
   char filename_vz[256];
   char filename_rho[256];

   sprintf(filename_vx,  "sproutinput_vx.txt");
   sprintf(filename_vz,  "sproutinput_vz.txt");
   sprintf(filename_rho, "sproutinput_rho.txt");

   vxInputFile = fopen(filename_vx,"r");
   printf("opened %s\n", filename_vx);
   vzInputFile = fopen(filename_vz,"r");
   printf("opened %s\n", filename_vz);
   rhoInputFile = fopen(filename_rho,"r");
   printf("opened %s\n", filename_rho);

   int i, j;
   for( i=0 ; i<XINPUT ; ++i ){
      for( j=0 ; j<ZINPUT ; ++j ){
         fscanf( vxInputFile,  "%lf", &vxInput[i][j]  );
         fscanf( vzInputFile,  "%lf", &vzInput[i][j]  );
         fscanf( rhoInputFile, "%lf", &rhoInput[i][j] );
      }
   }

   fclose(vxInputFile);
   fclose(vzInputFile);
   fclose(rhoInputFile);
   printf("Done reading input data, closed files\n");
}

void initial( double * prim , double * xi , double t , bool debug ){

   double x, y, z;
   double vx, vy, vz, vr, vcyl;
   int i, j, minIndex_i, minIndex_j;
   double dist2;
   double minDist2 = 1.0e100;
   double rhoRead;
   bool inEjecta;

   x = xi[0];
   y = xi[1];
   z = xi[2];

   // shift z values if we're doing a quadrant
   if( quadrant ) {
      z = z - Lz/2.0;
   }

   // determine fluid velocities
   vx = x/t0;
   vy = y/t0;
   vz = z/t0;
   vr = sqrt(vx*vx+vy*vy+vz*vz);
   vcyl = sqrt(vx*vx+vy*vy); // velocity in x-y plane

   // ignore initial data with v > vmax
   inEjecta = vr < vmax;

   // do a nearest neighbor search if we're within the ejecta
   if( inEjecta ) {
      for( i=0 ; i<XINPUT ; ++i ){
         for( j=0 ; j<ZINPUT ; ++j ){
            dist2 = ( vcyl - vxInput[i][j] ) * ( vcyl - vxInput[i][j] )
                  + ( vz   - vzInput[i][j] ) * ( vz   - vzInput[i][j] );
            if(dist2 < minDist2) {
               minDist2 = dist2;
               minIndex_i = i;
               minIndex_j = j;
            }
         }
      }
   }

   // set density to value of nearest neighbor
   rhoRead = rhoInput[minIndex_i][minIndex_j];

   // various debug messages
   if ( debug ) {
      printf("found neighbor with indices %d %d\n",minIndex_i,minIndex_j);
      printf("read rho %5.3e\n",rhoRead);
      printf("x y z minDist2 %5.3e %5.3e %5.3e %5.3e\n",x,y,z,minDist2);
   }

   // define primitives
   if( inEjecta ) {
      prim[RHO] = rhoRead;
      prim[UU1] = vx;
      prim[UU2] = vy;
      prim[UU3] = vz;
      prim[XXX] = 1.0; // tracks ejecta fraction
   } else {
      prim[RHO] = rhoISM;
      prim[UU1] = 0.0;
      prim[UU2] = 0.0;
      prim[UU3] = 0.0;
      prim[XXX] = 0.0;
   }

   prim[YYY] = 0.0; // what should we use this for?

   // set pressure to small fraction of ram pressure
   prim[PPP] = ramPressureFactor*vmax*vmax*prim[RHO];
}

