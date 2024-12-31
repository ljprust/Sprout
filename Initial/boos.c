
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
static double NInput[XINPUT][ZINPUT];
static double OInput[XINPUT][ZINPUT];
static double SiInput[XINPUT][ZINPUT];
static double FeInput[XINPUT][ZINPUT];

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
   vmax   = 5.0e9; // 2.53e9;
   ramPressureFactor = 1.0e-5;

   // CSM parameters
   vwind  = 10.0e5;
   Mdot   = 0.0; // 4.0e-5*Msun/yr;
   rhoISM = 6.31e-25; // 1.7e-24; (constant)

   // epoch of intial data
   tinput = 53.059; // m100_70
   //tinput = 66.74676; // m100_70_nosecdet

   // model a quadrant of the cube or just an octant
   quadrant = true;

   ////// READ DATA FROM TXT FILES //////

   FILE *vxInputFile;
   FILE *vzInputFile;
   FILE *rhoInputFile;
   FILE *NInputFile;
   FILE *OInputFile;
   FILE *SiInputFile;
   FILE *FeInputFile;

   char filename_vx[256];
   char filename_vz[256];
   char filename_rho[256];
   char filename_n[256];
   char filename_o[256];
   char filename_si[256];
   char filename_fe[256];

   sprintf(filename_vx,  "sproutinput_vx.txt");
   sprintf(filename_vz,  "sproutinput_vz.txt");
   sprintf(filename_rho, "sproutinput_rho.txt");
   sprintf(filename_n,   "sproutinput_n.txt");
   sprintf(filename_o,   "sproutinput_o.txt");
   sprintf(filename_si,  "sproutinput_si.txt");
   sprintf(filename_fe,  "sproutinput_fe.txt");

   vxInputFile  = fopen(filename_vx,"r");
   //printf("opened %s\n", filename_vx);
   vzInputFile  = fopen(filename_vz,"r");
   rhoInputFile = fopen(filename_rho,"r");
   NInputFile   = fopen(filename_n,"r");
   OInputFile   = fopen(filename_o,"r");
   SiInputFile  = fopen(filename_si,"r");
   FeInputFile  = fopen(filename_fe,"r");

   int i, j;
   for( i=0 ; i<XINPUT ; ++i ){
      for( j=0 ; j<ZINPUT ; ++j ){
         fscanf( vxInputFile,  "%lf", &vxInput[i][j]  );
         fscanf( vzInputFile,  "%lf", &vzInput[i][j]  );
         fscanf( rhoInputFile, "%lf", &rhoInput[i][j] );
         fscanf( NInputFile,   "%lf", &NInput[i][j]   );
         fscanf( OInputFile,   "%lf", &OInput[i][j]   );
         fscanf( SiInputFile,  "%lf", &SiInput[i][j]  );
         fscanf( FeInputFile,  "%lf", &FeInput[i][j]  );
      }
   }

   fclose(vxInputFile);
   fclose(vzInputFile);
   fclose(rhoInputFile);
   fclose(NInputFile);
   fclose(OInputFile);
   fclose(SiInputFile);
   fclose(FeInputFile);
   printf("Done reading input data, closed files\n");
}

void initial( double * prim , double * xi , double t , bool debug ){

   double x, y, z;
   double vx, vy, vz, vcyl;
   int i, j;
   int minIndex_i = 0;
   int minIndex_j = 0;
   double dist2;
   double minDist2 = 1.0e100;
   double rhoRead, NRead, ORead, SiRead, FeRead;
   //bool isAtmosphere = false;
   double scaleFactor;

   // scale rho to current time
   scaleFactor = tinput*tinput*tinput/t0/t0/t0;

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
   //vr = sqrt(vx*vx+vy*vy+vz*vz);
   vcyl = sqrt(vx*vx+vy*vy); // velocity in x-y plane

   // ignore initial data with v > vmax
   //isAtmosphere = vr > vmax;

   // do a nearest neighbor search if we're within the ejecta
   //if( !isAtmosphere ) {
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
   //}

   //if( !isAtmosphere ) {
      // set density to value of nearest neighbor
   rhoRead = rhoInput[minIndex_i][minIndex_j]*scaleFactor;
      //if( rhoRead < rhoISM ) isAtmosphere = true;
   //}

   // set abundances to those of nearest neighbor
   NRead  = NInput[minIndex_i][minIndex_j];
   ORead  = OInput[minIndex_i][minIndex_j];
   SiRead = SiInput[minIndex_i][minIndex_j];
   FeRead = FeInput[minIndex_i][minIndex_j];

   // various debug messages
   //if ( debug || false ) {
   //   printf("found neighbor with indices %d %d\n",minIndex_i,minIndex_j);
   //   printf("read rho %5.3e\n",rhoRead);
   //   printf("x y z minDist2 %5.3e %5.3e %5.3e %5.3e\n",x,y,z,minDist2);
   //}

   // define primitives
   if( rhoRead > rhoISM ) {
      prim[RHO] = rhoRead;
      prim[UU1] = vx;
      prim[UU2] = vy;
      prim[UU3] = vz;
      prim[XXX] = 1.0; // tracks ejecta fraction
      prim[TRACER_N]  = NRead;
      prim[TRACER_O]  = ORead;
      prim[TRACER_SI] = SiRead;
      prim[TRACER_FE] = FeRead;
   } else {
      prim[RHO] = rhoISM;
      prim[UU1] = 0.0;
      prim[UU2] = 0.0;
      prim[UU3] = 0.0;
      prim[XXX] = 0.0;
      prim[TRACER_N]  = 0.0;
      prim[TRACER_O]  = 0.0;
      prim[TRACER_SI] = 0.0;
      prim[TRACER_FE] = 0.0;
   }

   // set pressure to small fraction of ram pressure
   prim[PPP] = ramPressureFactor*vmax*vmax*prim[RHO];
}

