
#include "../defs.h"
#include <string.h>

static double Lx    = 0.0;
static double Ly    = 0.0;
static double d     = 0.0;
static double x_cen = 0.0;
static double y_cen = 0.0;
static double z_cen = 0.0;


static int NL = 0;
static double * rr     = NULL;
static double * rho    = NULL;
static double * xx_Si  = NULL;
static double * xx_S   = NULL;
static double * xx_Fe  = NULL;


// initial.dat has to be in the 'Initial' directory
// has to have 5 columns: r, rho, XX1, XX2 and XX3


int countlines(char * filename){
   FILE *pFile = fopen(filename, "r");
   int lines=0;
   char c;
   while ((c = fgetc(pFile)) != EOF){
      if (c == '\n') ++lines;
   }
   fclose(pFile);
   return(lines);
}


int getTable( void ){
   int nL = countlines("Initial/initial_ddt_2013_n100.dat"); printf("Nl = %i\n",nL);
   rr     = (double *) malloc( nL*sizeof(double) );
   rho    = (double *) malloc( nL*sizeof(double) );
   xx_Si  = (double *) malloc( nL*sizeof(double) );
   xx_S   = (double *) malloc( nL*sizeof(double) );
   xx_Fe  = (double *) malloc( nL*sizeof(double) );
   FILE * pFile = fopen("Initial/initial_ddt_2013_n100.dat","r");
   int l;
   for( l=0 ; l<nL ; ++l ){
      fscanf(pFile,"%lf %lf %lf %lf %lf\n",&(rr[l]),&(rho[l]),&(xx_Si[l]),&(xx_S[l]),&(xx_Fe[l]));
   }
   fclose(pFile);

   FILE * pFile2 = fopen("Initial/initial_write.dat","w");
   for( l=0 ; l<nL ; ++l ){
      fprintf(pFile2,"%e %e %e %e %e\n",(rr[l]),(rho[l]),(xx_Si[l]),(xx_S[l]),(xx_Fe[l]));
   }
   fclose(pFile2);


   return(nL);
}


void setICParams( struct domain * theDomain ){
   NL = getTable();

   Lx = theDomain->theParList.Lx;
   Ly = theDomain->theParList.Ly;
   x_cen = theDomain->theParList.MM_x0 * theDomain->theParList.Lx;
   y_cen = theDomain->theParList.MM_y0 * theDomain->theParList.Ly;
   z_cen = theDomain->theParList.MM_z0 * theDomain->theParList.Lz;
   if( theDomain->theParList.Num_x!=1 ) ++d;
   if( theDomain->theParList.Num_y!=1 ) ++d;
   if( theDomain->theParList.Num_z!=1 ) ++d;
}



void initial( double * prim , double * xi , double t ){

   double x   = xi[0] - x_cen;
   double y   = 0.;
   if(d>1.) y = xi[1] - y_cen;
   double z   = 0.;
   if(d>2.) z = xi[2] - z_cen;
   double r   = sqrt( x*x + y*y +z*z );
   //printf("r[0] = %e, r = %e, r[-1] = %e\n",rr[0],r,rr[NL-1]);

   int l=0;
   while( rr[l] < r && l < NL-1 ) ++l;

   if(l>NL-2){
      prim[UU1] = 0.;
      prim[UU2] = 0.;
      prim[UU3] = 0.;
      prim[RHO] = rho[NL-1];
      prim[PPP] = 1e-5*rho[NL-1];
      if( NUM_N > 0 ) prim[NUM_C] = 0.0;
      if( NUM_N > 1 ) prim[NUM_C+1] = 0.0;
      if( NUM_N > 2 ) prim[NUM_C+2] = 0.0;
      if( NUM_N > 3 ) prim[NUM_C+3] = 0.0;
   }else{
      double rp = rr[l];
      double rm = rr[l-1];
      double drm = fabs(r-rm);
      double drp = fabs(rp-r);

      double rh    = (rho[l-1]*drp + rho[l]*drm)/(drp+drm);
      double X_Si  = (xx_Si[l-1]*drp  + xx_Si[l]*drm )/(drp+drm);
      double X_S   = (xx_S[l-1] *drp  + xx_S[l]*drm  )/(drp+drm);
      double X_Fe  = (xx_Fe[l-1]*drp  + xx_Fe[l]*drm )/(drp+drm);

      prim[UU1] = x/t;
      prim[UU2] = y/t;
      prim[UU3] = z/t;
      prim[RHO] = rh;
      prim[PPP] = 1e-5*rh;
      if( NUM_N > 0 ) prim[NUM_C] = 1.0;
      if( NUM_N > 1 ) prim[NUM_C+1] = X_Si;
      if( NUM_N > 2 ) prim[NUM_C+2] = X_S;
      if( NUM_N > 3 ) prim[NUM_C+3] = X_Fe;
   }

   

   
}
