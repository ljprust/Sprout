
#include "../defs.h"

static double g_acc = 0.0;
static double x_len = 0.0;
static double y_len = 0.0;
static double z_len = 0.0;
static double K     = 0.0;


void setICParams( struct domain * theDomain ){
   x_len = theDomain->theParList.Lx;
   y_len = theDomain->theParList.Ly;
   z_len = theDomain->theParList.Lz;
   g_acc = -1.0 * theDomain->theParList.Central_Mass * theDomain->theParList.Grav_G;
   K     = 1.0;
}

double f( double x , double y ){
   return( ( 1.-cos(2.*M_PI*x/x_len) )*( 1. - cos(2.*M_PI*y/y_len) ) );
}

void initial( double * prim , double * xi , double t ){
   
   double x = xi[0];
   double y = xi[1];
   double z = xi[2];
   double rho,Pp;
   //rho = 1.0;
   //Pp  = 8. + prim[RHO] * x * g_acc;
   
   Pp  = pow( g_acc*pow(K,-3./4.)/4.*(10.-x) , 4.);
   rho = pow( P/K , 3./4. );

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   prim[UU3] = 0.0;

}
