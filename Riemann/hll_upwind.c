
#include "../defs.h"

void prim2cons( double * , double * , double * , double );
void flux( double * , double * , double * , double * );
void vel( double * , double * , double * , double * , double * , double * );

double get_wn( double * , double , double , double , double , int );

void get_flux_coefficients( int no_of_dims , int first_step , int last_step , double W , double dt , double * C_F , double * C_U ){

   if( first_step==1 ){
      if( no_of_dims==1 ){
         *C_F = 1.;
         *C_U = 1.;
      }
      if( no_of_dims==2 ){
         *C_F = 1. + W*dt/2.;
         *C_U = 1. + W*dt/2.;
      }
      if( no_of_dims==3 ){
         *C_F = 1. + W*dt + W*W*dt*dt/3.;
         *C_U = 1. + W*dt + W*W*dt*dt/3.;
      }
   }

   if( first_step==0 ){
      if( no_of_dims==1 ){
         *C_F = 1.;
         *C_U = 1.;
      }
      if( no_of_dims==2 ){
         *C_F = (1. + W*dt) /(1. + W*dt*2.);
         *C_U = (1. + W*dt)/(1. + W*dt*2.)/(1. + W*dt*2.);
      }
      if( no_of_dims==3 ){
         *C_F = (1. + W*dt*2. + W*W*dt*dt*4./3.)/(1. + W*dt*2.)/(1. + W*dt*2.);
         *C_U = (1. + W*dt*2. + W*W*dt*dt*4./3.)/(1. + W*dt*2.)/(1. + W*dt*2.)/(1. + W*dt*2.);
      }
   }

}

void riemann1D( struct cell * cL , struct cell * cR , double dx , double dy , double dz , double dt , double W , int no_of_dims , int theDIM , int first_step , int last_step ){

   double primL[NUM_Q];
   double primR[NUM_Q];

   double n[3] = {0.0};
   n[theDIM] = 1.0;
   double faceVelocity[3] = {0.0};

   double * xl , * xr;
   xl = cL->xi;
   xr = cR->xi;

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + .5 * (cL->gradx[q]*dx*n[0] + cL->grady[q]*dy*n[1] + cL->gradz[q]*dz*n[2]);
      primR[q] = cR->prim[q] - .5 * (cR->gradx[q]*dx*n[0] + cR->grady[q]*dy*n[1] + cR->gradz[q]*dz*n[2]);
   }

   double Sl,Sr,Ss;

   double Fl[NUM_Q];
   double Fr[NUM_Q];
   double Ul[NUM_Q];
   double Ur[NUM_Q];
   double Ul_unboosted[NUM_Q];
   double Ur_unboosted[NUM_Q];

   double fluxFace[NUM_Q];
   double uFace[NUM_Q];
   double fluxCorrection[NUM_Q];
   double Flux[NUM_Q];

   double wn = get_wn( xl , dx , dy , dz , W , theDIM );

   faceVelocity[theDIM] = wn;

   double C_F, C_U;
   get_flux_coefficients( no_of_dims , first_step , last_step , W , dt , &C_F , &C_U );

   prim2cons( primL , Ul_unboosted , xl , 1.0 );
   prim2cons( primR , Ur_unboosted , xr , 1.0 );

   // boost prims into face frame
   primL[UU1] -= wn*n[0];
   primL[UU2] -= wn*n[1];
   primL[UU3] -= wn*n[2];
   primR[UU1] -= wn*n[0];
   primR[UU2] -= wn*n[1];
   primR[UU3] -= wn*n[2];

   vel( primL , primR , &Sl , &Sr , &Ss , n );

   // convert to cons
   prim2cons( primL , Ul , xl , 1.0 );
   prim2cons( primR , Ur , xr , 1.0 );

   // compute flux function
   flux( primL , Fl , xl , n );
   flux( primR , Fr , xr , n );

   if( 0.0 < Sl ){
      for( q=0 ; q<NUM_Q ; ++q ){
         fluxFace[q] = Fl[q];
         uFace[q] = Ul_unboosted[q];
         //Flux[q] = C_F*Fl[q] - C_U*wn*Ul[q];
      }
   }else if( 0.0 > Sr ){
      for( q=0 ; q<NUM_Q ; ++q ){
         fluxFace[q] = Fr[q];
         uFace[q] = Ur_unboosted[q];
         //Flux[q] = C_F*Fr[q] - C_U*wn*Ur[q];
      }
   }else{

      double Fstar[NUM_Q];
      double Ustar[NUM_Q];

      double aL =  Sr;
      double aR = -Sl;

      for( q=0 ; q<NUM_Q ; ++q ){

         Fstar[q] = ( aL*Fl[q] + aR*Fr[q] + aL*aR*( Ul[q] - Ur[q] ) )/( aL + aR );
         Ustar[q] = ( aR*Ul[q] + aL*Ur[q] + Fl[q] - Fr[q] )/( aL + aR );

         fluxFace[q] = Fstar[q];
         //uFace[q] = Ustar;
      }
      double primFace[NUM_Q];
      cons2prim( Ustar, primFace, xl, 1.0 );
      primFace[UU1] += wn*n[0];
      primFace[UU2] += wn*n[1];
      primFace[UU3] += wn*n[2];
      prim2cons( primFace, uFace, xl, 1.0 );
      //Flux[q] = C_F*Fstar - C_U*wn*Ustar;
   }

   double Q2dotw = fluxFace[SS1] * faceVelocity[0]
                 + fluxFace[SS2] * faceVelocity[1]
                 + fluxFace[SS3] * faceVelocity[2];

   fluxCorrection[DEN] = (C_U-C_F) * wn * uFace[DEN];
   fluxCorrection[SS1] = -C_F * fluxFace[DEN] * faceVelocity[0]
                       + (C_U-C_F) * wn * uFace[SS1];
   fluxCorrection[SS2] = -C_F * fluxFace[DEN] * faceVelocity[1]
                       + (C_U-C_F) * wn * uFace[SS2];
   fluxCorrection[SS3] = -C_F * fluxFace[DEN] * faceVelocity[2]
                       + (C_U-C_F) * wn * uFace[SS3];
   fluxCorrection[TAU] = -C_F * ( Q2dotw + 0.5 * fluxFace[DEN] * wn * wn )
                       + (C_U-C_F) * wn * uFace[TAU];
   fluxCorrection[XXX] = (C_U-C_F) * wn * uFace[XXX];

   for( q=0 ; q<NUM_Q ; ++q ){
      Flux[q] = C_F*fluxFace[q] - fluxCorrection[q];
   }

   double dA = dy*dz*n[0] + dz*dx*n[1] + dx*dy*n[2];
   for( q=0 ; q<NUM_Q ; ++q ){
      cL->cons[q] -= Flux[q]*dt*dA;
      cR->cons[q] += Flux[q]*dt*dA;
   }

   //printf("%e %e\n", xl[0]+dx/2-.5,Flux[SS1]); //for file

   //printf("x = %e, Delta_F = %e\n", xl[0]+dx/2-.5, Flux[DEN]); //mass

   //printf("x = %e, Delta_F = %e\n", xl[0]+dx/2-.5, -Flux[SS1] + (xl[0]+dx/2-.5)*(xl[0]+dx/2-.5) + 30. - (xl[0]+dx/2) );
   //printf("x = %e, Delta_F = %e\n", xl[0]+dx/2-.5,-Flux[SS1]+primL[PPP]+(xl[0]+dx/2-.5)*(xl[0]+dx/2-.5)); //momentum

   //printf("x = %e, Delta_F = %e\n", xl[0]+dx/2-.5,-Flux[TAU]+0.5*pow((xl[0]+dx/2-.5),3.) + 2.5*(30.-(xl[0]+dx/2))*(xl[0]+dx/2-.5));  //energy
}
