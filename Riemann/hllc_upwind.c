#include "../defs.h"

void prim2cons( double * , double * , double * , double );
void flux( double * , double * , double * , double * );
void vel( double * , double * , double * , double * , double * , double * );
void getUstar( double * , double * , double * , double , double , double * );

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
         *C_U = 1./(1. + W*dt*2.);
      }
      if( no_of_dims==2 ){
         *C_F = (1. + W*dt)/(1. + W*dt*2.);
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
   double faceV2;

   double * xl , * xr;
   xl = cL->xi;
   xr = cR->xi;

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + .5 * (cL->gradx[q]*dx*n[0] + cL->grady[q]*dy*n[1] + cL->gradz[q]*dz*n[2]);
      primR[q] = cR->prim[q] - .5 * (cR->gradx[q]*dx*n[0] + cR->grady[q]*dy*n[1] + cR->gradz[q]*dz*n[2]);
   }

   double Sl,Sr,Ss;

   double Fk[NUM_Q];
   double Uk[NUM_Q];

   double fluxFace[NUM_Q];
   double uFace[NUM_Q];
   double fluxCorrection[NUM_Q];
   double Flux[NUM_Q];

   double wn = get_wn( xl , dx , dy , dz , W , theDIM );

   double C_F, C_U;
   get_flux_coefficients( no_of_dims , first_step , last_step , W , dt , &C_F , &C_U );

   // boost prims into face frame
   primL[UU1] -= wn*n[0];
   primL[UU2] -= wn*n[1];
   primL[UU3] -= wn*n[2];
   primR[UU1] -= wn*n[0];
   primR[UU2] -= wn*n[1];
   primR[UU3] -= wn*n[2];

   vel( primL , primR , &Sl , &Sr , &Ss , n );

   if( 0.0 < Sl ){
      flux( primL , Fk , xl , n );
      prim2cons( primL , Uk , xl , 1.0 );

      for( q=0 ; q<NUM_Q ; ++q ){
         //Flux[q] = C_F*Fk[q] - C_U*wn*Uk[q];
         fluxFace[q] = Fk[q];
         uFace[q] = Uk[q];
      }
   }else if( 0.0 > Sr ){
      flux( primR , Fk , xr , n );
      prim2cons( primR , Uk , xr , 1.0 );

      for( q=0 ; q<NUM_Q ; ++q ){
         //Flux[q] = C_F*Fk[q] - C_U*wn*Uk[q];
         fluxFace[q] = Fk[q];
         uFace[q] = Uk[q];
      }
   }else{
      double Ustar[NUM_Q];
      if( 0.0 < Ss ){
         prim2cons( primL , Uk , xl , 1.0 );
         getUstar( primL , Ustar , xl , Sl , Ss , n );
         flux( primL , Fk , xl , n );

         for( q=0 ; q<NUM_Q ; ++q )
            //Flux[q] = C_F*( Fk[q] + Sl*( Ustar[q] - Uk[q] ) ) - C_U*wn*Ustar[q];
            fluxFace[q] = Fk[q] + Sl*( Ustar[q] - Uk[q] );
            uFace[q] = Ustar[q];
      }else{
         prim2cons( primR , Uk , xr , 1.0 );
         getUstar( primR , Ustar , xr , Sr , Ss , n );
         flux( primR , Fk , xr , n );

         for( q=0 ; q<NUM_Q ; ++q )
            //Flux[q] = C_F*( Fk[q] + Sr*( Ustar[q] - Uk[q] ) ) - C_U*wn*Ustar[q];
            fluxFace[q] = Fk[q] + Sl*( Ustar[q] - Uk[q] );
            uFace[q] = Ustar[q];
      }
   }

   double Q2dotw = fluxFace[SS1] * faceVelocity[0]
                 + fluxFace[SS2] * faceVelocity[1]
                 + fluxFace[SS3] * faceVelocity[2];
   fluxCorrection[DEN] = (C_F-C_U) * fluxFace[DEN];
   fluxCorrection[SS1] = C_F * fluxFace[DEN] * faceVelocity[0]
                       + (C_F-C_U)*(fluxFace[SS1]+fluxFace[DEN]*faceVelocity[0]);
   fluxCorrection[SS2] = C_F * fluxFace[DEN] * faceVelocity[1]
                       + (C_F-C_U)*(fluxFace[SS2]+fluxFace[DEN]*faceVelocity[1]);
   fluxCorrection[SS3] = C_F * fluxFace[DEN] * faceVelocity[2]
                       + (C_F-C_U)*(fluxFace[SS3]+fluxFace[DEN]*faceVelocity[2]);
   fluxCorrection[TAU] = C_F * ( -Q2dotw - 0.5 * fluxFace[DEN] * faceV2 )
                       + (C_F-C_U)*(fluxFace[TAU]+Q2dotw+0.5*fluxFace[DEN]*faceV2);

   for( q=0 ; q<NUM_Q ; ++q ){
      Flux[q] = fluxFace[q] - fluxCorrection[q];
   }

   double dA = dy*dz*n[0] + dz*dx*n[1] + dx*dy*n[2];
   for( q=0 ; q<NUM_Q ; ++q ){
      cL->cons[q] -= Flux[q]*dt*dA;
      cR->cons[q] += Flux[q]*dt*dA;
   }

}
