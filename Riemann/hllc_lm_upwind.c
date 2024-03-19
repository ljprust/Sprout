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
   double Ul[NUM_Q];
   double Fr[NUM_Q];
   double Ur[NUM_Q];
   double Ul_unboosted[NUM_Q];
   double Ur_unboosted[NUM_Q];

   double fluxFace[NUM_Q];
   double uFace[NUM_Q];
   double fluxCorrection[NUM_Q];
   double Flux[NUM_Q];

   double wn = get_wn( xl , dx , dy , dz , W , theDIM );

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

   double Mach_L, Mach_R, Mach_local, Mach_limit;
   Mach_L = fabs( (primL[UU1+theDIM]-wn)/sqrt(5./3.*primL[PPP]/primL[RHO]) );
   Mach_R = fabs( (primR[UU1+theDIM]-wn)/sqrt(5./3.*primR[PPP]/primR[RHO]) );
   Mach_local = Mach_L;
   if(Mach_R>Mach_L) Mach_local = Mach_R;
   Mach_limit = 0.5;

   double phi = 1.0;
   if( Mach_local<Mach_limit ) phi = sin(Mach_local/Mach_limit*M_PI/2.);

   double Sl_LM = Sl*phi;
   double Sr_LM = Sr*phi;


   if( 0.0 < Sl_LM ){
      flux( primL , Fl , xl , n );
      //prim2cons( primL , Ul , xl , 1.0 );

      for( q=0 ; q<NUM_Q ; ++q ){
         //Flux[q] = C_F*Fl[q] - C_U*wn*Ul[q];
         fluxFace[q] = Fl[q];
         uFace[q] = Ul_unboosted[q];
      }
   }else if( 0.0 > Sr_LM ){
      flux( primR , Fr , xr , n );
      //prim2cons( primR , Ur , xr , 1.0 );

      for( q=0 ; q<NUM_Q ; ++q ){
         //Flux[q] = C_F*Fr[q] - C_U*wn*Ur[q];
         fluxFace[q] = Fr[q];
         uFace[q] = Ur_unboosted[q];
      }
   }else{
      double UstarL[NUM_Q];
      double UstarR[NUM_Q];
      double Ustar[NUM_Q];
      prim2cons( primL , Ul , xl , 1.0 );
      prim2cons( primR , Ur , xr , 1.0 );
      getUstar( primL , UstarL , xl , Sl , Ss , n );
      getUstar( primR , UstarR , xr , Sr , Ss , n );
      flux( primL , Fl , xl , n );
      flux( primR , Fr , xr , n );

      if( 0.0 < Ss ){
         for( q=0 ; q<NUM_Q ; ++q ) {
            //Flux[q] = C_F*0.5*( Fl[q] + Fr[q] + Sl_LM*(UstarL[q]-Ul[q]) + Sr_LM*(UstarR[q]-Ur[q]) + Ss*(UstarL[q]-UstarR[q]) ) - C_U*wn*UstarL[q];
            //Flux[q] = C_F*( Fl[q] + Sl*( UstarL[q] - Ul[q] ) ) - C_U*wn*UstarL[q];
            fluxFace[q] = 0.5*( Fl[q] + Fr[q] + Sl_LM*(UstarL[q]-Ul[q]) + Sr_LM*(UstarR[q]-Ur[q]) + Ss*(UstarL[q]-UstarR[q]) );
            Ustar[q] = UstarL[q];
         }
      }else{
         for( q=0 ; q<NUM_Q ; ++q ) {
            //Flux[q] = C_F*0.5*( Fl[q] + Fr[q] + Sl_LM*(UstarL[q]-Ul[q]) + Sr_LM*(UstarR[q]-Ur[q]) - Ss*(UstarL[q]-UstarR[q]) ) - C_U*wn*UstarR[q];
            //Flux[q] = C_F*( Fr[q] + Sr*( UstarR[q] - Ur[q] ) ) - C_U*wn*UstarR[q];
            fluxFace[q] = 0.5*( Fl[q] + Fr[q] + Sl_LM*(UstarL[q]-Ul[q]) + Sr_LM*(UstarR[q]-Ur[q]) - Ss*(UstarL[q]-UstarR[q]) );
            Ustar[q] = UstarR[q];
         }
      }
      double primFace[NUM_Q];
      cons2prim( Ustar, primFace, xl, 1.0 );
      primFace[UU1] += wn*n[0];
      primFace[UU2] += wn*n[1];
      primFace[UU3] += wn*n[2];
      prim2cons( primFace, uFace, xl, 1.0 );
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

   for( q=0 ; q<NUM_Q ; ++q ){
      Flux[q] = C_F*fluxFace[q] - fluxCorrection[q];
   }

   double dA = dy*dz*n[0] + dz*dx*n[1] + dx*dy*n[2];
   for( q=0 ; q<NUM_Q ; ++q ){
      cL->cons[q] -= Flux[q]*dt*dA;
      cR->cons[q] += Flux[q]*dt*dA;
   }

}
