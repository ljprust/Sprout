enum{RHO,PPP,UU1,UU2,UU3,XXX,XX2,XX3,XX4,XX5,XX6};
enum{DEN,TAU,SS1,SS2,SS3};

// XXX = ejecta fraction
// XX2 = He
// XX3 = N
// XX4 = O
// XX5 = Si
// XX6 = Fe

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#define NUM_C 5
#define NUM_N 6
#define NUM_Q (NUM_C+NUM_N)
#define NUM_G 2

struct param_list{
   double t_min, t_max;
   int NumRepts, NumSnaps, NumChecks;
   int Out_LogTime;

   int Num_x, Num_y, Num_z;
   double Lx , Ly , Lz;

   double CFL, PLM;
   double Density_Floor, Pressure_Floor;

   double W0, eta_on;
   double MM_x0, MM_y0, MM_z0;

   double Adiabatic_Index;
   int Gravity_Switch;
   double Grav_G , Central_Mass;
   int Nozzle_Switch;
   double Nozzle_x0, Nozzle_y0, Nozzle_z0;
   double Nozzle_pow, Nozzle_v;

   int restart_flag;
};

struct cell{
   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double gradx[NUM_Q];
   double grady[NUM_Q];
   double gradz[NUM_Q];
   double pblax[NUM_Q];
   double pblay[NUM_Q];
   double pblaz[NUM_Q];

   double xi[3];
   double RKxi[3];
};

struct domain{
   struct cell * theCells;
   int Ng,Nx,Ny,Nz;
   double dx,dy,dz;

   double W; 

   time_t Wallt_init;
   int count_steps;

   int rank,size;
   int dim_rank[3];
   int dim_size[3];
   int left_rank[3];
   int right_rank[3];
   MPI_Comm theComm;

   struct param_list theParList;

   double t;
   double t_init, t_fin;
   int nrpt;
   int N_rpt;
   int nsnp;
   int N_snp;
   int nchk;
   int N_chk;

   int final_step;
};
