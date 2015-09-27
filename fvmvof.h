#ifndef FVMVOF_H_
#define FVMVOF_H_

#define DELTA 1e-10
#define PI 2.0*acos(0.0)

double *r_x, *r_y, *r_y,
       *u_x, *u_y, *u_z,
       *ust_x, *ust_y, *ust_z,
       *rho,
       *p,
       *omega_x, *omega_y, *omega_z,
       *a_w,*a_e,*a_s,*a_n,*a_t,*a_b,
       *b;
double dx,dy,dz,dt;
double mu, CFL, p_ref, rho_ref, u_ref;
int N_cells, N_cells_x, N_cells_y, N_cells_z ;
int N_points, N_points_x, N_points_y, N_points_z ;
typedef enum{
  NONE,
  WALL,
  AMBIENT,
  INLET
} BCs;
BCs *bc;
#endif
