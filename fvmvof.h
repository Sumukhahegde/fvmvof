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
int N_cells, N_cells_x, N_cells_y, N_cells_z ,N_flux_x,N_flux_y;
int N_points, N_points_x, N_points_y, N_points_z ;
typedef enum{
  INSIDE,
  XMIN,
  XMAX,
  YMIN,
  YMAX,
  ZMIN,
  ZMAX,
  SOLID
} PATCH_type;
PATCH_type patch;

typedef enum{
  POISSON,
  HELMHOLTZ
}equation_type;
/*
typedef enum{
  NONE,
  WALL,
  AMBIENT,
  INLET
} BC_type;
*/
typedef enum{
  NONE,
  DIRICHLET,
  NEUMANN,
} BC_type;

BC_type *bc, *p_bc, *u_x_bc, *u_y_bc, *u_z_bc;
//BC_type p_bc[8], u_x_bc[8], u_y_bc[8], u_z_bc[8];
double p_bc_val[8],u_x_bc_val[8],u_y_bc_val[8],u_z_bc_val[8];

#endif
