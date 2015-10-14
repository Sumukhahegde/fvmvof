#ifndef FVMVOF_H_
#define FVMVOF_H_

#define DELTA 1e-10
#define PI 2.0*acos(0.0)

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
  CENTERED,
  X_STAGGERED,
  Y_STAGGERED
}data_location;

typedef enum{
  POISSON,
  HELMHOLTZ
}equation_type;

typedef enum{
  NONE,
  DIRICHLET,
  NEUMANN,
} BC_type;

typedef struct _Field Field;
struct _Field{
  data_location grid;  
  int N_x, N_y, N_z;
  int N;
  double bc_val[8];
  BC_type *bc;
  double *val;
};

typedef struct _Domain Domain;
struct _Domain{
  Field * p;
  Field * u_x;
  Field * u_y;
  Field * u_z;
  Field * phi;
};
/*
Field  *r_x, *r_y, *r_y,
       *u_x, *u_y, *u_z,
       *ust_x, *ust_y, *ust_z,
       *rho,
       *p,
       *phi,
       *temp,
       *omega_x, *omega_y, *omega_z,
       *acc_x, *acc_y,
       *divergence,
       *a_w,*a_e,*a_s,*a_n,*a_t,*a_b,
       *b;
       */
double dx,dy,dz,dt;
double mu, CFL, p_ref, rho_ref, u_ref;
/* Declare Functions */
void advection(Field * ,Field * ,Field * ,
    Field * ,Field * );
void diffusion(Field * ,double ,Field * );

void set_ghosts(Domain);
void set_bc(Field * phi);



#endif
