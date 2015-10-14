#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>
#include "fvmvof.h"

/* Finite Volume Method */
/* Code to solve NS. First step, solve poisson equation in 2D */
/* Define variables */

void Compute_AX(double * );
int solve_BiCGSTAB(void);
static void write_vtk(int);
//void set_ghosts(void);
//void set_bc(double * var , BC_type * bc, double * bc_value);

void allocate_variable(field_variable * phi);
int main(int argc, char *argv[])
{
  double l_x = 1.0;
  double l_y = 1.0;

  p   = malloc(sizeof(field_variable));
  phi   = malloc(sizeof(field_variable));
  u_x = malloc(sizeof(field_variable));
  u_y = malloc(sizeof(field_variable));
  u_z = malloc(sizeof(field_variable));
  acc_x = malloc(sizeof(field_variable));
  acc_y = malloc(sizeof(field_variable));
  divergence = malloc(sizeof(field_variable));
  // cell dimensions of each field variable
  int N_x = 250+2, N_y = 500+2;
  p->N_x = N_x; p->N_y = N_y; p->N_z = 1;
  phi->N_x = N_x; phi->N_y = N_y; phi->N_z = 1;
  u_x->N_x = N_x; u_x->N_y = N_y; u_x->N_z = 1;
  u_y->N_x = N_x; u_y->N_y = N_y; u_y->N_z = 1;
  acc_x->N_x = N_x; acc_x->N_y = N_y; acc_x->N_z = 1;
  acc_y->N_x = N_x; acc_y->N_y = N_y; acc_y->N_z = 1;
  divergence->N_x = N_x; divergence->N_y = N_y; divergence->N_z = 1;
//  u_x->N_x = (p->N_x -1) ; u_x->N_y =p->N_y; u_x->N_z = p->N_z;
//  u_y->N_x = p->N_x  ; u_y->N_y =p->N_y-1; u_y->N_z = p->N_z;
int N_cells_x = p->N_x;
int N_cells_y = p->N_y;
int N_cells_z = p->N_z;
int N_cells = N_cells_x *  N_cells_y;
  allocate_variable(p);
  allocate_variable(phi);
  allocate_variable(u_x);
  allocate_variable(acc_x);
  allocate_variable(u_y);
  allocate_variable(u_z);
  allocate_variable(divergence);
  allocate_variable(acc_y);
double * temp = malloc(N_cells*sizeof(double));
  dx = l_x / (N_cells_x-2);
  dy = l_y / (N_cells_y-2);
  dz = 1.0;
  int i , l, m; 
  //initial velocity field for advection
  for(i=0;i<N_cells;i++){
      l= i%N_cells_x;
      m =(int) i/N_cells_x;
    p->val[i] = 0.0;
    u_x->val[i] = 0.0;
    u_y->val[i] = 0.0;
    u_z->val[i] = 0.0;
    temp[i] = 0.0;
  }
  for(i=0;i<8;i++){
    p->bc_val[i] = 0.0;
    phi->bc_val[i] = 0.0;
    u_x->bc_val[i] = 0.0;
    u_y->bc_val[i] = 0.0;
  }
//  u_y->bc_val[YMAX] = 1.0;

  set_ghosts();
  for(i=0;i<N_cells;i++){
      l= i%N_cells_x;
      m =(int) i/N_cells_x;
    p->val[i] = 0.0;
    u_x->val[i] = 1.0;
    u_y->val[i] = 1.0;
    u_z->val[i] = 0.0;
    phi->val[i] = 0.0;
//    phi->bc[i] = NONE;
    if(l>10 && m >20 && l<50 && m<100)
      phi->val[i] = 10.0;
  }
  
  set_bc(phi); set_bc(u_x);
  set_bc(u_y);
  set_bc(u_z);
  dt = 0.0001;
  double nu = 0.005;

  double peclet = 1.0*1.0 *dx/nu;
  if(peclet>=2.0){
    printf("peclet number not less than 2.0 \n");
    exit(1);
  }
  if(dt>=0.5*dx/1.0){
    printf("t not within CFL criterion \n");
    exit(1);
  }
    //advection
  int qq;
  for(qq = 0;  qq<1000 ; qq++){
    for(i=0;i<N_cells;i++)
      temp[i] = 0.0;
    advection(phi, u_x, u_y, u_z, temp );
    diffusion(phi, nu, temp);
    for( i=0;i<N_cells;i++){
      phi->val[i] = phi->val[i] + dt*temp[i]/(dx*dy);
    }
 set_bc(phi);
    if(qq%10 == 0)  write_vtk(qq); 
  }
  /*
  // only for poisson ;  
  b = malloc(N_cells*sizeof(double));
  //initial and boundary conditions
  int i,l,m;
  for(i=0;i<N_cells;i++){
    p[i] = 0.0;
    u_x[i] = 0.0;
    u_y[i] = 0.0;
    u_z[i] = 0.0;
  }
  for(i=0;i<8;i++){
    u_x_bc_val[i] = 0.0;
    u_y_bc_val[i] = 0.0;
    u_z_bc_val[i] = 0.0;
    p_bc_val[i] = 0.0;
  }
  u_x_bc_val[4] = 1.0;
  //set_bc();
  set_bc(p, p_bc, p_bc_val);
  set_bc(u_x, u_x_bc, u_x_bc_val);
  set_bc(u_y, u_y_bc, u_y_bc_val);
  set_bc(u_z, u_z_bc, u_z_bc_val);
  //Solve convection, get u star 
  int test  = solve_BiCGSTAB();*/
  return 0;
}
/*
int solve_BiCGSTAB() 
{
  int i, j ;
  double *r0_star, *sj, *rj, *pj, *pstar, *sstar ; 
  double *Temp, *Uj, *Var  ;
  double rhoj_Minus, alphaj, omegaj, rhoj, betaj, H1, H2 ;
  double norm, BICGEPS = 1.0E-12; 
  int BICG_ITER ;
  bool STOP = false ;
  int N =  N_cells;

  r0_star = malloc(N*sizeof(double));
  sj = malloc(N*sizeof(double));
  rj = malloc(N*sizeof(double));
  pj = malloc(N*sizeof(double));
  pstar = malloc(N*sizeof(double));
  sstar = malloc(N*sizeof(double));
  Temp = malloc(N*sizeof(double));
  Uj = malloc(N*sizeof(double));
  Var = malloc(N*sizeof(double));

  // Start BICGSTAB iterations
  // set initial solution vector x_0 = (Uj, Vj) 
  Compute_AX(Temp) ;
  // Initial vector r_0 = b - Ax_0, and r0* = r_0
  for(i = 0 ; i < N ; i++) {
    Uj[i] = p[i] ;
    rj[i] = b[i] - Temp[i] ; 
    r0_star[i] = rj[i] ;
  }
  BICG_ITER = 0 ; norm = 0.0 ;
  do {
    // compute rhoj = (r0, r0*)
    rhoj = 0.0 ;
    for(i = 0 ; i < N ; i++) 
      rhoj += rj[i]*r0_star[i] ; 
    if( sqrt(rhoj/((double)N)) < BICGEPS ) STOP = true ;
    else {
      if( BICG_ITER == 0 ) {
        for(i = 0 ; i < N ; i++) 
          pj[i] = rj[i]; // p0 = r0 
      } else {
        betaj = (rhoj/rhoj_Minus)*(alphaj/omegaj) ;
        for(i = 0 ; i < N ; i++) 
          pj[i] = rj[i] + betaj*(pj[i] - omegaj*Var[i]);
      }
      // Solve for Upstar, Vpstar from Ku* = u...., where K is the preconditioning matrix
      // No preconditioning
      for(i = 0 ; i < N ; i++)
        pstar[i] = pj[i] ;
      // compute vj = A*pstar
      for(i = 0 ; i < N ; i++) 
        p[i] = pstar[i] ;
      Compute_AX(Temp) ;
      for(i = 0 ; i < N ; i++) 
        Var[i] = Temp[i] ;
      H1 = 0.0 ;
      for(i = 0 ; i < N ; i++) 
        H1 += Var[i]*r0_star[i] ;
      alphaj = rhoj/H1 ;
      // find sj
      for(i = 0 ; i < N ; i++) 
        sj[i] = rj[i] - alphaj*Var[i] ;
      // Solve for Upstar, Vpstar from Ku* = u..., where K is the preconditioning matrix
      // No preconditioning
      for(i = 0 ; i < N ; i++) 
        sstar[i] = sj[i] ;
      norm = 0.0 ;
      for(i = 0 ; i < N ; i++) 
        norm += sstar[i]*sstar[i] ;
      norm = sqrt(norm/((double)N)) ;
      if( norm < BICGEPS) {
        STOP = true ; //if ||s||_2 is small x_i = x_{i-1}+alphai*p_i
        for(i = 0 ; i < N ; i++) 
          Uj[i] += alphaj*pstar[i] ;
      } else {
        // compute t = As
        for(i = 0 ; i < N ; i++) 
          p[i] = sstar[i] ; 
        Compute_AX(Temp) ;
        H1 = H2 = 0.0 ;
        for(i =0 ; i < N ; i++) {
          H1 += Temp[i]*sj[i] ;
          H2 += Temp[i]*Temp[i] ;
        }	
        omegaj = H1/H2;
        // find xj 
        norm = 0.0 ;
        for(i = 0 ; i < N ; i++) {
          H1 = (alphaj*pstar[i] + omegaj*sstar[i]) ;
          Uj[i] += H1 ;
          norm += H1*H1 ;
        }
        norm = sqrt(norm/((double)N)) ;
        if(norm < BICGEPS) STOP = true ;
        // find rjplusone
        for(i = 0 ; i < N ; i++) 
          rj[i] = sj[i] - omegaj*Temp[i];
        rhoj_Minus = rhoj ;
      }
    }
    BICG_ITER++;
    if(BICG_ITER%100 == 0) printf("%d \t %lf \n", BICG_ITER, norm );
  }while( (BICG_ITER < 10000) && (!STOP) ) ;
  for(i = 0 ; i < N ; i++)
    p[i] = Uj[i] ;
  return 0;
}

void Compute_AX(double * Temp){
  int i,l,m;
  double phi_w, phi_e, phi_n, phi_s;
  for(i=0;i<N_cells;i++){
    if(p_bc[i] == DIRICHLET){
      Temp[i] = p[i];
    }else{
      l= i%N_cells_x;
      m =(int) i/N_cells_x;
      int south = (m-1)*N_cells_x + l, north =(m+1)*N_cells_x + l,
          west =m*N_cells_x + (l-1), east = m*N_cells_x + (l+1);
 //     if(patch[south] != NONE )
        if(p_bc[south] == DIRICHLET)
          phi_s = 2.0*p[south] - p[i];
        else if(p_bc[south] == NEUMANN)
          phi_s = p[i];
        else phi_s = p[south];

      if(p_bc[north] == DIRICHLET)
        phi_n = 2.0*p[north] - p[i];
      else if(p_bc[north] == NEUMANN)
        phi_n = p[i];
      else phi_n = p[north];

      if(p_bc[east] == DIRICHLET)
        phi_e = 2.0*p[east] - p[i];
      else if(p_bc[east] == NEUMANN)
        phi_e = p[i];
      else phi_e = p[east];

      if(p_bc[west] == DIRICHLET)
        phi_w = 2.0*p[west] - p[i];
      else if(p_bc[west] == NEUMANN)
        phi_w = p[i];
      else phi_w = p[west];

      Temp[i] = -2.0*(dy/dx + dx/dy)*p[i] + phi_w*dy/dx + phi_e * dy/dx + phi_s * dx/dy + phi_n * dx/dy ;    
    }
  }
  return;
}
*/
static void write_vtk(int q)
{
  char filename[30]; 
  sprintf(filename, "output_%05d.vtk",q);
  FILE *fp = fopen(filename, "w");

  int Nx = p->N_x+1;
  int Ny = p->N_y+1;
  int Nz = 2; //N_cells_z+1;
  int N_cells_y = p->N_y;
  int N_cells_x = p->N_x;
  int N_cells_z = p->N_z;
  int N_cells = N_cells_x * N_cells_y;
  fprintf(fp,"# vtk DataFile Version 3.0\n");     
  fprintf(fp,"particle point data\n");           
  fprintf(fp,"ASCII\n");                         
  fprintf(fp,"DATASET STRUCTURED_GRID\n");       
  fprintf(fp,"DIMENSIONS %d %d %d\n",Nx,Ny,Nz);  
  fprintf(fp,"POINTS %d double\n",Nx*Ny*Nz);
  int l,m,n, i;
  for(n = 0; n<Nz; n++){
    for(m = 0; m<Ny; m++){
      for( l = 0; l<Nx ; l ++){
        fprintf(fp,"%2.8lf %2.8lf %2.8lf\n",l*dx , m*dy, n*1.0);
      }
    }
  }
  fprintf(fp,"CELL_DATA %d\n SCALARS pressure double 1\n LOOKUP_TABLE default\n",N_cells);  
  for( l = 0; l<N_cells_y ; l ++){
    for(m = 0; m<N_cells_x; m++){
      fprintf(fp,"%2.8lf\n",p->val[l*N_cells_x + m]);
    }
  }
  fprintf(fp,"SCALARS phi double 1\n LOOKUP_TABLE default\n");  
  for( l = 0; l<N_cells_y ; l ++){
    for(m = 0; m<N_cells_x; m++){
      fprintf(fp,"%2.8lf\n",phi->val[l*N_cells_x + m]);
    }
  }
  fprintf(fp,"SCALARS boundary int 1\n LOOKUP_TABLE default\n");  
  for( l = 0; l<N_cells_y ; l ++){
    for(m = 0; m<N_cells_x; m++){
      fprintf(fp,"%d\n",p->bc[l*N_cells_x + m]);
    }
  }
  fprintf(fp,"VECTORS velocity double \n");  
  for( l = 0; l<u_x->N_y ; l ++){
    for(m = 0; m<u_x->N_x; m++){
      fprintf(fp,"%2.8lf %2.8lf %2.8lf \n",u_x->val[l*u_x->N_x + m], u_y->val[l*u_y->N_x + m], 0.0 );
    }
  }
  //  fprintf(fp,"%2.8lf %2.8lf %2.8lf\n", x.x, x.y, x.z);
  return;
}



void allocate_variable(field_variable * phi)
{
  phi->N = p->N_x * p->N_y * p->N_z;
  phi->val = malloc(phi->N * sizeof(double));
  phi->bc = malloc(phi->N * sizeof(BC_type));
  return ;
}
/*
void advection_2d(void )
{

}
*/
