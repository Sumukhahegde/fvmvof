#include<stdio.h>        
#include<stdlib.h>       
#include<math.h>         
#include<string.h>  
#include<stdbool.h>      
#include "fvmvof.h"

void diffusion_implicit( Field * phi, Constant constant, 
    double * tmp
    )
{
  double dx = constant.dx, dy = constant.dy, dz =constant.dz;
  double nu = constant.nu;
  double dt = constant.dt;
  int i,j, l, m;
  int N = phi->N;
  int N_x = phi->N_x;
  int N_y = phi->N_y;
  double phi_s, phi_n , phi_e, phi_w;
  double * tmp_1, * val;
  BC_type * bc = phi->bc;
  val  = phi->val;

  for(i = 0 ; i <N; i++)
    tmp[i] = val[i];

  for(j=2;j<N_y-2;j++){
    l = j*N_x;
    double * val_e, *val_w, * val_s, * val_n, *val_p;
    BC_type * bc_e, *bc_w, * bc_s, * bc_n, *bc_p;

    val_p = &val[l];
    val_e = &val[l+1];
    val_w = &val[l-1];
    val_n = &val[l+N_x];
    val_s = &val[l-N_x];
    bc_p = &bc[l];
    bc_e = &bc[l+1];
    bc_w = &bc[l-1];
    bc_n = &bc[l+N_x];
    bc_s = &bc[l-N_x];
    tmp_1= &tmp[l];
    for(i = 2 ;i<  N_x - 2 ;i++){
      tmp_1[i] = val_p[i] - dt* nu* ((
            ( val_e[i] )
            +( val_w[i] )
            -2.0*val_p[i])*dy/dx + (
              (val_n[i] )
              +(val_s[i])
              -2.0*val_p[i])*dx/dy)/(dx*dy) ; 
    }
    //  asm volatile("": "+m"(val_p),"+m"(val_e), "+m"(val_w),"+m"(val_n),"+m"(val_s),"+m"(bc_p),"+m"(bc_e),"+m"(bc_w),"+m"(bc_n),"+m"(bc_s), "+m"(tmp_1));
  }

  for(m=1;m<N_y-1;m++){
    for(l= 1; l <N_x-1; l++){
      if(l ==1 || m == 1 || l==N_x-2 || m == N_y -2){
        i = m*N_x + l;

        tmp[i] = phi->val[i] - dt* nu* ((                                      
              ((phi->bc[EAST]==NONE)?phi->val[EAST] : ((phi->bc[EAST]==DIRICHLET)?2.0*phi->val[EAST] - phi->val[i] : phi->val[i]))
              +( (phi->bc[WEST]==NONE)?phi->val[WEST] : ((phi->bc[WEST]==DIRICHLET)?2.0*phi->val[WEST] - phi->val[i] : phi->val[i]))
              -2.0*phi->val[i])*dy/dx + (                                            
              ((phi->bc[NORTH]==NONE)?phi->val[NORTH] : ((phi->bc[NORTH]==DIRICHLET)? 2.0*phi->val[NORTH] - phi->val[i]: phi->val[i]))
             +((phi->bc[SOUTH]==NONE)?phi->val[SOUTH] : ((phi->bc[SOUTH]==DIRICHLET)? 2.0*phi->val[SOUTH] - phi->val[i]: phi->val[i]))
                -2.0*phi->val[i])*dx/dy)/(dx*dy) ;                                     
      }                                                                                                                                                    
    }
  }
  return;
}

void laplacian( Field * phi, Constant constant, 
    double * tmp
    )
{
  double dx = constant.dx, dy = constant.dy, dz =constant.dz;
  double nu = constant.nu;
  double dt = constant.dt;
  int i,j, l, m;
  int N = phi->N;
  int N_x = phi->N_x;
  int N_y = phi->N_y;
  double phi_s, phi_n , phi_e, phi_w;
  double * tmp_1, * val;
  BC_type * bc = phi->bc;
  val  = phi->val;

  for(i = 0 ; i <N; i++)
    tmp[i] = 0.0;

  for(j=2;j<N_y-2;j++){
    l = j*N_x;
    double * val_e, *val_w, * val_s, * val_n, *val_p;
    BC_type * bc_e, *bc_w, * bc_s, * bc_n, *bc_p;

    val_p = &val[l];
    val_e = &val[l+1];
    val_w = &val[l-1];
    val_n = &val[l+N_x];
    val_s = &val[l-N_x];
    bc_p = &bc[l];
    bc_e = &bc[l+1];
    bc_w = &bc[l-1];
    bc_n = &bc[l+N_x];
    bc_s = &bc[l-N_x];
    tmp_1= &tmp[l];
    for(i = 2 ;i<  N_x - 2 ;i++){
      tmp_1[i] =  ((
             val_e[i] 
            + val_w[i] 
            -2.0*val_p[i])*dy/dx + (
              (val_n[i] )
              +(val_s[i])
              -2.0*val_p[i])*dx/dy) ; 
    }
    //  asm volatile("": "+m"(val_p),"+m"(val_e), "+m"(val_w),"+m"(val_n),"+m"(val_s),"+m"(bc_p),"+m"(bc_e),"+m"(bc_w),"+m"(bc_n),"+m"(bc_s), "+m"(tmp_1));
  }

  for(m=1;m<N_y-1;m++){
    for(l= 1; l <N_x-1; l++){
      if(l ==1 || m == 1 || l==N_x-2 || m == N_y -2){
        i = m*N_x + l;

        tmp[i] =  ((                                      
              ((phi->bc[EAST]==NONE)?phi->val[EAST] : ((phi->bc[EAST]==DIRICHLET)?2.0*phi->val[EAST] - phi->val[i] : phi->val[i]))
              +( (phi->bc[WEST]==NONE)?phi->val[WEST] : ((phi->bc[WEST]==DIRICHLET)?2.0*phi->val[WEST] - phi->val[i] : phi->val[i]))
              -2.0*phi->val[i])*dy/dx + (                                            
              ((phi->bc[NORTH]==NONE)?phi->val[NORTH] : ((phi->bc[NORTH]==DIRICHLET)? 2.0*phi->val[NORTH] - phi->val[i]: phi->val[i]))
             +((phi->bc[SOUTH]==NONE)?phi->val[SOUTH] : ((phi->bc[SOUTH]==DIRICHLET)? 2.0*phi->val[SOUTH] - phi->val[i]: phi->val[i]))
                -2.0*phi->val[i])*dx/dy) ;                                     
      }                                                                                                                                                    
    }
  }
  return;
}

/*
void laplacian( Field * phi, Constant constant,
    double * tmp
    )
{
  double dx = constant.dx, dy = constant.dy, dz =constant.dz;
  int i,j, l, m;
  int N = phi->N;
  int N_x = phi->N_x;

  int N_y = phi->N_y;
  double phi_s, phi_n , phi_e, phi_w;
  for(i = 0 ; i <N; i++)
    tmp[i] = 0.0;


  for(j=1;j<N_y-1;j++){
    l = j*N_x;
    for(i = l+1 ;i< l + N_x - 1 ;i++){

      tmp[i] =  ((
            ( (phi->bc[EAST]==NONE)?phi->val[EAST] : ((phi->bc[EAST]==DIRICHLET)?2.0*phi->val[EAST] - phi->val[i] : phi->val[i])  )
            +( (phi->bc[WEST]==NONE)?phi->val[WEST] : ((phi->bc[WEST]==DIRICHLET)?2.0*phi->val[WEST] - phi->val[i] : phi->val[i]))
            -2.0*phi->val[i])*dy/dx + (
              ((phi->bc[NORTH]==NONE)?phi->val[NORTH] : ((phi->bc[NORTH]==DIRICHLET)? 2.0*phi->val[NORTH] - phi->val[i]: phi->val[i]))
              +((phi->bc[SOUTH]==NONE)?phi->val[SOUTH] : ((phi->bc[SOUTH]==DIRICHLET)? 2.0*phi->val[SOUTH] - phi->val[i]: phi->val[i]))
              -2.0*phi->val[i])*dx/dy) ;
    }   
  }
  return;
}
*/
