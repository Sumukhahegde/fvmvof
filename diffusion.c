#include<stdio.h>        
#include<stdlib.h>       
#include<math.h>         
#include<string.h>  
#include<stdbool.h>      
#include "fvmvof.h"

void diffusion( Field * phi, double nu,
    Field * tmp,
    Constant constant
    )
{
  double dx = constant.dx, dy = constant.dy, dz =constant.dz;
  int i,j, l, m;
  int N = phi->N;
  int N_x = phi->N_x;
  int N_y = phi->N_y;
  for(j = 1;j<N_y-1;j++){
    l= j*N_x;
    for(i = l+1 ;i<N_x-1;i++){
// looking for a strategy for solid boundaries within domain
     /* if(phi->bc[east] != NONE )
        phi_e = 2.0*(phi->val[east]*abs(2-phi->bc[east]) + phi->val[i]*abs(1-phi->bc[east])) - phi->val[i];
      if(phi->bc[west] != NONE )
        phi_w = 2.0*(phi->val[west]*abs(2-phi->bc[west]) + phi->val[i]*abs(1-phi->bc[west])) - phi->val[i];
      if(phi->bc[north] != NONE )
        phi_n = 2.0*(phi->val[north]*abs(2-phi->bc[north]) + phi->val[i]*abs(1-phi->bc[north])) - phi->val[i];
      if(phi->bc[south] != NONE )
        phi_s = 2.0*(phi->val[south]*abs(2-phi->bc[south]) + phi->val[i]*abs(1-phi->bc[south])) - phi->val[i];
*/
//      phi_e, phi_w;
//      tmp->val[i] = nu*(((phi->bc[EAST]==NONE) ? phi->val[EAST]: (2.0*(phi->val[EAST]*(2-phi->bc[EAST]) + phi->val[P]*abs(phi->bc[EAST] - 1)) - phi->val[i])
      tmp->val[i] = nu*(
            ( (phi->bc[EAST]==NONE) ? phi->val[EAST] : ((phi->bc[EAST]==DIRICHLET) ? 2.0*phi->val[EAST]-phi->val[P] : phi->val[P])
            + (phi->bc[WEST]==NONE) ? phi->val[WEST] : ((phi->bc[WEST]==DIRICHLET) ? 2.0*phi->val[WEST]-phi->val[P] : phi->val[P])
            -2.0*phi->val[P])*dy/dx 
          + ( (phi->bc[NORTH]==NONE)? phi->val[NORTH]: ((phi->bc[NORTH]==DIRICHLET) ? 2.0*phi->val[NORTH]-phi->val[P]:phi->val[P])            
            + (phi->bc[SOUTH]==NONE)? phi->val[SOUTH]: ((phi->bc[SOUTH]==DIRICHLET) ? 2.0*phi->val[SOUTH]-phi->val[P]:phi->val[P]) 
            -2.0*phi->val[P])*dx/dy
          );
    }
  }
return;
}
/*
void diffusion_implicit( Field * phi, Constant constant,
        double * tmp
            )
{
  double dx = constant.dx, dy = constant.dy, dz =constant.dz;
  double nu = constant.nu;
  double dt = constant.dt;
  int i, l, m;
  int N = phi->N;
  int N_x = phi->N_x;
  double phi_s, phi_n , phi_e, phi_w;
  for(i = 0;i<N;i++){
    if(phi->bc[i] == NONE ){
      l= i%N_x;
      m =(int) i/N_x;
      int south = (m-1)*N_x + l, north =(m+1)*N_x + l,
          west =m*N_x + (l-1), east = m*N_x + (l+1);
      phi_e = phi->val[east];
      phi_w = phi->val[west];
      phi_s = phi->val[south];
      phi_n = phi->val[north];
      if(phi->bc[east] != NONE )
        phi_e = 2.0*(phi->val[east]*abs(2-phi->bc[east]) + phi->val[i]*abs(1-phi->bc[east])) - phi->val[i];
      if(phi->bc[west] != NONE )
        phi_w = 2.0*(phi->val[west]*abs(2-phi->bc[west]) + phi->val[i]*abs(1-phi->bc[west])) - phi->val[i];
      if(phi->bc[north] != NONE )
        phi_n = 2.0*(phi->val[north]*abs(2-phi->bc[north]) + phi->val[i]*abs(1-phi->bc[north])) - phi->val[i];
      if(phi->bc[south] != NONE )
        phi_s = 2.0*(phi->val[south]*abs(2-phi->bc[south]) + phi->val[i]*abs(1-phi->bc[south])) - phi->val[i];

      tmp[i] = phi->val[i] - dt* nu* ((phi_e+phi_w-2.0*phi->val[i])*dy/dx + (phi_n+phi_s-2.0*phi->val[i])*dx/dy)/(dx*dy) ;
    } else
      tmp[i] = phi->val[i];
  }
  return;
}
*/

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
  for(i=0;i<N;i++)
    tmp[i] = phi->val[i];

  for(j = 1;j<N_y-1;j++){
    l= j*N_x;
    for(i = l+1 ;i< l+ N_x-1;i++){
      tmp[i] = phi->val[P] - dt*nu*(
          ( (phi->bc[EAST]==NONE) ? phi->val[EAST] : ((phi->bc[EAST]==DIRICHLET) ? 2.0*phi->val[EAST]-phi->val[P] : phi->val[P])
            + (phi->bc[WEST]==NONE) ? phi->val[WEST] : ((phi->bc[WEST]==DIRICHLET) ? 2.0*phi->val[WEST]-phi->val[P] : phi->val[P])
            -2.0*phi->val[P])*dy/dx 
          + ( (phi->bc[NORTH]==NONE)? phi->val[NORTH]: ((phi->bc[NORTH]==DIRICHLET) ? 2.0*phi->val[NORTH]-phi->val[P]:phi->val[P])     
            + (phi->bc[SOUTH]==NONE)? phi->val[SOUTH]: ((phi->bc[SOUTH]==DIRICHLET) ? 2.0*phi->val[SOUTH]-phi->val[P]:phi->val[P]) 
            -2.0*phi->val[P])*dx/dy
          )/(dx*dy);
    } 

  }

  return;
}
void laplacian( Field * phi, Constant constant,
    double * tmp
    )
{
  double dx = constant.dx, dy = constant.dy, dz =constant.dz;
  int i,j, l, m;
  int N = phi->N;
  int N_x = phi->N_x;
  int N_y = phi->N_y;
  for(i=0;i<N;i++)
    tmp[i] = 0.0;

  for(j = 1;j<N_y-1;j++){
    l= j*N_x;
    for(i = l+1 ;i< l + N_x-1;i++){
      tmp[i] = ((phi->bc[EAST]==NONE) ? phi->val[EAST] : ((phi->bc[EAST]==DIRICHLET) ? 2.0*phi->val[EAST]-phi->val[P] : phi->val[P])
            + (phi->bc[WEST]==NONE) ? phi->val[WEST] : ((phi->bc[WEST]==DIRICHLET) ? 2.0*phi->val[WEST]-phi->val[P] : phi->val[P])
            -2.0*phi->val[P])*dy/dx 
          + ( (phi->bc[NORTH]==NONE)? phi->val[NORTH]: ((phi->bc[NORTH]==DIRICHLET) ? 2.0*phi->val[NORTH]-phi->val[P]:phi->val[P])  
            + (phi->bc[SOUTH]==NONE)? phi->val[SOUTH]: ((phi->bc[SOUTH]==DIRICHLET) ? 2.0*phi->val[SOUTH]-phi->val[P]:phi->val[P]) 
            -2.0*phi->val[P])*dx/dy
          ;
    }
  }

return;
}
