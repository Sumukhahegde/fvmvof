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

     tmp->val[i] += nu* ((phi_e+phi_w-2.0*phi->val[i])*dy/dx + (phi_n+phi_s-2.0*phi->val[i])*dx/dy) ;
    } else 
      tmp->val[i] = 0.0;
  }
return;
}
void laplacian( Field * phi, Constant constant,
    double * tmp
    )
{
  double dx = constant.dx, dy = constant.dy, dz =constant.dz;
  double nu = constant.nu;
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

     tmp[i] =  ((phi_e+phi_w-2.0*phi->val[i])*dy/dx + (phi_n+phi_s-2.0*phi->val[i])*dx/dy) ;
    } else 
      tmp[i] = 0.0;
  }
return;
}
