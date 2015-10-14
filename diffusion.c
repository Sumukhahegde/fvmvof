#include<stdio.h>                                                               
#include<stdlib.h>                                                              
#include<math.h>                                                                
#include<string.h>  
#include<stdbool.h>                                                             
#include "fvmvof.h"

void diffusion( field_variable * phi, double nu,
    double * tmp
    )
{

  int i, l, m;
  int N = u_x->N;
  int N_x = u_x->N_x;
  int N_y = u_x->N_y;

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

  /*    if(u_x->bc[east] != NONE )
        u_x_e = 2.0*(u_x->val[east]*abs(2-u_x->bc[east]) + u_x->val[i]*abs(1-u_x->bc[east])) - u_x->val[i];
      if(u_x->bc[west] != NONE )
        u_x_w = 2.0*(u_x->val[west]*abs(2-u_x->bc[west]) + u_x->val[i]*abs(1-u_x->bc[west])) - u_x->val[i];
      if(u_y->bc[north] != NONE )
        u_y_n = 2.0*(u_y->val[north]*abs(2-u_y->bc[north]) + u_y->val[i]*abs(1-u_y->bc[north])) - u_y->val[i];
      if(u_y->bc[south] != NONE )
        u_y_s = 2.0*(u_y->val[south]*abs(2-u_y->bc[south]) + u_y->val[i]*abs(1-u_y->bc[south])) - u_y->val[i];
*/
      if(phi->bc[east] != NONE )
        phi_e = 2.0*(phi->val[east]*abs(2-phi->bc[east]) + phi->val[i]*abs(1-phi->bc[east])) - phi->val[i];
      if(phi->bc[west] != NONE )
        phi_w = 2.0*(phi->val[west]*abs(2-phi->bc[west]) + phi->val[i]*abs(1-phi->bc[west])) - phi->val[i];
      if(phi->bc[north] != NONE )
        phi_n = 2.0*(phi->val[north]*abs(2-phi->bc[north]) + phi->val[i]*abs(1-phi->bc[north])) - phi->val[i];
      if(phi->bc[south] != NONE )
        phi_s = 2.0*(phi->val[south]*abs(2-phi->bc[south]) + phi->val[i]*abs(1-phi->bc[south])) - phi->val[i];

     tmp[i] += nu* ((phi_e+phi_w-2.0*phi->val[i])*dy/dx + (phi_n+phi_s-2.0*phi->val[i])*dx/dy) ;
    } else 
      tmp[i] = 0.0;

  }

//        printf("%lf %lf %lf %lf %lf\n", tmp[12*N_x + 12], phi->val[12*N_x + 12], dy, u_x->val[12*N_x + 12], u_y->val[12*N_x+12]);
return;
}
