#include<stdio.h>                                                               
#include<stdlib.h>                                                              
#include<math.h>                                                                
#include<string.h>  
#include<stdbool.h>                                                             
#include "fvmvof.h"

void advection( Field * phi, Field * u_x, 
    Field * u_y, Field * u_z,
    Field * tmp,
    Constant constant
    )
{
 
  double dx = constant.dx, dy = constant.dy, dz =constant.dz;
  int i, l, m;
  int N = phi->N;
  int N_x = phi->N_x;
  //int N_y = u_x->N_y;

  double u_x_e, u_x_w; //u_x_s, u_x_n, 
  double u_y_s, u_y_n; // u_y_e, u_y_w;
  double phi_s, phi_n , phi_e, phi_w;

  for(i = 0;i<N;i++){
    if(phi->bc[i] == NONE ){
      l= i%N_x;                                                           
      m =(int) i/N_x;

      int south = (m-1)*N_x + l, north =(m+1)*N_x + l,              
          west =m*N_x + (l-1), east = m*N_x + (l+1);

      u_x_e = u_x->val[east];
      u_x_w = u_x->val[west];
      u_y_s = u_y->val[south];
      u_y_n = u_y->val[north];
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
      if(u_x->bc[east] != NONE )
        u_x_e = 2.0*(u_x->val[east]*abs(2-u_x->bc[east]) + u_x->val[i]*abs(1-u_x->bc[east])) - u_x->val[i];
      if(u_x->bc[west] != NONE )
        u_x_w = 2.0*(u_x->val[west]*abs(2-u_x->bc[west]) + u_x->val[i]*abs(1-u_x->bc[west])) - u_x->val[i];
      if(u_y->bc[north] != NONE )
        u_y_n = 2.0*(u_y->val[north]*abs(2-u_y->bc[north]) + u_y->val[i]*abs(1-u_y->bc[north])) - u_y->val[i];
      if(u_y->bc[south] != NONE )
        u_y_s = 2.0*(u_y->val[south]*abs(2-u_y->bc[south]) + u_y->val[i]*abs(1-u_y->bc[south])) - u_y->val[i];

//      tmp[i] =( (phi->val[i]*u_x->val[i] - phi_w*u_x_w)/dx +(phi->val[i]*u_y->val[i] -phi_s*u_y_s)/dy);
      //tmp[i] =(u_x_w*(phi->val[i] - phi_w)/dx + u_y_s*(phi->val[i] -phi_s)/dy);
//      tmp[i] = dy*(u_x_e +u_x->var[i])*0.5*(6./8. * phi->val[i] + 3./8.*phi_e - 1./8.0 * phi_w) - (u_x_w+u_x->var[i])*0.5*()
  //    tmp[i] =(u_x_w*(phi->val[i] - phi_w)/dx + u_y_s*(phi->val[i] -phi_s)/dy);
     //tmp[i] =0.5* ((phi_e*u_x_e - phi_w*u_x_w)*dy +(phi_n*u_y_n -phi_s*u_y_s)*dx);
   //  tmp[i] = ((phi->val[i]*u_x_e - phi_w*u_x_w)*dy +(phi->val[i]*u_y_n -phi_s*u_y_s)*dx);
     tmp->val[i] += - 0.5* ((phi_e*u_x_e - phi_w*u_x_w)*dy +(phi_n*u_y_n -phi_s*u_y_s)*dx);
//      tmp[i] =0.5* (u_x->val[i]*(phi_e - phi_w)/dx +u_y->val[i]*(phi_n -phi_s)/dy);
       
      // must be plus equal to 
    } else 
      tmp->val[i] = 0.0;

  }

//        printf("%lf %lf %lf %lf %lf\n", tmp[12*N_x + 12], phi->val[12*N_x + 12], dy, u_x->val[12*N_x + 12], u_y->val[12*N_x+12]);
return;
}
