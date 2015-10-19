#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>
#include "../fvmvof.h"

static Field * allocate_field( int N_x, int  N_y, int N_z);
int main(){
  int i, j, l, m, n;

  int N_x = 1000; int N_y = 10000;
/*  Field * phi = allocate_field(N_x,N_y,1);
  Field * tmp = allocate_field(N_x,N_y,1);
  double nu = 0.1;
  double dx = 0.5; double dy =0.9;
  int N = N_x * N_y;
*/
  double nu = 0.1;
//  for(j = 1;j<N_y-1;j++){
  //  l= j*N_x;
   // double * __restrict__ tmp1 = &tmp->val[l+1];
   // double * __restrict__ tmp2 = &tmp->val[l+2];
    //for(i = l+1 ;i<l+N_x-2;i++){
double  *tmp1,  *tmp2;
tmp1 = malloc(N_x*sizeof(double));
tmp2 = malloc(N_x*sizeof(double));
    for(i = 0 ;i<N_x-1;i++){
  //    tmp->val[i] = nu* ((phi->val[EAST]+phi->val[WEST]-2.0*phi->val[P])*dy/dx + (phi->val[NORTH]+phi->val[SOUTH]-2.0*phi->val[P])*dx/dy) ;
      tmp1[i] = nu* tmp2[i] * 5.0*tmp2[i];
    //}
  }

  asm volatile("": "+m"(tmp1), "+m"(tmp2));
  return 0;
  //asm volatile("": "+m"(a), "+m"(b), "+m"(c)::"memory");
}
static Field * allocate_field( int N_x, int  N_y, int N_z)
{
  Field * phi;
  phi      = malloc(sizeof(Field));
  phi->grid = CENTERED;
  phi->N_x = N_x;
  phi->N_y = N_y;
  phi->N_z = N_z;
  phi->N   = N_x*N_y*N_z;
  phi->val = malloc(phi->N * sizeof(double));
  phi->bc  = malloc(phi->N * sizeof(BC_type));
  return phi;
}

