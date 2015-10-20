#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>
#include "../fvmvof.h"

//static Field * allocate_field( int N_x, int  N_y, int N_z);

int main(){
  int i, j, l, m, n;

  int N_x = 1000; int N_y = 10000;
  double * x, * y;

  x = malloc(N_x*sizeof(double));
  y = malloc(N_x*sizeof(double));
  double * tmp1 = x;
  double * tmp2 = y;

  for(i = 1 ;i<100;i++){
    tmp1[i] = 5.0*tmp2[i] +  tmp2[i+129] ;
//    x[i] = 5.0*y[i] +  y[i+129] ;
    
  }
  asm volatile("": "+m"(tmp1), "+m"(tmp2));
  return 0;
}

