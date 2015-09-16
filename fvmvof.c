#include<stdio.h>
#include "fvmvof.h"

#define num_fpc 6
#define num_npc 8

int main(int argc, char *argv[])
{
  int nx =NX+3, ny = Ny+3, nz = 2; nout = 40, it;
  
  struct _domain
  {
    double *u, *v, *w, *x, *y, *z, *p, *div;
  } domain;


