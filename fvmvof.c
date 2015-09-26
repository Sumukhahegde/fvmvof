#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "fvmvof.h"
/* Finite Volume Method */
/* Code to solve NS. First step, solve poisson equation in 2D */
/* Define variables */


int main(int argc, char *argv[])
{

/* define grid sizes  - cell centered
 * allocate arrays for pos, vel 
 * allocate RHS of poisson
 * set BCs
 *
 * write function for compute_Ax
 *
 * call bicgstab - which calls compute_Ax
 * we get p 
 *write p and the grid to file */
  double l_x = 1.0;
  double l_y = 1.0;
  N_cells_x = 10 + 2;
  N_cells_y = 10 + 2;
  N_cells_z = 1;
  N_cells = N_cells_x * N_cells_y * N_cells_z ;

  u_x = malloc(N_cells*sizeof(double));
  u_y = malloc(N_cells*sizeof(double));
  u_z = malloc(N_cells*sizeof(double));
/*r_x = malloc(N_cells*sizeof(double));
  r_y = malloc(N_cells*sizeof(double));
  r_z = malloc(N_cells*sizeof(double));*/
  p = malloc(N_cells*sizeof(double));
  rho = malloc(N_cells*sizeof(double));
  omega_z = malloc(N_cells*sizeof(double));
  omega_x = malloc(N_cells*sizeof(double));
  omega_y = malloc(N_cells*sizeof(double));

  dx = l_x / (N_cells_x-2);
  dy = l_y / (N_cells_y-2);
  dz = 1.0;
  // only for poisson ;
  
  b = malloc(N_cells*sizeof(double));
  
  int i;

  for(i=0;i<N_cells;i++){
    b[i] = 0.0;
    p[i] = 0.0;
  }
  //Boundary conditions are: xmin: 50 ymin: 0 xmax: 50 ymax: 100 
  double p_bc_W = 50.0, p_bc_E=50.0, p_bc_S=0.0, p_bc_N=100.0;
  // generate rhs
/*  for(i=0;i<N_cells;i++){
    p = i % N_cells_x;
    q = (int) i/N_cells_x;
  if(p==0)
    b[i] =    ;
  if(q==0)
    b[i] += - p[i]* a_s  ;
  if(p==0)
    b[i] += - p[i]* a_w  ;
  if(p==0)
    b[i] += - p[i]* a_w  ;
  }*/

  int test  = solve();


  
  return 0;
}
int solve() {
	int i, j ;
	double rhoj_Minus, alphaj, omegaj, rhoj, betaj, H1, H2 ;
        double norm, BICGEPS = 1.0E-12; 
  	int BICG_ITER ;
	bool STOP = false ;

	// Start BICGSTAB iterations
	// set initial solution vector x_0 = (Uj, Vj) 
	ComputeAX() ;
	// Initial vector r_0 = b - Ax_0, and r0* = r_0
        for(i = 0 ; i < Nx ; i++) {
                for(j = 0 ; j < Ny ; j++) {
			Uj[i][j] = u[i][j] ;
			rj[i][j] = f[i][j] - Temp[i][j] ; 
        	        r0_star[i][j] = rj[i][j] ;

		//	if( (i == ip) && (j == jp) ) { Uj[i][j] = 0.0 ; rj[i][j] = 0.0 ; r0_star[i][j] = 0.0 ; }
		}
	}

	BICG_ITER = 0 ; norm = 0.0 ;
	do {
                // compute rhoj = (r0, r0*)
                rhoj = 0.0 ;
                for(i = 0 ; i < Nx ; i++) {
                	for(j = 0 ; j < Ny ; j++) {
				rhoj += rj[i][j]*r0_star[i][j] ; 
			}
		}
		if( sqrt(rhoj/double((Nx)*(Ny))) < BICGEPS ) STOP = true ;
		else {
			if( BICG_ITER == 0 ) {
				for(i = 0 ; i < Nx ; i++) {
               				for(j = 0 ; j <= Ny ; j++) pj[i][j] = rj[i][j]; // p0 = r0 
				}
			} else {
				betaj = (rhoj/rhoj_Minus)*(alphaj/omegaj) ;
				for(i = 0 ; i < Nx ; i++) {
                			for(j = 0 ; j < Ny ; j++) {
                       				pj[i][j] = rj[i][j] + betaj*(pj[i][j] - omegaj*Var[i][j]);
					}
				}
			}
			// Solve for Upstar, Vpstar from Ku* = u...., where K is the preconditioning matrix
			for(i = 0 ; i < Nx ; i++) {
               			for(j = 0 ; j <= Ny ; j++) {
					// No preconditioning
					pstar[i][j] = pj[i][j] ;
				}
			}

			// compute vj = A*pstar
	        	for(i = 0 ; i < Nx ; i++) {
        	        	for(j = 0 ; j < Ny ; j++) {
			//		if( (i == ip) && (j == jp) ) u[i][j] = 0.0 ;
					 u[i][j] = pstar[i][j] ;
				}
			}
			ComputeAX() ;

			for(i = 0 ; i < Nx ; i++) {
        	        	for(j = 0 ; j < Ny ; j++) {
					Var[i][j] = Temp[i][j] ;
				}
			}
			H1 = 0.0 ;
                	for(i = 0 ; i < Nx ; i++) {
                		for(j = 0 ; j < Ny ; j++) {
					H1 += Var[i][j]*r0_star[i][j] ;
				}
			}
                	alphaj = rhoj/H1 ;

                	// find sj
			for(i = 0 ; i < Nx ; i++) {
                		for(j = 0 ; j < Ny ; j++) {
					sj[i][j] = rj[i][j] - alphaj*Var[i][j] ;
				}
			}
			// Solve for Upstar, Vpstar from Ku* = u...., where K is the preconditioning matrix
			for(i = 0 ; i < Nx ; i++) {
               			for(j = 0 ; j < Ny ; j++) {
					// No preconditioning
					sstar[i][j] = sj[i][j] ;
				}
			}
			norm = 0.0 ;
        	        for(i = 0 ; i < Nx ; i++) {
                		for(j = 0 ; j < Ny ; j++) {
					norm += sstar[i][j]*sstar[i][j] ;
				}
			}
			norm = sqrt(norm/double((Nx)*(Ny))) ;
			if( norm < BICGEPS) {
				STOP = true ; // if ||s||_2 is small x_i = x_{i-1} + alphai*p_i
				for(i = 0 ; i < Nx ; i++) {
                			for(j = 0 ; j < Ny ; j++) {
//						if( (i == ip) && (j == jp) ) Uj[i][j] = 0.0 ;
						 Uj[i][j] += alphaj*pstar[i][j] ;
					}
				}
			} else {

				// compute t = As
		        	for(i = 0 ; i < Nx ; i++) {
        		        	for(j = 0 ; j < Ny ; j++) {
//						if( (i == ip) && (j == jp) ) u[i][j] = 0.0 ;
						 u[i][j] = sstar[i][j] ;
					}
				}
				ComputeAX() ;
	 	               	H1 = H2 = 0.0 ;
        		        for(i =01 ; i < Nx ; i++) {
        		        	for(j = 0 ; j <= Ny ; j++) {
						H1 += Temp[i][j]*sj[i][j] ;
						H2 += Temp[i][j]*Temp[i][j] ;
					}
				}	
	                	omegaj = H1/H2;
		
	                	// find xj 
	        	        norm = 0.0 ;
        	        	for(i = 0 ; i < Nx ; i++) {
        		        	for(j = 0 ; j < Ny ; j++) {
						H1 = (alphaj*pstar[i][j] + omegaj*sstar[i][j]) ;
//						if( (i == ip) && (j == jp) ) Uj[i][j] = 0.0 ;
						 Uj[i][j] += H1 ;
	        	                	norm += H1*H1 ;
					}
				}
                		norm = sqrt(norm/double((Nx)*(Ny))) ;
	
				if(norm < BICGEPS) STOP = true ;
        		        // find rjplusone
                		for(i = 0 ; i < Nx ; i++) {
        		        	for(j = 0 ; j < Ny ; j++) {
						rj[i][j] = sj[i][j] - omegaj*Temp[i][j];
					}
				}
				rhoj_Minus = rhoj ;
			}
		}
                BICG_ITER++;
		if(BICG_ITER%100 == 0) cout << BICG_ITER << "\t" << norm <<  endl;
        }while( (BICG_ITER < 10000) && (!STOP) ) ;
//        if(BICG_ITER > 1000-3) cout << BICG_ITER << "\t" << norm << endl;
//	cout << BICG_ITER << "\t" << norm << endl;
	for(i = 0 ; i < Nx ; i++) {
		for(j = 0 ; j < Ny ; j++) {
			u[i][j] = Uj[i][j] ;	
		}
	}
}

int solve()
{



}

