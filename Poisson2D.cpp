//  Sample code for solving Poisson equation with homogeneous Neumann BC.
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<cstring>
#include<math.h>

using namespace std ;

double const PI = 4.0*atan(1.0); // Value of PI

double *xx, *yy ;
double **u,**f, Exact ; // u will be the array for final solution, 
		        // f will store the right hand side of Poisson's equation 
		        //
double L2Error = 0.0, MaxError = 0.0 ;
double TOL = 1E-12 ; // Error tolerance for stopping iterations
enum BCs{
  NEUMANN,
  DIRICHLET
};

BCs BC=NEUMANN;
int Nx = 100, Ny = 100, iter, ip, jp ;
double hxx, hyy ; // grid size in x and y


double **r0_star, **sj, **rj, **pj, **pstar, **sstar ; // Variables for BICGSTAB
double **Temp, **Uj, **Var  ;

// Allocate two-dimensional arrays
void Allocate_2D_R(double**& m, int d1, int d2) {
        m=new double* [d1];
        for (int i=0; i<d1; ++i) {
                m[i]=new double [d2];
                for (int j=0; j<d2; ++j)
                        m[i][j]=0.0;
        }
}

void ComputeAX() {
  int i, j ;
  double u_xx, u_yy;

  for(i = 0 ; i < Nx ; i++) {
    for(j = 0 ; j < Ny ; j++) {

      if(i==0){
        if(BC==NEUMANN)
          u_xx = (u[i+1][j] - 1.0*u[i][j]  )/(hxx*hxx) ;
        else if(BC==DIRICHLET)
          u_xx = (u[i+1][j] - 2.0*u[i][j] + 0.0 )/(hxx*hxx) ;
      }else if(i==Nx-1){
        if(BC==NEUMANN)
          u_xx = ( - 1.0*u[i][j] + u[i-1][j])/(hxx*hxx) ;
        else if(BC==DIRICHLET)
          u_xx = (0.0 - 2.0*u[i][j] + u[i-1][j])/(hxx*hxx) ;
      }else{
        u_xx = (u[i+1][j] - 2.0*u[i][j] + u[i-1][j])/(hxx*hxx) ;
      }
      
      if(j==0){
        if(BC==NEUMANN)
          u_yy = (u[i][j+1] - 1.0*u[i][j]  )/(hyy*hyy) ;
        else if(BC==DIRICHLET)
          u_yy = (u[i][j+1] - 2.0*u[i][j] + 0.0 )/(hyy*hyy) ;
      }else if(j==Ny-1){
        if(BC==NEUMANN)
          u_yy = ( - 1.0*u[i][j] + u[i][j-1])/(hyy*hyy) ;
        else if(BC==DIRICHLET)
          u_yy = (0.0 - 2.0*u[i][j] + u[i][j-1])/(hyy*hyy) ;
      }else{
        u_yy = (u[i][j+1] - 2.0*u[i][j] + u[i][j-1])/(hyy*hyy) ;
      }


      Temp[i][j] = u_xx + u_yy ; 
    }
  } 
}

double InnerProduct(double **x, double **y, int Nx, int Ny) {
	double temp = 0.0 ;
	for(int i = 0 ; i < Nx ; i++) {
		for(int j = 0 ; j < Ny ; j++) {
			temp += x[i][j]*y[i][j] ;
		}
	}
	return temp ;
}

void AllocateMemory() {
	xx = new double[Nx] ; yy = new double[Ny] ; 
	Allocate_2D_R(u,Nx,Ny) ; Allocate_2D_R(f,Nx,Ny) ;

	// Variables to be used in BICGSTAB

	Allocate_2D_R(r0_star,Nx,Ny) ; Allocate_2D_R(sj,Nx,Ny) ; Allocate_2D_R(rj,Nx,Ny) ; Allocate_2D_R(pj,Nx,Ny) ;
	Allocate_2D_R(pstar,Nx,Ny) ; Allocate_2D_R(sstar,Nx,Ny) ; Allocate_2D_R(Temp,Nx,Ny) ; Allocate_2D_R(Uj,Nx,Ny) ; Allocate_2D_R(Var,Nx,Ny) ;

}

void BICGSTABIter() {
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

int main() {
	int i, j ;

	cout.flags( ios::dec | ios::scientific );
	cout.precision(5);

	AllocateMemory() ;
    
	// set a uniform grid.
	hxx = (1.0 -0.0)/(Nx-1) ; hyy = (1.0-0.)/(Ny-1) ;
	for(i = 0 ; i < Nx ; i++) xx[i] = i*hxx ;	
	for(i = 0 ; i < Ny ; i++) yy[i] = i*hyy ;

	ip = 0 ; jp = 0 ; // index at which the value of u for Poisson Neumann equation is set equal to zero

	// set the right hand side in the interior
	for(i = 1 ; i < Nx ; i++) {
		for(j = 0 ; j < Ny ; j++) 
                  f[i][j] = -8.0*PI*PI*cos(2.0*PI*xx[i])*cos(2.0*PI*yy[j]) ;
	}

	// Initial guess....
	for(i = 0 ; i < Nx ; i++) {
		for(j = 0 ; j < Ny ; j++) u[i][j] = 0.0 ;
	}
//	u[ip][jp] = 0.0 ;
	
	BICGSTABIter() ;

	// Set the boundary values
        /*
	for(j = 0 ; j <= Np ; j++) {
		u[0][j] = u[1][j] ;
		u[Nr][j] = u[Nr-1][j] ;
	}*/

	// Write the final solution and compute errors.
	ofstream FileD("Neumann2drect.dat", ios::out) ;
	FileD.flags( ios::dec | ios::scientific );
	FileD.precision(10) ;
	if(!FileD) {cerr<< "Error: Output file Neumann.dat couldnot be opened.\n"; exit(1) ;}
	L2Error = MaxError = 0.0 ;

	FileD << "TITLE = Flow" << endl << "VARIABLES = x, y, u, Exact " << endl;
	FileD << "Zone T = Omega I = " << Ny << " J = " << Nx << endl ;

        for(i = 0 ; i < Nx ; i++) {
          for(j = 0 ; j < Ny ; j++) {

            Exact = cos(2.0*PI*xx[i])*cos(2.0*PI*yy[j]) ;

            L2Error += (u[i][j]-Exact)*(u[i][j]-Exact)/((Nx)*(Ny)) ;
            if(fabs(u[i][j]-Exact) > MaxError) MaxError = fabs(u[i][j]-Exact);

            FileD << xx[i] << "\t" << yy[j] << "\t" << u[i][j] << "\t" << Exact << endl ;
          }		
         // Exact = ( log(r[i]) - (1.0/(RLower + RUpper))*(r[i] - RLower*RUpper/r[i])  )*sin(phi[0]) ;
         // FileD << r[i]*cos(phi[0]) << "\t" << r[i]*sin(phi[0]) << "\t" << u[i][0] << "\t" << Exact << endl ;
        }
        FileD.close() ;
        L2Error = sqrt(L2Error) ; 
        cout << "Error for BICGSTAB method" << endl ;
        cout << " L2 Error : " << L2Error << "\t Max Error : " << MaxError << endl ;

        // -------------------------------------------------------------------------------------------------------------------------------------------
        // Free all the memory and exit
        delete[] xx ; delete [] yy ; 
        return 0;
}
