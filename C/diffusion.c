/* Assignment: 
 * Parallelize the code, using fftw-mpi
 * This amount to 
 *   - distribute the data contained in diffusivity, conc, dconc, aux1, aux2 in the way fftw-mpi expects
 *   - modify the fftw calls in fftw-mpi in p_fftw_wrapper
 * You will need to modify the files
 *   - diffusion.c   
 *   - derivative.c
 *   - fftw_wrapper.c
 * In these files you will find some HINTs, make good use of them :)
 *
 * Created by G.P. Brandino, I. Girotto, R. Gebauer
 * Last revision: March 2016
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "utilities.h"
// HINT: include mpi and fftw3-mpi 
//       http://www.fftw.org/doc/MPI-Files-and-Data-Types.html#MPI-Files-and-Data-Types

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

int main(){

    // Dimensions of the system
    double L1 = 10., L2 = 10., L3 = 20.;
    // Grid size
    int n1 = 48, n2 = 48, n3 = 96;
    // time step for time integration
    double dt = 2.e-3; 
    // number of time steps
    int nstep = 101; 
    // Radius of diffusion channel
    double rad_diff = 0.7;
    // Radius of starting concentration
    double rad_conc = 0.6;
    double start, end;
  
    double *diffusivity, *conc, *dconc, *aux1, *aux2;
   
    int i1, i2, i3, ipol, istep, index;
  
    double f1conc, f2conc, f3conc, f1diff, f2diff, f3diff, fac, ss;
    double x1, x2 , x3, rr, r2mean;
    fftw_mpi_handler fft_h;

    /* 
     * Initializzation of the MPI environment 
     *
     */

    /*
     * initialize the fftw system and local dimension
     * as the value returned from the parallel FFT grid initializzation 
     *
     */
    init_fftw( &fft_h, n1, n2, n3, MPI_COMM_WORLD);

  
    /*
     * Allocate distribute memory arrays
     * HINT: the arrays need to be distributed, so you have to set the correct sizes 
     *       Use fftw_mpi_local_size_3d to calculate sizes
     *       http://www.fftw.org/doc/MPI-Data-Distribution-Functions.html
     *
     * More hints are reported into the fft_wrapper.c file where the fft_init_mpi is defined
     *
     */
    diffusivity = (double*)malloc(n1*n2*n3*sizeof(double));
    conc = (double*)malloc(n1*n2*n3*sizeof(double));
    dconc = (double*)malloc(n1*n2*n3*sizeof(double));
    aux1 = (double*)malloc(n1*n2*n3*sizeof(double));
    aux2 = (double*)malloc(n1*n2*n3*sizeof(double));


    /* 
     * Define the diffusivity inside the system and 
     * the starting concentration
     *
     * ss is to integrate (and normalize) the concentration
     *
     */   

    ss = 0.0;

    for (i3 = 0; i3 < n3; ++i3)
      {  
	x3=L3*((double)i3)/n3;
	f3diff = exp( -pow((x3-0.5*L3)/rad_diff,2));
	f3conc = exp( -pow((x3-0.5*L3)/rad_conc,2));
	
	
        for (i2 = 0; i2 < n2; ++i2)
	  {
            x2=L2*((double)i2)/n2;
            f2diff = exp( -pow((x2-0.5*L2)/rad_diff,2));
            f2conc = exp( -pow((x2-0.5*L2)/rad_conc,2));
	    
	
	    for (i1 = 0; i1 < n1; ++i1)
	      {
		x1=L1*((double)i1)/n1;
		f1diff = exp( -pow((x1-0.5*L1)/rad_diff,2));
		f1conc = exp( -pow((x1-0.5*L1)/rad_conc,2));
		
		index = index_f(i1, i2, i3, n1, n2, n3);
		diffusivity[index]  = MAX( f1diff * f2diff, f2diff * f3diff);
		conc[index] = f1conc * f2conc * f3conc;
		ss += conc[index]; 
		
	      }   
	  }
      }
    
    /*
     * HINT: The parallel version of  the output routines is provided in the mpi_output_routines folder
     *
     */
    plot_data_2d("diffusivity", n1, n2, n3, 1, diffusivity);
    plot_data_2d("diffusivity", n1, n2, n3, 2, diffusivity);
    plot_data_2d("diffusivity", n1, n2, n3, 3, diffusivity);
    
    
    fac= L1*L2*L3/(n1*n2*n3);
  
    /*
     * Now normalize the concentration
     *
     * HINT:the global ss must be computed and propagated to all processes 
     *      ss = 1.0/(ss*fac);
     *
     */
    ss = 1.0/(ss*fac);
    for (i1=0; i1< n1*n2*n3; ++i1)
      conc[i1]*=ss;
      
    /*
     * Now everything is defined: system size, diffusivity inside the system, and
     * the starting concentration
     *
     * Start the dynamics
     *
     */
    start = seconds();
    for (istep = 1; istep <= nstep; ++istep)
      {
        for (i1=0; i1< n1*n2*n3; ++i1)
	  dconc[i1] = 0.0;
        for (ipol =1; ipol<=3; ++ipol )
	  {
	    derivative(&fft_h, n1, n2, n3, L1, L2, L3, ipol, conc, aux1);
	    for (i1=0; i1< n1*n2*n3; ++i1)
	      {
                aux1[i1] *= diffusivity[i1];
	      }
	    derivative(&fft_h, n1, n2, n3, L1, L2, L3, ipol, aux1, aux2);
            // summing up contributions from the three spatial directions
            for (i1=0; i1< n1*n2*n3; ++i1)
	      dconc[i1] += aux2[i1];
	  } 
        for (i1=0; i1< n1*n2*n3; ++i1)
	  conc[i1] += dt*dconc[i1];
	
        if (istep%30 == 1)
	  {
            // Check the normalization of conc
            ss = 0.;
            r2mean = 0.;
            // HINT: the conc array is distributed, so only a part of it is on each processor
            for (i3 = 0; i3 < n3; ++i3)
	      {
                x3=L3*((double)i3)/n3 - 0.5*L3;
                for (i2 = 0; i2 < n2; ++i2)
		  {
                    x2=L2*((double)i2)/n2 - 0.5*L2;
                    for (i1 = 0; i1 < n1; ++i1)
		      {
			x1=L1*((double)i1)/n1 - 0.5*L1;
			rr = pow( x1, 2)  + pow( x2, 2) + pow( x3, 2);
			index = index_f(i1, i2, i3, n1, n2, n3); 
			ss += conc[index]; 
			r2mean += conc[index]*rr;
		      }   
		  }
	      }

            /*
	     * HINT: global values of ss and r2mean must be globally computed and distributed to all processes
	     *
	     */
            ss *= fac;
            r2mean *= fac;
            end = seconds();
            printf(" %d %17.15f %17.15f Elapsed time per iteration %f \n ", istep, r2mean, ss, (end-start)/istep);
            // HINT: Use parallel version of output routines
            plot_data_2d("concentration", n1, n2, n3, 2, conc);
            plot_data_1d("1d_conc", n1, n2, n3, 3, conc);
	  }
	
      } 
    
    close_fftw(&fft_h);
    free(diffusivity);
    free(conc);
    free(dconc);
    free(aux1);
    free(aux2);

    /*
     * Finalize the MPI environment
     */

    return 0;
} 
