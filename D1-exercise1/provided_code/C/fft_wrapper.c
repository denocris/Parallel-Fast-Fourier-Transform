/* Assignement:
 * Here you have to modify the includes, the array sizes and the fftw calls, to use the fftw-mpi
 *
 * Regarding the fftw calls. here is the substitution 
 * fftw_plan_dft_3d -> fftw_mpi_plan_dft_3d
 * ftw_execute_dft  > fftw_mpi_execute_dft 
 * use fftw_mpi_local_size_3d for local size of the arrays
 * 
 * Created by G.P. Brandino, I. Girotto, R. Gebauer
 * Last revision: March 2016
 */ 

#include <complex.h>
#include <fftw3.h>
#include <stdbool.h>
#include <string.h>
#include "utilities.h"

double seconds(){

  /* 
   * Return the second elapsed since Epoch (00:00:00 UTC, January 1, 1970) 
   *
   */
  struct timeval tmp;
  double sec;

  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;

  return sec;
}

/* 
 * Index linearization is computed following row-major order.
 * For more informtion see FFTW documentation:
 * http://www.fftw.org/doc/Row_002dmajor-Format.html#Row_002dmajor-Format
 *
 */
int index_f ( int i1, int i2, int i3, int n1, int n2, int n3)
{
  return n3*n2*i1 + n3*i2 + i3; 
}


void init_fftw(fftw_mpi_handler *fft, int n1, int n2, int n3, MPI_Comm mpi_comm)
{
  /*
   * Call to fftw_mpi_init is needed to initialize a parallel enviroment for the fftw_mpi
   * See also: http://www.fftw.org/doc/MPI-Initialization.html
   */

  /*
   *  Allocate a distributed grid for complex FFT using aligned memory allocation
   *  See details here:
   *  http://www.fftw.org/fftw3_doc/Allocating-aligned-memory-in-Fortran.html#Allocating-aligned-memory-in-Fortran
   *  HINT: the allocation size is given by fftw_mpi_local_size_3d (see also http://www.fftw.org/doc/MPI-Data-Distribution-Functions.html)
   *
   */
  fft->fftw_data = ( fftw_complex* ) fftw_malloc( fft->local_size_grid * sizeof( fftw_complex ) );

  /*
   * Create an FFTW plan for a distributed FFT grid
   * Use fftw_mpi_plan_dft_3d: http://www.fftw.org/doc/MPI-Plan-Creation.html#MPI-Plan-Creation
   *
   */
  fft->fw_plan = NULL;
  fft->bw_plan = NULL;

}

void close_fftw(fftw_mpi_handler *fft)
{
    fftw_destroy_plan(fft->bw_plan);
    fftw_destroy_plan(fft->fw_plan);
    fftw_free(fft->fftw_data);
}


/* 
 * This subroutine uses fftw to calculate 3-dimensional discrete FFTs.
 * The data in direct space is assumed to be real-valued
 * The data in reciprocal space is complex. 
 * direct_to_reciprocal indicates in which direction the FFT is to be calculated
 * 
 * Note that for real data in direct space (like here), we have
 * F(N-j) = conj(F(j)) where F is the array in reciprocal space.
 * Here, we do not make use of this property.
 * Also, we do not use the special (time-saving) routines of FFTW which
 * allow one to save time and memory for such real-to-complex transforms.
 *
 * f: array in direct space
 * F: array in reciprocal space
 * 
 * F(k) = \sum_{l=0}^{N-1} exp(- 2 \pi I k*l/N) f(l)
 * f(l) = 1/N \sum_{k=0}^{N-1} exp(+ 2 \pi I k*l/N) F(k)
 * 
 */
void fft_3d(fftw_mpi_handler* fft, int n1, int n2, int n3, double *data_direct, fftw_complex* data_rec, bool direct_to_reciprocal)
{
    double fac;
    int i;
    
    // Now distinguish in which direction the FFT is performed
    if ( direct_to_reciprocal)
      {
	for(i = 0; i < n1*n2*n3; i++)
	  {
	    fft->fftw_data[i]  = data_direct[i] + 0.0 * I;
	  } 
	
	fftw_execute_dft(fft->fw_plan, fft->fftw_data, fft->fftw_data);

	memcpy(data_rec, fft->fftw_data, n1*n2*n3*sizeof(fftw_complex)); 
      }
    else
      {
	memcpy(fft->fftw_data, data_rec, n1*n2*n3*sizeof(fftw_complex));
	  
	fftw_execute_dft(fft->bw_plan, fft->fftw_data, fft->fftw_data);
	
	fac = 1.0 / ( n1 * n2 * n3 );
	
	for( i = 0; i < n1 * n2 * n3; ++i )
	  {
	    data_direct[i] = creal(fft->fftw_data[i])*fac;
	  }
      }
}

