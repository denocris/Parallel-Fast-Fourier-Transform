/* 
 * Created by G.P. Brandino, I. Girotto, R. Gebauer
 * Last revision: March 2016
 */

#ifndef _FFTW_UTLITIES_
#define _FFTW_UTLITIES_
#include <complex.h>
#include <fftw3.h>
#include <sys/time.h>
#include <stdbool.h>
#define pi 3.14159265358979323846

typedef struct 
{
    fftw_plan fw_plan; 
    fftw_plan bw_plan;
    fftw_complex *fftw_data;

} fftw_handler;



double seconds();
inline int index_f ( int i1, int i2, int i3, int n1, int n2, int n3);


void plot_data_1d( char* name, int n1, int n2, int n3, int dir, double* data);
void plot_data_2d( char* name, int n1, int n2, int n3, int dir, double* data);
void init_fftw(fftw_handler *fft, int n1, int n2, int n3);
void close_fftw(fftw_handler *fft);

void derivative(fftw_handler* fft,int n1, int n2, int n3, double L1, double L2, double L3, int ipol, double* data, double* deriv);
void fft_3d(fftw_handler* fft, int n1, int n2, int n3, double *data_direct, fftw_complex* data_rec, bool direct_to_reciprocal);
 


#endif
