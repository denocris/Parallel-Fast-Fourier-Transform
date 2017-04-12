F90           = mpif90 
LIBS          = -lfftw3_mpi -lfftw3		
F90FLAGS      = -O3 -I$(FFTW_INC)  
LDFLAGS       = -L$(FFTW_LIB)
OBJS         = p_errore.o p_plot_data_1d.o  p_plot_data_2d.o p_fft_wrapper.o p_derivative.o 


all : diffusion.x 

%.o : %.f90
	$(F90) $(F90FLAGS) -c $<

diffusion.x : $(OBJS) p_diffusion.o 
	$(F90) $(LDFLAGS) $^ -o $@ $(LIBS) 

flush : 
	/bin/rm *.dat *.png *.gif

clean :
	- /bin/rm -f *.x *.o *.mod *.a *.so
