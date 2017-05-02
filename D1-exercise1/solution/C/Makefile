CC=mpicc
CFLAGS = -O0 -g -Wall -I./ -I$(FFTW_INC) -I/u/shared/programs/x86_64/openmpi/1.8.3/gnu/4.9.2/torque/include/ 
LIBS = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3 -lm		
EXE = diffusion.x
SRC = diffusion.c plot_data.c fft_wrapper.c derivative.c 
HEAD = utilities.h
OBJ = $(SRC:.c=.o)


all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $^ $(LIBS) -o $@

%.o : %.c $(HEAD)
	$(CC) $(CFLAGS) -c $<

flush:
	rm -f *.dat

clean: 
	rm -f *.o *.x *~
