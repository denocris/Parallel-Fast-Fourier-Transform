CC=mpicc

#FFTW_INC =/opt/local/include # On Mac Os X
#FFTW_LIB =/opt/local/lib # On Mac Os X

FFTW_INC =/u/shared/programs/x86_64/fftw/3.3.4/gnu/4.9.2/include
FFTW_LIB =/u/shared/programs/x86_64/fftw/3.3.4/gnu/4.9.2/lib

CFLAGS = -O3 -Wall -I$(FFTW_INC)
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
	rm -f *.o *.x *.dat *~
