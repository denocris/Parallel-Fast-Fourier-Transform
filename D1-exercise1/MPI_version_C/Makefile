CC=mpicc
FFTW_INC =/opt/local/include
FFTW_LIB =/opt/local/lib

CFLAGS = -O0 -g -Wextra -Wall -I$(FFTW_INC)

LIBS = -L$(FFTW_LIB) -lfftw3 -lfftw3_mpi -lm
EXE = diffusion.x
SRC = diffusion.c plot_data.c fft_wrapper.c derivative.c
HEAD = utilities.h
OBJ = $(SRC:.c=.o)


all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $^ $(LIBS) -o $@

diffusion.o: plot_data.o

%.o : %.c $(HEAD)
	$(CC) $(CFLAGS) -c $<

flush:
	rm -f *.dat

clean:
	rm -f *.o *.x *~
