CC = gcc
RM = rm -f
#RM = del

CFLAGS = -Wall \
	 -O3 -g -funroll-loops -march=native -mtune=native \
         -std=c99
LIBS   = -lm
#LIBS   = -lm -lwinmm
         
CFLAGS += -DUSE_OPENMP -fopenmp

#CFLAGS += -DUSE_BLAS
#LIBS   += -lblas -lgfortran -llapack

all: p6_wave2d

basic.o: basic.c basic.h
	$(CC) -c $(CFLAGS) $< -o $@ 
	
miniblas.o: miniblas.c miniblas.h basic.h
	$(CC) -c $(CFLAGS) $< -o $@ 
	
gridfunction.o: gridfunction.c miniblas.h basic.h
	$(CC) -c $(CFLAGS) $< -o $@ 

iteration.o: iteration.c
	$(CC) -c $(CFLAGS) $< -o $@ 
	
p6_wave2d: p6_wave2d.o basic.o miniblas.o gridfunction.o iteration.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
	$(RM) *.o
	$(RM) p6_wave2d

