CC=g++

FFLAGS= -O3 -fopenmp -march=native

all: 3d2rewq.cpp
	$(CC) $(FFLAGS) 3d2rewq.cpp -o 3d2rewq

clean:
	rm -f 3d2rewq
