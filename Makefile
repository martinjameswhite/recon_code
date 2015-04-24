CC     = g++
OPT    = -funroll-loops -O # -fopenmp



recon: recon.o io.o multigrid.o grid.o shift.o
	$(CC) $(OPT) recon.o io.o multigrid.o grid.o shift.o -o recon -lm

recon.o: recon.cpp global.h lcdm.h multigrid.h
	$(CC) $(OPT) -c recon.cpp
io.o: io.cpp global.h lcdm.h
	$(CC) $(OPT) -c io.cpp
multigrid.o: multigrid.cpp global.h lcdm.h multigrid.h
	$(CC) $(OPT) -c multigrid.cpp
grid.o: grid.cpp global.h lcdm.h
	$(CC) $(OPT) -c grid.cpp
shift.o: shift.cpp global.h lcdm.h
	$(CC) $(OPT) -c shift.cpp


.(PHONY) clean:
	rm -f *.o
