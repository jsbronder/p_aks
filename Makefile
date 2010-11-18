#
# Makefile.par
#
# For use to compile the parallel version
# of the AKS Primality Test
#
PAR_CC			= mpicc
PAR_LD			= $(PAR_CC)
PAR_OFILES 		= p_aks.o poly.o gmp_mpi.o
CFLAGS			= -Wall -g 
LDFLAGS			= -lgmp


p_aks: $(PAR_OFILES)
	$(PAR_CC) $(CFLAGS) -DUSE_MPI $(LDFILES) -o p_aks $(PAR_OFILES) $(LIBS) $(LDFLAGS)

clean: 
	rm -f *.o

.c.o: aks.h
	$(PAR_CC) -c $(INC) -DUSE_MPI $(CFLAGS) $< 


