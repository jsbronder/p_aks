#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "gmp.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

int snd_cnt;
int recv_cnt;
typedef unsigned long int ulong;

struct polynomial {
    unsigned long int   degree;
    mpz_t *             coeff;
};

/*
 * Polynomial functions,
 * see poly.c
 *
 * All polynomial functions that
 * return poly_t create a new
 * polynomial inside, hence we also
 * have booleans for forcing the
 * freeing of memory assoiciated
 * with the arguments.
 *
 * For instance if computing
 * f(x) = g(x)*f(x), we would
 * call 
 * f = poly_mult( g, f, 0, 1 );
 * so that the original f is freed.
 */


/* 
 * Allocates a struct polynomial and
 * enough mpz_t * to store the coeffients
 */
struct polynomial * poly_alloc(unsigned long int degree);


/* 
 * Free a polynomial.
 */
void poly_free(struct polynomial *f);


/*
 * The following fuctions rely upon the global
 * variables n and r.  We are able to make assumptions
 * about the degree and binary bits used to store the
 * coeffients which allows for faster calculation.
 *
 * This means these functions will not work correctly
 * outside of our AKS implementation.
 */


/* 
 * Returns the product of a and b, uses single point evaluation
 * to compute.
 */
struct polynomial * poly_mult(
    struct polynomial *a,
    struct polynomial *b,
    short int freea,
    short int freeb,
    mpz_t n,
    unsigned long int r,
    unsigned long int log_b2_n);


/* Reduce f (mod X^r-1, n) */
struct polynomial * poly_mod(
    struct polynomial *f,
    short int free,
    mpz_t n,
    unsigned long int r);


/* Output:  a_0 + a_1*X + ... a_d*X^d */
void poly_print(struct polynomial *f);


/* Calculate f^n (mod X^r-1, n) */
struct polynomial * poly_expmod(
    struct polynomial *f,
    short int free,
    mpz_t n,
    unsigned long int r,
    unsigned long int log_b2_n);


/* 1 if f = g, 0 if not */
int poly_cmp(struct polynomial *f, struct polynomial *g);

#ifdef USE_MPI
/*
 * MPI Functions
 */


/*
 * z should be allocated and set before hand.
 * Return code is -1  on failure to
 * create a buffer to hold z.
 */
int mpi_mpz_send(mpz_t z, int dest, int tag, MPI_Comm comm );

int mpi_mpz_recv(mpz_t z, int src, int tag, MPI_Comm comm, MPI_Status *status);

int mpi_mpz_bcast(mpz_t z, int root, MPI_Comm comm, int id );
#endif
