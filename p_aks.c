#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include "aks.h"
#include "mpi.h"

#define INPUT_BASE 10
#define PACK_SIZE  20
mpz_t    TERM;            /* Terminal sent to the Children */

static void print_syntax(int do_exit) {
    printf("Syntax:  p_aks <options>\n");
    printf("Options:\n");
    printf("\t-f <filename>:  Read integer from ASCII file.\n");
    printf("\t-d <number>:    Integer to test is number, limited by unsigned long.\n");
    printf("\t-h :            Print this help.\n");

    if (do_exit)
        exit(0);
    else
        return;
}

static int read_n_from_file(mpz_t n, char *filename) {
    FILE *in_file;

    in_file = fopen(filename, "r");
    if (!in_file) {
        fprintf(stderr, "open %s: %s\n", filename, strerror(errno));
        return 0;
    }

    if (!mpz_inp_str(n, in_file, INPUT_BASE)) {
        printf("Error reading n from %s in base %d.\n",
            filename, INPUT_BASE);
        fclose(in_file);
        return 0;
    }

    fclose(in_file);
    return 1;
}


/*
 * Get the upper bound on the
 * log base 2 of n.  Simply
 * find the number of bits used
 * to store it.
 */
inline unsigned long int get_log_b2(mpz_t n) {
    return (unsigned long int)mpz_sizeinbase(n, 2);
}


/*
 * Find r such that the order of n (mod r)
 * is greater then log(n)*log(n).
 *
 * We'll computer n^j for 1 \leq j \leq log(n)^2
 * and reduce modulo q for q > log(n)^2 until
 * we finally find a q such that none of these
 * powers is equilavent to 1.
 *
 */
static int find_r(mpz_t n, mpz_t r, unsigned long int log_b2_n) {
    unsigned long int j, gcd, ulong_r;
    unsigned long int log_sqrd = (log_b2_n)*(log_b2_n);
    int found = 0;
    mpz_t tmp;


    mpz_init_set_ui(r, log_sqrd);
    mpz_init2(tmp, log_sqrd);

    /*
     * if( found ) then we have
     * a power of n that is \equiv 1 \mod q
     * so we need to increase q and keep going.
     *
     */
    while (true) {
        found = 0;
        for (j = 1; j <= log_sqrd; j++) {
            mpz_powm_ui(tmp, n,  j, r);

            /*
             * If the first bit is a 1 and the
             * number of bits needed is 1, then
             * tmp \equiv 1 \mod q
             *
             */
            if (mpz_sgn( tmp )== 0
                    || (1 == mpz_tstbit(tmp, 0)
                        && 1 == mpz_sizeinbase(tmp, 2))){
                found = 1;
                break;
            }
        }
        if (!found)
            break;
        mpz_add_ui(r, r, 1);
    }
    /*
     * At this point r is correctly set to the first
     * value found st o(n) \mod r > log(n)^2
     *
     * Now we need to verify that:
     * 1 < (a,n) < n for a = 1, ... ,r
     *
     */
    ulong_r = mpz_get_ui( r );

    for (j = 1; j <= ulong_r; j++) {
        gcd = mpz_gcd_ui(tmp, n, j);
        if (gcd > 1 && mpz_cmp(n, tmp) > 0) {
            gmp_printf("(%Zd, %lu) = %lu\n", n, j, gcd);
            gmp_printf("%Zd is not prime.\n", n);
            return 0;
        }
    }
    mpz_clear(tmp);

    if (mpz_cmp(n, r) <= 0) {
        gmp_printf("%Zd <= r = %Zd\n", n, r);
        gmp_printf("%Zd is prime.\n", n);
        return 0;
    }

    return 1;
}


int test_polys (int nprocs, mpz_t n, mpz_t gmp_r, unsigned long int log_b2_n) {
    mpz_t upper_bound;
    mpz_t to_send;
    mpz_t ret;
    MPI_Status status;
    short int failure = 0;
    int i;
    int missed_children = 0;

    mpz_init(upper_bound);
    mpz_sqrt(upper_bound, gmp_r);
    mpz_mul_ui(upper_bound, upper_bound, log_b2_n);

    mpi_mpz_bcast(n, 0, MPI_COMM_WORLD, 0);
    mpi_mpz_bcast(gmp_r, 0, MPI_COMM_WORLD, 0);
    mpi_mpz_bcast(upper_bound, 0, MPI_COMM_WORLD, 0);

    mpz_init(to_send);
    mpz_init(ret);
    mpz_set_ui(to_send, 1);

    /*
     * Send out the initial coefficients
     * for the polynomials.  Each node
     * tests PACK_SIZE polynomials at a time.
     */
    for (i = 1 ; i < nprocs; i++) {
        mpi_mpz_send(to_send, i, 0, MPI_COMM_WORLD);
        mpz_add_ui(to_send, to_send, PACK_SIZE);
        if (mpz_cmp(to_send, upper_bound) > 0) {
            missed_children = i + 1;
            break;
        }
    }

    if (missed_children != 0) {
        for (i = missed_children; i < nprocs; i++) {
            mpi_mpz_send(TERM, i, 0, MPI_COMM_WORLD);
        }
        nprocs = missed_children;
    }


    /*
     * Case 1:
     * We have more nodes then test coefficients.
     * So we're done, we only need to check for
     * failure.
     */
    if (mpz_cmp( to_send, upper_bound) > 0) {
        for (i = 1; i < nprocs; i++) {
            mpi_mpz_recv(ret, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            if (mpz_sgn(ret) == 0)
                failure = 1;
            mpi_mpz_send(TERM, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
        }
        mpz_clear(to_send);
        mpz_clear(upper_bound);
        mpz_clear(ret);
        return !failure;
    }

    /*
     * Case 2:
     * We have lots of polynomials to test.
     * So we'll continue to check the input
     * and send out new coefficients until
     * we hit the upper bound
     */
    while (!failure && mpz_cmp(to_send, upper_bound) <= 0) {
        mpi_mpz_recv(ret, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        failure = (mpz_sgn(ret) == 0);
        if (!failure)
            mpi_mpz_send(to_send, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
        else
            mpi_mpz_send(TERM, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
        mpz_add_ui(to_send, to_send, PACK_SIZE);
    }

    /*
     * So either we have a failure, or we've hit the
     * upper bound.  Test for failure first and kill
     * all processes.
     *
     * On failure we have nprocs-1 nodes left to kill
     */
    if (failure) {
        for (i = 1; i < nprocs-1; i++) {
            mpi_mpz_recv(ret, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            mpi_mpz_send(TERM, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
        }
        mpz_clear(to_send);
        mpz_clear(ret);
        mpz_clear(upper_bound);
        return 0;
    }

    /*
     * We didn't fail in the loop, so we've hit
     * the upper bound.  We'll just accept the
     * last round of packets and test for
     * failure
     */
    for (i = 1; i < nprocs; i++) {
        mpi_mpz_recv(ret, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        mpi_mpz_send(TERM, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
        if (!mpz_sgn(ret))
            failure = 1;
    }
    mpz_clear(to_send);
    mpz_clear(ret);
    mpz_clear(upper_bound);
    return !failure;
}

int child_process(int id) {
    MPI_Status status;
    struct polynomial *f;
    struct polynomial *f_pow;
    mpz_t gmp_r;
    mpz_t n;
    mpz_t NOT_TERM;
    mpz_t a;
    mpz_t upper_bound;
    unsigned long int log_b2_n;
    unsigned long int r;
    int cnt;
    int failure;

    mpz_init(n);
    mpz_init(gmp_r);
    mpz_init(upper_bound);

    mpi_mpz_bcast(n, 0, MPI_COMM_WORLD, id);
    mpi_mpz_bcast(gmp_r, 0, MPI_COMM_WORLD, id);
    mpi_mpz_bcast(upper_bound, 0, MPI_COMM_WORLD, id);

    if (!mpz_sgn(n))
        return 0;

    r = mpz_get_ui(gmp_r);

    if (!mpz_sgn(gmp_r))
        return 0;

    log_b2_n = get_log_b2(n);

    f = poly_alloc(1);

    mpz_init(NOT_TERM);
    mpz_init(a);
    mpz_set_ui(NOT_TERM, 1);

    mpi_mpz_recv(a, 0, 0, MPI_COMM_WORLD, &status);

    while (mpz_sgn(a)) {
        failure = 0;
        for (cnt = 0; (cnt < PACK_SIZE && failure != 1 ); cnt ++) {
            mpz_set(f->coeff[1], a);
            f_pow = poly_expmod(f, 0, n, r, log_b2_n);
            if (poly_cmp( f, f_pow ) != 0) {
                poly_print(f);
                poly_print(f_pow);
                printf("%d:  Not congruent --> Not Prime\n", id);
                gmp_printf("%d:  Happened with (X+%Zd)\n", id, a );
                failure = 1;
                break;
            }
            mpz_add_ui(a, a, 1);
            if (mpz_cmp(a, upper_bound) > 0)
                break;
        }
        mpi_mpz_send(failure == 0 ? NOT_TERM : TERM, 0, 0, MPI_COMM_WORLD);
        mpi_mpz_recv(a, 0, 0, MPI_COMM_WORLD, &status);
    }
    return 0;
}


void kill_children(int nprocs) {
    int i;
    for (i = 0; i < 3; i++)
        mpi_mpz_bcast(TERM, 0, MPI_COMM_WORLD, 0);

    for (i = 1; i < nprocs; i++)
        mpi_mpz_send(TERM, i, 0, MPI_COMM_WORLD);

    return;
}

int main(int argc, char **argv) {
    int id, nprocs;
    double t1, t2, t3;

    /*
     * Initalization:
     * It is highly recommended that
     * we call MPI_Init before doing
     * much of anything.  That being
     * said, the amount of time spent
     * on steps preceeding test_polys()
     * is inconsequential by comparison,
     * so we don't actually begin using
     * MPI until then
     */

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    mpz_init(TERM);
    mpz_set_ui(TERM, 0);

    snd_cnt = 0;
    recv_cnt = 0;

    if (id)
        child_process( id );
    else {
        int got_n = 0;
        int opts;
        mpz_t n;
        mpz_t gmp_r;
        unsigned long int log_b2_n;
        int rc;

        t1 = MPI_Wtime();

        printf("Processors:\t%d\n", nprocs);
        printf("Pack Size:\t%d\n", PACK_SIZE);
        mpz_init( n );

        while (true) {
            opts = getopt(argc, argv, "f:hd:");
            if (opts == -1)
                break;
            switch (opts) {
                case 'd':
                    rc = mpz_set_str(n, optarg, INPUT_BASE);
                    if (rc)
                        fprintf(stderr, "mpz_set_str(%s) failed.\n", optarg);
                    else
                        got_n = 1;
                    break;
                case 'f':
                    if (!read_n_from_file(n, optarg))
                        fprintf(stderr, "Failed to read n from %s.\n", optarg);
                    else
                        got_n = 1;
                    break;
                default :
                case 'h':
                    print_syntax(1);
                    break;

            }
        }


        if (!got_n) {
            printf("There was an error reading a number to test.\n");
            mpz_clear(n);
            kill_children(nprocs);
            goto end;
        }

        if (mpz_tstbit( n, 0 ) != 1) {
            printf("Either n=2 (prime), or it's even (not prime).\n");
            mpz_clear(n);
            kill_children(nprocs);
            goto end;
        }


        gmp_printf("Input base 10: %Zd\n", n);

        if (mpz_perfect_power_p(n) != 0) {
            gmp_printf("%Zd is not prime\n", n );
            kill_children(nprocs);
            goto end;
        }

        log_b2_n = get_log_b2(n);
        printf("Log(n) = %lu\n", log_b2_n);

        /*
         * find and set r such that the order
         * of n (mod r) is greater then log(n)^2
         *
         * Also test a <=r to make sure that
         * 1 < (a,n) < n.
         *
         */
        mpz_init(gmp_r);
        if (!find_r(n, gmp_r, log_b2_n)) {
            mpz_clear(n);
            mpz_clear(gmp_r);
            kill_children(nprocs);
            goto end;
        }
        gmp_printf("The order of %Zd (modulo %Zd) is greater then (%lu)^2\n",
                n, gmp_r, log_b2_n);



        t2 = MPI_Wtime();
        if (!test_polys(nprocs, n, gmp_r, log_b2_n)) {
            mpz_clear(gmp_r);
            t3 = MPI_Wtime();
            gmp_printf("%Zd is not prime\n", n);
            mpz_clear(n);
            goto end;
        }
        t3 = MPI_Wtime();

        gmp_printf("%Zd is prime.\n", n);
        mpz_clear(n);
        mpz_clear(gmp_r);
    }

end:
    if (!id) {
        printf(" Sequential timing: %f seconds\n", t2 - t1);
        printf(" Parallel timing: %f seconds\n", t3 - t2);
    }
    MPI_Finalize();
    return 0;
}






















