#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "aks.h"

struct polynomial *poly_alloc(unsigned long int degree) {
    unsigned long int cnt;
    struct polynomial *f;

    f = (struct polynomial *) malloc(sizeof(struct polynomial));
    if (!f) {
        fprintf(stderr, "poly_alloc:  %s\n", strerror(errno));
        return 0;
    }

    f->degree = degree;
    f->coeff = (mpz_t *) malloc((degree+1) * sizeof(mpz_t));
    if (!f->coeff) {
        fprintf(stderr, "poly_alloc:  %s\n", strerror(errno));
        return 0;
    }

    for (cnt = 0; cnt <= degree; cnt++)
        mpz_init(f->coeff[cnt]);

    return f;
}



void poly_free(struct polynomial *f) {
    unsigned long int cnt;

    for (cnt = 0; cnt <= f->degree; cnt++)
        mpz_clear(f->coeff[cnt]);

    free(f->coeff);
    free(f);
    return;
}



struct polynomial *poly_mult(
        struct polynomial *a,
        struct polynomial *b,
        short int freea,
        short int freeb,
        mpz_t n,
        unsigned long int r,
        unsigned long int log_b2_n ) {

    struct polynomial *res;
    unsigned long int cnt;
    mpz_t A, total, tmp, a_bar, b_bar, res_bar, tmp_pow;

    mpz_init(A);
    mpz_init(total);
    mpz_init(tmp);
    mpz_init(a_bar);
    mpz_init(b_bar);
    mpz_init(res_bar);
    mpz_init(tmp_pow);

    mpz_ui_pow_ui(A, 2, 2 * log_b2_n + 1);
    mpz_mul_ui(A, A, r);

    mpz_set(a_bar, a->coeff[0]);
    for (cnt = 1; cnt <= a->degree; cnt++) {
        mpz_pow_ui(tmp, A, cnt);
        mpz_mul(tmp, tmp, a->coeff[cnt]);
        mpz_add(a_bar, a_bar, tmp);
    }

    mpz_set(b_bar, b->coeff[0]);
    for (cnt = 1; cnt <= b->degree; cnt++) {
        mpz_pow_ui(tmp, A, cnt);
        mpz_mul(tmp, tmp, b->coeff[cnt]);
        mpz_add(b_bar, b_bar, tmp);
    }

    mpz_mul(res_bar, a_bar, b_bar);
    mpz_clear(a_bar);
    mpz_clear(b_bar);

    res = poly_alloc(a->degree + b->degree);

    mpz_set_ui(tmp_pow, 1);
    mpz_mod(res->coeff[0], res_bar, A);
    mpz_set(total, res->coeff[0]);

    for (cnt = 1; cnt <= res->degree; cnt++) {
        /*
         * tmp = res_bar - sum of lesser parts
         * of the polynomial res
         */
        mpz_sub(tmp, res_bar, total);
        mpz_mul(tmp_pow, tmp_pow, A);
        mpz_div(tmp, tmp, tmp_pow);
        mpz_mod(res->coeff[cnt], tmp, A);
        /* now tmp is c_i*(2A)^cnt */
        mpz_mul(tmp, res->coeff[cnt], tmp_pow);
        mpz_add(total, total, tmp);
    }

    if (freea)
        poly_free(a);
    if (freeb)
        poly_free(b);

    return res;
}



/*
 * Calculate f (mod x^r-1, n )
 *
 * f can have at most degree = 2r-2,
 * this makes sense for our AKS
 * implementation, but not in general.
 *
 */
struct polynomial *poly_mod(
        struct polynomial *f,
        short int free,
        mpz_t n,
        unsigned long int r) {

    unsigned long int cnt;
    struct polynomial * res;

    for (cnt = 0; cnt <= f->degree; cnt++)
        mpz_mod(f->coeff[cnt], f->coeff[cnt], n);

    res = poly_alloc(r - 1);
    /*
     * Copy lower degree terms
     */
    for (cnt = 0; cnt < r; cnt++)
        mpz_set(res->coeff[cnt], f->coeff[cnt]);

    for (cnt = r; cnt <= f->degree; cnt++) {
        if (mpz_sgn( f->coeff[cnt] ) != 0) {
            mpz_add(res->coeff[cnt-r], res->coeff[cnt-r], f->coeff[cnt]);
            mpz_mod(res->coeff[cnt-r], res->coeff[cnt-r], n);
        }
    }

    if (free)
        poly_free( f );

    return res;
}



void poly_print(struct polynomial *f) {
    unsigned long int cnt;

    gmp_printf("(%Zd) + ", f->coeff[0]);
    for (cnt = 1; cnt < f->degree; cnt++)
        gmp_printf("(%Zd)X[%lu] + ", f->coeff[cnt], cnt);

    gmp_printf("(%Zd)X[%lu]\n", f->coeff[f->degree], f->degree);
    return;
}



struct polynomial *poly_expmod(
        struct polynomial *f,
        short int free,
        mpz_t n,
        unsigned long int r,
        unsigned long int log_b2_n) {

    unsigned long int cnt;
    struct polynomial * res;
    struct polynomial * exp_tmp;
    struct polynomial * mult_tmp;
    int i;
    unsigned long int curr_exp = 0;

    /*
     * tmp will be at most the multiple of
     * two polynomials with degree < r.
     *
     * curr_exp represents the current
     * squared power of f that we're at.
     * ie exp_tmp = (f)^2^curr_exp
     *
     * exp_tmp starts at f
     *
     * As n is odd, we can assume that
     * we need f as part of our result
     */

    exp_tmp = poly_alloc(r - 1);
    for (cnt = 0; cnt <= f->degree; cnt++)
        mpz_set(exp_tmp->coeff[cnt], f->coeff[cnt]);

    mult_tmp = poly_alloc(r - 1);
    for (cnt = 0; cnt <= f->degree; cnt++)
        mpz_set(mult_tmp->coeff[cnt], f->coeff[cnt]);

    for (cnt = 1; cnt <= log_b2_n; cnt++){
        if (mpz_tstbit( n, cnt ) == 1){
            for (i = 0; i < (cnt - curr_exp); i++) {
                exp_tmp = poly_mult(exp_tmp, exp_tmp, 1, 0, n , r, log_b2_n);
                exp_tmp= poly_mod(exp_tmp, 1, n, r);
            }
            curr_exp = cnt;
            mult_tmp = poly_mult(mult_tmp, exp_tmp, 1, 0, n, r, log_b2_n);
            mult_tmp = poly_mod(mult_tmp, 1, n, r);
        }
    }
    res = poly_alloc(r - 1);

    for (cnt = 0; cnt < r; cnt++)
        mpz_set(res->coeff[cnt], mult_tmp->coeff[cnt]);

    /*
     * XXX: sanity check
     */

    for (cnt = r; cnt <= mult_tmp->degree; cnt++)
        if (mpz_sgn( mult_tmp->coeff[cnt] ) != 0)
            gmp_printf("Error at mult_tmp[%lu] = %Zd\n", cnt, mult_tmp->coeff[cnt]);

    poly_free(mult_tmp);
    poly_free(exp_tmp);
    if (free)
        poly_free(f);

    return res;
}



int poly_cmp(struct polynomial *f, struct polynomial *g) {
    unsigned long int cnt;

    if (f->degree != g->degree)
        return 0;

    for (cnt = 0; cnt <= f->degree; cnt++)
        if (mpz_cmp(f->coeff[cnt], g->coeff[cnt]) != 0)
            return 0;

    return 1;
}

