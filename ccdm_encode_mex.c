/* ccdm_encode_mex.c
 *
 * CCDM encoder (unranking) with GMP big integers.  Maps k input bits to a
 * constant-composition amplitude sequence of length n, following Bocherer,
 * "Probabilistic Amplitude Shaping" (2023), Section 2.5.
 *
 * The index (k bits) is a big integer; multinomial coefficients are exact
 * big integers.  Using GMP removes the 2^53 double-precision limit of a
 * pure-MATLAB implementation, so n = 200 (or larger) is exact.
 *
 * MATLAB call:
 *     a = ccdm_encode_mex(bits, comp, amps)
 *
 *   bits : (1 x k) double/logical row, values in {0,1}, MSB first
 *   comp : (1 x M) double row, integer composition [n_1..n_M], sum(comp)=n
 *   amps : (1 x M) double row, amplitude alphabet values
 *   a    : (1 x n) double row, amplitude sequence with exact composition comp
 *
 * Build (on your machine, with libgmp-dev installed):
 *     mex -lgmp ccdm_encode_mex.c
 *
 * Gabriel Cabrera -- PAS + HARQ thesis.
 */

#include "mex.h"
#include <gmp.h>
#include <string.h>

/* multinomial( counts, M ) = (sum counts)! / prod(counts_i!), exact, into res */
static void multinomial(mpz_t res, const long *counts, int M)
{
    long total = 0;
    for (int i = 0; i < M; ++i) total += counts[i];

    /* res = total! / (c_0! c_1! ... c_{M-1}!) built as a product of binomials:
     *   multinom = C(total, c_0) * C(total-c_0, c_1) * ...
     * which keeps intermediate values smaller than computing total! directly. */
    mpz_t binom;
    mpz_init(binom);
    mpz_set_ui(res, 1);

    long running = total;
    for (int i = 0; i < M; ++i) {
        /* binom = C(running, counts[i]) */
        mpz_bin_uiui(binom, (unsigned long)running, (unsigned long)counts[i]);
        mpz_mul(res, res, binom);
        running -= counts[i];
    }
    mpz_clear(binom);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3)
        mexErrMsgIdAndTxt("ccdm:encode:nargin",
                          "Usage: a = ccdm_encode_mex(bits, comp, amps)");

    const double *bits = mxGetPr(prhs[0]);
    const double *compd = mxGetPr(prhs[1]);
    const double *amps = mxGetPr(prhs[2]);

    mwSize k = mxGetNumberOfElements(prhs[0]);
    int    M = (int)mxGetNumberOfElements(prhs[1]);
    if ((int)mxGetNumberOfElements(prhs[2]) != M)
        mexErrMsgIdAndTxt("ccdm:encode:dim", "comp and amps must have equal length");

    /* remaining composition as long[] */
    long *rem = (long*)mxMalloc(M * sizeof(long));
    long n = 0;
    for (int i = 0; i < M; ++i) { rem[i] = (long)compd[i]; n += rem[i]; }

    /* idx = integer value of the k bits (MSB first) */
    mpz_t idx, acc, Ns;
    mpz_init(idx); mpz_init(acc); mpz_init(Ns);
    mpz_set_ui(idx, 0);
    for (mwSize b = 0; b < k; ++b) {
        mpz_mul_2exp(idx, idx, 1);                 /* idx <<= 1 */
        if (bits[b] != 0.0) mpz_add_ui(idx, idx, 1);
    }

    /* output */
    plhs[0] = mxCreateDoubleMatrix(1, (mwSize)n, mxREAL);
    double *a = mxGetPr(plhs[0]);

    long *remS = (long*)mxMalloc(M * sizeof(long));

    for (long pos = 0; pos < n; ++pos) {
        mpz_set_ui(acc, 0);
        int chosen = -1;
        for (int s = 0; s < M; ++s) {
            if (rem[s] == 0) continue;
            memcpy(remS, rem, M * sizeof(long));
            remS[s] -= 1;
            multinomial(Ns, remS, M);              /* # sequences starting with s */
            /* if idx < acc + Ns  -> choose s */
            mpz_add(acc, acc, Ns);                 /* acc now = previous acc + Ns */
            if (mpz_cmp(idx, acc) < 0) {
                /* idx -= (acc - Ns) = previous acc */
                mpz_sub(acc, acc, Ns);             /* restore acc to block start */
                mpz_sub(idx, idx, acc);
                chosen = s;
                break;
            }
        }
        if (chosen < 0)
            mexErrMsgIdAndTxt("ccdm:encode:range",
                              "index out of range: bits exceed |T^n(P)|");
        a[pos] = amps[chosen];
        rem[chosen] -= 1;
    }

    mpz_clear(idx); mpz_clear(acc); mpz_clear(Ns);
    mxFree(rem); mxFree(remS);
}
