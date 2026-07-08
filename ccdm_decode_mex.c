/* ccdm_decode_mex.c  (OPTIMIZED -- incremental multinomial caching)
 *
 * CCDM decoder (ranking), exact inverse of ccdm_encode_mex. Uses the same
 * incremental multinomial identity so the per-position work is O(M) big-int
 * mul/div instead of O(M^2) binomials.
 *
 *   At each position with remaining composition rem (total = sum rem) and
 *   cached N = N(rem): the rank contribution of symbol s is sum over t<s of
 *   N_t = N * rem[t] / total. After consuming the actual symbol s, update
 *   N <- N * rem[s] / total, total -= 1, rem[s] -= 1.
 *
 * MATLAB call:  bits = ccdm_decode_mex(a, comp, amps, k)
 * Build:        mex -lgmp ccdm_decode_mex.c
 * Gabriel Cabrera -- PAS + HARQ thesis.
 */

#include "mex.h"
#include <gmp.h>

static void multinomial_full(mpz_t res, const long *counts, int M)
{
    long total = 0;
    for (int i = 0; i < M; ++i) total += counts[i];
    mpz_t binom; mpz_init(binom);
    mpz_set_ui(res, 1);
    long running = total;
    for (int i = 0; i < M; ++i) {
        mpz_bin_uiui(binom, (unsigned long)running, (unsigned long)counts[i]);
        mpz_mul(res, res, binom);
        running -= counts[i];
    }
    mpz_clear(binom);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4)
        mexErrMsgIdAndTxt("ccdm:decode:nargin",
                          "Usage: bits = ccdm_decode_mex(a, comp, amps, k)");

    const double *a     = mxGetPr(prhs[0]);
    const double *compd = mxGetPr(prhs[1]);
    const double *amps  = mxGetPr(prhs[2]);
    long          k     = (long)mxGetScalar(prhs[3]);

    long n = (long)mxGetNumberOfElements(prhs[0]);
    int  M = (int)mxGetNumberOfElements(prhs[1]);

    long *rem = (long*)mxMalloc(M * sizeof(long));
    long total = 0;
    for (int i = 0; i < M; ++i) { rem[i] = (long)compd[i]; total += rem[i]; }

    /* map amplitude values -> symbol indices */
    int *symIdx = (int*)mxMalloc(n * sizeof(int));
    for (long pos = 0; pos < n; ++pos) {
        int found = -1;
        for (int s = 0; s < M; ++s) {
            if (a[pos] == amps[s]) { found = s; break; }
        }
        if (found < 0)
            mexErrMsgIdAndTxt("ccdm:decode:alpha",
                              "sequence value not in amplitude alphabet");
        symIdx[pos] = found;
    }

    mpz_t idx, N, Nt, tmp;
    mpz_init(idx); mpz_init(N); mpz_init(Nt); mpz_init(tmp);
    mpz_set_ui(idx, 0);

    multinomial_full(N, rem, M);           /* initial N(rem), ONCE */

    for (long pos = 0; pos < n; ++pos) {
        int s = symIdx[pos];
        /* add N_t = N*rem[t]/total for all t < s (still available) */
        for (int t = 0; t < s; ++t) {
            if (rem[t] == 0) continue;
            mpz_mul_ui(tmp, N, (unsigned long)rem[t]);
            mpz_divexact_ui(Nt, tmp, (unsigned long)total);
            mpz_add(idx, idx, Nt);
        }
        /* update N for next position: N <- N*rem[s]/total */
        mpz_mul_ui(tmp, N, (unsigned long)rem[s]);
        mpz_divexact_ui(N, tmp, (unsigned long)total);
        rem[s] -= 1;
        total  -= 1;
    }

    /* expand idx to k bits, MSB first: bits[b] = bit (k-1-b) of idx */
    plhs[0] = mxCreateDoubleMatrix(1, (mwSize)k, mxREAL);
    double *bits = mxGetPr(plhs[0]);
    for (long b = 0; b < k; ++b) {
        bits[b] = (double)mpz_tstbit(idx, (mp_bitcnt_t)(k - 1 - b));
    }

    mpz_clear(idx); mpz_clear(N); mpz_clear(Nt); mpz_clear(tmp);
    mxFree(rem); mxFree(symIdx);
}