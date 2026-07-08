/* ccdm_encode_mex.c  (OPTIMIZED -- incremental multinomial caching)
 *
 * CCDM encoder (unranking) with GMP big integers. Same exact, invertible map
 * as the original (Bocherer 2023, Sec. 2.5), but the per-position multinomial
 * is UPDATED in O(M) big-int mul/div instead of recomputed from scratch in
 * O(M^2) binomials. This removes the ~262 s that ccdm_encode_mex cost in the
 * 64-QAM sweep profiler (4.9M calls).
 *
 * Exact identity (integer arithmetic):
 *   N(rem) = (sum rem)! / prod(rem_i!).  Sequences starting with s:
 *   N_s = N(rem) * rem[s] / total,  total = sum rem.  After choosing s:
 *   N(rem with s decremented) = N_s, total -= 1.  Division is exact
 *   (mpz_divexact_ui) -> bit-identical to the scratch computation.
 *
 * MATLAB call:  a = ccdm_encode_mex(bits, comp, amps)
 * Build:        mex -lgmp ccdm_encode_mex.c
 * Gabriel Cabrera -- PAS + HARQ thesis.
 */

#include "mex.h"
#include <gmp.h>

/* full multinomial -- used ONCE for the initial N(rem) */
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
    if (nrhs != 3)
        mexErrMsgIdAndTxt("ccdm:encode:nargin",
                          "Usage: a = ccdm_encode_mex(bits, comp, amps)");

    const double *bits  = mxGetPr(prhs[0]);
    const double *compd = mxGetPr(prhs[1]);
    const double *amps  = mxGetPr(prhs[2]);

    mwSize k = mxGetNumberOfElements(prhs[0]);
    int    M = (int)mxGetNumberOfElements(prhs[1]);
    if ((int)mxGetNumberOfElements(prhs[2]) != M)
        mexErrMsgIdAndTxt("ccdm:encode:dim", "comp and amps must have equal length");

    long *rem = (long*)mxMalloc(M * sizeof(long));
    long total = 0;
    for (int i = 0; i < M; ++i) { rem[i] = (long)compd[i]; total += rem[i]; }
    long n = total;

    mpz_t idx, acc, Ns, N, tmp;
    mpz_init(idx); mpz_init(acc); mpz_init(Ns); mpz_init(N); mpz_init(tmp);
    mpz_set_ui(idx, 0);
    for (mwSize b = 0; b < k; ++b) {
        mpz_mul_2exp(idx, idx, 1);
        if (bits[b] != 0.0) mpz_add_ui(idx, idx, 1);
    }

    plhs[0] = mxCreateDoubleMatrix(1, (mwSize)n, mxREAL);
    double *a = mxGetPr(plhs[0]);

    multinomial_full(N, rem, M);           /* initial N(rem), ONCE */

    for (long pos = 0; pos < n; ++pos) {
        mpz_set_ui(acc, 0);
        int chosen = -1;
        for (int s = 0; s < M; ++s) {
            if (rem[s] == 0) continue;
            /* Ns = N * rem[s] / total   (exact) */
            mpz_mul_ui(tmp, N, (unsigned long)rem[s]);
            mpz_divexact_ui(Ns, tmp, (unsigned long)total);
            mpz_add(acc, acc, Ns);
            if (mpz_cmp(idx, acc) < 0) {
                mpz_sub(acc, acc, Ns);     /* acc = block start */
                mpz_sub(idx, idx, acc);
                mpz_set(N, Ns);            /* N for next position */
                rem[s] -= 1;
                total  -= 1;
                chosen = s;
                break;
            }
        }
        if (chosen < 0)
            mexErrMsgIdAndTxt("ccdm:encode:range",
                              "index out of range: bits exceed |T^n(P)|");
        a[pos] = amps[chosen];
    }

    mpz_clear(idx); mpz_clear(acc); mpz_clear(Ns);
    mpz_clear(N);   mpz_clear(tmp);
    mxFree(rem);
}