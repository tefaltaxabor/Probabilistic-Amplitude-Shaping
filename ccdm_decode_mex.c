/* ccdm_decode_mex.c
 *
 * CCDM decoder (ranking) with GMP big integers.  Exact inverse of
 * ccdm_encode_mex: maps a constant-composition amplitude sequence of length n
 * back to its k input bits.  Bocherer, "Probabilistic Amplitude Shaping"
 * (2023), Section 2.5.
 *
 * MATLAB call:
 *     bits = ccdm_decode_mex(a, comp, amps, k)
 *
 *   a    : (1 x n) double row, amplitude sequence from ccdm_encode_mex
 *   comp : (1 x M) double row, integer composition [n_1..n_M]
 *   amps : (1 x M) double row, amplitude alphabet
 *   k    : scalar, number of output bits (= floor(log2 |T^n(P)|))
 *   bits : (1 x k) double row, recovered input bits, MSB first
 *
 * Build:
 *     mex -lgmp ccdm_decode_mex.c
 *
 * Gabriel Cabrera -- PAS + HARQ thesis.
 */

#include "mex.h"
#include <gmp.h>
#include <string.h>

static void multinomial(mpz_t res, const long *counts, int M)
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
    for (int i = 0; i < M; ++i) rem[i] = (long)compd[i];

    /* map each amplitude value to its symbol index 0..M-1 */
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

    mpz_t idx, term;
    mpz_init(idx); mpz_init(term);
    mpz_set_ui(idx, 0);

    long *remT = (long*)mxMalloc(M * sizeof(long));

    for (long pos = 0; pos < n; ++pos) {
        int s = symIdx[pos];
        /* add blocks of all symbols t < s still available */
        for (int t = 0; t < s; ++t) {
            if (rem[t] == 0) continue;
            memcpy(remT, rem, M * sizeof(long));
            remT[t] -= 1;
            multinomial(term, remT, M);
            mpz_add(idx, idx, term);
        }
        rem[s] -= 1;
    }

    /* expand idx to k bits, MSB first */
    plhs[0] = mxCreateDoubleMatrix(1, (mwSize)k, mxREAL);
    double *bits = mxGetPr(plhs[0]);
    for (long b = k - 1; b >= 0; --b) {
        bits[b] = (double)mpz_tstbit(idx, (mp_bitcnt_t)(k - 1 - b));
    }
    /* note: tstbit(idx, j) gives bit j (LSB=0). We want MSB-first output:
     * bits[0] is the MSB = bit (k-1). The loop above writes bits[b] using
     * tstbit(k-1-b); for b=0 -> tstbit(k-1) = MSB. Correct. */

    mpz_clear(idx); mpz_clear(term);
    mxFree(rem); mxFree(symIdx); mxFree(remT);
}
