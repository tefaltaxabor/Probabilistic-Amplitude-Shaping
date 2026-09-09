function ccdm = ccdm_init(pA, amps, n)
% CCDM_INIT  Rate metadata of a real, invertible constant-composition
%   distribution matcher of block length n targeting pA.
%
%   ccdm = pro.ccdm_init(pA, amps, n)
%   Use pro.ccdm_encode / pro.ccdm_decode for the actual bit<->amplitude maps.
%
%   Rate fields
%   -----------
%     Rccdm       k/n, the rate the matcher delivers        [bits/amplitude]
%     Rloss       H(A) - k/n, the loss against the TARGET entropy. This is the
%                 figure to quote; it is the one the write-up defines.
%     Rloss_int   H(p^(n)) - k/n, the part caused by k = floor(log2|T^n(P)|)
%                 alone. Non-negative by construction. (This is what this
%                 function used to return as `Rloss`.)
%     Rloss_quant H(A) - H(p^(n)), the part caused by quantising the target to
%                 an n-type. NOT sign-definite: the n-type may carry MORE
%                 entropy than pA, e.g. -1.43e-3 at n=540, nu=0.05.
%     Dtype       D(p^(n) || pA), single-letter divergence of the composition.
%                 It does NOT approximate Rloss -- it decays as O(n^-2) while
%                 Rloss ~ (K-1)*log2(n)/(2n), so the two diverge with n.
%     Dseq        (1/n) D(P_A || pA^{(x)n}), normalised divergence between the
%                 matcher's distribution over SEQUENCES and the i.i.d. target.
%                 This is the quantity that does converge to Rloss:
%                     Dseq = H(p^(n)) + Dtype - k/n  ->  Rloss.

    pA   = pA(:).';
    amps = amps(:).';
    M    = numel(amps);
    assert(numel(pA) == M, 'pA and amps must have equal length');

    % --- 1. Quantize pA to an n-type (VD-optimal, Bocherer Alg. 2.5.4) ---
    comp   = pro.quantize_composition(n, pA);
    pQuant = comp / n;

    % --- 2. |T^n(P)| = n!/(n_1!...n_M!), in log2 via gammaln to avoid overflow
    logNperm = (gammaln(n+1) - sum(gammaln(comp+1))) / log(2);

    % --- 3. CCDM rate: k = floor(log2 |T^n(P)|) input bits (2.58) ---
    k     = floor(logNperm);
    Rccdm = k / n;

    % --- 4. Entropies, rate-loss decomposition and divergences ---
    p     = pQuant(pQuant > 0);
    Hbar  = -sum(p .* log2(p));            % H(Abar), entropy of the n-type
    q     = pA(pA > 0);
    HA    = -sum(q .* log2(q));            % H(A),    entropy of the target

    sup   = pQuant > 0;
    Dtype = sum(pQuant(sup) .* log2(pQuant(sup) ./ pA(sup)));

    ccdm.n           = n;
    ccdm.M           = M;
    ccdm.amps        = amps;
    ccdm.comp        = comp;
    ccdm.pQuant      = pQuant;
    ccdm.k           = k;
    ccdm.HA          = HA;
    ccdm.Hbar        = Hbar;
    ccdm.Rccdm       = Rccdm;
    ccdm.Rloss       = HA   - Rccdm;       % write-up definition
    ccdm.Rloss_int   = Hbar - Rccdm;       % integrality part
    ccdm.Rloss_quant = HA   - Hbar;        % quantisation part (signed)
    ccdm.Dtype       = Dtype;
    ccdm.Dseq        = Hbar + Dtype - Rccdm;
    ccdm.logNperm    = logNperm;
end
