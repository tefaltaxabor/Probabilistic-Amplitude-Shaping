function ccdm = ccdm_init(pA, amps, n)
% CCDM_INIT  Build a real, invertible Constant-Composition Distribution Matcher.
%
%   ccdm = pro.ccdm_init(pA, amps, n)
%
%   Implements the CCDM of Bocherer, "Probabilistic Amplitude Shaping" (2023),
%   Section 2.5: a one-to-one map from k input bits to length-n amplitude
%   sequences of FIXED composition (type) n_1,...,n_M. All output sequences
%   have exactly the same empirical distribution, hence "constant composition".
%
%   The bit<->sequence map is realized by arithmetic coding over the multiset
%   permutations of the fixed composition (ranking/unranking), which is exact,
%   invertible, and needs no lookup table.
%
%   Inputs
%   ------
%   pA   : (1 x M) target amplitude probabilities (Maxwell-Boltzmann), sum = 1
%   amps : (1 x M) amplitude alphabet values (e.g. [1 3 5 7] for 8-ASK)
%   n    : DM output length (amplitudes per DM call). Choose n | (nfec/m) so a
%          DM sequence is never shared across codewords (Bocherer Ex. 2.1).
%
%   Output
%   ------
%   ccdm : struct with fields
%          n        : output length
%          M        : alphabet size
%          amps     : alphabet
%          comp     : (1 x M) integer composition [n_1,...,n_M], sum = n
%          pQuant   : (1 x M) realized n-type PA' = comp/n
%          k        : number of input bits, k = floor(log2 |T^n(P)|)
%          Hbar     : entropy H(Abar) of the quantized type [bits/amp]
%          Rccdm    : CCDM rate k/n [bits/amp]
%          Rloss    : rate loss Hbar - Rccdm [bits/amp]
%          logNperm : log2 of the multinomial |T^n(P)| (number of sequences)
%
%   Use pro.ccdm_encode / pro.ccdm_decode for the actual bit<->amplitude maps.

    pA   = pA(:).';
    amps = amps(:).';
    M    = numel(amps);
    assert(numel(pA) == M, 'pA and amps must have equal length');

    % --- 1. Quantize pA to an n-type PA' (Bocherer Alg. 2.5.4, VD-optimal) ---
    %   Work in INTEGER counts to avoid floating-point round-off in the type.
    base = floor(n * pA);                  % (2.50) integer base counts
    L    = round(n - sum(base));           % (2.51) integer by construction
    err  = pA - base / n;                  % approximation error per symbol
    [~, order] = sort(err, 'descend');     % largest error first
    comp = base;                           % base counts (integers)
    comp(order(1:L)) = comp(order(1:L)) + 1;   % add 1 to the L largest-error
    assert(sum(comp) == n, 'composition must sum to n');
    pQuant = comp / n;

    % --- 2. Number of constant-composition sequences: multinomial coefficient ---
    % |T^n(P)| = n! / (n_1! ... n_M!). Work in log2 via gammaln to avoid overflow.
    logNperm = (gammaln(n+1) - sum(gammaln(comp+1))) / log(2);   % log2 of count

    % --- 3. CCDM rate: k = floor(log2 |T^n(P)|) input bits (2.58) ---
    k = floor(logNperm);

    % --- 4. Entropy of the quantized type and rate loss ---
    p = pQuant(pQuant > 0); %support vector
    Hbar  = -sum(p .* log2(p)); % H(Abar) [bits/amp]
    Rccdm = k / n; % (2.58)
    Rloss = Hbar - Rccdm; % rate loss (Remark 2.4)

    ccdm.n        = n;
    ccdm.M        = M;
    ccdm.amps     = amps;
    ccdm.comp     = comp;
    ccdm.pQuant   = pQuant;
    ccdm.k        = k;
    ccdm.Hbar     = Hbar;
    ccdm.Rccdm    = Rccdm;
    ccdm.Rloss    = Rloss;
    ccdm.logNperm = logNperm;
end
