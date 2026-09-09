function ccdm = ccdm_init(pA, amps, n)
%   invertible Constant-Composition Distribution Matcher.
%
%   ccdm = pro.ccdm_init(pA, amps, n)
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
