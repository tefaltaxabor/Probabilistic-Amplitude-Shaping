function [R, HX, HbGivenY] = bmd_rate(cstll, snr_dB, nSym, seed)
% BMD_RATE  Bit-metric decoding rate of a real ASK constellation over AWGN.
%
%   [R, HX, HbGivenY] = src.bmd_rate(cstll, snr_dB, nSym, seed)
%
%   Monte-Carlo estimate of the rate achievable by BIT-METRIC DECODING, i.e.
%   by the receiver this chain actually uses (soft demapper producing per-bit
%   LLRs, followed by a binary LDPC decoder):
%
%       R_BMD = [ H(X) - sum_j H(B_j|Y) ]^+        [bits / real dimension]
%
%   over the real channel Y = X + N, X drawn from cstll.alphabet with priors
%   cstll.px, N ~ N(0, sigma^2). The alphabet is rescaled internally so that
%   E[X^2] = 1, hence sigma^2 = 10^(-snr_dB/10) and snr_dB is the same SNR
%   used by channel.real_channel / channel.complex_channel (per real
%   dimension the complex chain has signal power Es/2 and noise sigma^2/2, so
%   the per-dimension SNR equals the 2D SNR).
%
%   R_BMD is the relevant upper bound for PAS: it sits below the constellation
%   -constrained MI (the gap is the cost of bit-metric instead of symbol-metric
%   decoding) and above what any binary code on this mapping can deliver.
%   Comparing it against the coded operating points splits the total gap to
%   capacity into a MAPPING penalty (C - R_BMD) and a CODING penalty
%   (R_BMD - achieved rate).
%
%   Inputs
%   ------
%   cstll  : ASK struct from pro.dig_mod_ASK, with .px set to the desired
%            priors (uniform, or the shaped px from pro.build_shaping)
%   snr_dB : (1,S) SNR grid [dB]
%   nSym   : symbols per SNR point (default 2e5)
%   seed   : optional RNG seed for reproducibility
%
%   Outputs
%   -------
%   R        : (1,S) BMD rate [bits/real dimension]
%   HX       : scalar, input entropy H(X) [bits]
%   HbGivenY : (S,m) per-bit-level conditional entropies H(B_j|Y)

    if nargin < 3 || isempty(nSym), nSym = 2e5; end
    if nargin >= 4 && ~isempty(seed), rng(seed); end

    px = cstll.px(:).';                     % (1,M) priors, row
    A  = cstll.alphabet(:).';               % (1,M) alphabet,  row
    m  = size(cstll.label, 2);

    assert(numel(px) == numel(A), 'px and alphabet must have equal length');

    % --- normalise to E[X^2] = 1 so that sigma^2 = 1/SNR ---
    A = A / sqrt(sum(px .* A.^2));

    % --- input entropy ---
    pp = px(px > 0);
    HX = -sum(pp .* log2(pp));

    % constellation copy used by the demapper (normalised alphabet)
    cc = cstll;  cc.alphabet = A;  cc.px = px;

    nS       = numel(snr_dB);
    R        = nan(1, nS);
    HbGivenY = nan(nS, m);

    edges = [0, cumsum(px)];  edges(end) = 1;      % for inverse-CDF sampling

    for i = 1:nS
        sigma2 = 10^(-snr_dB(i)/10);               % E[X^2] = 1

        % --- draw symbols from px ---
        [~, ~, k] = histcounts(rand(nSym,1), edges);
        k = min(max(k, 1), numel(A));              % guard against edge cases
        bits = double(cstll.label(k, :));          % (nSym, m), row k <-> alphabet(k)
        x    = A(k).';                             % (nSym, 1)

        % --- real AWGN + soft demapper with the SAME priors ---
        y = x + sqrt(sigma2) * randn(nSym, 1);
        L = pro.demap(y, cc, sigma2, 'SD');        % (nSym, m), LLR = log(P0/P1)

        % --- H(B_j|Y) = E[ log2(1 + exp(-(1-2b) L)) ] ---
        s  = 1 - 2*bits;                           % +1 if bit=0, -1 if bit=1
        Hj = mean(log2_1pexp(-s .* L), 1);         % (1,m)

        HbGivenY(i,:) = Hj;
        R(i)          = max(HX - sum(Hj), 0);
    end
end

function v = log2_1pexp(z)
    % log2(1 + exp(z)), overflow-safe for large |z|.
    v        = zeros(size(z));
    pos      = z > 0;
    v(pos)   = (z(pos) + log1p(exp(-z(pos)))) / log(2);
    v(~pos)  = log1p(exp(z(~pos))) / log(2);
end
