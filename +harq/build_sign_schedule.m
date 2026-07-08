function sch = build_sign_schedule(n, maxTx, frac1, order, amps_or_seed)
% BUILD_SIGN_SCHEDULE  Parity(sign)-release schedule for Approach B (sign-level IR).
%
%   sch = harq.build_sign_schedule(n, maxTx)
%   sch = harq.build_sign_schedule(n, maxTx, frac1, order, amps_or_seed)
%
%   Partitions the n parity (sign) indices {1..n} into maxTx disjoint subsets
%   signSets{1..maxTx}. Round 1 releases a fraction `frac1` of the signs; the
%   remaining signs are split evenly across rounds 2..maxTx. All amplitude bits
%   are always transmitted (handled in run_point_signIR), so only the sign
%   partition is defined here.
%
%   Inputs
%   ------
%   n      : symbols per dimension (= number of parity/sign bits)
%   maxTx  : number of HARQ rounds (= number of sign subsets)
%   frac1  : fraction of signs released in round 1 (default 0.5). Higher -> more
%            parity up front -> stronger round 1 but less incremental gain.
%   order  : 'random' (default) baseline
%            'weakfirst'  (H1) release signs on SMALLEST amplitudes first
%            'strongfirst'       release signs on LARGEST amplitudes first
%   amps_or_seed :
%            - for 'random'    : integer seed (default 12345)
%            - for 'weak/strong': (n x 1) amplitude magnitudes of THIS codeword
%              (realization-dependent ordering; must be passed per codeword)
%
%   Output
%   ------
%   sch.mode      = "SIGN-IR"
%   sch.maxTx     = maxTx
%   sch.order     = order
%   sch.frac1     = frac1
%   sch.n         = n
%   sch.signSets  = {1 x maxTx} cell of released sign indices per round
%   sch.cumSigns  = (1 x maxTx) cumulative released signs after round r
%   sch.roundSigns= (1 x maxTx) signs released in round r

    if nargin < 3 || isempty(frac1), frac1 = 0.5;     end
    if nargin < 4 || isempty(order), order = 'random'; end
    assert(maxTx >= 1, 'maxTx must be >= 1');
    assert(frac1 > 0 && frac1 < 1, 'frac1 must be in (0,1)');

    % --- Determine the ORDER in which signs are released (index list) ---
    switch lower(order)
        case 'random'
            if nargin < 5 || isempty(amps_or_seed), seed = 12345;
            else, seed = amps_or_seed; end
            rs   = RandStream('twister', 'Seed', seed);
            perm = randperm(rs, n);                 % release order

        case 'weakfirst'   % H1: smallest amplitude signs released first
            assert(nargin >= 5 && numel(amps_or_seed) == n, ...
                'weakfirst needs per-codeword amplitudes (n x 1)');
            [~, perm] = sort(amps_or_seed(:), 'ascend');

        case 'strongfirst' % largest amplitude signs released first
            assert(nargin >= 5 && numel(amps_or_seed) == n, ...
                'strongfirst needs per-codeword amplitudes (n x 1)');
            [~, perm] = sort(amps_or_seed(:), 'descend');

        otherwise
            error('unknown order "%s" (random|weakfirst|strongfirst)', order);
    end

    % --- Split the ordered indices into round-1 chunk + even rounds 2..maxTx ---
    nR1 = round(frac1 * n);
    nR1 = min(max(nR1, 1), n - (maxTx - 1));   % leave >=1 sign for each later round
    signSets = cell(1, maxTx);
    signSets{1} = sort(perm(1:nR1));

    rest = perm(nR1+1:end);
    nRest = numel(rest);
    if maxTx > 1
        edges = round(linspace(0, nRest, maxTx));   % maxTx-1 later rounds
        for k = 1:(maxTx-1)
            seg = rest(edges(k)+1 : edges(k+1));
            signSets{k+1} = sort(seg);
        end
    end

    roundSigns = cellfun(@numel, signSets);
    cumSigns   = cumsum(roundSigns);

    sch.mode       = "SIGN-IR";
    sch.maxTx      = maxTx;
    sch.order      = string(order);
    sch.frac1      = frac1;
    sch.n          = n;
    sch.signSets   = signSets;
    sch.roundSigns = roundSigns;
    sch.cumSigns   = cumSigns;
end