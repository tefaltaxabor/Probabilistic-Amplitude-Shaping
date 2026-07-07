function sch = build_schedule(n, maxTx, mode, punctureFrac, seed)
% BUILD_SCHEDULE  Per-round transmit schedule for HARQ over a PAS codeword.
%
%   sch = harq.build_schedule(n, maxTx, mode)
%   sch = harq.build_schedule(n, maxTx, mode, punctureFrac, seed)
%
%   Both HARQ types share ONE fixed-rate PAS mother code (n coded symbols per
%   real dimension, parity -> sign). The two strategies differ ONLY in which
%   symbols each (re)transmission carries -- the receive/decode loop is
%   identical (see harq.run_point). This function encodes that difference.
%
%   Round 1 is IDENTICAL for both modes: it transmits a subset S1 of the n
%   symbols (the rest are punctured), so the initial code rate is higher than
%   the mother rate. On retransmissions:
%     - 'IR' (incremental redundancy): each round reveals NEW, previously
%            punctured symbols. The effective code rate DROPS every round.
%     - 'CC' (Chase combining): each round repeats S1. The rate stays fixed;
%            the receiver combines repeated copies for an SNR gain.
%
%   Inputs
%   ------
%   n            : coded symbols per dimension (= cfg.n of the PAS FEC code)
%   maxTx        : maximum number of (re)transmissions (HARQ rounds)
%   mode         : 'IR' or 'CC'
%   punctureFrac : fraction of the n symbols HIDDEN at round 1 (default 0.25).
%                  Round-1 code rate Rc1 = Rc_mother * n / round((1-f)*n).
%                  Set 0 for classic full-codeword Chase combining.
%   seed         : seed for the (reproducible) puncturing permutation.
%                  Uses a LOCAL RandStream, so it never disturbs the global rng.
%
%   Output
%   ------
%   sch : struct with fields
%         mode, maxTx, punctureFrac, n, nTx1 (symbols in round 1),
%         txSets      : {1 x maxTx} cell, symbol indices transmitted each round
%         cumSyms     : (1 x maxTx) cumulative DISTINCT symbols seen after round r
%                       (for IR this grows; for CC it stays nTx1)
%         roundSyms   : (1 x maxTx) symbols transmitted in round r (channel uses)

    if nargin < 4 || isempty(punctureFrac), punctureFrac = 0.25; end
    if nargin < 5 || isempty(seed),          seed = 12345;       end

    mode = upper(string(mode));
    assert(maxTx >= 1, 'maxTx must be >= 1');
    assert(punctureFrac >= 0 && punctureFrac < 1, 'punctureFrac must be in [0,1)');

    % Reproducible puncturing pattern, isolated from the global rng stream so
    % that the schedule is identical on every parfor worker and every run.
    rs      = RandStream('twister', 'Seed', seed);
    perm    = randperm(rs, n);
    nHidden = round(punctureFrac * n);
    nHidden = min(max(nHidden, 0), n - 1);      % keep at least 1 symbol in S1
    nTx1    = n - nHidden;

    S1     = sort(perm(1:nTx1));                % round-1 symbols (both modes)
    hidden = perm(nTx1+1:end);                  % symbols kept back for IR

    txSets = cell(1, maxTx);

    switch mode
        case "CC"
            % Chase combining: every round repeats the round-1 packet.
            for r = 1:maxTx
                txSets{r} = S1;
            end

        case "IR"
            % Incremental redundancy: reveal new symbols over rounds 2..maxTx.
            txSets{1} = S1;
            nR = maxTx - 1;                      % number of reveal rounds
            if nR >= 1 && ~isempty(hidden)
                if numel(hidden) < nR
                    warning(['harq:build_schedule: only %d hidden symbols for ' ...
                             '%d reveal rounds; some IR rounds will be empty.'], ...
                             numel(hidden), nR);
                end
                edges = round(linspace(0, numel(hidden), nR + 1));
                for k = 1:nR
                    seg = hidden(edges(k)+1 : edges(k+1));
                    txSets{k+1} = sort(seg);    % may be empty if hidden is small
                end
            else
                for r = 2:maxTx, txSets{r} = []; end
                if punctureFrac == 0
                    warning(['harq:build_schedule: IR with punctureFrac=0 has no ' ...
                             'incremental symbols; behaves as a single-shot code.']);
                end
            end

        otherwise
            error('harq:build_schedule: unknown mode "%s" (use IR or CC)', mode);
    end

    % Bookkeeping: distinct symbols seen and channel uses per round.
    roundSyms = cellfun(@numel, txSets);
    cumSyms   = zeros(1, maxTx);
    seen      = false(1, n);
    for r = 1:maxTx
        seen(txSets{r}) = true;
        cumSyms(r) = nnz(seen);
    end

    sch.mode         = mode;
    sch.maxTx        = maxTx;
    sch.punctureFrac = punctureFrac;
    sch.n            = n;
    sch.nTx1         = nTx1;
    sch.txSets       = txSets;
    sch.roundSyms    = roundSyms;
    sch.cumSyms      = cumSyms;
    sch.seed         = seed;
end
