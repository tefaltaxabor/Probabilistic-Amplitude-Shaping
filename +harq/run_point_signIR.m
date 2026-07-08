function out = run_point_signIR(snr, cfg, cstll, pA, amp_label, sch, maxFrames, targetCwErr, maxIter)
% RUN_POINT_SIGNIR  Monte-Carlo of one SNR point for Approach B (sign-level IR).
%
%   out = harq.run_point_signIR(snr, cfg, cstll, pA, amp_label, sch, ...
%                               maxFrames, targetCwErr, maxIter)
%
%   APPROACH B (proposed): incremental redundancy realized by releasing PARITY
%   (sign) bits progressively, while ALL amplitude (information) bits are
%   delivered in every round. Contrast with Approach A (run_point.m), which
%   punctures whole symbols.
%
%   PHYSICAL MODEL (Model 2, "genie sign reveal"):
%   -----------------------------------------------
%   Every round transmits the full set of n symbols x = s.*a (true sign, true
%   amplitude). The receiver ALWAYS observes all n symbols, but the sign LLR of
%   a symbol is ADMITTED into the decoder buffer only once its index has been
%   "released" by the schedule. Sign LLRs of not-yet-released indices are forced
%   to 0 (erasure). Amplitude LLRs are always admitted. Thus:
%     - amplitudes: refined every round (Chase-combining style),
%     - signs (parity): revealed incrementally (IR style).
%
%   This isolates the question "how much parity must be revealed to decode?"
%   without a separate amplitude-retransmission confound. Channel uses per round
%   are constant = n (all symbols sent every round).
%
%   SCHEDULE (sch) fields required (from a sign-partition builder):
%     sch.maxTx       : number of rounds (= number of parity subsets)
%     sch.signSets    : {1 x maxTx} cell, PARITY (sign) indices released in round r
%                       Their union must be 1:n (all signs eventually revealed).
%
%   Column convention (from fec.encode): column 1 = sign (parity),
%   columns 2:m = amplitude (information).
%
%   Output struct `out`: same fields as run_point.m (bler, thr, succRate,
%   avgTx, nCw, nFail). Here thr counts channel uses = n per round.

    if nargin < 9, maxIter = 20; end
    maxTx    = sch.maxTx;
    signSets = sch.signSets;          % parity indices released per round
    n        = cfg.n;
    K        = cfg.K;

    % Cumulative released-sign mask per round: relMask{r} = signs admitted by round r
    relMask = cell(1, maxTx);
    acc = false(n, 1);
    for r = 1:maxTx
        acc(signSets{r}) = true;
        relMask{r} = acc;             % logical (n x 1)
    end
    % sanity: after last round all signs must be released
    assert(all(relMask{maxTx}), ...
        'run_point_signIR: signSets do not cover all n signs by the last round');

    succCount    = zeros(1, maxTx);
    nFail        = 0;
    nCw          = 0;
    sumDelivered = 0;
    sumUsed      = 0;
    sumTx        = 0;

    for f = 1:maxFrames
        % --- TX: shaped amplitude bits (fake CCDM) + PAS FEC (parity -> sign) ---
        amp  = pro.draw_amplitude_bits(n, pA, amp_label);   % (n, m-1)
        bits = fec.encode(amp, cfg);                        % (n, m) = [sign, amp]
        x    = pro.map(bits, cstll);                        % (n, 1) REAL ASK (true s.*a)

        % --- HARQ rounds: full n symbols each round; signs admitted incrementally ---
        Lbuf      = zeros(n, cfg.m);   % per-bit LLR buffer
        succRound = 0;
        usedSyms  = 0;

        for r = 1:maxTx
            % All n symbols are transmitted every round (constant channel use).
            [yr, sig2r] = channel.real_channel(x, snr, n);
            Lr          = pro.demap(yr, cstll, sig2r, 'SD');   % (n, m)

            % --- SELECTIVE SIGN ERASURE: admit sign LLR only for released idx ---
            notReleased = ~relMask{r};        % signs still withheld this round
            Lr(notReleased, 1) = 0;           % erase sign (parity) LLR; keep amplitudes

            % Accumulate. Amplitudes accumulate every round (Chase on info);
            % sign column accumulates only where already released (IR on parity).
            Lbuf = Lbuf + Lr;
            usedSyms = usedSyms + n;          % constant n per round

            amp_hat = fec.decode(Lbuf, cfg, maxIter);
            if sum(amp_hat(:) ~= amp(:)) == 0     % genie ACK (CRC in practice)
                succRound = r;
                break;
            end
        end

        % --- Bookkeeping ---
        nCw = nCw + 1;
        if succRound > 0
            succCount(succRound) = succCount(succRound) + 1;
            sumDelivered = sumDelivered + K;
            sumTx        = sumTx + succRound;
        else
            nFail = nFail + 1;
            sumTx = sumTx + maxTx;
        end
        sumUsed = sumUsed + usedSyms;

        nRetx = nCw - succCount(1);
        if nFail >= targetCwErr || nRetx >= 200, break; end
    end

    cumSucc = cumsum(succCount);
    out.bler     = 1 - cumSucc / nCw;
    out.thr      = sumDelivered / sumUsed;    % info bits per channel use (uses = n/round)
    out.succRate = sum(succCount) / nCw;
    out.avgTx    = sumTx / nCw;
    out.nCw      = nCw;
    out.nFail    = nFail;
end