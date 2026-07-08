function out = run_point(snr, cfg, cstll, pA, amp_label, sch, maxFrames, targetCwErr, maxIter)
% RUN_POINT  Monte-Carlo of one SNR point for a HARQ chain over PAS.
%
%   out = harq.run_point(snr, cfg, cstll, pA, amp_label, sch, ...
%                        maxFrames, targetCwErr, maxIter)   
%
%   Simulates one real-dimension PAS codeword per frame under HARQ. The SAME
%   loop serves BOTH strategies: 'IR' and 'CC' differ only through the transmit
%   schedule `sch` (from harq.build_schedule). Each round transmits the symbols
%   sch.txSets{r}, the receiver ACCUMULATES per-bit LLRs into a buffer and
%   re-decodes; the codeword stops as soon as it decodes (genie ACK).
%
%   The unit is one codeword = one real dimension (I or Q), consistent with the
%   rest of the framework. One real branch of an m-ASK / QAM constellation: the
%   real-dimension SNR in dB equals the complex Es/N0 in dB.
%
%   Output struct `out`:
%     bler      : (1 x maxTx) residual BLER after each round (P[not decoded yet])
%     thr       : throughput in info bits per channel use (transmitted symbol)
%     succRate  : fraction of codewords decoded within maxTx rounds
%     avgTx     : average number of transmissions per codeword
%     nCw       : codewords simulated
%     nFail     : codewords still in error after maxTx rounds

    if nargin < 9, maxIter = 20; end

    maxTx  = sch.maxTx;
    txSets = sch.txSets;
    n      = cfg.n;
    K      = cfg.K;                     

    succCount    = zeros(1, maxTx);     % codewords first decoded at round r
    nFail        = 0;
    nCw          = 0;
    sumDelivered = 0;                   % info bits successfully delivered
    sumUsed      = 0;                   % channel uses (symbols) spent
    sumTx        = 0;                   % transmissions spent (for avgTx)

    for f = 1:maxFrames
        % --- TX: shaped amplitude bits (fake CCDM) + PAS FEC (parity -> sign) ---
        amp  = pro.draw_amplitude_bits(n, pA, amp_label);   % (n, m-1)
        bits = fec.encode(amp, cfg);                        % (n, m) = [sign, amp]
        x    = pro.map(bits, cstll);                        % (n, 1) REAL ASK

        % --- HARQ rounds: accumulate LLRs until the codeword decodes ---
        Lbuf      = zeros(n, cfg.m);    % per-bit LLR buffer (0 = not yet / erased)
        succRound = 0;
        usedSyms  = 0;
        for r = 1:maxTx
            idx = txSets{r};
            if ~isempty(idx)
                xr          = x(idx);
                [yr, sig2r] = channel.real_channel(xr, snr, numel(idx));
                Lr          = pro.demap(yr, cstll, sig2r, 'SD');   % (|idx|, m)
                Lbuf(idx,:) = Lbuf(idx,:) + Lr;   % CC: repeats -> sum; IR: fills new
                usedSyms    = usedSyms + numel(idx);
            end

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

        % corta si juntaste suficientes fallos finales, O si ya viste
        % suficientes codewords que necesitaron retransmisión (estadística de rondas)
        nRetx = nCw - succCount(1);   % codewords que NO decodificaron en ronda 1
        if nFail >= targetCwErr || nRetx >= 200, break; end

    end

    % Residual BLER after each round = 1 - (cumulative first-successes)/nCw.
    cumSucc = cumsum(succCount);

    out.bler     = 1 - cumSucc / nCw;
    out.thr      = sumDelivered / sumUsed;      % info bits per transmitted symbol
    out.succRate = sum(succCount) / nCw;
    out.avgTx    = sumTx / nCw;
    out.nCw      = nCw;
    out.nFail    = nFail;
end
