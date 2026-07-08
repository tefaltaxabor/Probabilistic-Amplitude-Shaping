function out = run_point_ccdm(snr, cfg, cstll, ccdm, amp_label, sch, ...
                             maxFrames, targetCwErr, maxIter, lut)
% RUN_POINT_CCDM  Monte-Carlo of one SNR point for a HARQ chain over PAS with a
%                 REAL invertible CCDM (Approach A: whole-symbol puncturing).
%
%   out = harq.run_point_ccdm(snr, cfg, cstll, ccdm, amp_label, sch, ...
%                             maxFrames, targetCwErr, maxIter, lut)
%
%   TX:  uniform info bits --ccdm_encode_mex--> shaped amplitudes
%        --(LUT)--> amp_bits --fec.encode--> PAS codeword --map--> symbols
%   HARQ: Approach A. Round r transmits symbols sch.txSets{r}; the receiver
%        accumulates per-bit LLRs and re-decodes until the codeword decodes.
%   ACK: HONEST info criterion. After LDPC decode, recover amplitudes (LUT,
%        vectorized) and info bits (ccdm_decode_mex), compare INFO bits.
%        Shortcut: if decoded amp_bits == transmitted, info matches for sure,
%        so ccdm_decode is skipped (the expensive inverse runs only on frames
%        that actually differ).
%
%   Speed notes (from profiling):
%     - amp<->bit conversions use the precomputed LUT `lut` (no ismember).
%     - The LDPC decode dominates, as it should.
%
%   `lut` struct (build once, pass in): fields
%       wbin      : (m-1 x 1) binary weights
%       bits2amp  : (2^(m-1) x 1) label value -> amplitude
%       amp2idx   : (max(amps) x 1) amplitude value -> row index in amp_label
%
%   Output: bler (1 x maxTx), thr (INFO bits/channel use), succRate, avgTx,
%           nCw, nFail, infoBitsPerCw.

    if nargin < 9 || isempty(maxIter), maxIter = 25; end
    maxTx  = sch.maxTx;
    txSets = sch.txSets;
    n      = cfg.n;
    m      = cfg.m;

    nDM     = ccdm.n;
    assert(mod(n,nDM)==0, 'cfg.n (%d) must be a multiple of ccdm.n (%d)', n, nDM);
    nBlocks = n / nDM;
    kDM     = ccdm.k;
    comp    = ccdm.comp;
    amps    = ccdm.amps;
    infoBitsPerCw = nBlocks * kDM;

    wbin     = lut.wbin;
    bits2amp = lut.bits2amp;
    amp2idx  = lut.amp2idx;

    succCount    = zeros(1, maxTx);
    nFail        = 0;
    nCw          = 0;
    sumDelivered = 0;
    sumUsed      = 0;
    sumTx        = 0;

    for f = 1:maxFrames
        % ================= TX (REAL CCDM) =================
        infoBits = randi([0 1], nBlocks, kDM);
        ampSeq   = zeros(n, 1);
        for bIdx = 1:nBlocks
            aBlk = pro.ccdm_encode_mex(infoBits(bIdx,:), comp, amps);
            ampSeq((bIdx-1)*nDM + (1:nDM)) = aBlk(:);
        end
        idx      = amp2idx(ampSeq);              % vectorized amp -> row index
        amp_bits = amp_label(idx, :);            % (n x m-1)
        bits = fec.encode(amp_bits, cfg);        % (n x m) = [sign, amp]
        x    = pro.map(bits, cstll);             % (n x 1) real ASK

        % ================= HARQ rounds =================
        Lbuf      = zeros(n, m);
        succRound = 0;
        usedSyms  = 0;
        for r = 1:maxTx
            id = txSets{r};
            if ~isempty(id)
                xr          = x(id);
                [yr, sig2r] = channel.real_channel(xr, snr, numel(id));
                Lr          = pro.demap(yr, cstll, sig2r, 'SD');
                Lbuf(id,:)  = Lbuf(id,:) + Lr;
                usedSyms    = usedSyms + numel(id);
            end
            amp_hat = fec.decode(Lbuf, cfg, maxIter);      % (n x m-1)

            % ---- HONEST ACK with shortcut ----
            if isequal(amp_hat, amp_bits)
                succRound = r; break;              % amp bits perfect -> info ok
            end
            % amp bits differ: recover amplitudes (vectorized) and check info
            labelValHat = double(amp_hat) * wbin;
            ampHat = bits2amp(labelValHat + 1);
            ok = true;
            for bIdx = 1:nBlocks
                seg = ampHat((bIdx-1)*nDM + (1:nDM));
                good = true;
                for i = 1:numel(amps)
                    if sum(seg==amps(i)) ~= comp(i), good=false; break; end
                end
                if ~good, ok=false; break; end     % invalid composition
                ib = pro.ccdm_decode_mex(seg, comp, amps, kDM);
                if ~isequal(ib, infoBits(bIdx,:)), ok=false; break; end
            end
            if ok, succRound = r; break; end
        end

        % ================= Bookkeeping =================
        nCw = nCw + 1;
        if succRound > 0
            succCount(succRound) = succCount(succRound) + 1;
            sumDelivered = sumDelivered + infoBitsPerCw;
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
    out.bler          = 1 - cumSucc / nCw;
    out.thr           = sumDelivered / sumUsed;   % INFO bits / channel use
    out.succRate      = sum(succCount) / nCw;
    out.avgTx         = sumTx / nCw;
    out.nCw           = nCw;
    out.nFail         = nFail;
    out.infoBitsPerCw = infoBitsPerCw;
end