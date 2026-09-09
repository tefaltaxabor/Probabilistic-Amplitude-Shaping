function out = run_point_symbolwise(snr, cfg, cstll, ccdm, amp_label, sch, ...
                                    maxFrames, targetCwErr, maxIter, lut, useApriori)
% RUN_POINT_SYMBOLWISE  HARQ over PAS with SYMBOL-WISE puncturing on the label
%   sequence of shaped symbols, following Shen et al., "Symbol-Wise Puncturing
%   for HARQ Integration With PAS", IEEE WCL 2022.
%
%   out = harq.run_point_symbolwise(snr, cfg, cstll, ccdm, amp_label, sch, ...
%                                   maxFrames, targetCwErr, maxIter, lut, useApriori)
%
%   Two key ideas from the paper, adapted to a real 8-ASK dimension (1 sign +
%   (m-1) amplitude bits per symbol; the paper uses M^2-QAM with 2 sign bits):
%
%   1) SYMBOL-WISE PUNCTURING (Sec. III, eq. 1-2). Puncturing is done on whole
%      SYMBOL LABELS, not on random bits. A punctured symbol withholds its
%      entire label (sign + amplitude bits) together, which preserves the
%      shaping distribution of the transmitted symbols in every round. This is
%      what `sch.txSets` already encodes (whole-symbol schedule), so the round
%      loop is unchanged from Approach A.
%
%   2) A-PRIORI INFO ON UNTRANSMITTED AMPLITUDE BITS (Sec. III-B, eq. 9).
%      Instead of filling the LLR buffer of a NOT-yet-transmitted symbol with 0
%      (neutral erasure), we fill the AMPLITUDE-bit LLRs with the shaping
%      a-priori L_{A,j} = log( sum_{a:a_j=0} P(a) / sum_{a:a_j=1} P(a) ).
%      With uniform P this is 0 (recovers the uniform/erasure case); with a
%      shaped P it injects the shaping knowledge, which is the source of the
%      shaping gain the paper reports (~0.6 dB).
%
%      `useApriori` (logical): true -> PAS with a-priori (paper's scheme).
%                              false -> uniform baseline (L_A = 0, erasure).
%
%   `lut` : struct from harq_chain (wbin, bits2amp, amp2idx) for amp<->bit.
%
%   Output: same fields as run_point_ccdm (bler, thr, succRate, avgTx, nCw,
%           nFail, infoBitsPerCw).

    if nargin < 9  || isempty(maxIter),    maxIter = 25;   end
    if nargin < 11 || isempty(useApriori), useApriori = true; end

    maxTx  = sch.maxTx;
    txSets = sch.txSets;
    n      = cfg.n;
    m      = cfg.m;

    nDM     = ccdm.n;
    assert(mod(n,nDM)==0, 'cfg.n must be a multiple of ccdm.n');
    nBlocks = n / nDM;
    kDM     = ccdm.k;
    comp    = ccdm.comp;
    amps    = ccdm.amps;
    infoBitsPerCw = nBlocks * kDM;

    wbin     = lut.wbin;
    bits2amp = lut.bits2amp;
    amp2idx  = lut.amp2idx;

    % ---- Precompute the a-priori LLRs of the AMPLITUDE bits (eq. 9) ----
    % Amplitude bits are columns 2..m of the per-symbol label. For each such
    % bit position, L_A = log( P(bit=0) / P(bit=1) ) under the amplitude
    % distribution pQuant. With uniform amps this is 0.
    %   amp_label(i,:) is the (m-1)-bit label of amplitude amps(i).
    pAmp = ccdm.pQuant(:);                 % (M x 1) amplitude probabilities
    LA_amp = zeros(1, m-1);                % a-priori LLR per amplitude-bit column
    for c = 1:(m-1)
        bitcol = amp_label(:, c);          % (M x 1) value of this bit per amp
        p0 = sum(pAmp(bitcol==0));
        p1 = sum(pAmp(bitcol==1));
        if p0 <= 0, p0 = eps; end
        if p1 <= 0, p1 = eps; end
        LA_amp(c) = log(p0 / p1);
    end
    if ~useApriori
        LA_amp = zeros(1, m-1);            % uniform baseline: no a-priori
    end
    % The sign bit (column 1) carries parity -> no shaping a-priori (L=0).
    LA_full = [0, LA_amp];                 % (1 x m) a-priori row per symbol

    succCount    = zeros(1, maxTx);
    nFail = 0; nCw = 0;
    sumDelivered = 0; sumUsed = 0; sumTx = 0;

    for f = 1:maxFrames
        % ================= TX (REAL CCDM, shaped) =================
        infoBits = randi([0 1], nBlocks, kDM);
        ampSeq   = zeros(n, 1);
        for bIdx = 1:nBlocks
            aBlk = pro.ccdm_encode_mex(infoBits(bIdx,:), comp, amps);
            ampSeq((bIdx-1)*nDM + (1:nDM)) = aBlk(:);
        end
        idx      = amp2idx(ampSeq);
        amp_bits = amp_label(idx, :);
        bits = fec.encode(amp_bits, cfg);
        x    = pro.map(bits, cstll);

        % ================= HARQ rounds with symbol-wise puncturing =========
        % Initialize the LLR buffer with the A-PRIORI info on amplitude bits for
        % symbols not yet transmitted (eq. 9). Sign gets 0. Transmitted symbols
        % overwrite these with channel LLRs (accumulated).
        Lbuf      = repmat(LA_full, n, 1);   % a-priori everywhere initially
        seenSym   = false(n, 1);             % which symbols have been transmitted
        succRound = 0;
        usedSyms  = 0;

        for r = 1:maxTx
            id = txSets{r};
            if ~isempty(id)
                % For symbols transmitted THIS round for the FIRST time, remove
                % the a-priori seed before accumulating channel LLRs (so we do
                % not double-count a-priori + channel on the amplitude bits).
                firstTime = id(~seenSym(id));
                Lbuf(firstTime, :) = 0;
                seenSym(id) = true;

                xr          = x(id);
                [yr, sig2r] = channel.real_channel(xr, snr, numel(id));
                Lr          = pro.demap(yr, cstll, sig2r, 'SD');
                Lbuf(id,:)  = Lbuf(id,:) + Lr;
                usedSyms    = usedSyms + numel(id);
            end
            amp_hat = fec.decode(Lbuf, cfg, maxIter);

            % ---- honest ACK with shortcut ----
            if isequal(amp_hat, amp_bits)
                succRound = r; break;
            end
            labelValHat = double(amp_hat) * wbin;
            ampHat = bits2amp(labelValHat + 1);
            ok = true;
            for bIdx = 1:nBlocks
                seg = ampHat((bIdx-1)*nDM + (1:nDM));
                good = true;
                for i = 1:numel(amps)
                    if sum(seg==amps(i)) ~= comp(i), good=false; break; end
                end
                if ~good, ok=false; break; end
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
    out.thr           = sumDelivered / sumUsed;
    out.succRate      = sum(succCount) / nCw;
    out.avgTx         = sumTx / nCw;
    out.nCw           = nCw;
    out.nFail         = nFail;
    out.infoBitsPerCw = infoBitsPerCw;
    out.useApriori    = useApriori;
end