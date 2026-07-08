function out = run_point_ccdm(snr, cfg, cstll, ccdm, amp_label, sch, ...
                              maxFrames, targetCwErr, maxIter)
% RUN_POINT_CCDM  Monte-Carlo of one SNR point for a HARQ chain over PAS,
%                 using a REAL invertible CCDM (not the fake amplitude draw).
%
%   out = harq.run_point_ccdm(snr, cfg, cstll, ccdm, amp_label, sch, ...
%                             maxFrames, targetCwErr, maxIter)
%
%   Difference from run_point.m (Approach A, fake CCDM):
%     TX:  uniform info bits --ccdm_encode--> shaped amplitudes
%          --amp_to_bits--> amp_bits --fec.encode--> PAS codeword
%     RX:  LDPC decode --> amp_bits_hat --bits_to_amp--> amplitudes_hat
%          --ccdm_decode--> info_bits_hat
%     ACK: compare INFO BITS (ccdm input), not amplitudes. This is the honest
%          information word-error criterion (Bocherer eq. 2.1): errors confined
%          to parity (signs) do not count as info errors.
%
%   Throughput is measured in INFORMATION bits per channel use, where the
%   information is the CCDM input (k bits per CCDM call), i.e. the truly
%   delivered payload after shaping -- rigorous, unlike counting K systematic
%   bits with the fake matcher.
%
%   Inputs
%   ------
%   ccdm : struct from pro.ccdm_init (fields: n (DM length), M, amps, comp, k, ...)
%          Requires n_sym = cfg.n to be a multiple of ccdm.n (no DM sequence
%          shared across codewords; Bocherer Ex. 2.1).
%   sch  : symbol-puncturing schedule (harq.build_schedule), Approach A style.
%          (A sign-IR CCDM variant would use run_point_signIR logic instead.)
%
%   The MEX functions ccdm_encode_mex / ccdm_decode_mex must be on the path.
%
%   Output struct `out`: bler, thr, succRate, avgTx, nCw, nFail  (as before),
%   plus out.infoBitsPerCw = number of info bits delivered per successful cw.

    if nargin < 9, maxIter = 20; end
    maxTx  = sch.maxTx;
    txSets = sch.txSets;
    n      = cfg.n;                 % symbols (= amplitudes) per codeword
    m      = cfg.m;

    % --- CCDM blocking: how many DM calls fill one codeword ---
    nDM = ccdm.n;
    assert(mod(n, nDM) == 0, ...
        'cfg.n (%d) must be a multiple of ccdm.n (%d)', n, nDM);
    nBlocks = n / nDM;              % CCDM calls per codeword
    kDM     = ccdm.k;               % info bits per CCDM call
    infoBitsPerCw = nBlocks * kDM;  % total info bits carried by one codeword
    comp = ccdm.comp;
    amps = ccdm.amps;

    succCount    = zeros(1, maxTx);
    nFail        = 0;
    nCw          = 0;
    sumDelivered = 0;               % INFO bits successfully delivered
    sumUsed      = 0;               % channel uses (symbols)
    sumTx        = 0;

    for f = 1:maxFrames
        % ================= TX =================
        % 1. Uniform information bits -> shaped amplitudes via REAL CCDM.
        infoBits = randi([0 1], nBlocks, kDM);     % (nBlocks x kDM) payload
        ampSeq   = zeros(1, n);
        for bIdx = 1:nBlocks
            aBlk = pro.ccdm_encode_mex(infoBits(bIdx,:), comp, amps);  % (1 x nDM)
            ampSeq((bIdx-1)*nDM + (1:nDM)) = aBlk;
        end
        % 2. Amplitudes -> amplitude bits (systematic part), per PAS convention.
        amp_bits = local_amp_to_bits(ampSeq(:), amp_label, amps, m); % (n x m-1)
        % 3. PAS encode: parity -> sign.
        bits = fec.encode(amp_bits, cfg);           % (n x m) = [sign, amp]
        x    = pro.map(bits, cstll);                % (n x 1) real ASK

        % ================= HARQ rounds =================
        Lbuf      = zeros(n, m);
        succRound = 0;
        usedSyms  = 0;
        for r = 1:maxTx
            idx = txSets{r};
            if ~isempty(idx)
                xr          = x(idx);
                [yr, sig2r] = channel.real_channel(xr, snr, numel(idx));
                Lr          = pro.demap(yr, cstll, sig2r, 'SD');
                Lbuf(idx,:) = Lbuf(idx,:) + Lr;
                usedSyms    = usedSyms + numel(idx);
            end
            amp_hat = fec.decode(Lbuf, cfg, maxIter);   % (n x m-1) amp bits

            % ---- HONEST GENIE ACK: recover info bits through inverse CCDM ----
            ampSeq_hat = local_bits_to_amp(amp_hat, amp_label, amps, m); % (1 x n)
            ok = true;
            for bIdx = 1:nBlocks
                seg = ampSeq_hat((bIdx-1)*nDM + (1:nDM));
                % composition must match; if LDPC produced an invalid amplitude
                % multiset, ccdm_decode is undefined -> treat as failure.
                if ~local_valid_composition(seg, comp, amps)
                    ok = false; break;
                end
                ib_hat = pro.ccdm_decode_mex(seg, comp, amps, kDM);
                if ~isequal(ib_hat, infoBits(bIdx,:))
                    ok = false; break;
                end
            end
            if ok
                succRound = r;
                break;
            end
        end

        % ================= Bookkeeping =================
        nCw = nCw + 1;
        if succRound > 0
            succCount(succRound) = succCount(succRound) + 1;
            sumDelivered = sumDelivered + infoBitsPerCw;   % real info bits
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
    out.thr           = sumDelivered / sumUsed;   % INFO bits per channel use
    out.succRate      = sum(succCount) / nCw;
    out.avgTx         = sumTx / nCw;
    out.nCw           = nCw;
    out.nFail         = nFail;
    out.infoBitsPerCw = infoBitsPerCw;
end

% ======================================================================
%  Chain-specific converters. These depend on YOUR amp_label convention.
%  Wire them to the same mapping pro.draw_amplitude_bits / pro.map use.
%  The stubs below assume amp_label is an (M x (m-1)) bit matrix whose row i
%  is the (m-1)-bit label of amplitude amps(i) (the amplitude Gray labeling).
% ======================================================================

function amp_bits = local_amp_to_bits(ampSeq, amp_label, amps, m)
% amplitudes (n x 1) -> amplitude bits (n x m-1) using amp_label rows.
    n = numel(ampSeq);
    amp_bits = zeros(n, m-1);
    [tf, loc] = ismember(ampSeq, amps);
    assert(all(tf), 'amp_to_bits: amplitude not in alphabet');
    amp_bits = amp_label(loc, :);       % (n x m-1)
end

function ampSeq = local_bits_to_amp(amp_bits, amp_label, amps, m)
% amplitude bits (n x m-1) -> amplitudes (1 x n): inverse label lookup.
    n = size(amp_bits, 1);
    ampSeq = zeros(1, n);
    % match each row of amp_bits to a row of amp_label
    for j = 1:n
        [tf, loc] = ismember(amp_bits(j,:), amp_label, 'rows');
        if tf
            ampSeq(j) = amps(loc);
        else
            ampSeq(j) = NaN;   % invalid label (shouldn't happen post-decode)
        end
    end
end

function ok = local_valid_composition(seg, comp, amps)
% true if amplitude segment seg has exactly the CCDM composition comp.
    ok = true;
    if any(isnan(seg)), ok = false; return; end
    for i = 1:numel(amps)
        if sum(seg == amps(i)) ~= comp(i), ok = false; return; end
    end
end