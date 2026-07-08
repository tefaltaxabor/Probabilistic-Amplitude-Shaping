%% Gabriel Cabrera
%% 64-QAM PAS with FEC and REAL CCDM -- SNR sweep (no HARQ)
%
%  Adapted from the fake-CCDM sweep: draw_amplitude_bits is replaced by a REAL
%  invertible CCDM (ccdm_encode_mex) applied independently to the I and Q
%  dimensions (64-QAM = 2 x 8-ASK). Metrics kept from the original sweep:
%    - BER pre-FEC   (hard decision on all coded bits)
%    - BER post-FEC  (on the systematic amplitude bits)
%    - BLER          (codeword/FECFRAME error rate)
%  ADDED: honest information FER, measured through the inverse CCDM (an error
%  confined to parity/signs does NOT count as an info error; Bocherer eq. 2.1).
%
%  Prereqs: ccdm_encode_mex / ccdm_decode_mex compiled and on the path.
%
%  Checklist ref (Steiner): p.62 red curve, p.83 CCDM 8-ASK, p.123 blue curve.

clear; rng(7);

% ---------------- Parameters ----------------
m       = 3;
nu      = 0.05;                      % Maxwell-Boltzmann parameter
SNR_dB  = 6:0.3:15;

maxFrames   = 1000;
targetCwErr = 40;
maxLDPCIter = 20;

% ---------------- Constellation and shaping ----------------
cstll = pro.dig_mod_ASK(m, "gray");
[amp_label, amps] = pro.get_amplitude_label(cstll);
[pA, px, ~] = pro.build_shaping(nu, cstll, amps);
cstll.px = px;

% ---------------- FEC code ----------------
cfg = fec.pas_config(m, 'dvbs2-2/3');
n   = cfg.n;

% ---------------- REAL CCDM setup ----------------
nDM = 200;                           % must divide n (21600 = 200*108)
assert(mod(n,nDM)==0, 'nDM must divide n');
ccdm    = pro.ccdm_init(pA, amps, nDM);
nBlocks = n / nDM;
kDM     = ccdm.k;
comp    = ccdm.comp;
infoBitsPerDim = nBlocks * kDM;      % info bits carried per real dimension

fprintf('FEC: %s | N=%d K=%d Rc=%.4f | n=%d symbols/dim\n', ...
        cfg.code, cfg.N, cfg.K, cfg.Rc, n);
fprintf('CCDM: nDM=%d, k=%d, Rccdm=%.4f, Rloss=%.4f bits/amp\n', ...
        ccdm.n, ccdm.k, ccdm.Rccdm, ccdm.Rloss);
fprintf('Info bits/dim=%d | SE_info=%.4f bits/cx-sym (2 dims)\n', ...
        infoBitsPerDim, 2*infoBitsPerDim/n);

% ---------------- Sweep ----------------
nPts    = numel(SNR_dB);
berPre  = nan(1, nPts);
berPost = nan(1, nPts);
bler    = nan(1, nPts);
ferInfo = nan(1, nPts);              % NEW: honest information FER

pool = gcp('nocreate'); if isempty(pool), parpool(6); end

parfor p = 1:nPts
    snr = SNR_dB(p);
    bitErrPre=0; nBitsPre=0; bitErrPost=0; nBitsPost=0;
    cwErr=0; nCw=0; infoErr=0; nInfoCw=0;
    t0 = tic;

    for f = 1:maxFrames
        % --- TX: REAL CCDM per dimension -> amplitudes -> amp_bits ---
        [ampI_bits, infoI] = local_ccdm_tx(nBlocks, nDM, kDM, comp, amps, amp_label, m);
        [ampQ_bits, infoQ] = local_ccdm_tx(nBlocks, nDM, kDM, comp, amps, amp_label, m);

        % --- PAS FEC: systematic = amplitude, parity = sign ---
        bitsI = fec.encode(ampI_bits, cfg);                 % (n, m)
        bitsQ = fec.encode(ampQ_bits, cfg);

        % --- Mapping to 64-QAM ---
        xI = pro.map(bitsI, cstll);
        xQ = pro.map(bitsQ, cstll);
        x  = xI + 1j*xQ;

        % --- Complex AWGN channel ---
        [y, sigma2] = channel.complex_channel(x, snr, n);

        % --- RX: demap per dimension (noise sigma2/2 per real dim) ---
        llrI = pro.demap(real(y), cstll, sigma2/2, 'SD');   % (n, m)
        llrQ = pro.demap(imag(y), cstll, sigma2/2, 'SD');

        % --- Pre-FEC BER (hard decision over all bits) ---
        hdI = uint8(llrI < 0);  hdQ = uint8(llrQ < 0);
        bitErrPre = bitErrPre + sum(hdI(:) ~= uint8(bitsI(:))) ...
                              + sum(hdQ(:) ~= uint8(bitsQ(:)));
        nBitsPre  = nBitsPre + numel(bitsI) + numel(bitsQ);

        % --- FEC decode ---
        ampI_hat = fec.decode(llrI, cfg, maxLDPCIter);      % (n, m-1)
        ampQ_hat = fec.decode(llrQ, cfg, maxLDPCIter);

        % --- Post-FEC BER and BLER (systematic amplitude bits) ---
        eI = sum(ampI_hat(:) ~= ampI_bits(:));
        eQ = sum(ampQ_hat(:) ~= ampQ_bits(:));
        bitErrPost = bitErrPost + eI + eQ;
        nBitsPost  = nBitsPost + numel(ampI_bits) + numel(ampQ_bits);
        cwErr = cwErr + (eI > 0) + (eQ > 0);
        nCw   = nCw + 2;

        % --- HONEST info FER through inverse CCDM (per dimension) ---
        infoErr = infoErr + local_ccdm_info_err(ampI_hat, infoI, nBlocks, nDM, kDM, comp, amps, amp_label);
        infoErr = infoErr + local_ccdm_info_err(ampQ_hat, infoQ, nBlocks, nDM, kDM, comp, amps, amp_label);
        nInfoCw = nInfoCw + 2;

        if cwErr >= targetCwErr, break; end
    end

    berPre(p)  = bitErrPre / nBitsPre;
    berPost(p) = bitErrPost / nBitsPost;
    bler(p)    = cwErr / nCw;
    ferInfo(p) = infoErr / nInfoCw;
    fprintf(['SNR=%4.1f dB | BERpre=%.3e | BERpost=%.3e | BLER=%.3e | ' ...
             'FERinfo=%.3e (%d cw, %.1fs)\n'], ...
             snr, berPre(p), berPost(p), bler(p), ferInfo(p), nCw, toc(t0));
end

% ---------------- Results ----------------
T = table(SNR_dB.', berPre.', berPost.', bler.', ferInfo.', ...
          'VariableNames', {'SNR_dB','BER_pre','BER_post','BLER','FER_info'});
disp(T);

figure('Name','PAS 64-QAM real-CCDM sweep','Color','w');
bpp = berPost; bpp(bpp==0)=NaN;
blp = bler;    blp(blp==0)=NaN;
fip = ferInfo; fip(fip==0)=NaN;
semilogy(SNR_dB, berPre, '-o', SNR_dB, bpp, '-s', ...
         SNR_dB, blp, '-^', SNR_dB, fip, '-d', 'LineWidth', 1.3); grid on;
xlabel('SNR [dB]'); ylabel('error rate');
legend('BER pre-FEC','BER post-FEC','BLER','FER info (CCDM)','Location','southwest');
title(sprintf('PAS 64-QAM real CCDM, %s, \\nu=%.2g', cfg.code, nu));

save('results/pas64qam_realccdm.mat','SNR_dB','berPre','berPost','bler','ferInfo','cfg','nu','ccdm');
fprintf('Saved results/pas64qam_realccdm.mat\n');

% ======================================================================
%  Helpers (chain-specific: match your amp_label convention)
% ======================================================================
function [amp_bits, info] = local_ccdm_tx(nBlocks, nDM, kDM, comp, amps, amp_label, m)
% Generate one real-dimension frame via REAL CCDM. Returns amp_bits (n x m-1)
% and the info bits (nBlocks x kDM) for the genie ACK.
    info   = randi([0 1], nBlocks, kDM);
    ampSeq = zeros(nBlocks*nDM, 1);
    for b = 1:nBlocks
        aBlk = pro.ccdm_encode_mex(info(b,:), comp, amps);
        ampSeq((b-1)*nDM + (1:nDM)) = aBlk(:);
    end
    [tf, loc] = ismember(ampSeq, amps);
    assert(all(tf));
    amp_bits = amp_label(loc, :);            % (n x m-1)
end

function e = local_ccdm_info_err(amp_hat, info, nBlocks, nDM, kDM, comp, amps, amp_label)
% Returns 1 if ANY info block fails after inverse CCDM, else 0 (info FER).
    n = size(amp_hat,1);
    ampHat = zeros(1,n);
    for j = 1:n
        [tf, loc] = ismember(amp_hat(j,:), amp_label, 'rows');
        if tf, ampHat(j) = amps(loc); else, e = 1; return; end
    end
    for b = 1:nBlocks
        seg = ampHat((b-1)*nDM + (1:nDM));
        % invalid composition -> undecodable -> info error
        ok = true;
        for i = 1:numel(amps)
            if sum(seg==amps(i)) ~= comp(i), ok=false; break; end
        end
        if ~ok, e = 1; return; end
        ib = pro.ccdm_decode_mex(seg, comp, amps, kDM);
        if ~isequal(ib, info(b,:)), e = 1; return; end
    end
    e = 0;
end