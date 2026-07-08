%% Gabriel Cabrera
%% 64-QAM PAS with FEC and REAL CCDM -- SNR sweep (no HARQ) -- OPTIMIZED
%
%  Speed fixes vs. the previous version:
%   (1) NO ismember. amp<->bit conversions use precomputed lookup tables and
%       vectorized indexing (thousands of times faster than row-wise ismember).
%   (2) ACK shortcut: ccdm_decode is skipped whenever the decoded amplitude
%       bits equal the transmitted ones (info then matches for sure). The CCDM
%       inverse runs only on frames that actually differ.
%  Net effect: the LDPC decoder (fec.decode) becomes the dominant cost, as it
%  should be.
%
%  Prereqs: pro.ccdm_encode_mex / pro.ccdm_decode_mex compiled and on the path.

clear; rng(7);

% ---------------- Parameters ----------------
m       = 3;
nu      = 0.05;
SNR_dB  = 6:0.3:15;

maxFrames   = 2e4;
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
nDM = 200;
assert(mod(n,nDM)==0, 'nDM must divide n');
ccdm    = pro.ccdm_init(pA, amps, nDM);
nBlocks = n / nDM;
kDM     = ccdm.k;
comp    = ccdm.comp;

% ---------------- PRECOMPUTED LOOKUP TABLES (replaces ismember) ----------
% amp_label is (M x m-1): row i is the (m-1)-bit label of amps(i).
% Interpret each label as a binary number -> index into a LUT.
wbin     = (2.^(size(amp_label,2)-1:-1:0)).';    % binary weights, e.g. [2;1]
labelVal = double(amp_label) * wbin;                     % (M x 1) decimal value per label
bits2amp = zeros(2^(m-1), 1);
bits2amp(labelVal + 1) = amps;                   % LUT: label value -> amplitude
amp2idx  = zeros(max(amps), 1);
amp2idx(amps) = 1:numel(amps);                   % amplitude value -> 1..M

% sanity: LUT reproduces the row-wise mapping exactly
for i = 1:numel(amps)
    assert(bits2amp(double(amp_label(i,:))*wbin + 1) == amps(i), 'LUT mismatch');
end

fprintf('FEC: %s | N=%d K=%d Rc=%.4f | n=%d symbols/dim\n', ...
        cfg.code, cfg.N, cfg.K, cfg.Rc, n);
fprintf('CCDM: nDM=%d, k=%d, Rccdm=%.4f, Rloss=%.4f bits/amp\n', ...
        ccdm.n, ccdm.k, ccdm.Rccdm, ccdm.Rloss);
fprintf('Info bits/dim=%d | SE_info=%.4f bits/cx-sym (2 dims)\n', ...
        nBlocks*kDM, 2*nBlocks*kDM/n);

% ---------------- Sweep ----------------
nPts    = numel(SNR_dB);
berPre  = nan(1, nPts);
berPost = nan(1, nPts);
bler    = nan(1, nPts);
ferInfo = nan(1, nPts);

%pool = gcp('nocreate'); if isempty(pool), parpool(6); end

for p = 1:nPts
    snr = SNR_dB(p);
    bitErrPre=0; nBitsPre=0; bitErrPost=0; nBitsPost=0;
    cwErr=0; nCw=0; infoErr=0; nInfoCw=0;
    t0 = tic;

    for f = 1:maxFrames
        % --- TX: REAL CCDM per dimension (vectorized amp->bits) ---
        [ampI_bits, infoI] = tx_ccdm(nBlocks, nDM, kDM, comp, amps, amp2idx, amp_label);
        [ampQ_bits, infoQ] = tx_ccdm(nBlocks, nDM, kDM, comp, amps, amp2idx, amp_label);

        bitsI = fec.encode(ampI_bits, cfg);
        bitsQ = fec.encode(ampQ_bits, cfg);

        xI = pro.map(bitsI, cstll);
        xQ = pro.map(bitsQ, cstll);
        x  = xI + 1j*xQ;

        [y, sigma2] = channel.complex_channel(x, snr, n);

        llrI = pro.demap(real(y), cstll, sigma2/2, 'SD');
        llrQ = pro.demap(imag(y), cstll, sigma2/2, 'SD');

        hdI = uint8(llrI < 0);  hdQ = uint8(llrQ < 0);
        bitErrPre = bitErrPre + sum(hdI(:) ~= uint8(bitsI(:))) ...
                              + sum(hdQ(:) ~= uint8(bitsQ(:)));
        nBitsPre  = nBitsPre + numel(bitsI) + numel(bitsQ);

        ampI_hat = fec.decode(llrI, cfg, maxLDPCIter);   % <-- the real cost
        ampQ_hat = fec.decode(llrQ, cfg, maxLDPCIter);

        eI = sum(ampI_hat(:) ~= ampI_bits(:));
        eQ = sum(ampQ_hat(:) ~= ampQ_bits(:));
        bitErrPost = bitErrPost + eI + eQ;
        nBitsPost  = nBitsPost + numel(ampI_bits) + numel(ampQ_bits);
        cwErr = cwErr + (eI > 0) + (eQ > 0);
        nCw   = nCw + 2;

        % --- info FER with ACK shortcut (vectorized amp recovery) ---
        infoErr = infoErr + info_err_fast(ampI_hat, ampI_bits, infoI, eI, ...
                              nBlocks, nDM, kDM, comp, amps, wbin, bits2amp);
        infoErr = infoErr + info_err_fast(ampQ_hat, ampQ_bits, infoQ, eQ, ...
                              nBlocks, nDM, kDM, comp, amps, wbin, bits2amp);
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
%  Helpers -- vectorized, no ismember
% ======================================================================
function [amp_bits, info] = tx_ccdm(nBlocks, nDM, kDM, comp, amps, amp2idx, amp_label)
% One real-dimension frame via REAL CCDM. amp->bits via LUT (no ismember).
    info   = randi([0 1], nBlocks, kDM);
    ampSeq = zeros(nBlocks*nDM, 1);
    for b = 1:nBlocks
        aBlk = pro.ccdm_encode_mex(info(b,:), comp, amps);
        ampSeq((b-1)*nDM + (1:nDM)) = aBlk(:);
    end
    idx      = amp2idx(ampSeq);            % amplitude value -> row index (vectorized)
    amp_bits = amp_label(idx, :);          % (n x m-1)
end

function e = info_err_fast(amp_hat, amp_bits_tx, info, nBitErr, ...
                           nBlocks, nDM, kDM, comp, amps, wbin, bits2amp)
% Info FER with shortcut: if decoded amp bits == transmitted, info matches for
% sure (CCDM is deterministic & invertible) -> no ccdm_decode needed.
    if nBitErr == 0
        e = 0; return;                     % fast path: perfect amp bits
    end
    % Only reached when amp bits differ. Recover amplitudes (vectorized).
    labelValHat = double(amp_hat) * wbin;          % (n x 1) label value per symbol
    % guard: values must be valid label indices
    ampHat = bits2amp(labelValHat + 1);    % (n x 1) amplitudes
    for b = 1:nBlocks
        seg = ampHat((b-1)*nDM + (1:nDM));
        ok = true;
        for i = 1:numel(amps)
            if sum(seg==amps(i)) ~= comp(i), ok=false; break; end
        end
        if ~ok, e = 1; return; end         % invalid composition -> info error
        ib = pro.ccdm_decode_mex(seg, comp, amps, kDM);
        if ~isequal(ib, info(b,:)), e = 1; return; end
    end
    e = 0;
end