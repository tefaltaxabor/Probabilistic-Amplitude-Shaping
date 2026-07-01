%% Gabriel Cabrera
%% 64-QAM PAS with FEC -- SNR sweep

%% check list
% (hecho) hacer con ASK y hacerlo doble, (para mostrar equivalencia con 2XASK -QAM)

% BMD(hacerlo...) bit metric decoding rate,

% fabian steiner
% page 62 (red curve)
% page 83 (ccdm 8-ask) (hecho)
% page 123 (curva azul) 

%44 CCDM, rate...
%% 

% Shaping for Wireless
% HARQ (not wifi) and ARQ
% encontrar el benchmark para comparar

%%
´% BCM, information theoreticals comparison.
clear; rng(7);

% ---------------- Parameters ----------------
m       = 3;                         
nu      = 0.05;                      % Maxwell-Boltzmann parameter
SNR_dB  = 6:0.3:15;              

maxFrames     = 10000;                 
targetCwErr   = 40;                 
maxLDPCIter   = 50;                 
% ---------------- Constellation and shaping ----------------
cstll = pro.dig_mod_ASK(m, "gray");
[amp_label, amps] = pro.get_amplitude_label(cstll);

[pA,px,~]= pro.build_shaping(nu,cstll,amps);

cstll.px = px;                       % shaped priors for the demapper

% ---------------- FEC code ----------------
cfg = fec.pas_config(m, 'dvbs2-2/3');
n   = cfg.n;                         % QAM symbols per codeword and dimension
fprintf('FEC: %s  | N=%d K=%d Rc=%.4f | n=%d symbols/dim/codeword\n', ...
        cfg.code, cfg.N, cfg.K, cfg.Rc, n);

% ---------------- Sweep ----------------
nPts    = numel(SNR_dB);
berPre  = nan(1, nPts);
berPost = nan(1, nPts);
bler    = nan(1, nPts);
% --- example frame for shaping statistics (independent of SNR) ---
x_ex = pro.map(fec.encode(pro.draw_amplitude_bits(n,pA,amp_label), cfg), cstll) ...
     + 1j*pro.map(fec.encode(pro.draw_amplitude_bits(n,pA,amp_label), cfg), cstll);

%% Parallel pool
pool = gcp('nocreate');
if isempty(pool)
    parpool(6);
end
%%
parfor p = 1:nPts
    snr = SNR_dB(p);
    bitErrPre = 0; nBitsPre = 0;
    bitErrPost = 0; nBitsPost = 0;
    cwErr = 0; nCw = 0;
    t0 = tic;

    for f = 1:maxFrames
        % --- TX: shaped amplitude bits per dimension (fake CCDM) ---
        ampI = pro.draw_amplitude_bits(n, pA, amp_label);   % (n, m-1)
        ampQ = pro.draw_amplitude_bits(n, pA, amp_label);

        % --- PAS FEC: systematic = amplitude, parity = sign ---
        bitsI = fec.encode(ampI, cfg);                      % (n, m)
        bitsQ = fec.encode(ampQ, cfg);

        % --- Mapping to 64-QAM ---
        xI = pro.map(bitsI, cstll);                         % (n,1)
        xQ = pro.map(bitsQ, cstll);
        x  = xI + 1j*xQ;

        % --- Complex AWGN channel ---
        [y, sigma2] = channel.complex_channel(x, snr, n);

        % --- RX: demap per dimension (noise sigma2/2 per real dimension) ---
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

        % --- Post-FEC BER and BLER (over the systematic amplitude bits) ---
        eI = sum(ampI_hat(:) ~= ampI(:));
        eQ = sum(ampQ_hat(:) ~= ampQ(:));
        bitErrPost = bitErrPost + eI + eQ;
        nBitsPost  = nBitsPost + numel(ampI) + numel(ampQ);
        cwErr = cwErr + (eI > 0) + (eQ > 0);   % +1 per failed codeword (FECFRAME)
        nCw   = nCw + 2;                        % 2 codewords (I and Q) per frame

        if cwErr >= targetCwErr, break; end
    end

    bpre = bitErrPre / nBitsPre;  bpost = bitErrPost / nBitsPost;  bl = cwErr / nCw;
    berPre(p)  = bpre;
    berPost(p) = bpost;
    bler(p)    = bl;
    fprintf(['SNR=%4.1f dB | BERpre=%.3e | BERpost=%.3e | BLER=%.3e ' ...
             '(%d/%d codewords, %d frames, %.1fs)\n'], ...
             snr, bpre, bpost, bl, cwErr, nCw, nCw/2, toc(t0));
end
%%

% ---------------- Results ----------------
T = table(SNR_dB.', berPre.', berPost.', bler.', ...
          'VariableNames', {'SNR_dB','BER_pre','BER_post','BLER'});
disp(T);

figure('Name','PAS 64-QAM FEC sweep','Color','w');
berPostPlot = berPost; berPostPlot(berPostPlot==0) = NaN;   % avoid log(0)
blerPlot    = bler;    blerPlot(blerPlot==0)       = NaN;
semilogy(SNR_dB, berPre, '-o', SNR_dB, berPostPlot, '-s', ...
         SNR_dB, blerPlot, '-^', 'LineWidth', 1.3); grid on;
xlabel('SNR [dB]'); ylabel('error rate');
legend('BER pre-FEC','BER post-FEC','BLER','Location','southwest');
title(sprintf('PAS 64-QAM, %s, \\nu=%.2g', cfg.code, nu));

% ---------------- Shaping statistics (example frame) ----------------
if ~isempty(x_ex)
    S = src.stats(x_ex, cstll, pA, amps);
end


%%
%% 16-ASK PAS -- DVB-S2 3/4, comparacion Fig. 4.10b
clear; rng(7);

m       = 4;
nu      = 0.0145;                    % H(A)≈2.55
SNR_dB  = 15:0.25:18;
maxFrames   = 3000;
targetCwErr = 50;
maxLDPCIter = 50;

cstll = pro.dig_mod_ASK(m, "gray");
[amp_label, amps] = pro.get_amplitude_label(cstll);
[pA,px,~] = pro.build_shaping(nu, cstll, amps);
cstll.px  = px;
fprintf('H(A) = %.4f bpcu\n', -sum(pA(pA>0).*log2(pA(pA>0))));

cfg = fec.pas_config(m, 'dvbs2-3/4');
n   = cfg.n;
fprintf('FEC: %s | N=%d K=%d Rc=%.4f | n=%d\n', cfg.name, cfg.N, cfg.K, cfg.Rc, n);

nPts = numel(SNR_dB);
berPre = nan(1,nPts); berPost = nan(1,nPts); fer = nan(1,nPts);

pool = gcp('nocreate'); if isempty(pool), parpool(6); end

parfor p = 1:nPts
    snr = SNR_dB(p);
    bitErrPre=0; nBitsPre=0; bitErrPost=0; nBitsPost=0; cwErr=0; nCw=0;
    t0=tic;
    for f = 1:maxFrames
        amp  = pro.draw_amplitude_bits(n, pA, amp_label);   % (n,m-1)
        bits = fec.encode(amp, cfg);                        % (n,m)=[sign,amp]
        x    = pro.map(bits, cstll);                        % (n,1) REAL

        [y, sigma2] = channel.real_channel(x, snr, n);      % REAL, sigma2 completo
        llr = pro.demap(y, cstll, sigma2, 'SD');            % sigma2 COMPLETO

        hd = uint8(llr < 0);
        bitErrPre = bitErrPre + sum(hd(:) ~= uint8(bits(:)));
        nBitsPre  = nBitsPre + numel(bits);

        amp_hat = fec.decode(llr, cfg, maxLDPCIter);
        e = sum(amp_hat(:) ~= amp(:));
        bitErrPost = bitErrPost + e;
        nBitsPost  = nBitsPost + numel(amp);
        cwErr = cwErr + (e>0); nCw = nCw + 1;
        if cwErr >= targetCwErr, break; end
    end
    berPre(p)=bitErrPre/nBitsPre; berPost(p)=bitErrPost/nBitsPost; fer(p)=cwErr/nCw;
    fprintf('SNR=%4.1f | BERpre=%.2e | FER=%.2e (%d/%d, %.1fs)\n', ...
            snr, berPre(p), fer(p), cwErr, nCw, toc(t0));
end

figure('Color','w');
ferPlot=fer; ferPlot(ferPlot==0)=NaN;
semilogy(SNR_dB, ferPlot, '-o','LineWidth',1.3,'Color',[0 0.45 0.85]); grid on;
xlabel('SNR [dB]'); ylabel('FER'); ylim([1e-5 1]);
legend('16-ASK PAS, DVB-S2 3/4','Location','southwest');
title('16-ASK PAS, R_{tx}≈2.5 bpcu (regimen Fig. 4.10b)');