%% Gabriel Cabrera
%% Shaping gain under PAS 64-QAM (DVB-S2 LDPC 2/3): nu sweep.
%% nu=0 = uniform (baseline). The gain is measured on the SNR axis as the
%% horizontal gap to AWGN capacity, because at fixed Rc=2/3 each nu has a
%% different net rate R2D = 2*H(A) bits/2D symbol.
clear; rng(7);

m    = 3;
code = 'dvbs2-2/3';
nus  = [0 0.025 0.05 0.075 0.10];

% SNR ranges per nu (centered on each waterfall; refined with the pre-scan).
ranges = { 13.2:0.1:14.0, ...    % nu=0    (uniform, waterfall ~13.5)
           11.6:0.1:12.4, ...    % nu=0.025
           10.0:0.1:10.8, ...    % nu=0.05
            8.4:0.1:9.2,  ...    % nu=0.075
            7.1:0.1:7.9 };       % nu=0.10

maxFrames    = 1200;
targetCwErr  = 80;              % stop the point at this many codeword (FECFRAME) errors
maxLDPCIter  = 50;
blerTarget   = 1e-2;            % target BLER for measuring the threshold

% ---------------- Constellation and FEC ----------------
cstll = pro.dig_mod_ASK(m, "gray");
[amp_label, amps] = pro.get_amplitude_label(cstll);
cfg = fec.pas_config(m, code);
fprintf('FEC: %s | n=%d | nu sweep = %s\n', cfg.name, cfg.n, mat2str(nus));

if isempty(gcp('nocreate')), parpool(6); end

% ---------------- Sweep per nu ----------------
nN = numel(nus);
S = struct('nu',{},'HA',{},'R2D',{},'SNR',{},'bler',{},'berPost',{});
for i = 1:nN
    [pA, px, HA] = pro.build_shaping(nus(i), cstll, amps);
    cc = cstll;  cc.px = px;  R2D = 2*HA;
    snr = ranges{i};  np = numel(snr);
    fprintf('\n=== nu=%.3f | H(A)=%.3f bits/amp | R=%.2f bits/2D ===\n', nus(i), HA, R2D);

    bpost = nan(1,np);  bl = nan(1,np);
    parfor p = 1:np
        [~, b, d] = fec.run_point(snr(p), cfg, cc, pA, amp_label, ...
                          maxFrames, targetCwErr, maxLDPCIter);
        bpost(p) = b;  bl(p) = d;
        fprintf('  nu=%.3f SNR=%5.2f | BERpost=%.2e BLER=%.2e\n', nus(i), snr(p), b, d);
    end
    S(i) = struct('nu',nus(i),'HA',HA,'R2D',R2D, 'SNR',snr, ...
                  'bler',bl, 'berPost',bpost);
end

% ---------------- Thresholds and shaping gain @ target BLER ----------------
% Metric = gap to AWGN capacity on the SNR axis (rate-independent):
%   SNRmin(R) = 2^R - 1  (minimum Shannon SNR for R bits/2D)
%   gap(nu)   = SNR threshold - SNRmin
%   gain      = gap(uniform) - gap(nu)  (what shaping recovers toward capacity)
thSnr = nan(1,nN);  snrMin = nan(1,nN);  gapCap = nan(1,nN);
for i = 1:nN
    thSnr(i)  = interp_threshold(S(i).SNR, S(i).bler, blerTarget);
    snrMin(i) = 10*log10(2^S(i).R2D - 1);     % Shannon AWGN for R2D bits/2D
    gapCap(i) = thSnr(i) - snrMin(i);
end
shGain = gapCap(1) - gapCap;                  % shaping gain vs uniform (nu=0)

fprintf('\n nu     R2D   SNR@BLER  SNRmin  gapCap  shGain[dB]\n');
for i = 1:nN
    fprintf(' %.3f  %.2f   %7.2f   %6.2f  %6.2f   %+5.2f\n', ...
            nus(i), S(i).R2D, thSnr(i), snrMin(i), gapCap(i), shGain(i));
end

% ---------------- Plots ----------------
co = turbo(nN);
labels = arrayfun(@(s) sprintf('\\nu=%.3g (R=%.2f)', s.nu, s.R2D), S, 'uni', 0);

% Fig 1: BLER vs SNR (each nu shifts with its waterfall -> rate-confounded, NOT a fair gain)
figure('Name','PAS: BLER vs SNR (nu sweep)','Color','w'); hold on; grid on;
for i = 1:nN
    yb = S(i).bler; yb(yb==0) = NaN;
    semilogy(S(i).SNR, yb, '-o', 'Color', co(i,:), 'LineWidth', 1.4, 'DisplayName', labels{i});
end
set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('BLER');
legend('Location','southwest'); title('PAS 64-QAM, DVB-S2 2/3 — BLER vs SNR');

% Fig 2: throughput R vs SNR threshold against capacity (shaping gain = horizontal gap to capacity)
snrGrid = 0:0.05:24;
figure('Name','PAS: R vs SNR threshold vs capacity (shaping gain)','Color','w'); hold on; grid on;
plot(snrGrid, log2(1 + 10.^(snrGrid/10)), 'k--', 'LineWidth', 1.3, ...
     'DisplayName', 'AWGN capacity log_2(1+SNR)');
plot(thSnr, [S.R2D], '-o', 'LineWidth', 1.6, 'MarkerFaceColor', 'auto', ...
     'DisplayName', 'PAS (sweep \nu)');
for i = 1:nN
    text(thSnr(i)+0.15, S(i).R2D-0.12, sprintf('\\nu=%.3g', nus(i)), 'FontSize', 8);
end
xlabel('SNR needed for BLER = 10^{-2}  [dB]'); ylabel('throughput R  [bits / 2D symbol]');
legend('Location','northwest'); xlim([4 16]); ylim([0 7]);
title('PAS 64-QAM, DVB-S2 2/3 — shaping gain = horizontal gap to capacity');

% Fig 3: gap to capacity and shaping gain (correct) vs nu
figure('Name','PAS: gap to capacity and shaping gain vs nu','Color','w');
yyaxis left;  plot(nus, gapCap, '-o', 'LineWidth', 1.6); ylabel('gap to AWGN capacity [dB]');
yyaxis right; plot(nus, shGain, '-s', 'LineWidth', 1.6); ylabel('shaping gain vs uniform [dB]');
grid on; xlabel('\nu (shaping strength)');
title(sprintf('Gap to capacity and shaping gain @ BLER=%.0e (DVB-S2 2/3)', blerTarget));

%% ---------------- helper ----------------
function th = interp_threshold(snr, bler, target)
    ok = isfinite(bler) & bler > 0;  snr = snr(ok);  bler = bler(ok);
    [snr, idx] = sort(snr);  bler = bler(idx);
    th = NaN;
    if numel(snr) < 2 || min(bler) > target || max(bler) < target, return; end
    lbler = log10(bler);  lt = log10(target);
    for k = 1:numel(snr)-1
        if (bler(k)-target)*(bler(k+1)-target) <= 0
            th = snr(k) + (lt - lbler(k))*(snr(k+1)-snr(k))/(lbler(k+1)-lbler(k));
            return;
        end
    end
end
