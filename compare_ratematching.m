%% Gabriel Cabrera
% Rate matching for 6G via distribution matching.
% Compares two ways to adapt the information rate R [bits / 2D symbol] to the
% channel SNR, both measured as the SNR needed to reach a target BLER:
%
%   (1) DM-based rate matching  -> fixed modulation (64-QAM) and fixed FEC
%       (DVB-S2 2/3); the rate is tuned CONTINUOUSLY by the distribution
%       matcher via nu  =>  R2D = 2*H(A) in (0, 4].
%   (2) Conventional adaptive modulation+coding -> uniform constellation,
%       switch order 16/64/256-QAM (DVB-S2 1/2, 2/3, 3/4) => R2D in {2,4,6},
%       a coarse STAIRCASE (this is the 5G-style MCS behaviour).
%
% The fake CCDM is valid here: every quantity depends only on the amplitude
%%
clear; rng(7);

blerTarget   = 1e-2;     % target BLER for the SNR threshold
%maxFrames = 1200;
maxFrames    = 1e5;
targetCwErr = 100;
%targetCwErr  = 80;      
maxLDPCIter  = 50;

if isempty(gcp('nocreate')), parpool(6); end

%% ---- DM-based rate matching (64-QAM + DVB-S2 2/3, sweep nu) --------
dm.m      = 3;
dm.code   = 'dvbs2-2/3';
dm.nus    = [0 0.025 0.05 0.075 0.10];
dm.ranges = { 13.2:0.1:14.2, 11.5:0.1:12.6, 9.9:0.1:11.0, 8.3:0.1:9.4, 7.0:0.1:8.2 };

%% ---- Conventional adaptive mod+code (uniform, switch order) --------
cv.m      = [2 3 4];
cv.code   = {'dvbs2-1/2','dvbs2-2/3','dvbs2-3/4'};   % 16 / 64 / 256-QAM, uniform
cv.ranges = { 5.5:0.2:8.5, 13.2:0.1:14.2, 18.5:0.2:21.5 }; %aumentar el snr 

%% ---- Sweep -----------------------------------------------------

% el dm
nD = numel(dm.nus);
dmR = nan(1,nD); dmTh = nan(1,nD);
for i = 1:nD
    [dmR(i), dmTh(i)] = sweep_config(dm.m, dm.code, dm.nus(i), dm.ranges{i}, ...
                          blerTarget, maxFrames, targetCwErr, maxLDPCIter);
    fprintf('[DM ] nu=%.3f  R2D=%.3f  SNR@BLER=%.2f dB\n', dm.nus(i), dmR(i), dmTh(i));
end

%el normal
nC = numel(cv.m);
cvR = nan(1,nC); cvTh = nan(1,nC);
for j = 1:nC
    [cvR(j), cvTh(j)] = sweep_config(cv.m(j), cv.code{j}, 0, cv.ranges{j}, ...
                          blerTarget, maxFrames, targetCwErr, maxLDPCIter);
    fprintf('[CONV] %s (uniform)  R2D=%.3f  SNR@BLER=%.2f dB\n', cv.code{j}, cvR(j), cvTh(j));
end

%% ---- Plot: throughput R [bits/2D] vs SNR threshold, against capacity ---------
snrGrid = 0:0.05:24;
Cawgn   = log2(1 + 10.^(snrGrid/10));         % AWGN capacity per 2D symbol

fig = figure('Name','Rate matching: DM continuum vs conventional staircase','Color','w');
hold on; grid on;
plot(snrGrid, Cawgn, 'k--', 'LineWidth', 1.3, 'DisplayName', 'AWGN capacity log_2(1+SNR)');
plot(dmTh, dmR, '-o', 'LineWidth', 1.8, 'MarkerFaceColor','auto', ...
     'DisplayName', 'DM rate matching (64-QAM, DVB-S2 2/3, sweep \nu)');
stairs([cvTh cvTh(end)+3], [cvR cvR(end)], ':', 'LineWidth', 1.2, 'HandleVisibility','off');
plot(cvTh, cvR, 's', 'MarkerSize', 10, 'LineWidth', 1.8, 'MarkerFaceColor','w', ...
     'DisplayName', 'Conventional (uniform 16/64/256-QAM)');
for j = 1:nC
    text(cvTh(j)+0.2, cvR(j), sprintf(' %d-QAM', 2^(2*cv.m(j))), 'FontSize', 9);
end
for i = 1:nD
    text(dmTh(i)+0.15, dmR(i)-0.12, sprintf('\\nu=%.3g', dm.nus(i)), 'FontSize', 8);
end
xlabel('SNR needed for BLER = 10^{-2}  [dB]');
ylabel('throughput R  [bits / 2D symbol]');
title('Rate matching for 6G: distribution matching vs conventional MCS switching');
legend('Location','northwest'); xlim([4 24]); ylim([0 7]);

%% ---- Save figure (.fig) and full workspace (.mat) ---------------------------
if ~isfolder('results'), mkdir('results'); end
savefig(fig, fullfile('results','compare_ratematching.fig'));
save(fullfile('results','compare_ratematching.mat'));

%% ---- Table ------------------------------------------------------------------
fprintf('\n--- DM rate matching (one FEC, one modulation) ---\n');
fprintf('  nu      R2D     SNR@BLER[dB]\n');
for i = 1:nD, fprintf('  %.3f   %.3f   %7.2f\n', dm.nus(i), dmR(i), dmTh(i)); end
fprintf('\n--- Conventional (switch modulation+code) ---\n');
fprintf('  code          R2D     SNR@BLER[dB]\n');
for j = 1:nC, fprintf('  %-12s  %.3f   %7.2f\n', cv.code{j}, cvR(j), cvTh(j)); end

%% ============================== helpers ======================================
function [R2D, thSnr] = sweep_config(m, code, nu, snrRange, blerTarget, maxFrames, targetCwErr, maxIter)
    
    cstll = pro.dig_mod_ASK(m, "gray");
    [amp_label, amps] = pro.get_amplitude_label(cstll);
    cfg = fec.pas_config(m, code);
    [pA, px, HA] = pro.build_shaping(nu, cstll, amps);
    cc = cstll; cc.px = px;
    R2D = 2*HA;

    np = numel(snrRange);
    bler = nan(1,np);
    parfor p = 1:np
        [~,~,d] = fec.run_point(snrRange(p), cfg, cc, pA, amp_label, ...
                                maxFrames, targetCwErr, maxIter);
        bler(p) = d;
    end
    thSnr = interp_threshold(snrRange, bler, blerTarget);

end

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


