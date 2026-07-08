%% Gabriel Cabrera
%% HARQ chain for PAS with REAL CCDM -- Approach A (whole-symbol puncturing)
%
%  Same HARQ mechanics as harq_chain.m (IR / CC by schedule), but the TX uses a
%  REAL invertible CCDM and the ACK is the honest information criterion (errors
%  confined to parity/signs do NOT count). Throughput is in INFO bits/channel
%  use, shaping-accurate (rate loss of the CCDM included).
%
%  Sequential FOR (not parfor): the LDPC decoder MEX is multi-threaded, so a
%  parfor over SNR points oversubscribes the cores and runs SLOWER. FOR lets
%  each point use the multi-threaded decoder fully.
%
%  Prereqs: pro.ccdm_encode_mex / pro.ccdm_decode_mex compiled and on path.

clear; rng(7);

% ================= HARQ SELECTOR =================
harqMode = 'IR';            % 'IR' | 'CC'
% ================================================

maxTx        = 4;
punctureFrac = 0.25;

m       = 3;
nu      = 0.05;
SNR_dB  = 16:0.5:25;    % round-1 waterfall region

maxFrames   = 3e4;
targetCwErr = 50;
maxLDPCIter = 20;

% ---------------- Constellation and shaping ----------------
cstll = pro.dig_mod_ASK(m, "gray");
[amp_label, amps] = pro.get_amplitude_label(cstll);
[pA, px, ~]       = pro.build_shaping(nu, cstll, amps);
cstll.px = px;

% ---------------- FEC ----------------
cfg = fec.pas_config(m, 'dvbs2-2/3');
n   = cfg.n;

% ---------------- REAL CCDM setup ----------------
nDM  = 200;
assert(mod(n,nDM)==0, 'nDM must divide n');
ccdm = pro.ccdm_init(pA, amps, nDM);

% ---------------- LUT for amp<->bit (no ismember) ----------------
wbin     = (2.^(size(amp_label,2)-1:-1:0)).';
labelVal = double(amp_label) * wbin;
bits2amp = zeros(2^(m-1), 1);
bits2amp(labelVal + 1) = amps;
amp2idx  = zeros(max(amps), 1);
amp2idx(amps) = 1:numel(amps);
for i = 1:numel(amps)
    assert(bits2amp(double(amp_label(i,:))*wbin + 1) == amps(i), 'LUT mismatch');
end
lut = struct('wbin', wbin, 'bits2amp', bits2amp, 'amp2idx', amp2idx);

% ---------------- HARQ schedule (Approach A: whole symbols) ----------------
sch = harq.build_schedule(n, maxTx, harqMode, punctureFrac);

fprintf('HARQ mode: %s | maxTx=%d | mother %s (Rc=%.3f) | n=%d\n', ...
        sch.mode, maxTx, cfg.name, cfg.Rc, n);
fprintf('CCDM: nDM=%d, k=%d, Rccdm=%.4f, Rloss=%.4f | info bits/cw=%d\n', ...
        ccdm.n, ccdm.k, ccdm.Rccdm, ccdm.Rloss, (n/nDM)*ccdm.k);

% ---------------- Sweep (sequential FOR) ----------------
nPts     = numel(SNR_dB);
blerAll  = nan(nPts, maxTx);
thr      = nan(1, nPts);
succRate = nan(1, nPts);
avgTx    = nan(1, nPts);

for p = 1:nPts
    snr = SNR_dB(p);
    t0  = tic;
    out = harq.run_point_ccdm(snr, cfg, cstll, ccdm, amp_label, sch, ...
                              maxFrames, targetCwErr, maxLDPCIter, lut);
    blerAll(p, :) = out.bler;
    thr(p)        = out.thr;
    succRate(p)   = out.succRate;
    avgTx(p)      = out.avgTx;
    fprintf(['SNR=%4.1f dB | BLER1=%.2e BLER%d=%.2e | thr=%.3f info b/sym | ' ...
             'avgTx=%.2f | (%d cw, %.1fs)\n'], ...
             snr, out.bler(1), maxTx, out.bler(end), out.thr, out.avgTx, ...
             out.nCw, toc(t0));
end

% ---------------- Results table ----------------
T = table(SNR_dB.', blerAll(:,1), blerAll(:,end), thr.', avgTx.', ...
          'VariableNames', {'SNR_dB','BLER_round1','BLER_final','Throughput','avgTx'});
disp(T);

% ---------------- Plots ----------------
fig1 = figure('Name', sprintf('HARQ %s CCDM -- BLER per round', sch.mode), 'Color','w');
blerPlot = blerAll; blerPlot(blerPlot==0) = NaN;
semilogy(SNR_dB, blerPlot, '-o', 'LineWidth', 1.3); grid on;
xlabel('SNR [dB]'); ylabel('residual BLER'); ylim([1e-4 1]);
legend(compose('after round %d', 1:maxTx), 'Location','southwest');
title(sprintf('PAS 64-QAM HARQ (%s) real CCDM, \\nu=%.2g, maxTx=%d', sch.mode, nu, maxTx));

fig2 = figure('Name', sprintf('HARQ %s CCDM -- throughput', sch.mode), 'Color','w');
plot(SNR_dB, thr, '-s', 'LineWidth', 1.3, 'Color', [0 0.45 0.85]); grid on;
xlabel('SNR [dB]'); ylabel('throughput [info bits / symbol]');
title(sprintf('Normalized throughput (info), HARQ %s real CCDM', sch.mode));

% ---------------- Save ----------------
stamp   = sprintf('harq_%s_m%d_maxTx%d_ccdm', lower(char(sch.mode)), m, maxTx);
outfile = fullfile('results', stamp);
save([outfile '.mat'], 'SNR_dB','blerAll','thr','succRate','avgTx', ...
                       'sch','cfg','nu','harqMode','ccdm');
savefig(fig1, [outfile '_bler.fig']);
savefig(fig2, [outfile '_thr.fig']);
fprintf('Saved results to %s.{mat,fig}\n', outfile);