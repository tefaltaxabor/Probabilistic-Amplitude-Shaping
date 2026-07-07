%% Gabriel Cabrera
%% HARQ chain for PAS -- one switch selects the retransmission strategy
%
%  Both HARQ types run through the SAME transmit/receive/decode loop. The only
%  thing that changes is the per-round transmit schedule (harq.build_schedule):
%
%     harqMode = 'IR'  -> incremental redundancy (reveal new symbols per round)
%     harqMode = 'CC'  -> Chase combining        (repeat + LLR-combine each round)
%
%  Metrics per the thesis plan: BLER per HARQ round and normalized throughput.
clear; rng(7);

% ================= HARQ SELECTOR (change this line only) =================
harqMode = 'IR';            % 'IR' (incremental redundancy) | 'CC' (Chase)
% ========================================================================

% ---------------- HARQ parameters ----------------
maxTx        = 4;           % max transmissions (HARQ rounds)
punctureFrac = 0.25;        % symbols hidden at round 1 (higher initial rate).
                            % IR reveals them over rounds; CC repeats round 1.
                            % Set 0 for classic full-codeword Chase combining.

% ---------------- Link parameters ----------------
m       = 3;                % 8-ASK per dimension (64-QAM)
nu      = 0.05;             % Maxwell-Boltzmann shaping parameter
SNR_dB  = 4:0.5:14;

maxFrames   = 20000;
targetCwErr = 40;           % residual codeword errors (after maxTx) to accumulate
maxLDPCIter = 50;

% ---------------- Constellation and shaping ----------------
cstll = pro.dig_mod_ASK(m, "gray");
[amp_label, amps] = pro.get_amplitude_label(cstll);
[pA, px, ~]       = pro.build_shaping(nu, cstll, amps);
cstll.px = px;

% ---------------- FEC (fixed-rate PAS mother code) ----------------
cfg = fec.pas_config(m, 'dvbs2-2/3');
n   = cfg.n;

% ---------------- HARQ schedule ----------------
sch = harq.build_schedule(n, maxTx, harqMode, punctureFrac);

fprintf('HARQ mode: %s | maxTx=%d | mother %s (Rc=%.3f) | n=%d symbols/dim\n', ...
        sch.mode, maxTx, cfg.name, cfg.Rc, n);
fprintf('Per-round effective code rate and spectral efficiency:\n');
for r = 1:maxTx
    Rc_r = cfg.K / (sch.cumSyms(r) * cfg.m);          % effective code rate
    se_r = cfg.m * Rc_r;                              % info bits / symbol
    fprintf('  round %d: %+4d new sym | %4d distinct | Rc=%.3f | SE=%.3f b/sym\n', ...
            r, sch.roundSyms(r), sch.cumSyms(r), Rc_r, se_r);
end

% ---------------- Sweep ----------------
nPts     = numel(SNR_dB);
blerAll  = nan(nPts, maxTx);       % residual BLER after each round
thr      = nan(1, nPts);
succRate = nan(1, nPts);
avgTx    = nan(1, nPts);

pool = gcp('nocreate'); if isempty(pool), parpool(6); end

parfor p = 1:nPts
    snr = SNR_dB(p);
    t0  = tic;
    out = harq.run_point(snr, cfg, cstll, pA, amp_label, sch, ...
                         maxFrames, targetCwErr, maxLDPCIter);
    blerAll(p, :) = out.bler;
    thr(p)        = out.thr;
    succRate(p)   = out.succRate;
    avgTx(p)      = out.avgTx;
    fprintf(['SNR=%4.1f dB | BLER1=%.2e BLER%d=%.2e | thr=%.3f b/sym | ' ...
             'avgTx=%.2f | (%d cw, %.1fs)\n'], ...
             snr, out.bler(1), maxTx, out.bler(end), out.thr, out.avgTx, ...
             out.nCw, toc(t0));
end

% ---------------- Results table ----------------
T = table(SNR_dB.', blerAll(:,1), blerAll(:,end), thr.', avgTx.', ...
          'VariableNames', {'SNR_dB', 'BLER_round1', 'BLER_final', ...
                            'Throughput', 'avgTx'});
disp(T);

% ---------------- Plots ----------------
fig1 = figure('Name', sprintf('HARQ %s -- BLER per round', sch.mode), 'Color', 'w');
blerPlot = blerAll; blerPlot(blerPlot == 0) = NaN;   % avoid log(0)
semilogy(SNR_dB, blerPlot, '-o', 'LineWidth', 1.3); grid on;
xlabel('SNR [dB]'); ylabel('residual BLER'); ylim([1e-4 1]);
legend(compose('after round %d', 1:maxTx), 'Location', 'southwest');
title(sprintf('PAS 64-QAM HARQ (%s), \\nu=%.2g, maxTx=%d', sch.mode, nu, maxTx));

fig2 = figure('Name', sprintf('HARQ %s -- throughput', sch.mode), 'Color', 'w');
plot(SNR_dB, thr, '-s', 'LineWidth', 1.3, 'Color', [0 0.45 0.85]); grid on;
xlabel('SNR [dB]'); ylabel('throughput [info bits / symbol]');
title(sprintf('Normalized throughput, HARQ %s', sch.mode));

% ---------------- Save (results/) ----------------
stamp   = sprintf('harq_%s_m%d_maxTx%d', lower(char(sch.mode)), m, maxTx);
outfile = fullfile('results', stamp);
save([outfile '.mat'], 'SNR_dB', 'blerAll', 'thr', 'succRate', 'avgTx', ...
                       'sch', 'cfg', 'nu', 'harqMode');
savefig(fig1, [outfile '_bler.fig']);
savefig(fig2, [outfile '_thr.fig']);
fprintf('Saved results to %s.{mat,fig}\n', outfile);
