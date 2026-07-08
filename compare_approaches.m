%% Gabriel Cabrera
%% compare_approaches.m  --  Approach A (symbol puncturing) vs Approach B (sign IR)
%
%  Runs BOTH HARQ realizations over the SAME PAS chain and SNR grid, then
%  overlays throughput and avgTx. This is the "base experiment" to decide
%  whether Approach B (sign-level IR, random partition) is competitive with
%  Approach A before investing in shaping-informed ordering (H1/H2).
%
%  Keep the per-round effective rate comparable between A and B:
%    - A: punctureFrac = 0.25  -> round-1 sends 75% of symbols (16200/21600)
%    - B: frac1 chosen so round-1 releases a comparable amount of PARITY.
%  NOTE: A and B are NOT identical experiments (A withholds info+parity per
%  symbol; B withholds only parity). The comparison IS the result.

clear; rng(7);

% ---------------- Link parameters ----------------
m       = 3;                % 8-ASK per dimension (64-QAM)
nu      = 0.05;             % Maxwell-Boltzmann shaping parameter
SNR_dB  = 12:0.5:17;        % operating range from the Approach-A waterfall

maxTx       = 4;
maxFrames   = 3e3;          % base experiment: modest, just to see the shape
targetCwErr = 60;
maxLDPCIter = 25;

% ---------------- Constellation and shaping ----------------
cstll = pro.dig_mod_ASK(m, "gray");
[amp_label, amps] = pro.get_amplitude_label(cstll);
[pA, px, ~]       = pro.build_shaping(nu, cstll, amps);
cstll.px = px;

% ---------------- FEC ----------------
cfg = fec.pas_config(m, 'dvbs2-2/3');
n   = cfg.n;

% ---------------- Schedules ----------------
% Approach A: symbol puncturing (existing builder)
schA = harq.build_schedule(n, maxTx, 'IR', 0.25);
% Approach B: sign-level IR, random partition (baseline order)
schB = harq.build_sign_schedule(n, maxTx, 0.5, 'strongfirst',amps);

fprintf('Comparing Approach A (symbol punct.) vs Approach B (sign IR)\n');
fprintf('  A round-1 symbols: %d/%d  (Rc1=%.3f)\n', ...
        schA.cumSyms(1), n, cfg.K/(schA.cumSyms(1)*cfg.m));
fprintf('  B round-1 signs:   %d/%d  released\n', schB.cumSigns(1), n);

% ---------------- Sweep ----------------
nPts   = numel(SNR_dB);
thrA   = nan(1,nPts);  avgTxA = nan(1,nPts);  blerA = nan(nPts,maxTx);
thrB   = nan(1,nPts);  avgTxB = nan(1,nPts);  blerB = nan(nPts,maxTx);

pool = gcp('nocreate'); if isempty(pool), parpool(6); end

parfor p = 1:nPts
    snr = SNR_dB(p);
    oA = harq.run_point(snr, cfg, cstll, pA, amp_label, schA, ...
                        maxFrames, targetCwErr, maxLDPCIter);
    oB = harq.run_point_signIR(snr, cfg, cstll, pA, amp_label, schB, ...
                        maxFrames, targetCwErr, maxLDPCIter);
    thrA(p)=oA.thr; avgTxA(p)=oA.avgTx; blerA(p,:)=oA.bler;
    thrB(p)=oB.thr; avgTxB(p)=oB.avgTx; blerB(p,:)=oB.bler;
    fprintf('SNR=%4.1f | A: thr=%.3f avgTx=%.2f | B: thr=%.3f avgTx=%.2f\n', ...
            snr, oA.thr, oA.avgTx, oB.thr, oB.avgTx);
end

% ---------------- Compare table ----------------
T = table(SNR_dB.', thrA.', thrB.', avgTxA.', avgTxB.', ...
    'VariableNames', {'SNR_dB','thr_A','thr_B','avgTx_A','avgTx_B'});
disp(T);

% ---------------- Plots ----------------
figure('Color','w','Name','Throughput: A vs B');
plot(SNR_dB, thrA, '-o', 'LineWidth',1.4); hold on;
plot(SNR_dB, thrB, '-s', 'LineWidth',1.4); grid on;
xlabel('SNR [dB]'); ylabel('throughput [info bits / channel use]');
legend('A: symbol puncturing','B: sign-level IR','Location','southeast');
title('HARQ throughput: Approach A vs B');

figure('Color','w','Name','avgTx: A vs B');
plot(SNR_dB, avgTxA, '-o', 'LineWidth',1.4); hold on;
plot(SNR_dB, avgTxB, '-s', 'LineWidth',1.4); grid on;
xlabel('SNR [dB]'); ylabel('average transmissions');
legend('A: symbol puncturing','B: sign-level IR','Location','northeast');
title('HARQ average transmissions: Approach A vs B');

% ---------------- Save ----------------
save('results/compare_A_vs_B.mat', 'SNR_dB','thrA','thrB', ...
     'avgTxA','avgTxB','blerA','blerB','schA','schB','cfg','nu');
fprintf('Saved results/compare_A_vs_B.mat\n');