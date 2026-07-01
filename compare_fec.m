%% Gabriel Cabrera
%% FEC code comparison under PAS 64-QAM: DVB-S2 vs 5G NR LDPC (rate 2/3)
%% Same chain (shaping, 64-QAM, channel, demapper, parity->sign); only the LDPC changes.
clear; rng(7);

% ---------------- Parameters ----------------
m  = 3;                               % 64-QAM (2 x 8-ASK)
nu = 0.05;                            % Maxwell-Boltzmann

% Codes and SNR range specific to each one (centered on its waterfall, step 0.1 dB).
% Narrow ranges => no frames wasted on points with BLER=1 or BLER=0.
codes  = {'dvbs2-2/3',  'nr-bg1-2/3',  'nr-bg2-2/3'};
ranges = {9.7:0.1:10.6, 10.0:0.1:11.2, 10.4:0.1:11.7};

maxFrames    = 1200;                  % max frames/point (each frame = 2 codewords)
targetCwErr  = 100;                   % stop the point at this many codeword (FECFRAME) errors
maxLDPCIter  = 50;

% ---------------- Constellation and shaping (common to all) ----------------
cstll = pro.dig_mod_ASK(m, "gray");
[amp_label, amps] = pro.get_amplitude_label(cstll);
pA = exp(-nu * amps.^2);  pA = pA / sum(pA);
px = zeros(1, cstll.M);
for k = 1:cstll.M
    px(k) = 0.5 * pA(amps == abs(cstll.alphabet(k)));
end
cstll.px = px;

% ---------------- Parallel pool ----------------
if isempty(gcp('nocreate')), parpool(6); end

% ---------------- Sweep per code ----------------
nCodes = numel(codes);
R = struct('name',{},'SNR',{},'berPre',{},'berPost',{},'bler',{},'n',{});

for c = 1:nCodes
    cfg = fec.pas_config(m, codes{c});
    snr = ranges{c};  np = numel(snr);
    fprintf('\n=== %s | N=%d K=%d Rc=%.4f | n=%d symbols/dim ===\n', ...
            cfg.name, cfg.N, cfg.K, cfg.Rc, cfg.n);

    bpre = nan(1,np); bpost = nan(1,np); bl = nan(1,np);
    parfor p = 1:np
        [a, b, d] = fec.run_point(snr(p), cfg, cstll, ...
                          pA, amp_label, maxFrames, targetCwErr, maxLDPCIter);
        bpre(p) = a;  bpost(p) = b;  bl(p) = d;
        fprintf('  %-16s SNR=%5.2f | BERpre=%.2e BERpost=%.2e BLER=%.2e\n', ...
                cfg.name, snr(p), a, b, d);
    end
    R(c) = struct('name',cfg.name, 'SNR',snr, 'berPre',bpre, ...
                  'berPost',bpost, 'bler',bl, 'n',cfg.n);
end

% ---------------- Results ----------------
fprintf('\n');
for c = 1:nCodes
    T = table(R(c).SNR.', R(c).berPre.', R(c).berPost.', R(c).bler.', ...
              'VariableNames', {'SNR_dB','BER_pre','BER_post','BLER'});
    fprintf('--- %s (n=%d) ---\n', R(c).name, R(c).n); disp(T);
end

% ---------------- Plots ----------------
co = lines(nCodes);

figure('Name','PAS 64-QAM: BLER DVB-S2 vs 5G NR','Color','w'); hold on; grid on;
for c = 1:nCodes
    yb = R(c).bler; yb(yb==0) = NaN;            % avoid log(0)
    semilogy(R(c).SNR, yb, '-o', 'Color', co(c,:), 'LineWidth', 1.4, ...
             'DisplayName', sprintf('%s, n=%d', R(c).name, R(c).n));
end
set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('BLER');
legend('Location','southwest'); title(sprintf('PAS 64-QAM, rate 2/3, \\nu=%.2g', nu));

figure('Name','PAS 64-QAM: post-FEC BER DVB-S2 vs 5G NR','Color','w'); hold on; grid on;
for c = 1:nCodes
    yb = R(c).berPost; yb(yb==0) = NaN;
    semilogy(R(c).SNR, yb, '-s', 'Color', co(c,:), 'LineWidth', 1.4, ...
             'DisplayName', sprintf('%s, n=%d', R(c).name, R(c).n));
end
set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('post-FEC BER');
legend('Location','southwest'); title(sprintf('PAS 64-QAM, rate 2/3, \\nu=%.2g', nu));
