%% Gabriel Cabrera
%% Blocklength study under PAS 64-QAM, rate 2/3.
%% Separates two effects that the DVB-S2 vs NR comparison mixed together:
%%   (1) block length  -> NR BG1 family at several Zc (same code design)
%%   (2) code design    -> NR BG1 vs NR BG2 at EQUAL n (1760)
clear; rng(7);

m  = 3;  nu = 0.05;
maxFrames    = 1200;                  % max frames/point (each frame = 2 codewords)
targetCwErr  = 100;                   % stop the point at this many codeword (FECFRAME) errors
maxLDPCIter  = 50;

% ---------------- Constellation and shaping (common) ----------------
cstll = pro.dig_mod_ASK(m, "gray");
[amp_label, amps] = pro.get_amplitude_label(cstll);
pA = exp(-nu * amps.^2);  pA = pA / sum(pA);
px = zeros(1, cstll.M);
for k = 1:cstll.M
    px(k) = 0.5 * pA(amps == abs(cstll.alphabet(k)));
end
cstll.px = px;

if isempty(gcp('nocreate')), parpool(6); end
P = {cstll, pA, amp_label, maxFrames, targetCwErr, maxLDPCIter};

%% ---------- Experiment 1: blocklength family (NR BG1) + DVB-S2 ----------
% job = {code, Zc, SNR_range}
fam = { {'nr-bg1-2/3',  96, 10.5:0.15:13.2}, ...   % n = 1056
        {'nr-bg1-2/3', 192, 10.2:0.15:12.6}, ...   % n = 2112
        {'nr-bg1-2/3', 384, 10.0:0.15:11.9}, ...   % n = 4224
        {'dvbs2-2/3',  NaN, 9.6:0.15:10.65} };     % n = 21600 (reference)
Rfam = sweep_jobs(fam, m, P);

figure('Name','PAS 64-QAM: BLER vs blocklength (NR BG1)','Color','w'); hold on; grid on;
co = parula(numel(Rfam));
for c = 1:numel(Rfam)
    yb = Rfam(c).bler; yb(yb==0) = NaN;
    sty = '-o';  if startsWith(Rfam(c).name,'DVB'), sty = '-^k'; end
    semilogy(Rfam(c).SNR, yb, sty, 'Color', co(c,:), 'LineWidth', 1.5, ...
             'MarkerFaceColor', co(c,:), 'DisplayName', Rfam(c).label);
end
if startsWith(Rfam(end).name,'DVB')   % DVB-S2 reference in black
    h = findobj(gca,'DisplayName',Rfam(end).label); set(h,'Color','k');
end
set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('BLER');
legend('Location','southwest');
title(sprintf('Blocklength effect (NR BG1 vs DVB-S2), rate 2/3, \\nu=%.2g', nu));

%% ---------- Experiment 2: EQUAL n (1760), BG1 vs BG2 ----------
% Same n and rate => isolates the base graph design.
mn = { {'nr-bg1-2/3', 160, 10.4:0.15:12.9}, ...    % n = 11*160 = 1760
       {'nr-bg2-2/3', 352, 10.4:0.15:12.9} };      % n =  5*352 = 1760
Rmn = sweep_jobs(mn, m, P);

figure('Name','PAS 64-QAM: BG1 vs BG2 at equal n (1760)','Color','w'); hold on; grid on;
com = lines(numel(Rmn));
for c = 1:numel(Rmn)
    yb = Rmn(c).bler; yb(yb==0) = NaN;
    semilogy(Rmn(c).SNR, yb, '-o', 'Color', com(c,:), 'LineWidth', 1.5, ...
             'DisplayName', Rmn(c).label);
end
set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('BLER');
legend('Location','southwest');
title(sprintf('NR BG1 vs BG2 at n=1760 (same blocklength and rate), \\nu=%.2g', nu));

% ---------------- Tables ----------------
fprintf('\n');
for c = 1:numel(Rfam)
    fprintf('--- %s ---\n', Rfam(c).label);
    disp(table(Rfam(c).SNR.', Rfam(c).bler.', 'VariableNames', {'SNR_dB','BLER'}));
end
for c = 1:numel(Rmn)
    fprintf('--- %s ---\n', Rmn(c).label);
    disp(table(Rmn(c).SNR.', Rmn(c).bler.', 'VariableNames', {'SNR_dB','BLER'}));
end

%% ---------------- helper ----------------
function R = sweep_jobs(jobs, m, P)
    [cstll, pA, amp_label, maxFrames, targetCwErr, maxLDPCIter] = P{:};
    R = struct('name',{},'label',{},'SNR',{},'berPost',{},'bler',{},'n',{});
    for c = 1:numel(jobs)
        code = jobs{c}{1};  Zc = jobs{c}{2};  snr = jobs{c}{3};
        cfg = fec.pas_config(m, code, Zc);
        np  = numel(snr);
        fprintf('\n=== %s | Zc=%g N=%d K=%d | n=%d ===\n', cfg.name, Zc, cfg.N, cfg.K, cfg.n);
        bpost = nan(1,np); bl = nan(1,np);
        parfor p = 1:np
            [~, b, d] = fec.run_point(snr(p), cfg, cstll, ...
                              pA, amp_label, maxFrames, targetCwErr, maxLDPCIter);
            bpost(p) = b; bl(p) = d;
            fprintf('  %-16s SNR=%5.2f | BERpost=%.2e BLER=%.2e\n', cfg.name, snr(p), b, d);
        end
        R(c) = struct('name',cfg.name, 'label',sprintf('%s, n=%d', cfg.name, cfg.n), ...
                      'SNR',snr, 'berPost',bpost, 'bler',bl, 'n',cfg.n);
    end
end
