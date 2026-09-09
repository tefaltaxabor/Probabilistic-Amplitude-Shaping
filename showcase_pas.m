function showcase_pas(parts, opts)
% SHOWCASE_PAS  Figure set for the PAS chapter (Sec. "PAS System").
%
%   showcase_pas                 % run every panel with default settings
%   showcase_pas([1 2])          % run only the cheap analytical panels
%   showcase_pas(1:6, opts)      % override the Monte-Carlo budget
%
%   Produces, into results/, the six panels that document the PAS chain on
%   its own (before HARQ is introduced). Panels 1-2 are pure computation and
%   run in seconds; panels 3-6 are Monte Carlo and are meant for the cluster.
%
%     1  ccdm_rateloss    R_dm(n) and R_loss(n) of a real finite-length CCDM,
%                         decomposed into quantisation + integrality, against
%                         D(p_A^(n) || p_A).   [Sec. "Accounting for the CCDM
%                         rate loss", eqs. (22)-(24)]
%     2  ccdm_validation  What the transmitter actually emits vs the target
%                         p_A, and vs the i.i.d. sampler the text describes.
%     3  pas_waterfall    BER pre/post-FEC and BLER vs SNR at the operating
%                         point.
%     4  pas_nu_sweep     BLER vs SNR for a family of nu, thresholds at a
%                         target BLER, and the gap to capacity.
%     5  pas_throughput   Throughput vs SNR against C = log2(1+SNR), with the
%                         CCDM rate loss applied. Same axes as the HARQ
%                         chapter, so the two can be read together.
%     6  pas_bmd_rate     Bit-metric decoding rate vs SNR, shaped vs uniform,
%                         with the coded operating points marked on top.
%
%   Panel 5 consumes the data of panel 4 (from memory in the same run, or
%   from results/pas_nu_sweep.mat), and panel 6 marks panel 4's thresholds if
%   they are available.
%
%   All figures are created invisible so the script is safe under
%   `matlab -batch "showcase_pas"` on a headless node.

    if nargin < 1 || isempty(parts), parts = 1:6; end
    if nargin < 2, opts = struct(); end

    % ---------------- Defaults ----------------
    d.m           = 3;                 % 8-ASK per dimension -> 64-QAM
    d.code        = 'dvbs2-2/3';       % PAS needs Rc = (m-1)/m
    d.nu0         = 0.05;              % operating shaping parameter
    d.nus         = [0 0.025 0.05 0.075 0.10];
    d.nDM         = 540;               % CCDM block length; MUST divide cfg.n
    d.blerTarget  = 1e-2;
    d.maxFrames   = 2000;
    d.targetCwErr = 80;
    d.maxLDPCIter = 50;
    d.nSymBMD     = 4e5;               % symbols per SNR point for panel 6
    d.outdir      = 'results';
    d.seed        = 7;
    opts = merge_opts(d, opts);

    rng(opts.seed);
    if ~exist(opts.outdir, 'dir'), mkdir(opts.outdir); end

    % ---------------- Common constellation / FEC ----------------
    cstll = pro.dig_mod_ASK(opts.m, "gray");
    [amp_label, amps] = pro.get_amplitude_label(cstll);
    cfg = fec.pas_config(opts.m, opts.code);

    fprintf('=== showcase_pas ===\n');
    fprintf('FEC %s | N=%d K=%d Rc=%.4f | n=%d symbols/dim\n', ...
            cfg.name, cfg.N, cfg.K, cfg.Rc, cfg.n);

    if mod(cfg.n, opts.nDM) ~= 0
        warning('showcase_pas:nDM', ...
            ['nDM=%d does not divide n=%d (%.2f blocks per codeword): the ' ...
             'rate-loss accounting would not correspond to an integer number ' ...
             'of matcher blocks. Divisors near it: %s'], ...
             opts.nDM, cfg.n, cfg.n/opts.nDM, mat2str(nearby_divisors(cfg.n, opts.nDM)));
    else
        fprintf('CCDM block nDM=%d -> %d blocks per codeword\n', ...
                opts.nDM, cfg.n/opts.nDM);
    end

    ctx = struct('cstll',cstll, 'amps',amps, 'amp_label',amp_label, ...
                 'cfg',cfg, 'opts',opts);
    S = struct();   % results carried between panels

    if any(parts == 1), panel1_rateloss(ctx);          end
    if any(parts == 2), panel2_validation(ctx);        end
    if any(parts == 3), panel3_waterfall(ctx);         end
    if any(parts == 4), S.nu = panel4_nu_sweep(ctx);   end
    if any(parts == 5), panel5_throughput(ctx, S);     end
    if any(parts == 6), panel6_bmd(ctx, S);            end

    fprintf('=== done ===\n');
end

%% ------------------------------------------------------------------------
%  Panel 1: finite-length CCDM rate loss
%  ------------------------------------------------------------------------
function panel1_rateloss(ctx)
    fprintf('\n[1] CCDM rate loss vs block length\n');
    o = ctx.opts;

    nGrid = unique(round(logspace(log10(24), log10(ctx.cfg.n), 60)));
    nNu   = numel(o.nus);   nN = numel(nGrid);

    Rdm    = nan(nNu, nN);   % k/n, what a real CCDM delivers
    Ltot   = nan(nNu, nN);   % H(A) - k/n            (eq. 23, total loss)
    Lquant = nan(nNu, nN);   % H(A) - H(p^(n))       (quantisation of the type)
    Lint   = nan(nNu, nN);   % H(p^(n)) - k/n        (integrality of log2|T|)
    Dtype  = nan(nNu, nN);   % D(p^(n) || p_A)       (eq. 24)
    HA     = nan(1, nNu);

    for i = 1:nNu
        [pA, ~, HA(i)] = pro.build_shaping(o.nus(i), ctx.cstll, ctx.amps);
        for j = 1:nN
            cc = pro.ccdm_init(pA, ctx.amps, nGrid(j));
            Rdm(i,j)    = cc.Rccdm;
            Ltot(i,j)   = HA(i) - cc.Rccdm;
            Lquant(i,j) = HA(i) - cc.Hbar;
            Lint(i,j)   = cc.Rloss;              % = Hbar - Rccdm
            q           = cc.pQuant;  k = q > 0;
            Dtype(i,j)  = sum(q(k) .* log2(q(k) ./ pA(k)));
        end
        fprintf('  nu=%.3f  H(A)=%.4f | R_dm(n=%d)=%.4f  R_loss=%.4f\n', ...
                o.nus(i), HA(i), o.nDM, ...
                interp1(nGrid, Rdm(i,:), o.nDM), interp1(nGrid, Ltot(i,:), o.nDM));
    end

    fig = figure('Name','CCDM rate loss','Color','w','Visible','off', ...
                 'Position',[100 100 1250 400]);
    tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');
    co = lines(nNu);

    % (a) matcher rate approaching H(A)
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); set(ax,'XScale','log');
    for i = 1:nNu
        semilogx(ax, nGrid, Rdm(i,:), '-', 'Color', co(i,:), 'LineWidth', 1.4, ...
                 'DisplayName', sprintf('\\nu=%.3g', o.nus(i)));
        yline(ax, HA(i), '--', 'Color', co(i,:), 'HandleVisibility','off');
    end
    xline(ax, o.nDM, ':k', 'HandleVisibility','off');
    xlabel(ax,'CCDM block length n'); ylabel(ax,'R_{dm} = k/n  [bits/amplitude]');
    title(ax,'(a) matcher rate vs H(A) (dashed)');
    legend(ax,'Location','southeast');

    % (b) total rate loss
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on');
    set(ax,'XScale','log','YScale','log');
    for i = 1:nNu
        loglog(ax, nGrid, max(Ltot(i,:), eps), '-', 'Color', co(i,:), ...
               'LineWidth', 1.4, 'DisplayName', sprintf('\\nu=%.3g', o.nus(i)));
    end
    xline(ax, o.nDM, ':k', 'HandleVisibility','off');
    xlabel(ax,'CCDM block length n'); ylabel(ax,'R_{loss} = H(A) - R_{dm}');
    title(ax,'(b) finite-length rate loss');
    legend(ax,'Location','southwest');

    % (c) decomposition at the operating nu
    [~, i0] = min(abs(o.nus - o.nu0));
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on');
    set(ax,'XScale','log','YScale','log');
    loglog(ax, nGrid, max(Ltot(i0,:),eps),   '-k',  'LineWidth',1.6, 'DisplayName','total  H(A)-k/n');
    loglog(ax, nGrid, max(Lquant(i0,:),eps), '-o',  'LineWidth',1.1, 'MarkerSize',3, 'DisplayName','quantisation  H(A)-H(p^{(n)})');
    loglog(ax, nGrid, max(Lint(i0,:),eps),   '-s',  'LineWidth',1.1, 'MarkerSize',3, 'DisplayName','integrality  H(p^{(n)})-k/n');
    loglog(ax, nGrid, max(Dtype(i0,:),eps),  '--',  'LineWidth',1.4, 'DisplayName','D(p^{(n)} || p_A)   (eq. 24)');
    xline(ax, o.nDM, ':k', 'HandleVisibility','off');
    xlabel(ax,'CCDM block length n'); ylabel(ax,'bits/amplitude');
    title(ax, sprintf('(c) decomposition, \\nu=%.3g', o.nus(i0)));
    legend(ax,'Location','southwest');

    title(tl, sprintf('Finite-length CCDM: rate and rate loss (%d-ASK)', 2^ctx.opts.m));

    data = struct('nGrid',nGrid, 'nus',o.nus, 'HA',HA, 'Rdm',Rdm, ...
                  'Ltot',Ltot, 'Lquant',Lquant, 'Lint',Lint, 'Dtype',Dtype, ...
                  'nDM',o.nDM);
    src.save_result(fig, 'pas_ccdm_rateloss', o.outdir, data);
    close(fig);
end

%% ------------------------------------------------------------------------
%  Panel 2: what the transmitter really emits
%  ------------------------------------------------------------------------
function panel2_validation(ctx)
    fprintf('\n[2] Fake-CCDM validation\n');
    o = ctx.opts;  n = ctx.cfg.n;

    [pA, px, HA] = pro.build_shaping(o.nu0, ctx.cstll, ctx.amps);
    cs = ctx.cstll;  cs.px = px;
    cs.alphabet = cs.alphabet / sqrt(sum(px(:) .* (cs.alphabet(:).^2)));

    % (i) what pro.draw_amplitude_bits actually produces: a CONSTANT
    %     COMPOSITION (n-type), randomly permuted -> exact by construction.
    comp   = pro.build_composition(n, pA);
    pCC    = histcounts(comp, -0.5:1:(numel(pA)-0.5)) / n;

    % (ii) the VD-optimal type used by pro.ccdm_init (Bocherer Alg. 2.5.4)
    ccdm   = pro.ccdm_init(pA, ctx.amps, o.nDM);
    pVD    = ccdm.pQuant;

    % (iii) the i.i.d. sampler the LaTeX text describes (eq. 18)
    nIID   = 50 * n;
    edges  = [0, cumsum(pA(:).')];  edges(end) = 1;
    [~,~,kk] = histcounts(rand(nIID,1), edges);
    kk     = min(max(kk,1), numel(pA));
    pIID   = histcounts(kk, 0.5:1:(numel(pA)+0.5)) / nIID;

    dv = @(p) sum(p(p>0) .* log2(p(p>0) ./ pA(p>0)));
    fprintf('  H(A) = %.4f bits/amplitude\n', HA);
    fprintf('  D(p_CC^(n=%d)  || p_A) = %.3e   (constant composition, what the code does)\n', n, dv(pCC));
    fprintf('  D(p_VD^(n=%d)  || p_A) = %.3e   (VD-optimal type, what ccdm_init uses)\n', o.nDM, dv(pVD));
    fprintf('  D(p_iid        || p_A) = %.3e   (i.i.d. sampler, what the text claims)\n', dv(pIID));

    % --- a shaped frame, for the 2D picture and the energy gain ---
    xI = pro.map(fec.encode(pro.draw_amplitude_bits(n, pA, ctx.amp_label), ctx.cfg), cs);
    xQ = pro.map(fec.encode(pro.draw_amplitude_bits(n, pA, ctx.amp_label), ctx.cfg), cs);
    x  = xI + 1j*xQ;

    % Energy gain must be read at FIXED MINIMUM DISTANCE, i.e. on the raw ASK
    % levels {+-1,+-3,...}: comparing after normalising both constellations to
    % E[X^2]=1 is vacuous (it returns 0 dB by construction).
    A0      = ctx.cstll.alphabet(:).';
    Es      = 2 * sum(px(:).' .* A0.^2);      % shaped,  per 2D symbol
    Es_unif = 2 * mean(A0.^2);                % uniform, same levels
    gain_dB = 10*log10(Es_unif / Es);
    pSign   = mean([double(real(x) > 0); double(imag(x) > 0)]);
    fprintf('  E_s shaped=%.4f  uniform=%.4f  ->  energy reduction %.2f dB\n', ...
            Es, Es_unif, gain_dB);
    fprintf(['       (at fixed d_min, NOT the shaping gain: the rate also drops ' ...
             'from %.3f to %.3f bits/amplitude.\n' ...
             '        The shaping gain is the gap-to-capacity comparison in panel 4.)\n'], ...
             log2(numel(pA)), HA);
    fprintf('  P(sign>0) = %.4f  (LDPC parity fills the signs; should be ~0.5)\n', pSign);

    fig = figure('Name','Fake CCDM validation','Color','w','Visible','off', ...
                 'Position',[100 100 1150 420]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

    ax = nexttile(tl);
    bar(ax, ctx.amps, [pA(:), pCC(:), pIID(:)]);  grid(ax,'on');
    xlabel(ax,'amplitude a'); ylabel(ax,'probability');
    legend(ax, {'p_A  (MB target)', ...
                sprintf('constant composition, n=%d', n), ...
                'i.i.d. sampler'}, 'Location','northeast');
    title(ax, sprintf('(a) amplitude distribution, \\nu=%.3g', o.nu0));
    % The three bars are visually indistinguishable -- that IS the result;
    % the ordering only shows up in the divergences, so print them.
    subtitle(ax, sprintf(['D(\\cdot||p_A):  const.comp. %.1e   ' ...
                          'i.i.d. %.1e   VD-type(n=%d) %.1e'], ...
                          dv(pCC), dv(pIID), o.nDM, dv(pVD)));

    ax = nexttile(tl);
    ed = -(ctx.amps(end)+1):2:(ctx.amps(end)+1);
    ed = ed / sqrt(sum(px(:) .* (ctx.cstll.alphabet(:).^2)));
    histogram2(ax, real(x), imag(x), ed, ed, 'DisplayStyle','tile', ...
               'ShowEmptyBins','on', 'Normalization','probability');
    axis(ax,'square'); colorbar(ax);
    xlabel(ax,'I'); ylabel(ax,'Q');
    title(ax, sprintf('(b) P(x_I,x_Q)'));
    subtitle(ax, sprintf(['E_s %.1f vs %.1f uniform: %.2f dB less energy at fixed d_{min} ' ...
                          '(rate %.3f vs %.3f bit/amp)'], ...
                          Es, Es_unif, gain_dB, HA, log2(numel(pA))));

    title(tl, 'Shaped source: target vs transmitted');

    data = struct('pA',pA, 'pCC',pCC, 'pVD',pVD, 'pIID',pIID, 'amps',ctx.amps, ...
                  'HA',HA, 'nu',o.nu0, 'n',n, 'nDM',o.nDM, ...
                  'D_CC',dv(pCC), 'D_VD',dv(pVD), 'D_iid',dv(pIID), ...
                  'Es',Es, 'Es_unif',Es_unif, 'gain_dB',gain_dB, 'pSign',pSign);
    src.save_result(fig, 'pas_ccdm_validation', o.outdir, data);
    close(fig);
end

%% ------------------------------------------------------------------------
%  Panel 3: baseline waterfall
%  ------------------------------------------------------------------------
function panel3_waterfall(ctx)
    fprintf('\n[3] Waterfall at nu=%.3g\n', ctx.opts.nu0);
    o = ctx.opts;

    [pA, px, HA] = pro.build_shaping(o.nu0, ctx.cstll, ctx.amps);
    cs = ctx.cstll;  cs.px = px;
    cs.alphabet = cs.alphabet / sqrt(sum(px(:) .* (cs.alphabet(:).^2)));

    snr = 8.6:0.15:11.6;
    np  = numel(snr);
    bpre = nan(1,np); bpost = nan(1,np); bl = nan(1,np);

    for p = 1:np
        t0 = tic;
        [a, b, c] = fec.run_point(snr(p), ctx.cfg, cs, pA, ctx.amp_label, ...
                                  o.maxFrames, o.targetCwErr, o.maxLDPCIter);
        bpre(p) = a;  bpost(p) = b;  bl(p) = c;
        fprintf('  [%2d/%2d] SNR=%5.2f | BERpre=%.3e BERpost=%.3e BLER=%.3e (%.1fs)\n', ...
                p, np, snr(p), a, b, c, toc(t0));
    end

    fig = figure('Name','PAS waterfall','Color','w','Visible','off');
    semilogy(snr, nan_zeros(bpre), '-o', snr, nan_zeros(bpost), '-s', ...
             snr, nan_zeros(bl), '-^', 'LineWidth', 1.4);
    grid on; xlabel('SNR [dB]'); ylabel('error rate');
    legend('BER pre-FEC','BER post-FEC','BLER','Location','southwest');
    title(sprintf('PAS %d-QAM, %s, \\nu=%.3g (H(A)=%.3f, R=%.2f bit/2D)', ...
          4^ctx.opts.m, ctx.cfg.name, o.nu0, HA, 2*HA));

    data = struct('SNR_dB',snr, 'berPre',bpre, 'berPost',bpost, 'bler',bl, ...
                  'nu',o.nu0, 'HA',HA, 'code',ctx.cfg.code, 'n',ctx.cfg.n);
    src.save_result(fig, 'pas_waterfall', o.outdir, data);
    close(fig);
end

%% ------------------------------------------------------------------------
%  Panel 4: nu sweep, thresholds and gap to capacity
%  ------------------------------------------------------------------------
function S = panel4_nu_sweep(ctx)
    fprintf('\n[4] nu sweep\n');
    o = ctx.opts;

    % SNR windows per nu, wide enough to cover the full BLER transition
    % (needed by panel 5, which integrates throughput across the waterfall).
    ranges = { 12.8:0.15:14.6, 11.2:0.15:13.0, 9.6:0.15:11.4, ...
                8.0:0.15:9.8,   6.7:0.15:8.5 };
    assert(numel(ranges) == numel(o.nus), 'ranges and nus must match');

    nN = numel(o.nus);
    S  = struct('nu',{},'HA',{},'R2D',{},'Rdm',{},'R2D_real',{}, ...
                'SNR',{},'bler',{},'berPost',{});

    for i = 1:nN
        [pA, px, HA] = pro.build_shaping(o.nus(i), ctx.cstll, ctx.amps);
        cs = ctx.cstll;  cs.px = px;
        cs.alphabet = cs.alphabet / sqrt(sum(px(:) .* (cs.alphabet(:).^2)));

        cc  = pro.ccdm_init(pA, ctx.amps, o.nDM);   % real matcher at this nu
        snr = ranges{i};  np = numel(snr);
        fprintf('  nu=%.3g | H(A)=%.4f | R_dm=%.4f | R_2D: ideal %.3f real %.3f\n', ...
                o.nus(i), HA, cc.Rccdm, 2*HA, 2*cc.Rccdm);

        bpost = nan(1,np);  bl = nan(1,np);
        for p = 1:np
            t0 = tic;
            [~, b, c] = fec.run_point(snr(p), ctx.cfg, cs, pA, ctx.amp_label, ...
                                      o.maxFrames, o.targetCwErr, o.maxLDPCIter);
            bpost(p) = b;  bl(p) = c;
            fprintf('    [%2d/%2d] nu=%.3g SNR=%5.2f | BERpost=%.2e BLER=%.2e (%.1fs)\n', ...
                    p, np, o.nus(i), snr(p), b, c, toc(t0));
        end

        S(i) = struct('nu',o.nus(i), 'HA',HA, 'R2D',2*HA, 'Rdm',cc.Rccdm, ...
                      'R2D_real',2*cc.Rccdm, 'SNR',snr, 'bler',bl, 'berPost',bpost);
    end

    % --- thresholds and gap to capacity on the SNR axis ---
    th = nan(1,nN); snrMin = nan(1,nN); gap = nan(1,nN);
    for i = 1:nN
        th(i)     = interp_threshold(S(i).SNR, S(i).bler, o.blerTarget);
        snrMin(i) = 10*log10(2^S(i).R2D_real - 1);   % Shannon SNR for that rate
        gap(i)    = th(i) - snrMin(i);
    end
    gain = gap(1) - gap;                              % vs the uniform baseline

    fig = figure('Name','PAS nu sweep','Color','w','Visible','off', ...
                 'Position',[100 100 1150 430]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');
    co = lines(nN);

    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); set(ax,'YScale','log');
    for i = 1:nN
        semilogy(ax, S(i).SNR, nan_zeros(S(i).bler), '-o', 'Color', co(i,:), ...
                 'LineWidth',1.4, 'MarkerSize',4, ...
                 'DisplayName', sprintf('\\nu=%.3g (R=%.2f bit/2D)', S(i).nu, S(i).R2D_real));
    end
    yline(ax, o.blerTarget, ':k', 'HandleVisibility','off');
    xlabel(ax,'SNR [dB]'); ylabel(ax,'BLER');
    title(ax,'(a) BLER vs SNR'); legend(ax,'Location','southwest');

    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on');
    yyaxis(ax,'left');
    plot(ax, o.nus, gap, '-o', 'LineWidth',1.5);
    ylabel(ax,'gap to capacity [dB]');
    yyaxis(ax,'right');
    plot(ax, o.nus, gain, '-s', 'LineWidth',1.5);
    ylabel(ax,'shaping gain vs \nu=0 [dB]');
    xlabel(ax,'\nu');
    title(ax, sprintf('(b) gap and gain @ BLER=%g', o.blerTarget));

    title(tl, sprintf('Shaping sweep, %s, %d-QAM', ctx.cfg.name, 4^ctx.opts.m));

    data = struct('S',S, 'nus',o.nus, 'thr',th, 'snrMin',snrMin, ...
                  'gap',gap, 'gain',gain, 'blerTarget',o.blerTarget, ...
                  'nDM',o.nDM, 'code',ctx.cfg.code);
    src.save_result(fig, 'pas_nu_sweep', o.outdir, data);
    close(fig);

    for i = 1:nN
        fprintf('  nu=%.3g | R=%.3f bit/2D | SNR@BLER=%.2f dB | gap=%.2f dB | gain=%.2f dB\n', ...
                S(i).nu, S(i).R2D_real, th(i), gap(i), gain(i));
    end
end

%% ------------------------------------------------------------------------
%  Panel 5: throughput vs SNR against capacity
%  ------------------------------------------------------------------------
function panel5_throughput(ctx, Sin)
    fprintf('\n[5] Throughput vs SNR\n');
    o = ctx.opts;

    S = get_nu_data(ctx, Sin);
    if isempty(S), return; end
    nN = numel(S);

    fig = figure('Name','PAS throughput','Color','w','Visible','off', ...
                 'Position',[100 100 1150 430]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');
    co = lines(nN);

    snrGrid = 0:0.05:20;
    Cawgn   = log2(1 + 10.^(snrGrid/10));

    % (a) per-nu throughput, ideal (loss-free) vs real matcher
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on');
    plot(ax, snrGrid, Cawgn, '-k', 'LineWidth',1.8, 'DisplayName','C = log_2(1+SNR)');
    for i = 1:nN
        etaReal  = S(i).R2D_real .* (1 - S(i).bler);
        etaIdeal = S(i).R2D      .* (1 - S(i).bler);
        plot(ax, S(i).SNR, etaReal, '-o', 'Color', co(i,:), 'LineWidth',1.5, ...
             'MarkerSize',4, 'DisplayName', sprintf('\\nu=%.3g, real CCDM', S(i).nu));
        plot(ax, S(i).SNR, etaIdeal, '--', 'Color', co(i,:), 'LineWidth',1.0, ...
             'HandleVisibility','off');
    end
    xlabel(ax,'SNR [dB]'); ylabel(ax,'throughput [bits / 2D symbol]');
    xlim(ax, [min(cellfun(@(s) min(s), {S.SNR}))-0.5, ...
              max(cellfun(@(s) max(s), {S.SNR}))+0.5]);
    title(ax, sprintf('(a) per-\\nu throughput (dashed = loss-free, n_{DM}=%d)', o.nDM));
    legend(ax,'Location','northwest');

    % (b) envelope: what DM-based rate matching delivers as a continuum
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on');
    plot(ax, snrGrid, Cawgn, '-k', 'LineWidth',1.8, 'DisplayName','C = log_2(1+SNR)');
    allSNR = unique(round(cell2mat({S.SNR})*100)/100);
    env = nan(size(allSNR));
    for j = 1:numel(allSNR)
        best = 0;
        for i = 1:nN
            if allSNR(j) >= min(S(i).SNR) && allSNR(j) <= max(S(i).SNR)
                b = interp1(S(i).SNR, S(i).bler, allSNR(j), 'linear');
                best = max(best, S(i).R2D_real * (1 - b));
            end
        end
        env(j) = best;
    end
    plot(ax, allSNR, env, '-o', 'LineWidth',1.8, 'MarkerSize',4, ...
         'DisplayName','PAS envelope over \nu');
    xlabel(ax,'SNR [dB]'); ylabel(ax,'throughput [bits / 2D symbol]');
    xlim(ax, [min(allSNR)-0.5, max(allSNR)+0.5]);
    title(ax,'(b) rate-matching envelope');
    legend(ax,'Location','northwest');

    title(tl, 'PAS throughput with the finite-length CCDM rate loss applied');

    data = struct('S',S, 'snrGrid',snrGrid, 'Cawgn',Cawgn, ...
                  'envSNR',allSNR, 'env',env, 'nDM',o.nDM);
    src.save_result(fig, 'pas_throughput', o.outdir, data);
    close(fig);
end

%% ------------------------------------------------------------------------
%  Panel 6: bit-metric decoding rate
%  ------------------------------------------------------------------------
function panel6_bmd(ctx, Sin)
    fprintf('\n[6] BMD rate\n');
    o = ctx.opts;

    snr = 0:0.5:20;
    nus = o.nus;
    Rb  = nan(numel(nus), numel(snr));

    for i = 1:numel(nus)
        [~, px, ~] = pro.build_shaping(nus(i), ctx.cstll, ctx.amps);
        cs = ctx.cstll;  cs.px = px;
        Rb(i,:) = src.bmd_rate(cs, snr, o.nSymBMD, o.seed + i);
        fprintf('  nu=%.3g done\n', nus(i));
    end

    Cawgn = log2(1 + 10.^(snr/10));            % bits / 2D symbol
    Rb2D  = 2 * Rb;                            % two real dimensions

    fig = figure('Name','BMD rate','Color','w','Visible','off', ...
                 'Position',[100 100 700 480]);
    ax = axes(fig); hold(ax,'on'); grid(ax,'on');
    plot(ax, snr, Cawgn, '-k', 'LineWidth',1.8, 'DisplayName','C = log_2(1+SNR)');
    co = lines(numel(nus));
    for i = 1:numel(nus)
        nm = sprintf('\\nu=%.3g', nus(i));
        if nus(i) == 0, nm = [nm ' (uniform)']; end
        plot(ax, snr, Rb2D(i,:), '-', 'Color', co(i,:), 'LineWidth',1.5, ...
             'DisplayName', ['R_{BMD}, ' nm]);
    end

    % coded operating points from panel 4, if available
    S = get_nu_data(ctx, Sin, true);
    if ~isempty(S)
        th = nan(1,numel(S)); rr = nan(1,numel(S));
        for i = 1:numel(S)
            th(i) = interp_threshold(S(i).SNR, S(i).bler, o.blerTarget);
            rr(i) = S(i).R2D_real;
        end
        ok = isfinite(th);
        plot(ax, th(ok), rr(ok), 'kp', 'MarkerSize',12, 'MarkerFaceColor','y', ...
             'DisplayName', sprintf('coded, BLER=%g', o.blerTarget));
    end

    xlabel(ax,'SNR [dB]'); ylabel(ax,'rate [bits / 2D symbol]');
    ylim(ax, [0 2*ctx.opts.m]);
    title(ax, sprintf('Achievable BMD rate, %d-QAM (2 \\times %d-ASK)', ...
          4^ctx.opts.m, 2^ctx.opts.m));
    legend(ax,'Location','northwest');

    data = struct('snr',snr, 'nus',nus, 'Rbmd_dim',Rb, 'Rbmd_2D',Rb2D, ...
                  'Cawgn',Cawgn, 'nSym',o.nSymBMD);
    src.save_result(fig, 'pas_bmd_rate', o.outdir, data);
    close(fig);
end

%% ------------------------------------------------------------------------
%  helpers
%  ------------------------------------------------------------------------
function S = get_nu_data(ctx, Sin, quiet)
    if nargin < 3, quiet = false; end
    S = [];
    if isfield(Sin, 'nu') && ~isempty(Sin.nu)
        S = Sin.nu;  return;
    end
    f = fullfile(ctx.opts.outdir, 'pas_nu_sweep.mat');
    if exist(f, 'file')
        L = load(f, 'S');  S = L.S;
        if ~quiet, fprintf('  (loaded %s)\n', f); end
    elseif ~quiet
        warning('showcase_pas:noNuData', ...
            'Panel 4 data not found (%s); run showcase_pas(4) first.', f);
    end
end

function t = interp_threshold(snr, bler, target)
% First crossing of `target`, interpolated in log10(BLER) vs SNR.
    snr = snr(:).';  bler = bler(:).';
    ok = isfinite(bler) & bler > 0;
    snr = snr(ok);  bler = bler(ok);
    if numel(snr) < 2 || min(bler) > target || max(bler) < target
        t = NaN;  return;
    end
    lb = log10(bler);  lt = log10(target);
    k  = find(lb(1:end-1) >= lt & lb(2:end) <= lt, 1, 'first');
    if isempty(k), t = NaN; return; end
    t = snr(k) + (lt - lb(k)) * (snr(k+1) - snr(k)) / (lb(k+1) - lb(k));
end

function y = nan_zeros(x)
    y = x;  y(y == 0) = NaN;      % so semilogy does not drop the point silently
end

function d = nearby_divisors(n, target)
    all_d = divisors_of(n);
    [~, ix] = sort(abs(all_d - target));
    d = sort(all_d(ix(1:min(4, numel(ix)))));
end

function d = divisors_of(n)
    k = 1:floor(sqrt(n));
    k = k(mod(n, k) == 0);
    d = unique([k, n ./ k]);
end

function o = merge_opts(d, u)
    o = d;
    f = fieldnames(u);
    for i = 1:numel(f), o.(f{i}) = u.(f{i}); end
end
