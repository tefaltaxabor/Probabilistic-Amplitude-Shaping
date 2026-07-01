%% Gabriel Cabrera
%% Shaping gain (Richardson, slide 3): achievable rate of M-ASK vs SNR.
%
% Reproduces the two-panel "Motivation: Shaping Gain" figure:
%   left  panel -> UNIFORM input over the ASK constellation,
%   right panel -> MAXWELL-BOLTZMANN input  p(a) ~ exp(-nu a^2),
% for M = 4, 8, 16, 32 ASK, against the per-dimension AWGN capacity
%   C(SNR) = 0.5*log2(1 + SNR)   [bits / channel use].
%
% The achievable rate is the constellation-constrained mutual information
% I(X;Y) of the real AWGN channel Y = X + N, N ~ N(0,1), with the points
% scaled so that E[X^2] = SNR (i.e. SNR = E[X^2]/sigma^2, sigma = 1).
% I(X;Y) is evaluated by Gauss-Hermite quadrature. For the MB curves nu is
% optimized at every SNR to maximize I (nu=0 = uniform is in the grid, so
% the MB curve is always >= the uniform one -> shaping gain).
clear;

Ms       = [4 8 16 32];                 % ASK orders (per dimension)
snr_dB   = 0:0.5:40;                    % SNR axis [dB]
snr_lin  = 10.^(snr_dB/10);
nGH      = 60;                          % Gauss-Hermite nodes
nu_grid  = [0, logspace(-3, 0, 80)];    % MB shape search (nu=0 -> uniform)
R_ref    = 0.8;                         % reference rate (fraction of log2 M) for the gain readout

[xi, w]  = gauss_hermite(nGH);          % quadrature nodes / weights (weight e^{-x^2})

nM = numel(Ms);  nS = numel(snr_lin);
I_unif = nan(nM, nS);                   % uniform mutual information
I_mb   = nan(nM, nS);                   % MB-optimized mutual information
nu_opt = nan(nM, nS);                   % winning nu per SNR (diagnostics)

for im = 1:nM
    M = Ms(im);
    b = -(M-1):2:(M-1);                 % base ASK alphabet {±1,±3,...,±(M-1)} (row)
    fprintf('M=%2d-ASK ...\n', M);

    for is = 1:nS
        P = snr_lin(is);                % target power (sigma^2 = 1)

        % ---- uniform input ----
        pu  = ones(1,M)/M;
        cu  = sqrt(P / sum(pu .* b.^2)); % scale to E[X^2] = P
        I_unif(im,is) = ask_mi(cu*b, pu, 1, xi, w);

        % ---- MB input: maximize I over nu ----
        bestI = -inf;  bestNu = 0;
        for nu = nu_grid
            pm = exp(-nu * b.^2);  pm = pm / sum(pm);
            cm = sqrt(P / sum(pm .* b.^2));
            Im = ask_mi(cm*b, pm, 1, xi, w);
            if Im > bestI, bestI = Im; bestNu = nu; end
        end
        I_mb(im,is)   = bestI;
        nu_opt(im,is) = bestNu;
    end
end

Cawgn = 0.5*log2(1 + snr_lin);          % per-dimension AWGN capacity

%% ---- Shaping gain @ reference rate (horizontal SNR gap, in dB) ----------------
gain_dB = nan(1,nM);
for im = 1:nM
    R0 = R_ref * log2(Ms(im));
    su = interp_snr(I_unif(im,:), snr_dB, R0);
    sm = interp_snr(I_mb(im,:),   snr_dB, R0);
    gain_dB(im) = su - sm;
end
fprintf('\n Shaping gain @ R = %.2f*log2(M)  (horizontal SNR gap):\n', R_ref);
for im = 1:nM
    fprintf('   %2d-ASK : R0=%.2f bits -> %+.3f dB\n', Ms(im), R_ref*log2(Ms(im)), gain_dB(im));
end

%% ---- Plot: two panels like the slide ----------------------------------------
reds  = [0.96 0.62 0.62; 0.90 0.36 0.36; 0.78 0.15 0.15; 0.50 0.00 0.00];
blues = [0.60 0.80 0.96; 0.32 0.56 0.86; 0.10 0.36 0.70; 0.00 0.13 0.45];

fig = figure('Name','Shaping gain: uniform vs MB (Richardson slide 3)','Color','w', ...
             'Position',[100 100 1100 460]);

% --- left: uniform ---
ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
plot(ax1, snr_dB, Cawgn, 'k-', 'LineWidth', 1.8, 'DisplayName','AWGN capacity');
for im = 1:nM
    plot(ax1, snr_dB, I_unif(im,:), '-', 'Color', reds(im,:), 'LineWidth', 1.8, ...
         'DisplayName', sprintf('Uniform %d-ASK', Ms(im)));
end
xlabel(ax1,'SNR (dB)'); ylabel(ax1,'bits / channel use');
xlim(ax1,[0 40]); ylim(ax1,[0 6]); title(ax1,'Uniform input');
legend(ax1,'Location','southeast');

% --- right: MB ---
ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
plot(ax2, snr_dB, Cawgn, 'k-', 'LineWidth', 1.8, 'DisplayName','AWGN capacity');
for im = 1:nM
    plot(ax2, snr_dB, I_mb(im,:), '-', 'Color', blues(im,:), 'LineWidth', 1.8, ...
         'DisplayName', sprintf('MB %d-ASK', Ms(im)));
end
xlabel(ax2,'SNR (dB)'); ylabel(ax2,'bits / channel use');
xlim(ax2,[0 40]); ylim(ax2,[0 6]); title(ax2,'Maxwell-Boltzmann input');
legend(ax2,'Location','southeast');

% annotate the 32-ASK shaping gain on the MB panel (matches the slide)
R0 = R_ref*log2(Ms(end));
sm = interp_snr(I_mb(end,:), snr_dB, R0);
text(ax2, sm-1, R0+0.5, sprintf('\\approx %.3f dB shaping gain\nover uniform %d-ASK', ...
     gain_dB(end), Ms(end)), 'Color',[0 0.45 0.74], 'FontSize',9, 'HorizontalAlignment','right');
plot(ax2, sm, R0, 'o', 'Color',[0 0.45 0.74], 'MarkerFaceColor',[0 0.45 0.74], 'MarkerSize',7, ...
     'HandleVisibility','off');

title('Shaping gain over the AWGN channel (M-ASK)');

%% ---- Save .fig + .mat into results/ -----------------------------------------
if ~isfolder('results'), mkdir('results'); end
savefig(fig, fullfile('results','capacity_mb_vs_uniform.fig'));
save(fullfile('results','capacity_mb_vs_uniform.mat'), ...
     'Ms','snr_dB','snr_lin','I_unif','I_mb','nu_opt','Cawgn','gain_dB','R_ref');
fprintf('\nSaved results/capacity_mb_vs_uniform.{fig,mat}\n');

%% ============================== helpers ======================================
function I = ask_mi(a, p, sigma, xi, w)
% Mutual information I(X;Y) [bits] of a 1D AWGN channel Y = X+N, N~N(0,sigma^2),
% with discrete input X over points a with probabilities p. Gauss-Hermite over N.
%   I = sum_i p_i E_N[ log2( g(N) / sum_j p_j g(N + a_i - a_j) ) ],  g = N(0,sigma^2)
% Evaluated with a numerically stable log-sum-exp on the denominator.
    a = a(:).';  p = p(:).';  xi = xi(:).';  w = w(:).';
    M = numel(a);
    z = sqrt(2)*sigma*xi;                          % (1,K) noise sample points
    I = 0;
    for i = 1:M
        shift = a(i) - a(:);                       % (M,1)
        arg   = -((z + shift).^2)/(2*sigma^2) + log(p(:));   % (M,K)
        m     = max(arg,[],1);
        logD  = m + log(sum(exp(arg - m),1));      % (1,K) = ln sum_j p_j g(...)
        integrand = (-xi.^2 - logD)/log(2);        % log2( exp(-xi^2) / D )
        I = I + p(i) * sum(w .* integrand)/sqrt(pi);
    end
end

function [xi, w] = gauss_hermite(n)
% Nodes/weights for ∫ e^{-x^2} f(x) dx ≈ sum_k w_k f(xi_k) (Golub-Welsch).
    k = 1:n-1;
    J = diag(sqrt(k/2),1) + diag(sqrt(k/2),-1);
    [V, D] = eig(J);
    [xi, idx] = sort(diag(D).');
    w = sqrt(pi) * (V(1,idx).^2);
end

function s = interp_snr(I_curve, snr_dB, R0)
% SNR [dB] at which the (monotone) rate curve I_curve reaches level R0.
    s = NaN;
    if R0 <= I_curve(1) || R0 >= I_curve(end), return; end
    j = find(I_curve >= R0, 1, 'first');
    s = snr_dB(j-1) + (R0 - I_curve(j-1)) * ...
        (snr_dB(j)-snr_dB(j-1)) / (I_curve(j)-I_curve(j-1));
end
