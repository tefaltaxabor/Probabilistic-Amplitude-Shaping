function S = stats(x, cstll, pA, amps)
% STATS  Statistical analysis of the PAS chain shaping (fake CCDM).
%
%   S = src.stats(x, cstll, pA, amps)
%
%   Inputs
%   ------
%   x     : (n,1) transmitted complex QAM symbols (x = xI + j*xQ)
%   cstll : ASK constellation struct (from pro.dig_mod_ASK)
%   pA    : (1,K) target Maxwell-Boltzmann distribution over amplitudes
%   amps  : (1,K) positive amplitudes [1 3 5 7 ...]
%
%   Output
%   ------
%   S : struct with all the metrics (also prints and plots them).

    x = x(:);
    n = numel(x);
    pA = pA(:).';  amps = amps(:).';          % ensure rows

    % ---------- (a) Empirical amplitude distribution ----------
    a_all  = [abs(real(x)); abs(imag(x))];     % I and Q amplitudes together
    edges  = [amps - 1, amps(end) + 1];        % edges centered on amps (step 2)
    pA_emp = histcounts(a_all, edges) / numel(a_all);

    % ---------- (b) Divergence D(P_emp || P_MB) ----------
    idx   = pA_emp > 0;
    D_emp = sum(pA_emp(idx) .* log2(pA_emp(idx) ./ pA(idx)));

    % ---------- (c) n-type composition and its divergence (~rate loss) ----------
    ni      = round(n * pA);  ni(end) = n - sum(ni(1:end-1));
    pA_type = ni / n;
    jdx     = pA_type > 0;
    D_type  = sum(pA_type(jdx) .* log2(pA_type(jdx) ./ pA(jdx)));

    % ---------- (d) Entropies and rate loss of the real CCDM ----------
    H_A    = -sum(pA .* log2(pA));             % target entropy (bits/amplitude)
    k      = floor( log2_multinomial(n, ni) ); % bits the real CCDM would map
    R_dm   = k / n;
    R_loss = H_A - R_dm;

    % ---------- (e) Power and energy gain (shaped vs uniform) ----------
    Es      = mean(abs(x).^2);                  % mean shaped energy (2D)
    Es_unif = 2 * mean(cstll.alphabet.^2);      % uniform, same constellation
    gain_dB = 10 * log10(Es_unif / Es);

    % ---------- Printing ----------
    fprintf('\n===== Shaping statistics (n = %d QAM symbols) =====\n', n);
    fprintf(' amplitude | P_MB tgt | P empirical\n');
    fprintf(' ----------|----------|------------\n');
    for i = 1:numel(amps)
        fprintf('   %4d    |  %.4f  |  %.4f\n', amps(i), pA(i), pA_emp(i));
    end
    fprintf('\n D(P_emp || P_MB)   = %.3e bits/amplitude\n', D_emp);
    fprintf(' D(n-type || P_MB)  = %.3e bits/amplitude  (= D_emp for fake CCDM)\n', D_type);
    fprintf(' H(A)               = %.4f bits/amplitude\n', H_A);
    fprintf(' Rate CCDM k/n      = %.4f bits/amplitude\n', R_dm);
    fprintf(' Rate loss          = %.4f bits/amplitude\n', R_loss);
    fprintf(' Es shaped          = %.4f\n', Es);
    fprintf(' Es uniform         = %.4f\n', Es_unif);
    fprintf(' Energy gain        = %.2f dB\n\n', gain_dB);

    % ---------- (f) Visualization ----------
    figure('Name','Shaping stats','Color','w');
    subplot(1,2,1);
    bar(amps, [pA(:), pA_emp(:)]);
    legend('P_{MB} target','P empirical','Location','northeast');
    xlabel('amplitude |x|'); ylabel('probability');
    title('Amplitude distribution'); grid on;

    subplot(1,2,2);
    edges = -(amps(end)+1):2:(amps(end)+1);    % edges centered on the ASK levels
    histogram2(real(x), imag(x), edges, edges, ...
               'DisplayStyle','bar', 'ShowEmptyBins','on', ...
               'Normalization','probability', 'FaceColor','flat');
    colorbar; view(3);
    xlabel('I'); ylabel('Q'); zlabel('probability');
    title('P(x_I, x_Q) empirical');

    % ---------- Pack output ----------
    S = struct('pA_emp',pA_emp, 'D_emp',D_emp, 'D_type',D_type, ...
               'H_A',H_A, 'R_dm',R_dm, 'R_loss',R_loss, ...
               'Es',Es, 'Es_unif',Es_unif, 'gain_dB',gain_dB);
end

function L = log2_multinomial(n, ni)
    % log2 of the multinomial coefficient n! / prod(ni!), stable via gammaln.
    L = (gammaln(n+1) - sum(gammaln(ni+1))) / log(2);
end
