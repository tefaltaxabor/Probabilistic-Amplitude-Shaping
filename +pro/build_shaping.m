function [pA, px, HA] = build_shaping(nu, cstll, amps)
% BUILD_SHAPING  Maxwell-Boltzmann distribution over amplitudes and 2D prior.
%
%   [pA, px, HA] = pro.build_shaping(nu, cstll, amps)
%
%   nu = 0  -> uniform amplitudes (baseline without shaping), H(A) = log2(K).
%   nu > 0  -> MB pA(a) ~ exp(-nu a^2).
%
%   Outputs
%   -------
%   pA : (1,K) target distribution over positive amplitudes
%   px : (1,M) per-symbol prior for the demapper (0.5*pA, symmetrized in sign)
%   HA : entropy H(A) in bits/amplitude (net information rate per dimension
%        in this PAS: shaped amplitude + sign = parity). R per 2D symbol = 2*HA.
    
    %larger nu penalizes large amplitudes and favors small ones
    pA = exp(-nu * amps.^2);  pA = pA / sum(pA); 
    
    px = zeros(1, cstll.M);
    for k = 1:cstll.M
        px(k) = 0.5 * pA(amps == abs(cstll.alphabet(k)));
    end
    HA = -sum(pA .* log2(pA));
end
