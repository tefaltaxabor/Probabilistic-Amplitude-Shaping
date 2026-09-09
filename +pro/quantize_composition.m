function comp = quantize_composition(n, pA)
% QUANTIZE_COMPOSITION  Variational-distance-optimal integer n-type of pA.
%
%   comp = pro.quantize_composition(n, pA)
%
%   Returns integer counts comp = [n_1,...,n_K] with sum(comp) = n, realising
%   the n-type p_A^(n) = comp/n closest to pA in variational distance
%   (Bocherer, "Probabilistic Amplitude Shaping" 2023, Alg. 2.5.4, eqs.
%   (2.50)-(2.51)): floor n*pA, then hand the L remaining counts to the
%   amplitudes with the largest approximation error.
%
%   This is the SINGLE definition of the composition used across the project.
%   The transmitter (pro.build_composition -> pro.draw_amplitude_bits), the
%   CCDM rate accounting (pro.ccdm_init) and the HARQ rate metadata
%   (harq.make_rate_info) all call it, so the composition that goes on the
%   channel is exactly the one whose rate loss is reported.
%
%   A plain round(n*pA) with the remainder dumped on the last amplitude is NOT
%   equivalent: it disagreed with this quantiser in 4 of 10 tested (n,nu)
%   pairs, and it forces the rounding residue onto the last amplitude, which
%   in a Maxwell-Boltzmann alphabet is the least probable and highest-energy
%   one -- the worst place to put it.

    pA   = pA(:).';
    base = floor(n * pA);                  % (2.50)
    L    = round(n - sum(base));           % (2.51), integer by construction
    err  = pA - base / n;                  % approximation error per amplitude
    [~, order] = sort(err, 'descend');     % largest error first
    comp = base;
    comp(order(1:L)) = comp(order(1:L)) + 1;
    assert(sum(comp) == n, 'composition must sum to n');
end
