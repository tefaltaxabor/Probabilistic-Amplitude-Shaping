function amp_bits = draw_amplitude_bits(n, pA, amp_label)
% DRAW_AMPLITUDE_BITS  Shaped amplitude bits for ONE dimension/codeword.
%
%   amp_bits = pro.draw_amplitude_bits(n, pA, amp_label)
%
%   Fake CCDM: takes a fixed-size composition (n-type composition of pA)
%   and permutes it randomly -> gives the correct amplitude statistics,
%   but is not invertible. Returns only the shaped part (amplitude bits),
%   WITHOUT sign bits: in PAS the signs are set by the FEC code parity.
%
%   Inputs
%   ------
%   n         : number of symbols (= cfg.n of the FEC code)
%   pA        : (1,K) Maxwell-Boltzmann distribution over amplitudes
%   amp_label : (K, m-1) amplitude bit pattern per magnitude
%
%   Output
%   ------
%   amp_bits  : (n, m-1) amplitude bits (uint8)

    comp     = pro.build_composition(n, pA);   % n-type composition, size n
    comp     = comp(randperm(n));              % random permutation
    amp_bits = amp_label(comp + 1, :);         % (n, m-1), +1 = base-1
end

