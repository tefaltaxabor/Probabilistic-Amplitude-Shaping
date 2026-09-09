function comp = build_composition(n, pA)
% BUILD_COMPOSITION  Length-n amplitude-index sequence with the n-type of pA.
%
%   comp = pro.build_composition(n, pA)
%
%   Returns a (1,n) row of 0-based amplitude indices whose composition is the
%   VD-optimal n-type of pA, in sorted order; pro.draw_amplitude_bits permutes
%   it. Sharing pro.quantize_composition with pro.ccdm_init is what makes the
%   transmitted composition equal to the one the rate accounting describes.

    ni   = pro.quantize_composition(n, pA);
    comp = repelem(0:numel(pA)-1, ni);
end
