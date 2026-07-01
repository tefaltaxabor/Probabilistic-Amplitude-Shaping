function [bits_sym, cw] = encode(amp_bits, cfg)
% ENCODE  PAS encoding: amplitude bits (systematic) -> parity (signs).
%
%   [bits_sym, cw] = fec.encode(amp_bits, cfg)
%
%   Inputs
%   ------
%   amp_bits : (n, m-1) shaped amplitude bits (from pro.draw_amplitude_bits)
%   cfg      : config from fec.pas_config
%
%   Output
%   ------
%   bits_sym : (n, m) per-symbol bits [sign, amp_1, ..., amp_{m-1}].
%              Column 1 (sign) are the LDPC parity bits.
%   cw       : (N,1) full codeword [systematic; parity]
%
%   Systematic order: info((j-1)*(m-1)+l) = amp_bits(j,l)  (per symbol).

    info = reshape(double(amp_bits).', [], 1);   % (K,1)
    cw   = ldpcEncode(info, cfg.enc);            % (N,1) = [info; parity]

    parity   = cw(cfg.K+1:end);                  % (n,1) -> sign bits
    bits_sym = [parity, double(amp_bits)];       % (n, m): [sign, amplitude]
end
