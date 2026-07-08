function [amp_hat, nIter, ok] = decode(llrs_sym, cfg, maxIter)
% DECODE  PAS LDPC decoding from the demapper's per-bit LLRs.
%
%   [amp_hat, nIter, ok] = fec.decode(llrs_sym, cfg, maxIter)
%
%   Inputs
%   ------
%   llrs_sym : (n, m) per-bit LLRs [sign, amp_1, ...], convention
%              LLR = log(P(bit=0)/P(bit=1))  (same as pro.demap and ldpcDecode)
%   cfg      : config from fec.pas_config
%   maxIter  : max BP iterations (default 50)
%
%   Output
%   ------
%   amp_hat  : (n, m-1) decoded amplitude bits (uint8)
%   nIter    : BP iterations performed
%   ok       : true if the syndrome is zero (all parity checks satisfied)

    if nargin < 3, maxIter = 50; end

    sign_llr = llrs_sym(:, 1);                 % (n,1) -> parity bits
    amp_llr  = llrs_sym(:, 2:end);             % (n, m-1) -> systematic
    info_llr = reshape(amp_llr.', [], 1);      % (K,1) same order as encode
    cw_llr   = [info_llr; sign_llr];           % (N,1) = [info; parity]

    [dec, nIter, pcheck] = ldpcDecode(cw_llr, cfg.dec, maxIter,"DecisionType","hard");

    amp_hat = uint8(reshape(dec, cfg.m-1, cfg.n).');   % (n, m-1)
    ok      = all(pcheck == 0);
end
