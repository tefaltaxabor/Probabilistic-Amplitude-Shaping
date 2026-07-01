function [berPre, berPost, bler, nCw] = run_point(snr, cfg, cstll, pA, amp_label, maxFrames, targetCwErr, maxIter)
% RUN_POINT  Monte-Carlo of one SNR point for a given PAS FEC config.
%
%   [berPre, berPost, bler, nCw] = fec.run_point(snr, cfg, cstll, pA, ...
%                                        amp_label, maxFrames, targetCwErr, maxIter)
%
%   Accumulates up to targetCwErr codeword (FECFRAME) errors (or maxFrames
%   frames). Each frame = 2 codewords (I and Q). Returns pre-FEC BER, post-FEC
%   BER (over the systematic amplitude bits), BLER and the number of codewords
%   processed. BLER (block error rate) is the error probability of one decoded
%   LDPC codeword (a.k.a. FER in DVB-S2); the unit is one codeword = one real
%   dimension (I or Q), NOT the I+Q frame.

    n = cfg.n;
    bitErrPre = 0; nBitsPre = 0;
    bitErrPost = 0; nBitsPost = 0;
    cwErr = 0; nCw = 0;

    for f = 1:maxFrames
        % --- TX: shaped amplitude bits (fake CCDM) + PAS FEC (parity->sign) ---
        ampI = pro.draw_amplitude_bits(n, pA, amp_label);
        ampQ = pro.draw_amplitude_bits(n, pA, amp_label);
        bitsI = fec.encode(ampI, cfg);
        bitsQ = fec.encode(ampQ, cfg);

        % --- Mapping to 64-QAM and complex AWGN channel ---
        x = pro.map(bitsI, cstll) + 1j*pro.map(bitsQ, cstll);
        [y, sigma2] = channel.complex_channel(x, snr, n);

        % --- RX: demap per dimension ---
        llrI = pro.demap(real(y), cstll, sigma2/2, 'SD');
        llrQ = pro.demap(imag(y), cstll, sigma2/2, 'SD');

        % --- Pre-FEC BER (hard decision over all bits) ---
        hdI = uint8(llrI < 0);  hdQ = uint8(llrQ < 0);
        bitErrPre = bitErrPre + sum(hdI(:) ~= uint8(bitsI(:))) ...
                              + sum(hdQ(:) ~= uint8(bitsQ(:)));
        nBitsPre  = nBitsPre + numel(bitsI) + numel(bitsQ);

        % --- FEC decode + post-FEC BER and BLER ---
        ampI_hat = fec.decode(llrI, cfg, maxIter);
        ampQ_hat = fec.decode(llrQ, cfg, maxIter);
        eI = sum(ampI_hat(:) ~= ampI(:));
        eQ = sum(ampQ_hat(:) ~= ampQ(:));
        bitErrPost = bitErrPost + eI + eQ;
        nBitsPost  = nBitsPost + numel(ampI) + numel(ampQ);
        cwErr = cwErr + (eI > 0) + (eQ > 0);   % +1 per failed codeword (FECFRAME)
        nCw   = nCw + 2;                        % 2 codewords (I and Q) per frame

        if cwErr >= targetCwErr, break; end
    end

    berPre  = bitErrPre  / nBitsPre;
    berPost = bitErrPost / nBitsPost;
    bler    = cwErr      / nCw;
end
