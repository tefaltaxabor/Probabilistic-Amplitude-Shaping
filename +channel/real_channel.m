function [y,sigma2] = real_channel(x,SNR_dB,n)
    % 1 sample/symbol real baseband AWGN. With Ps = mean(x.^2) the actual
    % transmitted power, SNR = Ps/sigma2, i.e. SNR_dB equals Es/N0_dB here.
    x      = x(:);
    Ps     = mean(x.^2);
    sigma2 = Ps / 10^(SNR_dB/10);
    noise  = sqrt(sigma2) * randn(n,1);     % full sigma2 (1 real dimension)
    y = x + noise;
end