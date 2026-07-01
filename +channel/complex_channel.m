function [y,sigma2] = complex_channel(x,SNR_dB,n)
    % 1 sample/symbol complex baseband AWGN. With Ps = mean(|x|^2) the actual
    % transmitted power, SNR = Ps/sigma2, i.e. SNR_dB equals Es/N0_dB here.
    x     = x(:);
    Ps    = mean(abs(x).^2);
    sigma2 = Ps / 10^(SNR_dB/10);
    noise  = sqrt(sigma2/2) * (randn(n,1) + 1j*randn(n,1));
    y = x + noise;

end

