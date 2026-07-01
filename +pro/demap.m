function L = demap(y,obj,noise_power,decision)
    if nargin < 3
        noise_power = 1.0;
    end
    if nargin < 4
        decision = 'SD';
    end

    sigma = sqrt(noise_power);
    [M, m] = size(obj.label);
    n = numel(y);
    y = y(:);       % ensure column (n,1)

    llrs = zeros(n, m);

    % if sigma is a per-symbol vector, reshape to (n,1)
    if ~isscalar(sigma) && numel(sigma) > 1
        sigma = sigma(:);   % (n,1)
    end

    for level = 1:m
        idx_0 = find(obj.label(:, level) == 0);            % indices with bit=0
        a0    = reshape(obj.alphabet(idx_0), 1, []);       % (1,k0) row
        px0   = reshape(obj.px(idx_0),       1, []);       % (1,k0) row
        tmp0  = y - a0;                                     % (n,k0) broadcasting
        p0    = sum(gauss(tmp0, sigma) .* px0, 2);         % (n,1)

        idx_1 = find(obj.label(:, level) == 1);            % indices with bit=1
        a1    = reshape(obj.alphabet(idx_1), 1, []);       % (1,k1) row
        px1   = reshape(obj.px(idx_1),       1, []);       % (1,k1) row
        tmp1  = y - a1;                                     % (n,k1) broadcasting
        p1    = sum(gauss(tmp1, sigma) .* px1, 2);         % (n,1)

        llrs(:, level) = log(p0) - log(p1);
    end

    if strcmp(decision, 'SD')
        L = llrs;
    elseif strcmp(decision, 'HD')
        L = 1 - 2 * (llrs < 0);
    else
        error('decision = %s not supported, options are SD and HD', decision);
    end
end

function p = gauss(x, sigma)
    % Gaussian pdf N(0,sigma^2) without the Statistics Toolbox.
    % sigma scalar or (n,1); x is (n,k). The constant 1/(sigma*sqrt(2*pi))
    % cancels in log(p0)-log(p1), but is included so that p is a real pdf.
    p = exp(-x.^2 ./ (2 * sigma.^2)) ./ (sigma * sqrt(2*pi));
end
