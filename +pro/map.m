function symbols = map(bits,obj)
    % bits: (n, m) matrix
    % symbols: (n, 1) vector
    indices = bi2de(bits, 'left-msb') + 1;
    symbols = obj.alphabet(obj.label2idx(indices) + 1).';   % .' forces column (n,1)
end
