function bits = ccdm(amplitude_composition, amplitude_label)
    n = numel(amplitude_composition);
    comp = amplitude_composition(randperm(n));          % random permutation
    amplitude_bits = amplitude_label(comp + 1, :);      % +1: base-1 indices
    sign_bits = randi([0 1], n, 1, 'uint8');
    bits = [sign_bits, amplitude_bits];                 % (n, m)
    %MUX
    bits = pro.mux(bits');
end