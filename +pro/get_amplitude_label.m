function [amp_label, amps] = get_amplitude_label(obj)
    half      = obj.M/2;
    amps      = obj.alphabet(half+1 : end);        % [1 3 5 7]
    amp_label = obj.label(half+1 : end, 2:end);    % amplitude bits (without the sign)
end