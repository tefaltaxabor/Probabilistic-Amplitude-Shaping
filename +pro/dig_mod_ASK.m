function [ASK_obj] = dig_mod_ASK(bits_per_symbol,labeling,px)
    ASK_obj = struct();
    ASK_obj.M = 2^bits_per_symbol;
    ASK_obj.bits_per_symbol = bits_per_symbol;
    ASK_obj.alphabet = -(ASK_obj.M-1) : 2 : (ASK_obj.M-1);

    if nargin < 3
        ASK_obj.px = ones(1,ASK_obj.M);
        ASK_obj.px = ASK_obj.px/sum(ASK_obj.px);
    else
        ASK_obj.px = px;
    end

    ASK_obj.label = get_label(bits_per_symbol,labeling);
    ASK_obj.label2idx = zeros(1,ASK_obj.M,'int32');
    %careful with MSB and LSB
    %'left-msb or right-msb'
    ASK_obj.label2idx(bi2de(ASK_obj.label,'left-msb') + 1) = 0:ASK_obj.M-1;


end

function label = get_label(m,labeling)
    if strcmp(labeling,"gray")
        if m == 1
            label = uint8([0;1]);
        else
            label_n = get_label(m-1,labeling);
            tmp = 2^(m-1);
            first_half = [zeros(tmp,1,'uint8'),label_n];
            second_half = [ones(tmp,1,'uint8'),flipud(label_n)];
            label = [first_half;second_half]; 
        end 
    elseif strcmp(labeling,"natural")
        d = (0:2^m - 1)';
        power = 2.^(0:m-1);
        label = floor(mod(d,2*power)./power);
        label = uint8(fliplr(label));
    else
        error('whichlabel = %s not supported, options are gray and natural', labeling);
    end

end

