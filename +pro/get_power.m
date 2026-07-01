function P_s = get_power(obj)
   
    P_s = sum(obj.px.*(obj.alphabet.^2));

end