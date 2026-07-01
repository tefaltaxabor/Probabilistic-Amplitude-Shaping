function comp = build_composition(n, pA)
    ni      = round(n * pA);
    ni(end) = n - sum(ni(1:end-1));  
    comp    = repelem(0 : numel(pA)-1, ni);  
end