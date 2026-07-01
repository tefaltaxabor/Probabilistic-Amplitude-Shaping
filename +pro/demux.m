function parallel_bits = demux(bits,m)
    n = numel(bits)/m;
    parallel_bits = reshape(bits,m,n)';

end
% 
% par = demux(ans,m);