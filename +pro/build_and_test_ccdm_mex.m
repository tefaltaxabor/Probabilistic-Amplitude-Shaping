% build_and_test_ccdm_mex.m
%
% Compiles the GMP-based CCDM MEX files and verifies exact invertibility for
% large n (200, 500, 1000) -- the regime where a pure-double implementation
% fails. Run this ONCE on your machine after installing libgmp-dev.
%
%   Ubuntu:  sudo apt-get install libgmp-dev
%   then in MATLAB:  build_and_test_ccdm_mex
%
% Gabriel Cabrera -- PAS + HARQ thesis.

%% 1. Compile
fprintf('Compiling CCDM MEX files (needs libgmp-dev)...\n');
try
    mex -lgmp ccdm_encode_mex.c
    mex -lgmp ccdm_decode_mex.c
    fprintf('  compilation OK\n');
catch ME
    fprintf(2, '  compilation FAILED: %s\n', ME.message);
    fprintf(2, ['  If gmp.h is missing: sudo apt-get install libgmp-dev\n' ...
                '  If mex is not set up: run "mex -setup C" first.\n']);
    return;
end

%% 2. Build a Maxwell-Boltzmann composition (same as the thesis chain)
amps = [1 3 5 7];
nu   = 0.05;
pmf  = exp(-nu*amps.^2); pmf = pmf/sum(pmf);

make_comp = @(pmf,n) local_make_comp(pmf,n);

%% 3. Invertibility test for several n
for n = [200 500 1000]
    comp = make_comp(pmf, n);
    % k = floor(log2 |T^n(P)|) via gammaln (only to SIZE the bit vector)
    logN = (gammaln(n+1) - sum(gammaln(comp+1)))/log(2);
    k    = floor(logN);

    nTest = 200; fails = 0;
    rng(3);
    for t = 1:nTest
        b  = randi([0 1], 1, k);
        a  = ccdm_encode_mex(b, comp, amps);
        % composition check
        okComp = all(arrayfun(@(i) sum(a==amps(i))==comp(i), 1:numel(amps)));
        b2 = ccdm_decode_mex(a, comp, amps, k);
        if ~isequal(b, b2) || ~okComp, fails = fails + 1; end
    end
    fprintf('n=%4d | k=%5d bits | invertibility: %d/%d exact\n', ...
            n, k, nTest-fails, nTest);
end

fprintf('\nIf all rows show N/N exact, the MEX CCDM is ready for n=200+.\n');

function comp = local_make_comp(pmf, n)
    base = floor(n*pmf);
    L    = round(n - sum(base));
    err  = pmf - base/n;
    [~, order] = sort(err, 'descend');
    comp = base;
    comp(order(1:L)) = comp(order(1:L)) + 1;
end
