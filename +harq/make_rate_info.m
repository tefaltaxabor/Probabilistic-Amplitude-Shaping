function rateInfo = make_rate_info(pA, amps, nDM)
% MAKE_RATE_INFO  Real CCDM rate of a chosen DM block length nDM, WITHOUT
%   running any CCDM. Used to report shaping-accurate throughput while the
%   physical chain uses the fake amplitude draw.
%
%   rateInfo = harq.make_rate_info(pA, amps, nDM)
%
%   Computes, for a constant-composition DM of block length nDM targeting the
%   distribution pA over alphabet amps:
%     .nDM   : block length
%     .comp  : quantized composition [n_1..n_M] (n-type)
%     .k     : floor(log2 |T^nDM(P)|)  info bits per block
%     .Rdm   : k/nDM                   CCDM rate [bits/amplitude]
%     .Hbar  : entropy of the quantized type [bits/amplitude]
%     .Rloss : H(A) - Rdm              rate loss against the TARGET entropy
%     .Rloss_int : Hbar - Rdm          integrality part only
%
%   Same quantization and rate formulas as pro.ccdm_init (Bocherer 2023,
%   Sec. 2.5), but returns only the rate metadata -- no encode/decode.

    pA   = pA(:).';
    amps = amps(:).';
    M    = numel(amps);
    assert(numel(pA)==M, 'pA and amps must have equal length');

    comp   = pro.quantize_composition(nDM, pA);   % shared VD-optimal quantiser
    pQuant = comp / nDM;

    logNperm = (gammaln(nDM+1) - sum(gammaln(comp+1))) / log(2);
    k = floor(logNperm);

    p = pQuant(pQuant>0);
    Hbar = -sum(p .* log2(p));
    q = pA(pA>0);
    HA = -sum(q .* log2(q));

    rateInfo.nDM   = nDM;
    rateInfo.comp  = comp;
    rateInfo.k     = k;
    rateInfo.Rdm   = k / nDM;
    rateInfo.Hbar  = Hbar;
    rateInfo.HA    = HA;
    rateInfo.Rloss     = HA   - k/nDM;   % write-up definition
    rateInfo.Rloss_int = Hbar - k/nDM;   % integrality part only
end