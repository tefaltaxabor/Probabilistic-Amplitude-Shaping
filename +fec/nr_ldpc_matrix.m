function H = nr_ldpc_matrix(bgn, Zc, Kb, M)
% NR_LDPC_MATRIX  Quasi-cyclic LDPC matrix from the 5G NR base graph.
%
%   H = fec.nr_ldpc_matrix(bgn, Zc, Kb, M)
%
%   Takes the KERNEL of base graph BGn (TS 38.212): the Kb systematic
%   columns and the first M parity columns, lifted with lifting size Zc.
%   Resulting rate Rc = Kb/(Kb+M). It does NOT apply NR's puncturing of the
%   first 2*Zc systematic bits: H is used as a generic systematic code
%   (ldpcEncode/ldpcDecode), just like DVB-S2, for an apples-to-apples
%   comparison of the code design under PAS.
%
%   Inputs
%   ------
%   bgn : 1 or 2 (base graph)
%   Zc  : valid NR lifting size (see ZSets)
%   Kb  : systematic columns of the base graph (22 for BG1, 10 for BG2)
%   M   : parity columns of the kernel (Rc = Kb/(Kb+M))
%
%   Output
%   ------
%   H   : (M*Zc) x ((Kb+M)*Zc) sparse logical

    persistent BGS
    if isempty(BGS)
        f = fullfile(matlabroot,'toolbox','5g','5g', ...
                     '+nr5g','+internal','+ldpc','baseGraph.mat');
        BGS = load(f);
    end

    % NR lifting size sets (TS 38.212 Table 5.3.2-1) -> set index
    ZSets = {[2 4 8 16 32 64 128 256], [3 6 12 24 48 96 192 384], ...
             [5 10 20 40 80 160 320],  [7 14 28 56 112 224], ...
             [9 18 36 72 144 288],     [11 22 44 88 176 352], ...
             [13 26 52 104 208],       [15 30 60 120 240]};
    setIdx = find(cellfun(@(z) any(Zc==z), ZSets), 1);
    assert(~isempty(setIdx), 'Zc=%d is not a valid NR lifting size', Zc);

    V  = BGS.(sprintf('BG%dS%d', bgn, setIdx));   % shift-value matrix
    Vk = V(1:M, 1:(Kb+M));                          % rate Kb/(Kb+M) kernel

    % calcShiftValues: P = mod(V,Zc) for V>=0, -1 for empty blocks
    P = Vk;
    nz = Vk ~= -1;
    P(nz) = mod(Vk(nz), Zc);

    H = ldpcQuasiCyclicMatrix(Zc, P);
end
