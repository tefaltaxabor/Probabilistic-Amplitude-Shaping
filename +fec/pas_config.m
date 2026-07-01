function cfg = pas_config(m, code, Zc)
% PAS_CONFIG  FEC block configuration for PAS (parity -> sign).
%
%   cfg = fec.pas_config(m)            % default DVB-S2 LDPC 2/3
%   cfg = fec.pas_config(m, code)
%   cfg = fec.pas_config(m, code, Zc) % Zc (lifting) only for NR codes
%
%   Builds a systematic LDPC code for PAS where the (m-1) amplitude bits per
%   symbol are the systematic part and the (N-K) parity bits occupy the signs.
%   This requires rate Rc = (m-1)/m, i.e. N-K = K/(m-1).
%
%   Inputs
%   ------
%   m    : bits per dimension (8-ASK -> m = 3)
%   code : 'dvbs2-2/3' (default) | 'nr-bg1-2/3' | 'nr-bg2-2/3'
%   Zc   : NR lifting size (default 384 = maximum). Controls the blocklength:
%          n = M*Zc with M = Kb/(m-1) (BG1: M=11, BG2: M=5 for m=3).
%
%   Output
%   ------
%   cfg  : struct with enc/dec configs, dimensions (N,K,n,m,Rc,Zc) and 'name'.

    if nargin < 2, code = 'dvbs2-2/3'; end
    if nargin < 3, Zc = 384; end

    switch code
        case 'dvbs2-1/2'                       % m=2 (16-QAM), Rc=(m-1) /m=1/2
            %H    = dvbs2ldpc(1/2);
            H = dvbsLDPCPCM('1/2');
            name = 'DVB-S2 (1/2)';
            Zc   = NaN;

        case 'dvbs2-2/3'
            H    = dvbsLDPCPCM('2/3');            % (N-K) x N, sparse logical
            name = 'DVB-S2 (2/3)';
            Zc   = NaN;                        % not applicable

        case 'dvbs2-3/4'                       % m=4 (256-QAM), Rc=(m-1)/m=3/4
            H    = dvbsLDPCPCM('3/4');
            name = 'DVB-S2 (3/4)';
            Zc   = NaN;
        
        case 'dvbs2-8/9'
            H    = dvbsLDPCPCM('8/9');
            name = 'DVB-S2 (8/9)';
            Zc   = NaN;   
        case 'nr-bg1-2/3'
            Kb = 22;  M = Kb/(m-1);            % BG1 kernel
            H    = fec.nr_ldpc_matrix(1, Zc, Kb, M);
            name = '5G NR BG1 (2/3)';

        case 'nr-bg2-2/3'
            Kb = 10;  M = Kb/(m-1);            % BG2 kernel
            H    = fec.nr_ldpc_matrix(2, Zc, Kb, M);
            name = '5G NR BG2 (2/3)';

        

        otherwise
            error('unsupported FEC code: %s', code);
    end

    [nmk, N] = size(H);                    % nmk = N-K (parity rows)
    K = N - nmk;
    n = K / (m - 1);                       % symbols per codeword and dimension

    assert(mod(n,1) == 0 && nmk == n, ...
        ['PAS parity->sign requires N-K == K/(m-1). ' ...
         'With m=%d you need Rc = (m-1)/m = %.4f, but the code gives Rc = %.4f.'], ...
         m, (m-1)/m, K/N);

    cfg.code = code;
    cfg.name = name;
    cfg.H    = H;
    cfg.enc  = ldpcEncoderConfig(H);
    cfg.dec  = ldpcDecoderConfig(cfg.enc);
    cfg.N    = N;
    cfg.K    = K;
    cfg.n    = n;
    cfg.m    = m;
    cfg.Rc   = K / N;
    cfg.Zc   = Zc;
end
