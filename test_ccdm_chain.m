%% test_ccdm_chain.m
%  Incremental validation of the REAL-CCDM PAS-HARQ chain (run_point_ccdm).
%
%  Follows the same debugging philosophy we used for Approach A: isolate each
%  stage, confirm it at high SNR first, THEN run the full point. Run the cells
%  in order and stop at the first one that fails -- that localizes any bug.
%
%  Prereqs:
%    - ccdm_encode_mex / ccdm_decode_mex compiled and on the path
%      (run build_and_test_ccdm_mex first).
%    - Your +pro, +fec, +channel, +harq packages on the path.
%
%  Gabriel Cabrera -- PAS + HARQ thesis.

%% ---------- Setup (shared by all tests) ----------
clear; rng(1);
m    = 3;                         % 8-ASK
nu   = 0.05;
amps_alphabet = [1 3 5 7];        % adjust if your alphabet differs

cstll = pro.dig_mod_ASK(m, "gray");
[amp_label, amps] = pro.get_amplitude_label(cstll);
[pA, px, ~]       = pro.build_shaping(nu, cstll, amps);
cstll.px = px;

cfg = fec.pas_config(m, 'dvbs2-2/3');
n   = cfg.n;                      % 21600 symbols

% CCDM setup. n must divide cfg.n. 21600 = 200*108 = 216*100 = 108*200 ...
nDM  = 200;                       % try 200 (Bocherer). Must divide n=21600.
assert(mod(n,nDM)==0, 'nDM must divide cfg.n');
ccdm = pro.ccdm_init(pA, amps, nDM);
fprintf('CCDM: nDM=%d, k=%d bits, Rccdm=%.4f, Rloss=%.4f bits/amp\n', ...
        ccdm.n, ccdm.k, ccdm.Rccdm, ccdm.Rloss);
fprintf('Blocks per codeword: %d | info bits per codeword: %d\n', ...
        n/nDM, (n/nDM)*ccdm.k);

%% ---------- TEST 1: MEX round-trip (bits -> amps -> bits) ----------
fprintf('\n[TEST 1] CCDM MEX invertibility on one block...\n');
b   = randi([0 1], 1, ccdm.k);
aa  = pro.ccdm_encode_mex(b, ccdm.comp, amps);
bb  = pro.ccdm_decode_mex(aa, ccdm.comp, amps, ccdm.k);
okComp = all(arrayfun(@(i) sum(aa==amps(i))==ccdm.comp(i), 1:numel(amps)));
fprintf('  invertible: %d | constant-composition: %d\n', isequal(b,bb), okComp);
assert(isequal(b,bb) && okComp, 'TEST 1 FAILED: MEX not invertible');

%% ---------- TEST 2: amp<->bits converters match your chain ----------
% CRITICAL: these must agree with pro.draw_amplitude_bits / pro.map labeling.
% We check that amp -> bits -> amp is identity on the alphabet.
fprintf('\n[TEST 2] amplitude <-> amp_bits label consistency...\n');
% Build one CCDM block, convert to bits and back, compare.
ampSeq = pro.ccdm_encode_mex(b, ccdm.comp, amps);          % (1 x nDM)
% --- forward: amp -> amp_bits (replicate run_point_ccdm converter) ---
[tf, loc] = ismember(ampSeq(:), amps);
assert(all(tf));
amp_bits  = amp_label(loc, :);                          % (nDM x m-1)
% --- backward: amp_bits -> amp ---
ampBack = zeros(1, nDM);
for j = 1:nDM
    [tf2, loc2] = ismember(amp_bits(j,:), amp_label, 'rows');
    ampBack(j) = amps(loc2*tf2 + (~tf2));   % loc2 if matched
end
fprintf('  amp->bits->amp identity: %d\n', isequal(ampBack, ampSeq));
assert(isequal(ampBack, ampSeq), ...
    ['TEST 2 FAILED: amp_label convention mismatch. Adjust ' ...
     'local_amp_to_bits/local_bits_to_amp in run_point_ccdm.m']);

%% ---------- TEST 3: full TX->encode->map, then noiseless decode ----------
% Confirms fec.encode + map + demap + decode round-trips with CCDM amplitudes.
fprintf('\n[TEST 3] noiseless full chain at very high SNR (30 dB)...\n');
nBlocks  = n/nDM;
infoBits = randi([0 1], nBlocks, ccdm.k);
ampSeq   = zeros(1,n);
for bi = 1:nBlocks
    ampSeq((bi-1)*nDM+(1:nDM)) = pro.ccdm_encode_mex(infoBits(bi,:), ccdm.comp, amps);
end
[tf, loc] = ismember(ampSeq(:), amps);
amp_bits  = amp_label(loc, :);
bits      = fec.encode(amp_bits, cfg);
x         = pro.map(bits, cstll);
[y, sig2] = channel.real_channel(x, 30, n);
L         = pro.demap(y, cstll, sig2, 'SD');
amp_hat   = fec.decode(L, cfg, 50);
nerr      = sum(amp_hat(:) ~= amp_bits(:));
fprintf('  amp-bit errors @30dB full codeword: %d / %d\n', nerr, numel(amp_bits));
assert(nerr==0, 'TEST 3 FAILED: base chain does not decode even at 30 dB');

%% ---------- TEST 4: honest genie ACK recovers info bits ----------
fprintf('\n[TEST 4] info-bit recovery through inverse CCDM (30 dB)...\n');
% Reuse amp_hat from TEST 3: convert to amplitudes, ccdm_decode, compare.
ampHatSeq = zeros(1,n);
for j = 1:n
    [tf3, loc3] = ismember(amp_hat(j,:), amp_label, 'rows');
    ampHatSeq(j) = amps(loc3*tf3 + (~tf3));
end
allok = true;
for bi = 1:nBlocks
    seg = ampHatSeq((bi-1)*nDM+(1:nDM));
    ibh = pro.ccdm_decode_mex(seg, ccdm.comp, amps, ccdm.k);
    if ~isequal(ibh, infoBits(bi,:)), allok = false; break; end
end
fprintf('  all %d blocks recover info bits: %d\n', nBlocks, allok);
assert(allok, 'TEST 4 FAILED: info bits not recovered through CCDM inverse');

%% ---------- TEST 5: run_point_ccdm at a few SNRs ----------
fprintf('\n[TEST 5] run_point_ccdm short sweep...\n');
sch = harq.build_schedule(n, 4, 'IR', 0.25);
for snr = [14 16 18 30]
    out = harq.run_point_ccdm(snr, cfg, cstll, ccdm, amp_label, sch, ...
                              200, 40, 25);
    fprintf('  SNR=%4.1f | bler1=%.3f bler4=%.3f thr=%.3f avgTx=%.2f nCw=%d\n', ...
        snr, out.bler(1), out.bler(4), out.thr, out.avgTx, out.nCw);
end
fprintf('\nIf TEST 5 shows bler4=0 and sensible thr at high SNR, the real-CCDM\n');
fprintf('chain works. thr is now in INFO bits/channel use (shaping-accurate).\n');