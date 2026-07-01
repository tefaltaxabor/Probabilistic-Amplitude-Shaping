cd('/home/gabriel-cabrera/Desktop/BA');
configs = {2,'dvbs2-1/2'; 3,'dvbs2-2/3'; 4,'dvbs2-3/4'};
for i = 1:size(configs,1)
    m = configs{i,1}; code = configs{i,2};
    cstll = pro.dig_mod_ASK(m, "gray");
    [amp_label, amps] = pro.get_amplitude_label(cstll);
    cfg = fec.pas_config(m, code);
    [pA, px, HA] = pro.build_shaping(0, cstll, amps);
    cc = cstll; cc.px = px;
    [bpre, bpost, bler, nCw] = fec.run_point(15, cfg, cc, pA, amp_label, 4, 8, 20);
    fprintf('OK m=%d %-12s N=%d K=%d n=%d HA=%.3f R2D=%.3f BERpre=%.2e BLER=%.2f nCw=%d\n', ...
            m, code, cfg.N, cfg.K, cfg.n, HA, 2*HA, bpre, bler, nCw);
end
disp('SMOKE PASS');
