# Probabilistic Amplitude Shaping for 6G — with HARQ Extension

MATLAB link-level simulation framework for **probabilistic amplitude shaping (PAS)** with 5G NR LDPC codes, developed as part of my Bachelor's thesis at the Institute for Communications Engineering (LNT), Technical University of Munich.

PAS combines constant composition distribution matching (CCDM) with systematic FEC to approach the Shannon capacity of the AWGN channel, closing most of the shaping gap (up to 1.53 dB) left by uniform signaling. This framework implements the full transmit/receive chain and reproduces reference results from the literature, and is being extended to study **HARQ retransmission strategies for shaped constellations in 6G**.

## Simulation chain

```
info bits ──► CCDM (distribution matcher) ──► amplitude bits ─┐
                                                              ├──► 5G NR LDPC encoder ──► M-ASK mapper ──► AWGN
sign bits ◄── systematic parity ◄─────────────────────────────┘
              ──► soft demapper (bit-metric decoding) ──► LDPC decoder ──► CCDM⁻¹ ──► info bits
```

- **Distribution matching** (`+pro/`): constant composition distribution matching (CCDM) with Maxwell–Boltzmann target distributions, amplitude labeling, bit mux/demux
- **FEC** (`+fec/`): 5G NR LDPC construction (BG1/BG2), systematic encoding, belief-propagation decoding, per-SNR-point simulation
- **Channel** (`+channel/`): real and complex AWGN models
- **Metrics**: post-FEC BER, BLER, BMD (bit-metric decoding) rate, rate loss, shaping gain, gap to capacity

## Experiments

| Script | Study |
|--------|-------|
| `capacity_mb_vs_uniform.m` | AWGN capacity: Maxwell–Boltzmann shaped vs. uniform ASK |
| `compare_shaping.m` | Shaped vs. uniform signaling — BLER waterfalls and shaping gain |
| `compare_fec.m` | 5G NR LDPC vs. DVB-S2 LDPC under PAS (64-QAM) |
| `compare_blocklength.m` | BLER vs. blocklength (NR BG1) |
| `compare_ratematching.m` | Rate matching interaction with shaping |
| `main.m` | Full-chain single-configuration run |

Selected result plots are in [`ran-results/`](ran-results/) — including BLER/BER waterfalls for DVB-S2 vs. 5G NR codes, BG1 vs. BG2 comparisons, shaping gain vs. the shaping parameter ν, and gap-to-capacity analysis.

## Ongoing work: HARQ for shaped constellations (6G)

The current research direction extends the chain with **hybrid ARQ**, investigating how retransmission strategies interact with probabilistic shaping — a key open question for 6G link adaptation:

- **Strategy A — incremental redundancy**: retransmit additional parity bits only, keeping the original shaped systematic part
- **Strategy B — re-shaping**: apply fresh distribution matching across retransmissions

Performance is evaluated in terms of **BLER per HARQ round** and **normalized throughput**, quantifying whether the shaping gain is preserved (or degraded) under retransmissions and which strategy wins in different SNR regimes.

## Usage

MATLAB R2023a+ with the Parallel Computing Toolbox (Monte Carlo loops are `parfor`-parallelized). Run any of the experiment scripts directly; configuration (constellation order, shaping parameter, code rate, blocklength, SNR grid) is set at the top of each script. Results are saved to `results/` as `.mat`/`.fig`.

## References

The implementation follows the PAS architecture of Böcherer, Steiner & Schulte ("Bandwidth Efficient and Rate-Matched Low-Density Parity-Check Coded Modulation", IEEE Trans. Commun., 2015) and the CCDM construction of Schulte & Böcherer (IEEE Trans. Inf. Theory, 2016).

## Author

Gabriel Cabrera — Bachelor's thesis, Electrical and Computer Engineering
Institute for Communications Engineering (LNT/ICE), Technical University of Munich
