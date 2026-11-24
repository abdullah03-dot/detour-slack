# DES validation -- summary

Environment: Python 3.12.15, SeQUeNCe 1.3.0, Windows-11-10.0.26200-SP0; workers 23.

## Automatic checks

| status | check | detail |
|---|---|---|
| PASS | PF: DES Mode F equals the closed form on every clean path (BBPSSW + twirl, ideal gates) | max |F_DES - F_closed| = 2.11e-15 |
| PASS | PF: same with gate noise (tolerance 1e-9: the paper's code rounds F to 1e-10 per round) | max |F_DES - F_closed| = 2.86e-10; 312 paths in total |
| PASS | PF: Mode F is deterministic given success | max sd over seeds = 3.42e-16 |
| PASS | PF: DES regime call s equals the closed-form s in every (configuration, channel, gate error) | 24 cells |
| NOTE | PF: regime calls under DEJMPS (ideal gates) | 12/12 cells change; min DEJMPS slack 8 |
| NOTE | PR: Mode R calls at R = 100000 Hz vs closed form at the same R | 12/12 agree; F_DES - F_closed between -0.0025 and -0.0008 |
| NOTE | PR: Mode R calls at R = 10000 Hz vs closed form at the same R | 7/12 agree; F_DES - F_closed between -0.0195 and -0.0051 |
| PASS | T: every DES threshold bracketed by the T2 grid (0.5x-16x of the closed form) | 48/48 bracketed |
| PASS | T: bootstrap CI available for every bracketed threshold (>= 90% of replicates cross) | 0 without CI |
| NOTE | T: T|r3|R1000|depol2|k0.002|s1: DES between optimistic and pessimistic accounting? | optimistic 0.591 s, DES 1.09 s, pessimistic 1.67 s -> inside |
| NOTE | T: T|r3|R1000|depol2|k0|s1: DES between optimistic and pessimistic accounting? | optimistic 0.269 s, DES 0.744 s, pessimistic 1.34 s -> inside |
| NOTE | T: T|r3|R100|depol2|k0.002|s1: DES between optimistic and pessimistic accounting? | optimistic 3.01 s, DES 8.00 s, pessimistic 13.8 s -> inside |
| NOTE | T: T|r3|R100|depol2|k0|s1: DES between optimistic and pessimistic accounting? | optimistic 2.69 s, DES 7.49 s, pessimistic 13.4 s -> inside |
| NOTE | T: T|r4|R1000|depol2|k0.002|s1: DES between optimistic and pessimistic accounting? | optimistic 0.472 s, DES 0.971 s, pessimistic 1.67 s -> inside |
| NOTE | T: T|r4|R1000|depol2|k0|s1: DES between optimistic and pessimistic accounting? | optimistic 0.300 s, DES 0.818 s, pessimistic 1.50 s -> inside |
| NOTE | T: T|r4|R100|depol2|k0.002|s1: DES between optimistic and pessimistic accounting? | optimistic 3.17 s, DES 7.98 s, pessimistic 15.2 s -> inside |
| NOTE | T: T|r4|R100|depol2|k0|s1: DES between optimistic and pessimistic accounting? | optimistic 3.00 s, DES 7.87 s, pessimistic 15.0 s -> inside |
| PASS | T: every DES threshold bracketed by the T2 grid (0.5x-16x of the closed form) | 12/12 bracketed |
| PASS | T: bootstrap CI available for every bracketed threshold (>= 90% of replicates cross) | 0 without CI |
| PASS | T: every DES threshold bracketed by the T2 grid (0.5x-16x of the closed form) | 8/8 bracketed |
| PASS | T: bootstrap CI available for every bracketed threshold (>= 90% of replicates cross) | 0 without CI |
| NOTE | T: T6|r4|R1000|depol2|k0.002|s1: DES between optimistic and pessimistic accounting? | optimistic 0.672 s, DES 1.51 s, pessimistic 2.64 s -> inside |
| NOTE | T: T6|r4|R100|depol2|k0.002|s1: DES between optimistic and pessimistic accounting? | optimistic 4.22 s, DES 12.8 s, pessimistic 23.9 s -> inside |
| NOTE | T: T6|r5|R1000|depol2|k0.002|s1: DES between optimistic and pessimistic accounting? | optimistic 0.722 s, DES 1.77 s, pessimistic 3.41 s -> inside |
| NOTE | T: T6|r5|R100|depol2|k0.002|s1: DES between optimistic and pessimistic accounting? | optimistic 5.56 s, DES 16.4 s, pessimistic 32.4 s -> inside |
| PASS | T: every DES threshold bracketed by the T2 grid (0.5x-16x of the closed form) | 4/4 bracketed |
| PASS | T: bootstrap CI available for every bracketed threshold (>= 90% of replicates cross) | 0 without CI |
| NOTE | T: T7|r5|R1000|depol2|k0.002|s1: DES between optimistic and pessimistic accounting? | optimistic 0.890 s, DES 2.39 s, pessimistic 4.69 s -> inside |
| NOTE | T: T7|r5|R100|depol2|k0.002|s1: DES between optimistic and pessimistic accounting? | optimistic 6.59 s, DES 22.1 s, pessimistic 44.6 s -> inside |
| PASS | S1: DES success equals qunet_tqe success on every trial and policy | 0 mismatches over 25000 (trial, policy) pairs |
| PASS | S1: DES fidelity equals qunet_tqe on every path without a swap failure | max |dF| = 2.22e-15 |
| PASS | S1: DES success equals the paper's stored E1 trials | 0 mismatches |
| PASS | S2 half|k0.002|R1000 s0: DES clean-path profile confirms s = 0 at T2 = 0.2426 s | F*(0) = 0.8182+-0.0007 (P>=.8 0.90); F*(1) = 0.7858+-0.0007 (P>=.8 0.14) |
| NOTE | S2 half|k0.002|R1000 s0: headroom oracle - BFS | +3.03 pp, exact McNemar p = 8.08e-28 |
| PASS | S2 half|k0.002|R1000 s1: DES clean-path profile confirms s = 1 at T2 = 0.4947 s | F*(0) = 0.8445+-0.0004 (P>=.8 1.00); F*(1) = 0.8147+-0.0004 (P>=.8 0.95); F*(2) = 0.7858+-0.0004 (P>=.8 0.00) |
| NOTE | S2 half|k0.002|R1000 s1: headroom oracle - BFS | +23.40 pp, exact McNemar p = 9.51e-212 |
| PASS | S2 deph2|k0.002|R1000 s0: DES clean-path profile confirms s = 0 at T2 = 0.5233 s | F*(0) = 0.8190+-0.0007 (P>=.8 0.90); F*(1) = 0.7854+-0.0007 (P>=.8 0.12) |
| NOTE | S2 deph2|k0.002|R1000 s0: headroom oracle - BFS | +2.13 pp, exact McNemar p = 1.08e-19 |
| PASS | S2 deph2|k0.002|R1000 s1: DES clean-path profile confirms s = 1 at T2 = 1.097 s | F*(0) = 0.8455+-0.0004 (P>=.8 1.00); F*(1) = 0.8154+-0.0004 (P>=.8 0.96); F*(2) = 0.7861+-0.0004 (P>=.8 0.00) |
| NOTE | S2 deph2|k0.002|R1000 s1: headroom oracle - BFS | +23.80 pp, exact McNemar p = 2.32e-215 |
| PASS | S2 depol2|k0.002|R1000 s0: DES clean-path profile confirms s = 0 at T2 = 0.8007 s | F*(0) = 0.8195+-0.0007 (P>=.8 0.91); F*(1) = 0.7860+-0.0007 (P>=.8 0.12) |
| NOTE | S2 depol2|k0.002|R1000 s0: headroom oracle - BFS | +2.67 pp, exact McNemar p = 1.65e-24 |
| PASS | S2 depol2|k0.002|R1000 s1: DES clean-path profile confirms s = 1 at T2 = 1.632 s | F*(0) = 0.8456+-0.0003 (P>=.8 1.00); F*(1) = 0.8148+-0.0004 (P>=.8 0.96); F*(2) = 0.7855+-0.0004 (P>=.8 0.00) |
| NOTE | S2 depol2|k0.002|R1000 s1: headroom oracle - BFS | +23.83 pp, exact McNemar p = 1.16e-215 |

## PF -- Mode F clean-path profiles (the paper's assumptions; deterministic given success)

Mode F is an implementation cross-check: it can disagree with the closed form only through an error in the closed-form algebra or in the simulator. Mode R (PR, T, S2) is the physics test.

| config | channel | s closed | margin delta | s DES Mode F | s DES DEJMPS | s DES CNOT, no twirl | F*(0) BBPSSW | F*(0) DEJMPS | F*(0) CNOT, no twirl |
|---|---|---|---|---|---|---|---|---|---|
| g5_r3 | deph2 | 0 | 0.0022 | 0 | >=8 | -1 | 0.8267 | 0.9340 | 0.5519 |
| g5_r3 | depol2 | 0 | 0.0031 | 0 | >=8 | -1 | 0.8259 | 0.9329 | 0.5362 |
| g5_r3 | half | 1 | 0.0039 | 1 | >=8 | -1 | 0.8313 | 0.9344 | 0.5523 |
| g5_r4 | deph2 | 2 | 0.0080 | 2 | >=8 | -1 | 0.8593 | 0.9359 | 0.5033 |
| g5_r4 | depol2 | 2 | 0.0093 | 2 | >=8 | -1 | 0.8585 | 0.9349 | 0.4857 |
| g5_r4 | half | 2 | 0.0009 | 2 | >=8 | -1 | 0.8627 | 0.9363 | 0.5037 |
| g6_r3 | deph2 | -1 | 0.0103 | -1 | >=8 | -1 | 0.7897 | 0.9179 | 0.5272 |
| g6_r3 | depol2 | -1 | 0.0116 | -1 | >=8 | -1 | 0.7884 | 0.9163 | 0.5067 |
| g6_r3 | half | -1 | 0.0034 | -1 | >=8 | -1 | 0.7966 | 0.9186 | 0.5278 |
| g6_r4 | deph2 | 1 | 0.0059 | 1 | >=8 | -1 | 0.8280 | 0.9202 | 0.4975 |
| g6_r4 | depol2 | 1 | 0.0045 | 1 | >=8 | -1 | 0.8267 | 0.9187 | 0.4755 |
| g6_r4 | half | 1 | 0.0074 | 1 | >=8 | -1 | 0.8331 | 0.9208 | 0.4981 |

'>=k': the slack reached the enumeration depth (h* + 8). CNOT without twirl = bilateral CNOT with neither twirl nor DEJMPS rotation; it detects bit flips only (an ablation of the twirl).

Gate noise (two-qubit depolarizing eps after each bilateral CNOT; half memory, T2 = 0.1 s):

| config | eps | s closed (BBPSSW) | s DES BBPSSW | s DES DEJMPS | F*(0) BBPSSW | F*(0) DEJMPS |
|---|---|---|---|---|---|---|
| g5_r3 | 0.002 | 0 | 0 | >=8 | 0.8206 | 0.9243 |
| g5_r3 | 0.005 | 0 | 0 | 7 | 0.8046 | 0.9093 |
| g5_r3 | 0.01 | -1 | -1 | 4 | 0.7785 | 0.8842 |
| g5_r4 | 0.002 | 2 | 2 | >=8 | 0.8498 | 0.9272 |
| g5_r4 | 0.005 | 1 | 1 | >=8 | 0.8308 | 0.9135 |
| g5_r4 | 0.01 | -1 | -1 | 5 | 0.7994 | 0.8904 |
| g6_r3 | 0.002 | -1 | -1 | >=8 | 0.7842 | 0.9065 |
| g6_r3 | 0.005 | -1 | -1 | 6 | 0.7659 | 0.8884 |
| g6_r3 | 0.01 | -1 | -1 | 3 | 0.7362 | 0.8586 |
| g6_r4 | 0.002 | 0 | 0 | >=8 | 0.8181 | 0.9099 |
| g6_r4 | 0.005 | -1 | -1 | 7 | 0.7959 | 0.8935 |
| g6_r4 | 0.01 | -1 | -1 | 4 | 0.7599 | 0.8659 |


## PR -- Mode R at the paper's T2, fast generation (per-memory storage resolved, kappa kept)

| config | channel | R (Hz) | s closed (same R) | s DES (mean F) | mean of F_DES - F_closed | best-class F*(k) DES +- SE (closed form, P(F>=0.8)) |
|---|---|---|---|---|---|---|
| g5_r3 | deph2 | 10000 | 0 | 0 | -0.0061 | k=0: 0.8177+-0.0002 (cf 0.8229, P>=.8 1.00); k=1: 0.7872+-0.0002 (cf 0.7942, P>=.8 0.00) |
| g5_r3 | deph2 | 100000 | 0 | 0 | -0.0009 | k=0: 0.8255+-0.0000 (cf 0.8264, P>=.8 1.00); k=1: 0.7964+-0.0000 (cf 0.7975, P>=.8 0.00) |
| g5_r3 | depol2 | 10000 | 0 | 0 | -0.0062 | k=0: 0.8167+-0.0002 (cf 0.8219, P>=.8 0.99); k=1: 0.7860+-0.0002 (cf 0.7931, P>=.8 0.00) |
| g5_r3 | depol2 | 100000 | 0 | 0 | -0.0009 | k=0: 0.8247+-0.0000 (cf 0.8255, P>=.8 1.00); k=1: 0.7954+-0.0000 (cf 0.7965, P>=.8 0.00) |
| g5_r3 | half | 10000 | 1 | 0 | -0.0067 | k=0: 0.8229+-0.0002 (cf 0.8279, P>=.8 1.00); k=1: 0.7940+-0.0002 (cf 0.8008, P>=.8 0.00); k=2: 0.7668+-0.0002 (cf 0.7752, P>=.8 0.00) |
| g5_r3 | half | 100000 | 1 | 1 | -0.0010 | k=0: 0.8302+-0.0000 (cf 0.8310, P>=.8 1.00); k=1: 0.8026+-0.0000 (cf 0.8036, P>=.8 1.00); k=2: 0.7766+-0.0000 (cf 0.7778, P>=.8 0.00) |
| g5_r4 | deph2 | 10000 | 2 | 1 | -0.0146 | k=0: 0.8406+-0.0004 (cf 0.8509, P>=.8 1.00); k=1: 0.8143+-0.0003 (cf 0.8280, P>=.8 0.96); k=2: 0.7901+-0.0003 (cf 0.8059, P>=.8 0.01); k=3: 0.7659+-0.0003 (cf 0.7847, P>=.8 0.00) |
| g5_r4 | deph2 | 100000 | 2 | 2 | -0.0019 | k=0: 0.8572+-0.0000 (cf 0.8585, P>=.8 1.00); k=1: 0.8334+-0.0000 (cf 0.8352, P>=.8 1.00); k=2: 0.8107+-0.0000 (cf 0.8128, P>=.8 1.00); k=3: 0.7889+-0.0000 (cf 0.7913, P>=.8 0.00) |
| g5_r4 | depol2 | 10000 | 2 | 1 | -0.0145 | k=0: 0.8403+-0.0003 (cf 0.8497, P>=.8 1.00); k=1: 0.8141+-0.0003 (cf 0.8266, P>=.8 0.96); k=2: 0.7876+-0.0004 (cf 0.8043, P>=.8 0.00); k=3: 0.7634+-0.0004 (cf 0.7829, P>=.8 0.00) |
| g5_r4 | depol2 | 100000 | 2 | 2 | -0.0019 | k=0: 0.8563+-0.0000 (cf 0.8576, P>=.8 1.00); k=1: 0.8324+-0.0000 (cf 0.8341, P>=.8 1.00); k=2: 0.8095+-0.0000 (cf 0.8116, P>=.8 1.00); k=3: 0.7875+-0.0000 (cf 0.7899, P>=.8 0.00) |
| g5_r4 | half | 10000 | 2 | 1 | -0.0141 | k=0: 0.8459+-0.0003 (cf 0.8549, P>=.8 1.00); k=1: 0.8201+-0.0003 (cf 0.8331, P>=.8 0.98); k=2: 0.7965+-0.0003 (cf 0.8124, P>=.8 0.34); k=3: 0.7741+-0.0003 (cf 0.7926, P>=.8 0.00) |
| g5_r4 | half | 100000 | 2 | 2 | -0.0018 | k=0: 0.8606+-0.0000 (cf 0.8619, P>=.8 1.00); k=1: 0.8381+-0.0000 (cf 0.8397, P>=.8 1.00); k=2: 0.8166+-0.0000 (cf 0.8186, P>=.8 1.00); k=3: 0.7962+-0.0000 (cf 0.7984, P>=.8 0.00) |
| g6_r3 | deph2 | 10000 | -1 | -1 | -0.0076 | k=0: 0.7796+-0.0002 (cf 0.7861, P>=.8 0.00); k=1: 0.7506+-0.0002 (cf 0.7593, P>=.8 0.00) |
| g6_r3 | deph2 | 100000 | -1 | -1 | -0.0012 | k=0: 0.7882+-0.0000 (cf 0.7893, P>=.8 0.00); k=1: 0.7610+-0.0000 (cf 0.7623, P>=.8 0.00) |
| g6_r3 | depol2 | 10000 | -1 | -1 | -0.0080 | k=0: 0.7774+-0.0002 (cf 0.7847, P>=.8 0.00); k=1: 0.7491+-0.0002 (cf 0.7578, P>=.8 0.00) |
| g6_r3 | depol2 | 100000 | -1 | -1 | -0.0012 | k=0: 0.7870+-0.0000 (cf 0.7881, P>=.8 0.00); k=1: 0.7597+-0.0000 (cf 0.7609, P>=.8 0.00) |
| g6_r3 | half | 10000 | -1 | -1 | -0.0076 | k=0: 0.7868+-0.0002 (cf 0.7935, P>=.8 0.00); k=1: 0.7600+-0.0002 (cf 0.7685, P>=.8 0.00) |
| g6_r3 | half | 100000 | -1 | -1 | -0.0011 | k=0: 0.7953+-0.0000 (cf 0.7963, P>=.8 0.00); k=1: 0.7699+-0.0000 (cf 0.7711, P>=.8 0.00) |
| g6_r4 | deph2 | 10000 | 0 | 0 | -0.0153 | k=0: 0.8066+-0.0004 (cf 0.8201, P>=.8 0.85); k=1: 0.7812+-0.0004 (cf 0.7984, P>=.8 0.00) |
| g6_r4 | deph2 | 100000 | 1 | 1 | -0.0021 | k=0: 0.8255+-0.0000 (cf 0.8272, P>=.8 1.00); k=1: 0.8031+-0.0000 (cf 0.8052, P>=.8 1.00); k=2: 0.7816+-0.0000 (cf 0.7840, P>=.8 0.00) |
| g6_r4 | depol2 | 10000 | 0 | 0 | -0.0148 | k=0: 0.8053+-0.0004 (cf 0.8184, P>=.8 0.81); k=1: 0.7801+-0.0004 (cf 0.7965, P>=.8 0.00) |
| g6_r4 | depol2 | 100000 | 1 | 1 | -0.0021 | k=0: 0.8242+-0.0000 (cf 0.8259, P>=.8 1.00); k=1: 0.8015+-0.0000 (cf 0.8037, P>=.8 0.95); k=2: 0.7798+-0.0000 (cf 0.7823, P>=.8 0.00) |
| g6_r4 | half | 10000 | 1 | 0 | -0.0152 | k=0: 0.8137+-0.0003 (cf 0.8260, P>=.8 0.96); k=1: 0.7904+-0.0003 (cf 0.8057, P>=.8 0.01); k=2: 0.7682+-0.0003 (cf 0.7863, P>=.8 0.00) |
| g6_r4 | half | 100000 | 1 | 1 | -0.0020 | k=0: 0.8307+-0.0000 (cf 0.8324, P>=.8 1.00); k=1: 0.8097+-0.0000 (cf 0.8117, P>=.8 1.00); k=2: 0.7897+-0.0000 (cf 0.7919, P>=.8 0.00) |

## T -- design chart check (5x5): smallest T2 with s >= 0 / s >= 1, Mode R

Mean rule: mean F*(k) >= F_th. Median rule: P(F*(k) >= F_th) >= 1/2. CI: 95% bootstrap over deliveries (replicates that do not cross inside the grid are counted as 'dropped').

| r | R (Hz) | channel | kappa (s/km) | target | T2 closed form, optimistic (s) | T2 pessimistic, depol2 (s) | T2 DES, mean rule (s) | 95% CI | T2 DES, median rule (s) | DES / closed | boot dropped |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 1000 | deph2 | 0.002 | s>=0 | 0.231 | -- | 0.375 | [0.367, 0.382] | 0.355 | 1.62 | 0 |
| 3 | 1000 | deph2 | 0.002 | s>=1 | 0.386 | -- | 0.731 | [0.718, 0.746] | 0.688 | 1.89 | 0 |
| 3 | 1000 | deph2 | 0 | s>=0 | 0.111 | -- | 0.252 | [0.243, 0.261] | 0.223 | 2.28 | 0 |
| 3 | 1000 | deph2 | 0 | s>=1 | 0.176 | -- | 0.492 | [0.478, 0.510] | 0.463 | 2.80 | 0 |
| 3 | 1000 | depol2 | 0.002 | s>=0 | 0.358 | -- | 0.589 | [0.576, 0.602] | 0.561 | 1.65 | 0 |
| 3 | 1000 | depol2 | 0.002 | s>=1 | 0.591 | 1.67 | 1.09 | [1.07, 1.11] | 1.03 | 1.84 | 0 |
| 3 | 1000 | depol2 | 0 | s>=0 | 0.171 | -- | 0.402 | [0.388, 0.415] | 0.365 | 2.34 | 0 |
| 3 | 1000 | depol2 | 0 | s>=1 | 0.269 | 1.34 | 0.744 | [0.727, 0.768] | 0.693 | 2.77 | 0 |
| 3 | 1000 | half | 0.002 | s>=0 | 0.102 | -- | 0.178 | [0.175, 0.182] | 0.169 | 1.75 | 0 |
| 3 | 1000 | half | 0.002 | s>=1 | 0.166 | -- | 0.330 | [0.324, 0.340] | 0.308 | 1.99 | 0 |
| 3 | 1000 | half | 0 | s>=0 | 0.0490 | -- | 0.122 | [0.119, 0.126] | 0.113 | 2.49 | 0 |
| 3 | 1000 | half | 0 | s>=1 | 0.0756 | -- | 0.243 | [0.236, 0.250] | 0.222 | 3.21 | 0 |
| 3 | 100 | deph2 | 0.002 | s>=0 | 1.23 | -- | 2.74 | [2.65, 2.81] | 2.49 | 2.23 | 0 |
| 3 | 100 | deph2 | 0.002 | s>=1 | 1.97 | -- | 5.21 | [5.09, 5.35] | 4.84 | 2.65 | 0 |
| 3 | 100 | deph2 | 0 | s>=0 | 1.11 | -- | 2.60 | [2.54, 2.67] | 2.46 | 2.35 | 0 |
| 3 | 100 | deph2 | 0 | s>=1 | 1.76 | -- | 4.95 | [4.83, 5.12] | 4.62 | 2.82 | 0 |
| 3 | 100 | depol2 | 0.002 | s>=0 | 1.90 | -- | 4.15 | [3.99, 4.30] | 3.72 | 2.18 | 0 |
| 3 | 100 | depol2 | 0.002 | s>=1 | 3.01 | 13.8 | 8.00 | [7.78, 8.25] | 7.47 | 2.66 | 0 |
| 3 | 100 | depol2 | 0 | s>=0 | 1.71 | -- | 3.96 | [3.82, 4.08] | 3.61 | 2.31 | 0 |
| 3 | 100 | depol2 | 0 | s>=1 | 2.69 | 13.4 | 7.49 | [7.33, 7.74] | 6.89 | 2.79 | 0 |
| 3 | 100 | half | 0.002 | s>=0 | 0.543 | -- | 1.27 | [1.23, 1.31] | 1.16 | 2.34 | 0 |
| 3 | 100 | half | 0.002 | s>=1 | 0.846 | -- | 2.46 | [2.38, 2.54] | 2.31 | 2.91 | 0 |
| 3 | 100 | half | 0 | s>=0 | 0.490 | -- | 1.25 | [1.22, 1.29] | 1.16 | 2.55 | 0 |
| 3 | 100 | half | 0 | s>=1 | 0.756 | -- | 2.40 | [2.31, 2.47] | 2.12 | 3.18 | 0 |
| 4 | 1000 | deph2 | 0.002 | s>=0 | 0.238 | -- | 0.436 | [0.427, 0.446] | 0.406 | 1.83 | 0 |
| 4 | 1000 | deph2 | 0.002 | s>=1 | 0.303 | -- | 0.639 | [0.622, 0.653] | 0.610 | 2.11 | 0 |
| 4 | 1000 | deph2 | 0 | s>=0 | 0.156 | -- | 0.347 | [0.336, 0.356] | 0.322 | 2.22 | 0 |
| 4 | 1000 | deph2 | 0 | s>=1 | 0.192 | -- | 0.520 | [0.510, 0.533] | 0.490 | 2.70 | 0 |
| 4 | 1000 | depol2 | 0.002 | s>=0 | 0.375 | -- | 0.670 | [0.657, 0.686] | 0.628 | 1.79 | 0 |
| 4 | 1000 | depol2 | 0.002 | s>=1 | 0.472 | 1.67 | 0.971 | [0.945, 0.997] | 0.918 | 2.06 | 0 |
| 4 | 1000 | depol2 | 0 | s>=0 | 0.246 | -- | 0.560 | [0.547, 0.573] | 0.530 | 2.27 | 0 |
| 4 | 1000 | depol2 | 0 | s>=1 | 0.300 | 1.50 | 0.818 | [0.798, 0.842] | 0.770 | 2.73 | 0 |
| 4 | 1000 | half | 0.002 | s>=0 | 0.109 | -- | 0.205 | [0.202, 0.210] | 0.193 | 1.88 | 0 |
| 4 | 1000 | half | 0.002 | s>=1 | 0.136 | -- | 0.301 | [0.294, 0.308] | 0.286 | 2.21 | 0 |
| 4 | 1000 | half | 0 | s>=0 | 0.0718 | -- | 0.168 | [0.163, 0.173] | 0.153 | 2.34 | 0 |
| 4 | 1000 | half | 0 | s>=1 | 0.0865 | -- | 0.242 | [0.237, 0.249] | 0.224 | 2.79 | 0 |
| 4 | 100 | deph2 | 0.002 | s>=0 | 1.65 | -- | 3.63 | [3.51, 3.73] | 3.38 | 2.21 | 0 |
| 4 | 100 | deph2 | 0.002 | s>=1 | 2.03 | -- | 5.22 | [5.11, 5.33] | 4.87 | 2.57 | 0 |
| 4 | 100 | deph2 | 0 | s>=0 | 1.56 | -- | 3.52 | [3.42, 3.61] | 3.27 | 2.25 | 0 |
| 4 | 100 | deph2 | 0 | s>=1 | 1.92 | -- | 5.22 | [5.08, 5.37] | 4.87 | 2.71 | 0 |
| 4 | 100 | depol2 | 0.002 | s>=0 | 2.59 | -- | 5.45 | [5.21, 5.64] | 4.94 | 2.10 | 0 |
| 4 | 100 | depol2 | 0.002 | s>=1 | 3.17 | 15.2 | 7.98 | [7.80, 8.15] | 7.43 | 2.51 | 0 |
| 4 | 100 | depol2 | 0 | s>=0 | 2.46 | -- | 5.33 | [5.14, 5.49] | 4.85 | 2.16 | 0 |
| 4 | 100 | depol2 | 0 | s>=1 | 3.00 | 15.0 | 7.87 | [7.71, 8.03] | 7.39 | 2.62 | 0 |
| 4 | 100 | half | 0.002 | s>=0 | 0.755 | -- | 1.73 | [1.69, 1.77] | 1.63 | 2.30 | 0 |
| 4 | 100 | half | 0.002 | s>=1 | 0.914 | -- | 2.56 | [2.49, 2.64] | 2.39 | 2.80 | 0 |
| 4 | 100 | half | 0 | s>=0 | 0.718 | -- | 1.70 | [1.66, 1.75] | 1.57 | 2.37 | 0 |
| 4 | 100 | half | 0 | s>=1 | 0.865 | -- | 2.61 | [2.52, 2.68] | 2.37 | 3.01 | 0 |

## TD -- design chart check (5x5): smallest T2 with s >= 0 / s >= 1, Mode R (DEJMPS)

Mean rule: mean F*(k) >= F_th. Median rule: P(F*(k) >= F_th) >= 1/2. CI: 95% bootstrap over deliveries (replicates that do not cross inside the grid are counted as 'dropped').

| r | R (Hz) | channel | kappa (s/km) | target | T2 closed form, optimistic (s) | T2 pessimistic, depol2 (s) | T2 DES, mean rule (s) | 95% CI | T2 DES, median rule (s) | DES / closed | boot dropped |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 1000 | deph2 | 0.002 | s>=1 | 0.0925 | -- | 0.130 | [0.128, 0.134] | 0.124 | 1.41 | 0 |
| 2 | 1000 | depol2 | 0.002 | s>=1 | 0.150 | -- | 0.241 | [0.236, 0.245] | 0.231 | 1.61 | 0 |
| 2 | 1000 | half | 0.002 | s>=1 | 0.0444 | -- | 0.0727 | [0.0712, 0.0740] | 0.0704 | 1.64 | 0 |
| 2 | 100 | deph2 | 0.002 | s>=1 | 0.334 | -- | 0.695 | [0.666, 0.721] | 0.644 | 2.08 | 0 |
| 2 | 100 | depol2 | 0.002 | s>=1 | 0.539 | -- | 1.43 | [1.40, 1.48] | 1.34 | 2.66 | 0 |
| 2 | 100 | half | 0.002 | s>=1 | 0.160 | -- | 0.449 | [0.436, 0.467] | 0.423 | 2.80 | 0 |
| 3 | 1000 | deph2 | 0.002 | s>=1 | 0.0975 | -- | 0.130 | [0.127, 0.134] | 0.122 | 1.33 | 0 |
| 3 | 1000 | depol2 | 0.002 | s>=1 | 0.160 | -- | 0.273 | [0.269, 0.278] | 0.259 | 1.71 | 0 |
| 3 | 1000 | half | 0.002 | s>=1 | 0.0481 | -- | 0.0859 | [0.0843, 0.0874] | 0.0819 | 1.78 | 0 |
| 3 | 100 | deph2 | 0.002 | s>=1 | 0.496 | -- | 0.862 | [0.835, 0.891] | 0.807 | 1.74 | 0 |
| 3 | 100 | depol2 | 0.002 | s>=1 | 0.811 | -- | 1.99 | [1.94, 2.04] | 1.86 | 2.45 | 0 |
| 3 | 100 | half | 0.002 | s>=1 | 0.244 | -- | 0.621 | [0.608, 0.634] | 0.584 | 2.54 | 0 |

## T6 -- design chart check (6x6): smallest T2 with s >= 0 / s >= 1, Mode R

Mean rule: mean F*(k) >= F_th. Median rule: P(F*(k) >= F_th) >= 1/2. CI: 95% bootstrap over deliveries (replicates that do not cross inside the grid are counted as 'dropped').

| r | R (Hz) | channel | kappa (s/km) | target | T2 closed form, optimistic (s) | T2 pessimistic, depol2 (s) | T2 DES, mean rule (s) | 95% CI | T2 DES, median rule (s) | DES / closed | boot dropped |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 1000 | deph2 | 0.002 | s>=1 | 0.435 | -- | 0.997 | [0.977, 1.01] | 0.951 | 2.29 | 0 |
| 4 | 1000 | depol2 | 0.002 | s>=1 | 0.672 | 2.64 | 1.51 | [1.48, 1.54] | 1.45 | 2.25 | 0 |
| 4 | 100 | deph2 | 0.002 | s>=1 | 2.73 | -- | 8.31 | [8.08, 8.52] | 7.73 | 3.04 | 0 |
| 4 | 100 | depol2 | 0.002 | s>=1 | 4.22 | 23.9 | 12.8 | [12.3, 13.1] | 11.9 | 3.02 | 0 |
| 5 | 1000 | deph2 | 0.002 | s>=1 | 0.460 | -- | 1.17 | [1.15, 1.19] | 1.13 | 2.55 | 0 |
| 5 | 1000 | depol2 | 0.002 | s>=1 | 0.722 | 3.41 | 1.77 | [1.74, 1.80] | 1.69 | 2.45 | 0 |
| 5 | 100 | deph2 | 0.002 | s>=1 | 3.54 | -- | 10.4 | [10.1, 10.7] | 9.84 | 2.95 | 0 |
| 5 | 100 | depol2 | 0.002 | s>=1 | 5.56 | 32.4 | 16.4 | [15.9, 16.8] | 15.5 | 2.95 | 0 |

## T7 -- design chart check (7x7): smallest T2 with s >= 0 / s >= 1, Mode R

Mean rule: mean F*(k) >= F_th. Median rule: P(F*(k) >= F_th) >= 1/2. CI: 95% bootstrap over deliveries (replicates that do not cross inside the grid are counted as 'dropped').

| r | R (Hz) | channel | kappa (s/km) | target | T2 closed form, optimistic (s) | T2 pessimistic, depol2 (s) | T2 DES, mean rule (s) | 95% CI | T2 DES, median rule (s) | DES / closed | boot dropped |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 1000 | deph2 | 0.002 | s>=1 | 0.571 | -- | 1.57 | [1.55, 1.60] | 1.50 | 2.76 | 0 |
| 5 | 1000 | depol2 | 0.002 | s>=1 | 0.890 | 4.69 | 2.39 | [2.35, 2.43] | 2.28 | 2.68 | 0 |
| 5 | 100 | deph2 | 0.002 | s>=1 | 4.23 | -- | 14.3 | [14.0, 14.6] | 13.7 | 3.38 | 0 |
| 5 | 100 | depol2 | 0.002 | s>=1 | 6.59 | 44.6 | 22.1 | [21.6, 22.5] | 21.1 | 3.35 | 0 |

## S1 -- Mode F routing, 5x5, r = 3, T2 = 0.1 s, half (the paper's E1 cell, same seeds)

n = 5000 trials (seeds 20260922 + t, identical to the paper's E1).

| policy | success DES (%) | success qunet_tqe (%) | success paper E1 (%) | DES vs E1 mismatching trials |
|---|---|---|---|---|
| shortest | 73.16 | 73.16 | 73.16 | 0 |
| baseline | 93.24 | 93.24 | 93.24 | 0 |
| full | 96.62 | 96.62 | 96.62 | 0 |
| dijkstra | 98.24 | 98.24 | 98.24 | 0 |
| oracle | 98.24 | 98.24 | 98.24 | 0 |

## S2 -- Mode R routing, 5x5, r = 3, R = 1 kHz: headroom, latency, delivery rate

'oracle' = the paper model's exact oracle path (model at the same R), evaluated in the DES; it need not be the DES-optimal path. Latency = t_used = establishment + kappa*D storage. Delivery rate = successes / total attempt time (attempts back to back).


### half|k0.002|R1000: DES thresholds T2(s>=0) = 0.1785 s, T2(s>=1) = 0.3298 s

**s0** (T2 = 0.2426 s; intended DES regime s = 0; verified: True), n = 3000

DES clean-path profile at this T2: F*(0) = 0.8182 +- 0.0007, P(F>=0.8) = 0.90; F*(1) = 0.7858 +- 0.0007, P(F>=0.8) = 0.14

| policy | success % [95% CI] | mean hops | latency mean / median (ms) | delivery rate (pairs/s) | diagnostic: rate over clean-path attempts only | raw pairs per attempt |
|---|---|---|---|---|---|---|
| shortest | 66.97 [65.3, 68.6] | 4.00 | 28.63 / 26.17 | 23.39 | 32.99 | 42.1 |
| baseline | 68.97 [67.3, 70.6] | 4.29 | 28.36 / 26.83 | 24.32 | 25.10 | 45.8 |
| full | 69.47 [67.8, 71.1] | 4.28 | 27.92 / 26.77 | 24.88 | 24.89 | 45.3 |
| dijkstra | 69.97 [68.3, 71.6] | 4.26 | 27.83 / 26.70 | 25.14 | 25.15 | 45.0 |
| oracle | 70.00 [68.3, 71.6] | 4.25 | 27.81 / 26.64 | 25.17 | 25.63 | 44.9 |

Headroom oracle - BFS = +3.03 pp, exact 95% conditional on the discordant trials [+2.79, +3.03] (its upper limit equals the estimate when no trial favours BFS), exact McNemar p = 8.08e-28 (oracle-only 91, BFS-only 0; 91 of the oracle-only successes used a detour). PG - BFS = +2.50 pp, p = 5.29e-23. P(success | oracle path is a clean one-hop detour, swaps ok) = 0.122.

Establishment time (without kappa*D) of clean paths by hop count (ms, mean +- SE, n): h=4: 15.93 +- 0.11 (2236); h=5: 17.15 +- 0.15 (1420); h=6: 17.30 +- 0.47 (128); h=7: 18.02 +- 1.69 (12); h=8: 20.72 +- 0.00 (1)

Latency incl. kappa*D of clean paths by hop count (ms, mean +- SE, n): h=4: 27.24 +- 0.11 (2236); h=5: 29.63 +- 0.15 (1420); h=6: 31.61 +- 0.47 (128); h=7: 34.09 +- 1.88 (12); h=8: 41.69 +- 0.00 (1); paths with a noisy link: 32.77 ms (780); BFS paths with a noisy link: 32.70 ms (764).

**s1** (T2 = 0.4947 s; intended DES regime s = 1; verified: True), n = 3000

DES clean-path profile at this T2: F*(0) = 0.8445 +- 0.0004, P(F>=0.8) = 1.00; F*(1) = 0.8147 +- 0.0004, P(F>=0.8) = 0.95; F*(2) = 0.7858 +- 0.0004, P(F>=0.8) = 0.00

| policy | success % [95% CI] | mean hops | latency mean / median (ms) | delivery rate (pairs/s) | diagnostic: rate over clean-path attempts only | raw pairs per attempt |
|---|---|---|---|---|---|---|
| shortest | 73.37 [71.8, 74.9] | 4.00 | 28.60 / 26.18 | 25.65 | 36.13 | 41.9 |
| baseline | 91.87 [90.8, 92.8] | 4.29 | 28.14 / 26.86 | 32.64 | 33.64 | 45.5 |
| full | 94.90 [94.1, 95.6] | 4.28 | 27.83 / 26.76 | 34.10 | 34.11 | 45.1 |
| dijkstra | 96.83 [96.1, 97.4] | 4.26 | 27.80 / 26.80 | 34.83 | 34.84 | 44.9 |
| oracle | 96.77 [96.1, 97.3] | 4.26 | 27.80 / 26.73 | 34.80 | 35.37 | 44.9 |

Headroom oracle - BFS = +23.40 pp, exact 95% conditional on the discordant trials [+23.15, +23.40] (its upper limit equals the estimate when no trial favours BFS), exact McNemar p = 9.51e-212 (oracle-only 702, BFS-only 0; 702 of the oracle-only successes used a detour). PG - BFS = +21.53 pp, p = 6.85e-195. P(success | oracle path is a clean one-hop detour, swaps ok) = 0.944.

Establishment time (without kappa*D) of clean paths by hop count (ms, mean +- SE, n): h=4: 15.93 +- 0.12 (2236); h=5: 16.74 +- 0.14 (1436); h=6: 17.65 +- 0.43 (127); h=7: 17.75 +- 1.27 (12); h=8: 14.40 +- 0.00 (1)

Latency incl. kappa*D of clean paths by hop count (ms, mean +- SE, n): h=4: 27.24 +- 0.12 (2236); h=5: 29.22 +- 0.14 (1436); h=6: 31.97 +- 0.43 (127); h=7: 33.82 +- 1.37 (12); h=8: 35.37 +- 0.00 (1); paths with a noisy link: 32.67 ms (775); BFS paths with a noisy link: 32.58 ms (764).


### deph2|k0.002|R1000: DES thresholds T2(s>=0) = 0.3746 s, T2(s>=1) = 0.7312 s

**s0** (T2 = 0.5233 s; intended DES regime s = 0; verified: True), n = 3000

DES clean-path profile at this T2: F*(0) = 0.8190 +- 0.0007, P(F>=0.8) = 0.90; F*(1) = 0.7854 +- 0.0007, P(F>=0.8) = 0.12

| policy | success % [95% CI] | mean hops | latency mean / median (ms) | delivery rate (pairs/s) | diagnostic: rate over clean-path attempts only | raw pairs per attempt |
|---|---|---|---|---|---|---|
| shortest | 66.60 [64.9, 68.3] | 4.00 | 28.92 / 26.26 | 23.03 | 32.47 | 42.1 |
| baseline | 68.47 [66.8, 70.1] | 4.29 | 28.49 / 26.84 | 24.03 | 24.81 | 45.7 |
| full | 68.70 [67.0, 70.3] | 4.28 | 28.05 / 26.79 | 24.49 | 24.50 | 45.2 |
| dijkstra | 68.80 [67.1, 70.4] | 4.26 | 28.00 / 26.74 | 24.57 | 24.58 | 45.0 |
| oracle | 68.73 [67.1, 70.4] | 4.25 | 28.01 / 26.77 | 24.54 | 24.98 | 45.0 |

Headroom oracle - BFS = +2.13 pp, exact 95% conditional on the discordant trials [+1.89, +2.13] (its upper limit equals the estimate when no trial favours BFS), exact McNemar p = 1.08e-19 (oracle-only 64, BFS-only 0; 64 of the oracle-only successes used a detour). PG - BFS = +2.10 pp, p = 2.17e-19. P(success | oracle path is a clean one-hop detour, swaps ok) = 0.086.

Establishment time (without kappa*D) of clean paths by hop count (ms, mean +- SE, n): h=4: 16.20 +- 0.12 (2236); h=5: 16.72 +- 0.15 (1430); h=6: 17.70 +- 0.64 (130); h=7: 20.12 +- 1.55 (12); h=8: 14.29 +- 0.00 (1)

Latency incl. kappa*D of clean paths by hop count (ms, mean +- SE, n): h=4: 27.52 +- 0.12 (2236); h=5: 29.21 +- 0.15 (1430); h=6: 32.01 +- 0.65 (130); h=7: 36.19 +- 1.53 (12); h=8: 35.26 +- 0.00 (1); paths with a noisy link: 33.24 ms (780); BFS paths with a noisy link: 33.02 ms (764).

**s1** (T2 = 1.097 s; intended DES regime s = 1; verified: True), n = 3000

DES clean-path profile at this T2: F*(0) = 0.8455 +- 0.0004, P(F>=0.8) = 1.00; F*(1) = 0.8154 +- 0.0004, P(F>=0.8) = 0.96; F*(2) = 0.7861 +- 0.0004, P(F>=0.8) = 0.00

| policy | success % [95% CI] | mean hops | latency mean / median (ms) | delivery rate (pairs/s) | diagnostic: rate over clean-path attempts only | raw pairs per attempt |
|---|---|---|---|---|---|---|
| shortest | 73.37 [71.8, 74.9] | 4.00 | 29.19 / 26.20 | 25.14 | 36.30 | 42.3 |
| baseline | 92.00 [91.0, 92.9] | 4.29 | 28.23 / 26.72 | 32.59 | 33.83 | 45.4 |
| full | 95.37 [94.6, 96.1] | 4.28 | 27.70 / 26.67 | 34.43 | 34.45 | 44.8 |
| dijkstra | 97.13 [96.5, 97.7] | 4.26 | 27.68 / 26.65 | 35.09 | 35.10 | 44.6 |
| oracle | 97.17 [96.5, 97.7] | 4.26 | 27.63 / 26.62 | 35.17 | 35.74 | 44.5 |

Headroom oracle - BFS = +23.80 pp, exact 95% conditional on the discordant trials [+23.55, +23.80] (its upper limit equals the estimate when no trial favours BFS), exact McNemar p = 2.32e-215 (oracle-only 714, BFS-only 0; 714 of the oracle-only successes used a detour). PG - BFS = +22.00 pp, p = 4.18e-199. P(success | oracle path is a clean one-hop detour, swaps ok) = 0.960.

Establishment time (without kappa*D) of clean paths by hop count (ms, mean +- SE, n): h=4: 15.81 +- 0.11 (2236); h=5: 16.65 +- 0.14 (1440); h=6: 17.31 +- 0.50 (128); h=7: 16.66 +- 1.44 (12); h=8: 14.86 +- 0.00 (1)

Latency incl. kappa*D of clean paths by hop count (ms, mean +- SE, n): h=4: 27.12 +- 0.11 (2236); h=5: 29.14 +- 0.14 (1440); h=6: 31.63 +- 0.51 (128); h=7: 32.74 +- 1.63 (12); h=8: 35.83 +- 0.00 (1); paths with a noisy link: 35.42 ms (775); BFS paths with a noisy link: 35.24 ms (764).


### depol2|k0.002|R1000: DES thresholds T2(s>=0) = 0.5891 s, T2(s>=1) = 1.088 s

**s0** (T2 = 0.8007 s; intended DES regime s = 0; verified: True), n = 3000

DES clean-path profile at this T2: F*(0) = 0.8195 +- 0.0007, P(F>=0.8) = 0.91; F*(1) = 0.7860 +- 0.0007, P(F>=0.8) = 0.12

| policy | success % [95% CI] | mean hops | latency mean / median (ms) | delivery rate (pairs/s) | diagnostic: rate over clean-path attempts only | raw pairs per attempt |
|---|---|---|---|---|---|---|
| shortest | 67.03 [65.3, 68.7] | 4.00 | 28.64 / 26.20 | 23.40 | 32.84 | 42.1 |
| baseline | 69.17 [67.5, 70.8] | 4.29 | 28.36 / 26.85 | 24.39 | 25.21 | 45.8 |
| full | 69.40 [67.7, 71.0] | 4.28 | 27.90 / 26.80 | 24.87 | 24.88 | 45.2 |
| dijkstra | 69.73 [68.1, 71.4] | 4.26 | 27.87 / 26.78 | 25.02 | 25.03 | 45.0 |
| oracle | 69.70 [68.0, 71.3] | 4.25 | 27.85 / 26.77 | 25.02 | 25.50 | 45.0 |

Headroom oracle - BFS = +2.67 pp, exact 95% conditional on the discordant trials [+2.43, +2.67] (its upper limit equals the estimate when no trial favours BFS), exact McNemar p = 1.65e-24 (oracle-only 80, BFS-only 0; 80 of the oracle-only successes used a detour). PG - BFS = +2.37 pp, p = 8.47e-22. P(success | oracle path is a clean one-hop detour, swaps ok) = 0.108.

Establishment time (without kappa*D) of clean paths by hop count (ms, mean +- SE, n): h=4: 16.07 +- 0.12 (2236); h=5: 16.68 +- 0.14 (1440); h=6: 16.70 +- 0.49 (127); h=7: 17.83 +- 1.68 (12); h=8: 17.93 +- 0.00 (1)

Latency incl. kappa*D of clean paths by hop count (ms, mean +- SE, n): h=4: 27.39 +- 0.12 (2236); h=5: 29.16 +- 0.14 (1440); h=6: 31.02 +- 0.49 (127); h=7: 33.90 +- 1.78 (12); h=8: 38.90 +- 0.00 (1); paths with a noisy link: 32.58 ms (780); BFS paths with a noisy link: 32.31 ms (764).

**s1** (T2 = 1.632 s; intended DES regime s = 1; verified: True), n = 3000

DES clean-path profile at this T2: F*(0) = 0.8456 +- 0.0003, P(F>=0.8) = 1.00; F*(1) = 0.8148 +- 0.0004, P(F>=0.8) = 0.96; F*(2) = 0.7855 +- 0.0004, P(F>=0.8) = 0.00

| policy | success % [95% CI] | mean hops | latency mean / median (ms) | delivery rate (pairs/s) | diagnostic: rate over clean-path attempts only | raw pairs per attempt |
|---|---|---|---|---|---|---|
| shortest | 73.40 [71.8, 75.0] | 4.00 | 28.60 / 26.17 | 25.66 | 36.24 | 41.9 |
| baseline | 92.23 [91.2, 93.1] | 4.29 | 28.13 / 26.76 | 32.78 | 33.86 | 45.4 |
| full | 95.43 [94.6, 96.1] | 4.28 | 27.72 / 26.70 | 34.43 | 34.44 | 44.9 |
| dijkstra | 97.20 [96.5, 97.7] | 4.26 | 27.65 / 26.58 | 35.16 | 35.17 | 44.7 |
| oracle | 97.23 [96.6, 97.8] | 4.26 | 27.66 / 26.59 | 35.16 | 35.73 | 44.7 |

Headroom oracle - BFS = +23.83 pp, exact 95% conditional on the discordant trials [+23.59, +23.83] (its upper limit equals the estimate when no trial favours BFS), exact McNemar p = 1.16e-215 (oracle-only 715, BFS-only 0; 715 of the oracle-only successes used a detour). PG - BFS = +22.03 pp, p = 2.09e-199. P(success | oracle path is a clean one-hop detour, swaps ok) = 0.961.

Establishment time (without kappa*D) of clean paths by hop count (ms, mean +- SE, n): h=4: 15.86 +- 0.12 (2236); h=5: 16.62 +- 0.14 (1454); h=6: 17.46 +- 0.44 (133); h=7: 16.43 +- 2.32 (12); h=8: 11.34 +- 0.00 (1)

Latency incl. kappa*D of clean paths by hop count (ms, mean +- SE, n): h=4: 27.17 +- 0.12 (2236); h=5: 29.11 +- 0.14 (1454); h=6: 31.75 +- 0.44 (133); h=7: 32.50 +- 2.53 (12); h=8: 32.31 +- 0.00 (1); paths with a noisy link: 32.89 ms (775); BFS paths with a noisy link: 32.77 ms (764).

