# Discrete-event validation (SeQUeNCe kernel, density matrices)

The Monte Carlo simulator in `../code` and the closed form in the paper implement the same scalar
physics, so their agreement is a consistency check, not independent evidence. This package
re-evaluates the paper's regime calls in a discrete-event, density-matrix simulation built on the
kernel of [SeQUeNCe](https://github.com/sequence-toolbox/SeQUeNCe) 1.3:
- purification and swapping are gate circuits with Born-rule measurement outcomes;
- every stored pair decoheres for the time it is actually stored;
- the delivered fidelity is read from the density matrix.

It reads `../code` and `../results_tqe`, never modifies them, and writes only to `results_des/`.
That folder holds the stored full run used in the paper (Section VII-L, Tables IX–XI, Fig. 10 and Supplementary Tables S2–S6).

---

## 1. Run

SeQUeNCe 1.3 needs **Python ≥ 3.12**, so this package has its own environment.

```bash
cd des_validation
conda env create -f environment.yml        # creates "qrf-des"
conda activate qrf-des

python smoke_test.py          # ~10 s;  must end with "SMOKE TEST PASSED"
python test_des_physics.py    # ~1 min; must end with "18/18 passed"
python run_des.py --quick     # ~5-10 min pipeline check with small n (NOT for the paper)
python run_des.py             # blocks PF, PR, T, S1, S2 (see timing below)
python run_des.py --only S3    # softness sweep across the boundary (needs T; about 1 h on 2 workers)
python run_des.py --only TD,T6a,T6b,T7a,T7b   # the further blocks, run separately
python analyze.py             # summary.md, tables, macros, figures in results_des/
python paper_extra.py         # further macros and the softness table (des_extra.tex, tab_des_soft.tex)
```

- `run_des.py` uses all cores but one; change this with `--workers N` (capped at 61 on Windows).
  It uses the `spawn` start method on every OS, so launch it from a terminal, not an interactive
  window.
- Blocks can be run separately with `--only PF,PR,T,S1,S2,TD,T6a,T6b,T7a,T7b` (T6 and T7 run as two halves each, which `analyze.py` merges). **S2 needs T**, because it takes its
  coherence times from T's thresholds.
- A `--quick` run overwrites `results_des/`. Use `analyze.py --results DIR` to analyse another
  folder. `analyze.py` prints a WARNING when a folder holds `--quick` data.
- The simulation is CPU-only; a GPU is not used.
- Tested versions are in `requirements-tested.txt`. The smoke test warns if pip resolves a
  different SeQUeNCe version.

**Timing.** The full run is about 10–12 CPU-core-hours. The run used in the paper took 27 min on
23 worker processes.

| Block | What | n | Share of time |
|---|---|---|---|
| PF | Mode F clean-path profiles: 3 protocols, gate noise | deterministic | ~3 % |
| PR | Mode R at the paper's T₂, R = 10⁵ and 10⁴ Hz | 400 per path | ~20 % |
| T  | Mode R design-chart thresholds (5×5) | 400 per point, 528 points | ~55 % |
| S1 | Mode F routing on the paper's 5×5, r = 3 trials (same seeds) | 5000 trials | ~2 % |
| S2 | Mode R routing: 3 memory channels × 2 regimes | 3000 trials each | ~20 % |
| TD | Mode R design-chart thresholds with DEJMPS (5×5, r = 2, 3) | 400 per point, 132 points | run separately: 33 min on 2 workers (Linux, Python 3.13) |
| T6 | Mode R design-chart thresholds on a second grid size (6×6, BBPSSW r = 4, 5, per-qubit channels, 100 Hz and 1 kHz) | 400 per point | run separately as T6a + T6b: about 1.75 h on 2 workers |
| T7 | the same on 7×7 (BBPSSW r = 5) | 400 per point | run separately as T7a + T7b: about 2.1 h on 2 workers |

---

## 2. What is independent, and what is shared

| Component | Source |
|---|---|
| Event scheduling | SeQUeNCe `Timeline` / `Event` / `Process` |
| States, gates, measurement | SeQUeNCe `QuantumManagerDensity`, `Circuit` (Born-rule sampling) |
| Purification | gate circuits: bilateral CNOT (+ DEJMPS rotations), Z measurement, keep if outcomes agree |
| Swapping | CNOT + H + two measurements, Pauli correction X^m2 Z^m1 |
| Memory noise, twirl, gate noise | Kraus channels in `des_sim.py`, unit-tested against their definitions |
| **Not used** | any closed form (BBPSSW recurrence, Werner composition, decay law). Also not SeQUeNCe's `Memory.fidelity` bookkeeping: its swapping multiplies fidelities (F₁·F₂·degradation), which is only a lower bound on the Werner result. |
| Shared with the paper's model | inputs and algorithms, not results: raw-pair fidelity, network sampling, the routers, the threshold 0.8 |

The simulation runs in two modes:

- **Mode F** uses the model's assumptions: twirled BBPSSW, instantaneous purification, and storage
  κD charged to the delivered pair only. It can disagree with the closed form only through an
  error in the algebra or its implementation. It is an implementation cross-check, not
  independent evidence.
- **Mode R** resolves the storage time:
  - raw pairs are heralded as a Poisson process at rate R;
  - each link runs the nested recurrence, so input pairs wait for partners;
  - a node swaps as soon as both of its links are ready;
  - classical messages travel at 2×10⁵ km/s;
  - every stored pair decoheres for exactly its own storage time;
  - κD is still charged to the delivered pair.

---

## 3. Results of the full run (`results_des/summary.md`: 0 FAIL)

- **PF (Mode F).**
  - Every clean-path fidelity matches the closed form to 2×10⁻¹⁵.
  - All 12 regime calls (4 configurations × 3 memory channels) match, as do the 12 calls with
    gate error.
  - DEJMPS raises F*(0) on 5×5, r = 3 from 0.831 to 0.934 and gives slack ≥ 8 in every
    configuration, including 6×6, r = 3, which is infeasible under BBPSSW.
  - At 1 % gate error every BBPSSW configuration is infeasible, while DEJMPS keeps a slack of 3–5.
  - A bilateral CNOT without twirl or rotation is infeasible everywhere.
- **PR (Mode R, the paper's T₂).**
  - At R = 10⁵ Hz all 12 calls hold; the delivered fidelity is 0.0008–0.0025 lower.
  - At 10⁴ Hz it is 0.005–0.020 lower, and 5 of 12 calls drop by one.
  - The delivery-to-delivery standard deviation of the fidelity is 0.003–0.008 at 10⁴ Hz.
- **T (design chart).**
  - Routing needs 1.8–2.7× the closed-form coherence time under per-qubit channels (5×5,
    100 Hz–1 kHz, κ = 2 ms/km).
  - The requirement is 2.6–2.8× without the κD term (2.6–3.2× including the baseline map).
  - It lies between the closed form and the pessimistic all-links-wait accounting of `code/design_law.py` (every channel, κ = 2 ms/km) in 12 of 12 cases on 5×5, and 12 of 12 on 6×6 and 7×7 (T6, T7; factor 2.3–3.4). `paper_extra.py` computes these counts.
- **TD (DEJMPS, 5×5, r = 2, 3).** Routing needs 1.3–2.7× the closed-form coherence time under per-qubit channels, and, at three rounds with the same rate and channel, a coherence time 3.8–6.0× shorter than BBPSSW in the same simulation.
- **S1.** 0 mismatches over 25,000 trial–policy pairs, both against `qunet_tqe` and against the
  stored `results_tqe/trials/E1.npz`.
- **S2 (Mode R routing, 1 kHz).**
  - The oracle's headroom over shortest-path routing is +23.4 to +23.8 pp where the simulation
    gives s = 1.
  - Where s = 0 it is +2.1 to +3.0 pp (exact McNemar P < 10⁻¹⁸). All extra successes are one-hop
    detours whose mean fidelity is below threshold but which pass in 9–12 % of attempts. Near the
    boundary, the regime call is soft.
- **S3 (softness sweep).** The S2 setting with per-qubit depolarizing memories at six
  coherence times across the s = 0 → 1 boundary (2000 paired trials each): see `soft_sweep.json` and
  `tab_des_soft.tex` (Table XI of the paper).

## 4. Outputs (`results_des/`)

| File | Contents |
|---|---|
| `summary.md` | every number and every automatic check (PASS / FAIL / NOTE) |
| `tab_des_*.tex`, `des_numbers.tex`, `des_extra.tex` | the paper's discrete-event tables and the macros used in its text (`des_extra.tex` from `paper_extra.py`) |
| `fig_des_threshold.pdf` | F*(1) versus T2 on 5×5, r = 3 at 100 Hz and 1 kHz (not in the paper) |
| `fig_des_profiles.pdf` | F*(k) by protocol on 5×5, r = 3 (Fig. 10 of the paper) |
| `des_paragraph_draft.md` | auto-generated text; every qualitative sentence is conditional on the numbers |
| `PF.json`, `PR.json`, `T.json`, `S1.json`, `S2.json`, `TD.json`, `T6a.json`, `T6b.json`, `T7a.json`, `T7b.json` | raw results (T and TD keep every delivery's fidelity for the bootstrap) |
| `design_points.json`, `design_points_dejmps.json`, `design_points_six.json` | design-chart thresholds used by the paper's tables and `paper_extra.py` |

## 5. Model choices (stated in the paper's limitations)

- **Protocol stack in Mode R:** Poisson heralding, nested (symmetric) recurrence, swap as soon as
  possible, no cutoff, unlimited memories. Other schedules would move the thresholds.
- **Links with raw F ≤ ½ are not purified**, which is the paper's own rule. Noisy links with F in
  (½, 0.55] are purified, slowly.
- **The baseline memory map** F → ½ + (F − ½)e^(−t/T₂) is not a per-qubit channel. In Mode R it is
  a two-qubit replacement channel toward Werner(½), with the Pauli frame corrected at the swap.
  Physical conclusions rest on the per-qubit dephasing and depolarizing channels.
- **"Oracle" in Mode R** is the exact oracle of the scalar model, evaluated in the simulation. It
  need not be the best path there: Dijkstra beats it by at most 0.07 pp in 3 of 6 settings.
- **The regime call in Mode R** uses the mean delivered fidelity. The median rule is reported
  alongside it.

## 6. Files

| File | Role |
|---|---|
| `des_sim.py` | the simulation (one delivery attempt along a path); Mode F and Mode R are parameter sets |
| `qunet_bridge.py` | the only coupling to `../code`: routers, network sampling, closed-form references |
| `test_des_physics.py` | 18 tests (see below) |
| `run_des.py` | the study (blocks PF, PR, T, S1, S2, TD, T6, T7) |
| `analyze.py` | summary, tables, macros, figures |
| `smoke_test.py` | environment and process-pool check |
| `environment.yml`, `requirements-tested.txt` | environment |

The 18 tests in `test_des_physics.py` cover:
- the channels: definition, CPTP, semigroup, agreement with the paper;
- the twirl;
- BBPSSW against the recurrence and its success probability;
- DEJMPS against Deutsch et al.;
- swap corrections for all 16 Bell inputs × 4 outcomes;
- Werner composition;
- Mode F paths against the closed form and `qunet_tqe`;
- exact Mode R storage accounting;
- the faulty-link rule.
