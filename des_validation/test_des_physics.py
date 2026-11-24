"""
test_des_physics.py -- checks of every physical primitive of des_sim.py, then of
whole delivery attempts, against the paper's closed forms (qunet_tqe.py,
design_law.py) and against exact identities.

Run:  python test_des_physics.py          (about 1-2 min; exits non-zero on failure)
  or  python -m pytest -q test_des_physics.py
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

import des_sim as D
import qunet_bridge as B

Q = B.Q
RNG = np.random.default_rng(12345)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def rand_state(n_qubits=2):
    d = 2 ** n_qubits
    A = RNG.normal(size=(d, d)) + 1j * RNG.normal(size=(d, d))
    rho = A @ A.conj().T
    return rho / np.trace(rho)


BELL = {
    "phi+": np.array([1, 0, 0, 1]) / math.sqrt(2),
    "phi-": np.array([1, 0, 0, -1]) / math.sqrt(2),
    "psi+": np.array([0, 1, 1, 0]) / math.sqrt(2),
    "psi-": np.array([0, 1, -1, 0]) / math.sqrt(2),
}


def bell_diag(coeffs: dict):
    return sum(c * np.outer(BELL[k], BELL[k].conj()) for k, c in coeffs.items()).astype(complex)


def bell_coeffs(rho):
    return {k: float(np.real(v.conj() @ rho @ v)) for k, v in BELL.items()}


def partial_trace_qubit(rho, q):
    t = rho.reshape(2, 2, 2, 2)
    return np.trace(t, axis1=0, axis2=2) if q == 0 else np.trace(t, axis1=1, axis2=3)


def choi(channel):
    """Choi matrix of a map on two qubits (4x4 -> 4x4)."""
    C = np.zeros((16, 16), dtype=complex)
    for i in range(4):
        for j in range(4):
            E = np.zeros((4, 4), dtype=complex)
            E[i, j] = 1.0
            C += np.kron(E, channel(E))
    return C


def new_tl():
    from sequence.kernel.timeline import Timeline
    return Timeline(formalism=D.DENSITY_MATRIX_FORMALISM)


# ---------------------------------------------------------------------------
# 1. channels
# ---------------------------------------------------------------------------
def test_channel_definitions():
    for _ in range(5):
        rho = rand_state()
        for q in (0, 1):
            c = float(RNG.uniform(0.1, 0.99))
            # depolarizing: c rho + (1-c) Tr_q rho (x) I/2
            red = partial_trace_qubit(rho, q)
            ref = c * rho + (1 - c) * (np.kron(np.eye(2) / 2, red) if q == 0 else np.kron(red, np.eye(2) / 2))
            assert np.allclose(D.depolarize_qubit(rho, q, c), ref, atol=1e-13)
            # dephasing: coherences of qubit q scaled by c
            out = D.dephase_qubit(rho, q, c).reshape(2, 2, 2, 2)
            r4 = rho.reshape(2, 2, 2, 2)
            for a in range(2):
                for b in range(2):
                    for a2 in range(2):
                        for b2 in range(2):
                            same = (a == a2) if q == 0 else (b == b2)
                            f = 1.0 if same else c
                            assert abs(out[a, b, a2, b2] - f * r4[a, b, a2, b2]) < 1e-13


def test_channels_cptp_and_semigroup():
    for mem in ("half", "deph2", "depol2"):
        ch = lambda E, mem=mem: D.memory_channel(E, mem, 0.03, 0.1)
        C = choi(ch)
        assert np.min(np.linalg.eigvalsh((C + C.conj().T) / 2)) > -1e-12, mem          # CP
        for _ in range(3):
            rho = rand_state()
            assert abs(np.trace(ch(rho)) - 1) < 1e-13                                      # TP
            a, b = float(RNG.uniform(0, 0.2)), float(RNG.uniform(0, 0.2))
            two = D.memory_channel(D.memory_channel(rho, mem, a, 0.1), mem, b, 0.1)
            one = D.memory_channel(rho, mem, a + b, 0.1)
            assert np.allclose(two, one, atol=1e-13), mem                                    # semigroup


def test_memory_channels_vs_paper_decay():
    for mem in ("half", "deph2", "depol2"):
        ph = Q.Physics(memory=mem, T2=0.1)
        for F in (0.6, 0.8, 0.9, 0.97):
            for t in (0.0, 0.003, 0.02, 0.1, 0.5):
                got = D.fidelity(D.memory_channel(D.werner(F), mem, t, 0.1))
                assert abs(got - ph.decay(F, t)) < 1e-13, (mem, F, t, got, ph.decay(F, t))


# ---------------------------------------------------------------------------
# 2. twirl
# ---------------------------------------------------------------------------
def test_twirl():
    for _ in range(5):
        rho = rand_state()
        tw = D.twirl(rho)
        assert np.allclose(tw, D.werner(D.fidelity(rho)), atol=1e-13)
        assert np.allclose(D.twirl(tw), tw, atol=1e-13)
    assert np.allclose(D.twirl(D.werner(0.83)), D.werner(0.83), atol=1e-13)


# ---------------------------------------------------------------------------
# 3. purification circuits through SeQUeNCe
# ---------------------------------------------------------------------------
def purify_once(rho1, rho2, scheme, samp):
    """Run the DES purification circuit on (kept = rho1, meas = rho2); returns
    (outcomes equal?, kept pair state)."""
    tl = new_tl()
    qm = tl.quantum_manager
    a1, b1, a2, b2 = qm.new(), qm.new(), qm.new(), qm.new()
    qm.set([a1, b1], rho1)
    qm.set([a2, b2], rho2)
    out = qm.run_circuit(D.PURIF_CIRCUITS[scheme], [a1, a2, b1, b2], samp)
    qm.discard(a2); qm.discard(b2)
    st = qm.states[a1]
    rho = np.array(st.state)
    if list(st.keys) != [a1, b1]:
        rho = D.SWAP2 @ rho @ D.SWAP2
    return out[a2] == out[b2], rho


def test_bbpssw_twirled_vs_recurrence():
    for F in (0.55, 0.7, 0.85, 0.95):
        Fp_ref, P_ref = Q.bbpssw_round(F)
        assert abs(Fp_ref - Q.noisy_purify_round(F, 0.0)) < 1e-12
        # both success branches (outcome 00 and 11) give the recurrence value
        for samp in (1e-9, 1 - 1e-9):
            ok, rho = purify_once(D.werner(F), D.werner(F), "bbpssw_twirl", samp)
            assert ok
            assert abs(D.fidelity(D.twirl(rho)) - Fp_ref) < 1e-12, (F, samp)
        # success probability (Monte Carlo, 4 sigma)
        n, k = 1500, 0
        for s in RNG.random(n):
            k += purify_once(D.werner(F), D.werner(F), "bbpssw_twirl", float(s))[0]
        se = math.sqrt(P_ref * (1 - P_ref) / n)
        assert abs(k / n - P_ref) < 4 * se + 1e-3, (F, k / n, P_ref)


def test_dejmps_vs_deutsch_formula():
    """DEJMPS (Deutsch et al. 1996): with A, B, C, D the weights of
    Phi+, Psi-, Psi+, Phi-:  A' = (A^2+B^2)/N, B' = 2CD/N, C' = (C^2+D^2)/N,
    D' = 2AB/N, N = (A+B)^2 + (C+D)^2."""
    for _ in range(4):
        w = RNG.dirichlet([8, 1, 1, 1])
        A, Bc, C, Dd = w
        rho = bell_diag({"phi+": A, "psi-": Bc, "psi+": C, "phi-": Dd})
        N = (A + Bc) ** 2 + (C + Dd) ** 2
        ref = {"phi+": (A * A + Bc * Bc) / N, "psi-": 2 * C * Dd / N,
               "psi+": (C * C + Dd * Dd) / N, "phi-": 2 * A * Bc / N}
        for samp in (1e-9, 1 - 1e-9):
            ok, out = purify_once(rho, rho, "dejmps", samp)
            assert ok
            got = bell_coeffs(out)
            for kk in ref:
                assert abs(got[kk] - ref[kk]) < 1e-12, (kk, got, ref)


def test_cnot_protocol_without_twirl():
    """Bilateral CNOT without twirl/rotation on Bell-diagonal input:
    Phi+' = (A^2 + D^2)/N, Phi-' = 2AD/N, N = (A+D)^2 + (B+C)^2 (A: Phi+, D: Phi-,
    B: Psi-, C: Psi+). Detects bit flips only."""
    w = RNG.dirichlet([8, 1, 1, 1])
    A, Bc, C, Dd = w
    rho = bell_diag({"phi+": A, "psi-": Bc, "psi+": C, "phi-": Dd})
    N = (A + Dd) ** 2 + (Bc + C) ** 2
    for samp in (1e-9, 1 - 1e-9):
        ok, out = purify_once(rho, rho, "bbpssw", samp)
        got = bell_coeffs(out)
        assert ok
        assert abs(got["phi+"] - (A * A + Dd * Dd) / N) < 1e-12
        assert abs(got["phi-"] - 2 * A * Dd / N) < 1e-12


def test_gate_noise_round_vs_paper():
    """One twirled BBPSSW round with eps_gate, run inside the DES (single link,
    r = 1), equals qunet_tqe.noisy_purify_round (the paper's density-matrix round)."""
    for eps in (0.002, 0.01, 0.03):
        for F in (0.7, 0.9):
            P = D.mode_F(rounds=1, eps_gate=eps, kappa=0.0)
            res = D.run_path([F], [1.0], P, seed=3)
            assert abs(res.F - Q.noisy_purify_round(F, eps)) < 1e-12, (eps, F, res.F)
    # whole clean paths: closed form with eps_gate (clean links have F >= 0.6, where qunet uses the noisy round)
    T = Q.Topology.grid(5)
    for (h, b), p in B.grid_paths_by_composition(5, 1).items():
        L, _ = B.path_links(T, p)
        ph = Q.Physics(rounds=3, eps_gate=0.005)
        res = D.run_path([ph.clean_F(x) for x in L], L, D.mode_F(rounds=3, eps_gate=0.005), seed=1)
        # qunet_tqe.Physics.purify rounds F to 10 decimals before each cached noisy round,
        # so agreement is ~1e-10, not machine precision
        assert abs(res.F - B.closed_form_path(ph, L)) < 1e-9


def test_faulty_link_not_purified():
    res = D.run_path([0.45, 0.95], [1.0, 1.0], D.mode_F(rounds=3), seed=0)
    assert res.raw_pairs[0] == 1 and res.raw_pairs[1] >= 8
    res = D.run_path([0.45], [1.0], D.mode_F(rounds=3, skip_purify_below_half=False), seed=0)
    assert res.raw_pairs[0] >= 8


def test_unknown_parameter_rejected():
    try:
        D.mode_F(roudns=3)
    except TypeError:
        return
    raise AssertionError("misspelled parameter accepted")


# ---------------------------------------------------------------------------
# 4. swapping through SeQUeNCe
# ---------------------------------------------------------------------------
def swap_once(rho_ab, rho_cd, samp, lam=1.0):
    """Swap at the node holding b and c; returns (outcomes, state of (a, d))."""
    tl = new_tl()
    qm = tl.quantum_manager
    a, b, c, d = qm.new(), qm.new(), qm.new(), qm.new()
    qm.set([a, b], D.depolarize_qubit(rho_ab, 1, lam))
    qm.set([c, d], rho_cd)
    out = qm.run_circuit(D.BSM, [b, c], samp)
    mz, mx = out[b], out[c]
    qm.discard(b); qm.discard(c)
    if mz or mx:
        qm.run_circuit(D.CORR[(mz, mx)], [d])
    st = qm.states[a]
    rho = np.array(st.state)
    if list(st.keys) != [a, d]:
        rho = D.SWAP2 @ rho @ D.SWAP2
    return (mz, mx), rho


def test_swap_pauli_corrections_all_outcomes():
    seen = set()
    for samp in (0.1, 0.35, 0.6, 0.85):
        for k1 in BELL:
            for k2 in BELL:
                r1 = np.outer(BELL[k1], BELL[k1].conj())
                r2 = np.outer(BELL[k2], BELL[k2].conj())
                outc, rho = swap_once(r1, r2, samp)
                seen.add(outc)
                # Bell states compose as Pauli errors: Phi+ (I), Phi- (Z), Psi+ (X), Psi- (XZ)
                err = {"phi+": (0, 0), "phi-": (1, 0), "psi+": (0, 1), "psi-": (1, 1)}
                ez = err[k1][0] ^ err[k2][0]
                ex = err[k1][1] ^ err[k2][1]
                target = {v: k for k, v in err.items()}[(ez, ex)]
                assert abs(bell_coeffs(rho)[target] - 1) < 1e-12, (k1, k2, samp, bell_coeffs(rho))
    assert seen == {(0, 0), (0, 1), (1, 0), (1, 1)}


def test_swap_werner_composition():
    for F1, F2 in ((0.9, 0.8), (0.97, 0.95), (0.7, 0.99)):
        for samp in (0.1, 0.4, 0.6, 0.9):
            _, rho = swap_once(D.werner(F1), D.werner(F2), samp, lam=0.995)
            w = 0.995 * Q.werner_w(F1) * Q.werner_w(F2)
            assert np.allclose(rho, D.werner((1 + 3 * w) / 4), atol=1e-12)


# ---------------------------------------------------------------------------
# 5. whole delivery attempts, Mode F, vs the paper's closed forms
# ---------------------------------------------------------------------------
CONFIGS = [(5, 3), (6, 3), (5, 4), (6, 4)]
CHANNELS = [("half", 0.1), ("deph2", 0.2), ("depol2", 0.3)]


def test_modeF_clean_profiles_match_closed_form():
    worst = 0.0
    for N, r in CONFIGS:
        T = Q.Topology.grid(N)
        reps = B.grid_paths_by_composition(N, 3)
        for mem, T2 in CHANNELS:
            ph = B.physics(r, mem, T2)
            prof = Q.clean_profile(T, ph, 0, N * N - 1, extra=3)
            best = {}
            for (h, b), p in reps.items():
                L, _ = B.path_links(T, p)
                res = D.run_path([ph.clean_F(x) for x in L], L, D.mode_F(rounds=r, memory=mem, T2=T2), seed=h * 31 + b)
                ref = B.closed_form_path(ph, L)
                worst = max(worst, abs(res.F - ref))
                best[h] = max(best.get(h, 0.0), res.F)
            for h, Fh in prof.items():
                worst = max(worst, abs(best[h] - Fh))
                k = h - (N - 1)
                worst = max(worst, abs(best[h] - B.closed_form_fstar(N, k, r, mem, T2)))
    assert worst < 1e-10, worst


def test_modeF_noisy_paths_match_qunet():
    """Realisations with noisy links (E1 seeds), paths of several routers."""
    T = Q.Topology.grid(5)
    ph = Q.Physics(rounds=3)
    worst, n_cmp, n_fail, n_junk = 0.0, 0, 0, 0
    for t in range(40):
        R = Q.Realisation(T, ph, B.SEED0 + t)
        rt = Q.Router(R, 0, 24)
        paths = {tuple(rt.route(p)) for p in ("shortest", "full", "dijkstra", "highest_fidelity")}
        paths.add(tuple(rt.oracle()[0]))
        for p in paths:
            if len(p) < 2 or p[-1] != 24:
                continue
            L, E = B.path_links(T, list(p))
            h = len(L)
            sf = [bool(R.swap_u[j] < ph.p_swap_fail) for j in range(h - 1)]
            res = D.run_path([R.raw[e] for e in E], L, D.mode_F(rounds=3), seed=t, swap_fail=sf)
            Fq, okq, _ = R.evaluate(list(p))
            assert res.success == okq, (t, p, res.F, Fq)
            if res.swap_failed:
                n_fail += 1
            else:
                # includes paths with faulty links (raw F <= 1/2, not purified, as in qunet_tqe)
                worst = max(worst, abs(res.F - Fq))
                n_junk += int(min(R.raw[e] for e in E) <= 0.5)
            n_cmp += 1
    assert n_cmp > 60 and n_junk > 0
    assert worst < 1e-10, worst


# ---------------------------------------------------------------------------
# 6. timing and storage accounting, Mode R
# ---------------------------------------------------------------------------
def test_modeR_storage_accounting_exact():
    """r = 0 and depolarizing memories keep every pair Werner, so
    w_e2e = lambda^(h-1) prod w_i exp(- sum_q t_q / T2), where t_q is the total
    storage time of qubit q, read from the DES storage log."""
    lens = [1.0, math.sqrt(2), 1.0, math.sqrt(2), 1.0]
    Fs = [0.95, 0.9, 0.93, 0.88, 0.96]
    for seed in range(5):
        P = D.mode_R(1000.0, 0.002, rounds=0, memory="depol2", T2=0.05)
        res = D.run_path(Fs, lens, P, seed=seed, log_storage=True)
        assert len(res.storage) == 2 * len(lens)
        tot = 0.0
        for k, (t0, t1) in res.storage.items():
            assert t1 is not None and t1 >= t0 - 1e-15
            tot += t1 - t0
        w = 0.995 ** (len(lens) - 1) * np.prod([Q.werner_w(f) for f in Fs]) * math.exp(-tot / 0.05)
        assert abs(res.F - (1 + 3 * w) / 4) < 1e-9, (res.F, (1 + 3 * w) / 4)
        assert res.t_delivery >= max(res.t_link_ready)
        assert abs(res.t_used - (res.t_delivery + 0.002 * sum(lens))) < 1e-12


def test_modeR_without_decoherence_equals_modeF_without_storage():
    T = Q.Topology.grid(5)
    reps = B.grid_paths_by_composition(5, 1)
    ph = Q.Physics(rounds=3, kappa=0.0)
    for (h, b), p in reps.items():
        L, _ = B.path_links(T, p)
        Fl = [ph.clean_F(x) for x in L]
        ref = D.run_path(Fl, L, D.mode_F(rounds=3, kappa=0.0), seed=1).F
        for s in range(3):
            got = D.run_path(Fl, L, D.mode_R(1000.0, 0.0, rounds=3, T2=math.inf), seed=s).F
            assert abs(got - ref) < 1e-12


def test_raw_pair_cost_matches_recurrence():
    """Mean raw pairs per purified link pair (Mode F) = N_r = 2 N_{r-1} / P_r."""
    ph = Q.Physics(rounds=3)
    F = ph.clean_F(math.sqrt(2))
    Nr, _ = Q.purify_cost(F, 3)
    n = 400
    counts = [D.run_path([F], [math.sqrt(2)], D.mode_F(rounds=3), seed=s).raw_pairs[0] for s in range(n)]
    m, se = float(np.mean(counts)), float(np.std(counts) / math.sqrt(n))
    assert abs(m - Nr) < 4 * se, (m, Nr, se)


def test_pessimistic_reference_copy():
    p = os.path.join(B.RESULTS_TQE, "design_law.json")
    if not os.path.exists(p):
        print("      (skipped: ../results_tqe/design_law.json not found)")
        return
    ref = json.load(open(p))["latency_all_links_depol2"]
    for key in ("1000|5", "100|5"):
        R_, N = key.split("|")
        best = min(B.pessimistic_min_T2_depol2(int(N), r, float(R_)) for r in range(1, 7))
        assert abs(best - ref[key]["best_T2"]) / ref[key]["best_T2"] < 1e-9, (key, best, ref[key])


# ---------------------------------------------------------------------------
TESTS = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]

if __name__ == "__main__":
    t0 = time.time()
    bad = 0
    for f in TESTS:
        t = time.time()
        try:
            f()
            print(f"PASS  {f.__name__}  ({time.time() - t:.1f}s)", flush=True)
        except Exception as e:  # noqa: BLE001
            bad += 1
            print(f"FAIL  {f.__name__}: {type(e).__name__}: {e}", flush=True)
    print(f"{len(TESTS) - bad}/{len(TESTS)} passed in {time.time() - t0:.0f}s")
    sys.exit(1 if bad else 0)
