"""
run_des.py -- the discrete-event validation study (SeQUeNCe density-matrix DES).

Blocks (run all by default, or pick with --only PF,PR,T,S1,S2):

  PF  Mode F clean-path profiles F*(k) for the four configurations, three memory
      channels, three purification schemes (BBPSSW with twirl = the paper's
      model; DEJMPS; bilateral CNOT without twirl).  Regime call s from the DES.
      Plus a gate-noise sub-study (eps = 0.002, 0.005, 0.01; BBPSSW vs DEJMPS).
      Mode F is an implementation cross-check of the closed form, not a test of
      the physics (it uses the paper's assumptions).
  PR  Mode R clean-path profiles at the paper's coherence times with fast
      generation (R = 1e5, 1e4 Hz), kappa kept: per-memory storage resolved,
      everything else as in the paper.  Compared with the closed form at the
      same R (the paper's optimistic accounting t = kappa D + C_r / R).
  T   Mode R design-chart check: on 5x5, the smallest T2 at which the DES puts
      s >= 0 and s >= 1, for r = 3, 4, R = 1 kHz, 100 Hz, three channels,
      with and without kappa.  Compared with design_law.py (optimistic) and,
      for depol2, the pessimistic all-links-wait bound.
  S1  Mode F routing study on 5x5, r = 3 (the paper's E1 cell, same seeds):
      every router's path evaluated by the DES, trial by trial.
  S2  Mode R routing study on 5x5, r = 3, R = 1 kHz, kappa = 2 ms/km, three
      memory channels, at two coherence times taken from T (one where the DES
      puts s = 0, one where it puts s = 1; each verified by a DES clean-path
      profile before use): success, routing headroom, latency, delivery rate.

Usage:
  python run_des.py --quick              # ~10-20 min on 8 cores: every block, small n
  python run_des.py                      # full run (see README for timing)
  python run_des.py --only S2 --workers 16
Outputs: results_des/<block>.json (+ results_des/run_log.txt)
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
import zlib
import multiprocessing as mp
from multiprocessing import cpu_count
from typing import Dict, List

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results_des")

CONFIGS = [(5, 3), (6, 3), (5, 4), (6, 4)]                 # (N, rounds)
CHANNELS = [("half", 0.1), ("deph2", 0.2), ("depol2", 0.3)]  # paper's T2 per channel
SCHEMES = ["bbpssw_twirl", "dejmps", "bbpssw"]
EXTRA = {"bbpssw_twirl": 4, "dejmps": 8, "bbpssw": 2}       # enumeration depth h* + EXTRA
POLICIES = ["shortest", "baseline", "full", "dijkstra", "oracle"]
KAPPA = 0.002


def seed_of(*parts) -> int:
    """Platform-independent deterministic seed from a tuple of labels/ints."""
    s = "|".join(str(p) for p in parts).encode()
    return int(np.random.SeedSequence([zlib.crc32(s), len(s)]).generate_state(1)[0])


# ---------------------------------------------------------------------------
# worker functions (top level: required by Windows 'spawn' multiprocessing)
# ---------------------------------------------------------------------------
def make_params(pk: dict):
    import des_sim as D
    kw = dict(pk)
    mode = kw.pop("mode")
    if mode == "F":
        return D.mode_F(**kw)
    R = kw.pop("gen_rate")
    kap = kw.pop("kappa")
    return D.mode_R(R, kap, **kw)


def job_paths(job):
    """Deliveries on one fixed path, one per seed."""
    import des_sim as D
    P = make_params(job["P"])
    F, td, raw = [], [], []
    for s in job["seeds"]:
        r = D.run_path(job["link_F"], job["link_len"], P, seed=s)
        F.append(r.F)
        td.append(r.t_used)
        raw.append(int(sum(r.raw_pairs)))
    return job["key"], job["i0"], F, td, raw


def job_trials(job):
    """Routing trials: sample the network (qunet_tqe, E1 seeds), route with each
    policy (qunet_tqe routers), evaluate every distinct path in the DES."""
    import des_sim as D
    import qunet_bridge as B
    Q = B.Q
    N, rounds = job["N"], job["rounds"]
    T = Q.Topology.grid(N)
    # the paper's model at the same parameters (incl. gen_rate in Mode R): used by the routers
    # and the oracle, and as the model reference value F_model / ok_model
    gen = job["P"].get("gen_rate", math.inf) if job["P"]["mode"] == "R" else math.inf
    ph = Q.Physics(rounds=rounds, memory=job["memory"], T2=job["T2"], kappa=job["P"]["kappa"], gen_rate=gen)
    P = make_params(job["P"])
    dst = N * N - 1
    pols = job["policies"]
    rows = []
    for t in range(job["t0"], job["t1"]):
        R = Q.Realisation(T, ph, B.SEED0 + t)
        rt = Q.Router(R, 0, dst)
        paths = {}
        for p in pols:
            path = rt.oracle()[0] if p == "oracle" else rt.route(p)
            paths[p] = tuple(path) if (path and path[-1] == dst) else None
        cache = {}
        des_seed = seed_of(job["tag"], t)
        for p, path in paths.items():
            if path is None or path in cache:
                continue
            L, E = B.path_links(T, list(path))
            sf = [bool(R.swap_u[j] < ph.p_swap_fail) for j in range(len(L) - 1)]
            r = D.run_path([R.raw[e] for e in E], L, P, seed=des_seed, swap_fail=sf)
            Fq, okq, _ = R.evaluate(list(path))
            # noisy links carry a value drawn from U(lo, hi); clean links equal clean_F(d) exactly
            n_noisy = sum(1 for e, l in zip(E, L) if R.raw[e] != ph.clean_F(l))
            cache[path] = (r.F, bool(r.success), r.t_used, int(sum(r.raw_pairs)), len(L),
                           float(Fq), bool(okq), bool(r.swap_failed), n_noisy, r.t_delivery)
        row = {}
        for p, path in paths.items():
            row[p] = cache[path] if path is not None else (0.0, False, float("nan"), 0, 0, 0.0, False, False, -1,
                                                           float("nan"))
        rows.append((t, row))
    return job["tag"], rows


def pool_map(fn, jobs, workers, label):
    t0 = time.time()
    out = []
    n = len(jobs)
    if workers <= 1:
        for i, j in enumerate(jobs):
            out.append(fn(j))
            _progress(label, i + 1, n, t0)
        return out
    # 'spawn' everywhere: the Windows default, so Linux/macOS runs exercise the same code path
    with mp.get_context("spawn").Pool(workers) as pool:
        for i, r in enumerate(pool.imap_unordered(fn, jobs, chunksize=1)):
            out.append(r)
            _progress(label, i + 1, n, t0)
    return out


def _progress(label, i, n, t0):
    if i == n or i % max(1, n // 20) == 0:
        el = time.time() - t0
        eta = el / i * (n - i)
        print(f"  [{label}] {i}/{n} jobs  {el:.0f}s elapsed, ~{eta:.0f}s left", flush=True)


def split_seeds(tag: str, n: int, chunk: int):
    seeds = [seed_of(tag, i) for i in range(n)]
    return [(i0, seeds[i0:i0 + chunk]) for i0 in range(0, n, chunk)]


def run_path_block(specs: List[dict], n: int, chunk: int, workers: int, label: str) -> Dict[str, dict]:
    """specs: {key, link_F, link_len, P}; returns key -> {F, t, raw} arrays (length n)."""
    jobs = []
    for sp in specs:
        for i0, seeds in split_seeds(sp["key"], n, chunk):
            jobs.append({"key": sp["key"], "i0": i0, "seeds": seeds, "link_F": sp["link_F"],
                         "link_len": sp["link_len"], "P": sp["P"]})
    res = {sp["key"]: {"F": [None] * n, "t": [None] * n, "raw": [None] * n} for sp in specs}
    for key, i0, F, td, raw in pool_map(job_paths, jobs, workers, label):
        res[key]["F"][i0:i0 + len(F)] = F
        res[key]["t"][i0:i0 + len(F)] = td
        res[key]["raw"][i0:i0 + len(F)] = raw
    return res


def summarise_F(F: List[float]) -> dict:
    a = np.asarray(F, float)
    n = len(a)
    return {"n": n, "mean": float(a.mean()), "sd": float(a.std(ddof=1)) if n > 1 else 0.0,
            "se": float(a.std(ddof=1) / math.sqrt(n)) if n > 1 else 0.0,
            "min": float(a.min()), "max": float(a.max()),
            "p_pass": float(np.mean(a >= 0.8))}


# ---------------------------------------------------------------------------
# blocks
# ---------------------------------------------------------------------------
GATE_EPS = [0.002, 0.005, 0.01]     # gate-noise sub-study (half memory, T2 = 0.1 s), as in the paper's E8


def block_PF(args):
    import qunet_bridge as B
    Q = B.Q
    specs, meta = [], {}
    grid = [(mem, T2, sch, 0.0) for mem, T2 in CHANNELS for sch in SCHEMES]
    grid += [("half", 0.1, sch, e) for e in GATE_EPS for sch in ("bbpssw_twirl", "dejmps")]
    for N, r in CONFIGS:
        T = Q.Topology.grid(N)
        reps_all = {e: B.grid_paths_by_composition(N, e) for e in sorted(set(EXTRA.values()))}
        for mem, T2, sch, eps in grid:
            ph = Q.Physics(rounds=r, memory=mem, T2=T2, eps_gate=eps)
            for (h, b), p in sorted(reps_all[EXTRA[sch]].items()):
                L, _ = B.path_links(T, p)
                key = f"PF|g{N}_r{r}|{mem}|{sch}|e{eps:g}|h{h}_b{b}"
                specs.append({"key": key, "link_F": [ph.clean_F(x) for x in L], "link_len": L,
                              "P": {"mode": "F", "rounds": r, "memory": mem, "T2": T2, "scheme": sch,
                                    "eps_gate": eps}})
                meta[key] = {"N": N, "r": r, "memory": mem, "T2": T2, "scheme": sch, "eps": eps, "h": h, "b": b,
                             "h_star": N - 1, "path": p,
                             "closed_form": B.closed_form_path(ph, L) if sch == "bbpssw_twirl" else None}
    n = 3   # Mode F is deterministic given success; 3 seeds confirm it
    res = run_path_block(specs, n, n, args.workers, "PF")
    out = {k: {**meta[k], **summarise_F(v["F"])} for k, v in res.items()}
    # closed-form references per (config, channel, eps)
    ref = {}
    for N, r in CONFIGS:
        for mem, T2, sch, eps in grid:
            if sch != "bbpssw_twirl":
                continue
            ph = Q.Physics(rounds=r, memory=mem, T2=T2, eps_gate=eps)
            s = B.DL.slack(ph, N)
            wc, wd = B.DL.wcd(ph)
            ref[f"g{N}_r{r}|{mem}|e{eps:g}"] = {"slack": s, "Fstar": [B.DL.fstar(ph, N, k, wc, wd)
                                                                       for k in range(0, max(s, 0) + 2)]}
    return {"paths": out, "closed_form": ref, "n_per_path": n, "gate_eps": GATE_EPS}


def block_PR(args):
    import qunet_bridge as B
    Q = B.Q
    n = args.n_pr
    specs, meta = [], {}
    rates = [1e5, 1e4]
    for N, r in CONFIGS:
        T = Q.Topology.grid(N)
        reps = B.grid_paths_by_composition(N, 3)
        for mem, T2 in CHANNELS:
            for R in rates:
                ph = B.physics(r, mem, T2, gen_rate=R)
                s_cf = B.closed_form_slack(N, r, mem, T2, gen_rate=R)
                hmax = (N - 1) + max(s_cf, 0) + 1
                for (h, b), p in sorted(reps.items()):
                    if h > hmax:
                        continue
                    L, _ = B.path_links(T, p)
                    key = f"PR|g{N}_r{r}|{mem}|R{int(R)}|h{h}_b{b}"
                    specs.append({"key": key, "link_F": [ph.clean_F(x) for x in L], "link_len": L,
                                  "P": {"mode": "R", "gen_rate": R, "kappa": KAPPA, "rounds": r,
                                        "memory": mem, "T2": T2}})
                    meta[key] = {"N": N, "r": r, "memory": mem, "T2": T2, "R": R, "h": h, "b": b,
                                 "h_star": N - 1, "closed_form": B.closed_form_path(ph, L),
                                 "closed_form_slack": s_cf}
    res = run_path_block(specs, n, max(5, n // 8), args.workers, "PR")
    out = {k: {**meta[k], **summarise_F(v["F"]), "t_mean": float(np.mean(v["t"]))} for k, v in res.items()}
    return {"paths": out, "n_per_path": n, "rates": rates}


def t_grid(center: float):
    return [center * 2 ** (i / 2) for i in range(-2, 9)]       # 0.5x .. 16x, sqrt(2) steps


def block_T(args):
    import qunet_bridge as B
    Q = B.Q
    N = 5
    T = Q.Topology.grid(N)
    reps = B.grid_paths_by_composition(N, 1)
    n = args.n_t
    specs, curves = [], {}
    for r in (3, 4):
        for R in (1000.0, 100.0):
            for mem in ("half", "deph2", "depol2"):
                for kap in (KAPPA, 0.0):
                    for target in (0, 1):
                        cf = B.closed_form_min_T2(N, r, target, mem, R, kappa=kap)
                        pess = B.pessimistic_min_T2_depol2(N, r, R, target=target, kappa=kap) \
                            if (mem == "depol2" and target == 1) else None
                        h = N - 1 + target
                        classes = [(hh, b) for (hh, b) in reps if hh == h]
                        ckey = f"T|r{r}|R{int(R)}|{mem}|k{kap:g}|s{target}"
                        grid = t_grid(cf if math.isfinite(cf) else 1.0)
                        curves[ckey] = {"r": r, "R": R, "memory": mem, "kappa": kap, "target": target,
                                        "closed_form_T2": cf, "pessimistic_T2": pess, "grid": grid,
                                        "classes": classes, "points": {}}
                        for gi, T2 in enumerate(grid):
                            ph = B.physics(r, mem, T2, kappa=kap, gen_rate=R)
                            for (hh, b) in classes:
                                L, _ = B.path_links(T, reps[(hh, b)])
                                key = f"{ckey}|g{gi}|h{hh}_b{b}"
                                specs.append({"key": key, "link_F": [ph.clean_F(x) for x in L], "link_len": L,
                                              "P": {"mode": "R", "gen_rate": R, "kappa": kap, "rounds": r,
                                                    "memory": mem, "T2": T2}})
                                curves[ckey]["points"][key] = {"gi": gi, "T2": T2, "h": hh, "b": b,
                                                               "closed_form": B.closed_form_path(ph, L)}
    res = run_path_block(specs, n, max(5, n // 6), args.workers, "T")
    for ckey, c in curves.items():
        for key, pt in c["points"].items():
            pt.update(summarise_F(res[key]["F"]))
            pt["F"] = [float(x) for x in res[key]["F"]]           # kept for bootstrap CIs
            pt["t_mean"] = float(np.mean(res[key]["t"]))
    return {"curves": curves, "n_per_point": n, "N": N}


def block_T6(args, rounds=(4, 5), N=6, rates=(1000.0, 100.0), tag="T6"):
    """Design-chart check on further grid sizes: 6x6 (BBPSSW r = 4, 5) and 7x7 (r = 5), per-qubit channels, s >= 1."""
    import qunet_bridge as B
    Q = B.Q
    T = Q.Topology.grid(N)
    reps = B.grid_paths_by_composition(N, 1)
    n = args.n_t
    specs, curves = [], {}
    for r in rounds:
        for R in rates:
            for mem in ("deph2", "depol2"):
                kap, target = KAPPA, 1
                cf = B.closed_form_min_T2(N, r, target, mem, R, kappa=kap)
                pess = B.pessimistic_min_T2_depol2(N, r, R, target=target, kappa=kap) if mem == "depol2" else None
                h = N - 1 + target
                classes = [(hh, b) for (hh, b) in reps if hh == h]
                ckey = f"{tag}|r{r}|R{int(R)}|{mem}|k{kap:g}|s{target}"
                grid = t_grid(cf if math.isfinite(cf) else 1.0)
                curves[ckey] = {"r": r, "R": R, "memory": mem, "kappa": kap, "target": target, "N": N,
                                "closed_form_T2": cf, "pessimistic_T2": pess, "grid": grid,
                                "classes": classes, "points": {}}
                for gi, T2 in enumerate(grid):
                    ph = B.physics(r, mem, T2, kappa=kap, gen_rate=R)
                    for (hh, b) in classes:
                        L, _ = B.path_links(T, reps[(hh, b)])
                        key = f"{ckey}|g{gi}|h{hh}_b{b}"
                        specs.append({"key": key, "link_F": [ph.clean_F(x) for x in L], "link_len": L,
                                      "P": {"mode": "R", "gen_rate": R, "kappa": kap, "rounds": r,
                                            "memory": mem, "T2": T2}})
                        curves[ckey]["points"][key] = {"gi": gi, "T2": T2, "h": hh, "b": b,
                                                       "closed_form": B.closed_form_path(ph, L)}
    res = run_path_block(specs, n, max(5, n // 6), args.workers, tag)
    for ckey, c in curves.items():
        for key, pt in c["points"].items():
            pt.update(summarise_F(res[key]["F"]))
            pt["F"] = [float(x) for x in res[key]["F"]]
            pt["t_mean"] = float(np.mean(res[key]["t"]))
    return {"curves": curves, "n_per_point": n, "N": N}


def block_TD(args):
    """Design-chart check with DEJMPS (Mode R): smallest T2 with s >= 1 on 5x5, r = 2, 3."""
    import qunet_bridge as B
    Q = B.Q
    import design_law as DLM
    import qunet_bd as QB
    N = 5
    T = Q.Topology.grid(N)
    reps = B.grid_paths_by_composition(N, 1)
    n = args.n_t
    specs, curves = [], {}
    for r in (2, 3):
        for R in (1000.0, 100.0):
            for mem in ("half", "deph2", "depol2"):
                kap, target = KAPPA, 1
                cf = DLM.min_T2_s(N, r, target, "dejmps", memory=mem, gen_rate=R, kappa=kap)
                h = N - 1 + target
                classes = [(hh, b) for (hh, b) in reps if hh == h]
                ckey = f"TD|r{r}|R{int(R)}|{mem}|k{kap:g}|s{target}"
                grid = t_grid(cf if math.isfinite(cf) else 1.0)
                curves[ckey] = {"r": r, "R": R, "memory": mem, "kappa": kap, "target": target, "scheme": "dejmps",
                                "closed_form_T2": cf, "pessimistic_T2": None, "grid": grid,
                                "classes": classes, "points": {}}
                for gi, T2 in enumerate(grid):
                    ph = QB.PhysicsBD(scheme="dejmps", rounds=r, memory=mem, T2=T2, kappa=kap, gen_rate=R)
                    for (hh, b) in classes:
                        L, _ = B.path_links(T, reps[(hh, b)])
                        S = [0.0, 0.0, 0.0]
                        for x in L:
                            c = ph.link(ph.clean_F(x))[0]
                            S = [S[j] + math.log(c[j]) for j in range(3)]
                        key = f"{ckey}|g{gi}|h{hh}_b{b}"
                        specs.append({"key": key, "link_F": [ph.clean_F(x) for x in L], "link_len": L,
                                      "P": {"mode": "R", "gen_rate": R, "kappa": kap, "rounds": r,
                                            "memory": mem, "T2": T2, "scheme": "dejmps"}})
                        curves[ckey]["points"][key] = {"gi": gi, "T2": T2, "h": hh, "b": b,
                                                       "closed_form": ph.e2e_vec(S, len(L), sum(L))}
    res = run_path_block(specs, n, max(5, n // 6), args.workers, "TD")
    for ckey, c in curves.items():
        for key, pt in c["points"].items():
            pt.update(summarise_F(res[key]["F"]))
            pt["F"] = [float(x) for x in res[key]["F"]]
            pt["t_mean"] = float(np.mean(res[key]["t"]))
    return {"curves": curves, "n_per_point": n, "N": N}


def run_trials(tag, N, rounds, memory, T2, P, n, workers, chunk):
    jobs = [{"tag": tag, "N": N, "rounds": rounds, "memory": memory, "T2": T2, "P": P,
             "policies": POLICIES, "t0": a, "t1": min(a + chunk, n)} for a in range(0, n, chunk)]
    rows = {}
    for tg, rr in pool_map(job_trials, jobs, workers, tag):
        for t, row in rr:
            rows[t] = row
    # t_est = end-to-end pair established (corrections known); t_used = t_est + kappa * D (pair used)
    cols = ["F_des", "ok_des", "t_used", "raw", "hops", "F_model", "ok_model", "swap_failed", "n_noisy", "t_est"]
    arr = {p: {c: [rows[t][p][i] for t in range(n)] for i, c in enumerate(cols)} for p in POLICIES}
    return arr


def block_S1(args):
    n = args.n_s1
    P = {"mode": "F", "rounds": 3, "memory": "half", "T2": 0.1}
    arr = run_trials("S1", 5, 3, "half", 0.1, {**P, "kappa": KAPPA}, n, args.workers, max(20, n // (args.workers * 6)))
    return {"policies": POLICIES, "n": n, "params": P, "trials": arr, "seed0": 20260922}


def crossing_T2(points: List[dict]):
    """Smallest grid interval on which the mean F crosses F_TH upward; log-linear interpolation."""
    pts = sorted(points, key=lambda p: p["T2"])
    xs = [math.log(p["T2"]) for p in pts]
    ys = [p["mean"] - 0.8 for p in pts]
    if ys[0] >= 0:
        return {"T2": None, "bracket": [None, pts[0]["T2"]], "status": "below_grid"}
    for i in range(len(pts) - 1):
        if ys[i] < 0 <= ys[i + 1]:
            x = xs[i] + (xs[i + 1] - xs[i]) * (-ys[i]) / (ys[i + 1] - ys[i])
            return {"T2": math.exp(x), "bracket": [pts[i]["T2"], pts[i + 1]["T2"]], "status": "ok"}
    return {"T2": None, "bracket": [pts[-1]["T2"], None], "status": "above_grid"}


def best_class_points(curve: dict) -> List[dict]:
    """Per grid point, the class (h, b) with the highest mean F (here one class per h)."""
    best = {}
    for pt in curve["points"].values():
        g = pt["gi"]
        if g not in best or pt["mean"] > best[g]["mean"]:
            best[g] = pt
    return [best[g] for g in sorted(best)]


def des_profile(N, r, mem, T2, R, kap, kmax, n, workers, tag):
    """DES Mode R clean-path profile on the N x N grid: best-class mean F*(k), k = 0..kmax."""
    import qunet_bridge as B
    Q = B.Q
    T = Q.Topology.grid(N)
    ph = B.physics(r, mem, T2, kappa=kap, gen_rate=R)
    reps = B.grid_paths_by_composition(N, kmax)
    specs = []
    for (h, b), p in sorted(reps.items()):
        L, _ = B.path_links(T, p)
        specs.append({"key": f"{tag}|h{h}_b{b}", "link_F": [ph.clean_F(x) for x in L], "link_len": L,
                      "P": {"mode": "R", "gen_rate": R, "kappa": kap, "rounds": r, "memory": mem, "T2": T2}})
    res = run_path_block(specs, n, max(5, n // 6), workers, tag)
    prof = {}
    for sp in specs:
        h = int(sp["key"].split("|h")[-1].split("_")[0])
        k = h - (N - 1)
        st = summarise_F(res[sp["key"]]["F"])
        if k not in prof or st["mean"] > prof[k]["mean"]:
            prof[k] = {**st, "key": sp["key"]}
    return prof


def block_S2(args):
    path = os.path.join(OUT, "T.json")
    if not os.path.exists(path):
        raise SystemExit("S2 needs results_des/T.json (run block T first)")
    Tres = json.load(open(path, encoding="utf-8"))
    out = {"policies": POLICIES, "n": args.n_s2, "variants": {}}
    R, kap = 1000.0, KAPPA
    for mem in ("half", "deph2", "depol2"):
        c0 = Tres["curves"][f"T|r3|R{int(R)}|{mem}|k{kap:g}|s0"]
        c1 = Tres["curves"][f"T|r3|R{int(R)}|{mem}|k{kap:g}|s1"]
        x0, x1 = crossing_T2(best_class_points(c0)), crossing_T2(best_class_points(c1))
        if x0["status"] != "ok" or x1["status"] != "ok":
            print(f"  [S2] skipping {mem}: DES threshold not bracketed ({x0['status']}, {x1['status']})")
            continue
        vkey = f"{mem}|k{kap:g}|R{int(R)}"
        out["variants"][vkey] = {"memory": mem, "kappa": kap, "R": R, "T2_des_s0": x0["T2"],
                                 "T2_des_s1": x1["T2"], "runs": {}}
        # candidate T2 for each DES regime; each is VERIFIED with a DES clean-path profile
        g = math.sqrt(x0["T2"] * x1["T2"])
        cands = {"s0": [g, g * (x1["T2"] / g) ** -0.5, g * (x1["T2"] / g) ** 0.5],  # F*(0) >= F_th > F*(1)
                 "s1": [f * x1["T2"] for f in (1.5, 1.3, 1.15)]}                 # F*(1) >= F_th > F*(2)
        for label, lst in cands.items():
            want = 0 if label == "s0" else 1
            chosen, prof, ok = None, None, False
            for T2 in lst:
                prof = des_profile(5, 3, mem, T2, R, kap, want + 1, args.n_t, args.workers,
                                   f"S2v|{vkey}|{label}|{T2:.6g}")
                # the two deciding means must be > 3 SE from F_th (a call at 0.8005 +- 0.003 does not count)
                ok = all(prof[k]["mean"] >= 0.8 for k in range(want + 1)) and prof[want + 1]["mean"] < 0.8 \
                    and all(abs(prof[k]["mean"] - 0.8) > 3 * prof[k]["se"] for k in (want, want + 1))
                chosen = T2
                if ok:
                    break
            if not ok:
                print(f"  [S2] WARNING: {vkey} {label}: no candidate T2 verified the DES regime s = {want}")
            P = {"mode": "R", "gen_rate": R, "kappa": kap, "rounds": 3, "memory": mem, "T2": chosen}
            tag = f"S2|{vkey}|{label}"
            arr = run_trials(tag, 5, 3, mem, chosen, P, args.n_s2, args.workers, max(10, args.n_s2 // (args.workers * 6)))
            out["variants"][vkey]["runs"][label] = {
                "T2": chosen, "params": P, "trials": arr, "label_verified": ok,
                "profile": {str(k): {kk: v[kk] for kk in ("mean", "se", "p_pass", "n")} for k, v in prof.items()}}
    return out


def block_S3(args):
    """Softness sweep: the routing headroom of the S2 setting (5x5, r = 3, R = 1 kHz,
    kappa = 2 ms/km, per-qubit depolarizing memories) at six coherence times across the DES s = 0 -> 1
    boundary, each with its DES clean-path profile, so that the residual headroom can be reported as
    a function of the deciding margin F*(1) - F_th rather than at one point."""
    path = os.path.join(OUT, "T.json")
    Tres = json.load(open(path, encoding="utf-8"))
    R, kap, mem = 1000.0, KAPPA, "depol2"
    c1 = Tres["curves"][f"T|r3|R{int(R)}|{mem}|k{kap:g}|s1"]
    x1 = crossing_T2(best_class_points(c1))["T2"]
    out = {"policies": POLICIES, "n": args.n_s3, "memory": mem, "R": R, "kappa": kap, "T2_des_s1": x1, "runs": {}}
    for f in (0.55, 0.70, 0.85, 0.95, 1.05, 1.20):
        T2 = f * x1
        prof = des_profile(5, 3, mem, T2, R, kap, 2, args.n_t, args.workers, f"S3v|{T2:.6g}")
        P = {"mode": "R", "gen_rate": R, "kappa": kap, "rounds": 3, "memory": mem, "T2": T2}
        tag = f"S3|{mem}|{f:g}"
        arr = run_trials(tag, 5, 3, mem, T2, P, args.n_s3, args.workers, max(10, args.n_s3 // (args.workers * 6)))
        out["runs"][f"{f:g}"] = {"T2": T2, "factor": f, "params": P, "trials": arr,
                                 "profile": {str(k): {kk: v[kk] for kk in ("mean", "se", "p_pass", "n")}
                                             for k, v in prof.items()}}
        print(f"  [S3] T2 = {T2:.4g} s done", flush=True)
    return out


BLOCKS = {"PF": block_PF, "PR": block_PR, "T": block_T, "S1": block_S1, "S2": block_S2, "S3": block_S3, "TD": block_TD, "T6": block_T6,
          # T6 in two halves (one per round count), merged by analyze.py, for runs with a wall-clock limit
          "T6a": lambda a: block_T6(a, (4,)), "T6b": lambda a: block_T6(a, (5,)),
          # 7x7, r = 5, in two halves by generation rate (merged by analyze.py)
          "T7a": lambda a: block_T6(a, (5,), N=7, rates=(1000.0,), tag="T7"),
          "T7b": lambda a: block_T6(a, (5,), N=7, rates=(100.0,), tag="T7")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # Windows' multiprocessing cannot wait on more than 61 worker processes
    ap.add_argument("--workers", type=int, default=max(1, min(61, cpu_count() - 1)))
    ap.add_argument("--quick", action="store_true", help="small n everywhere (pipeline check)")
    ap.add_argument("--only", default="PF,PR,T,S1,S2")
    a = ap.parse_args()
    a.workers = max(1, min(a.workers, 61))
    if a.quick:
        a.n_pr, a.n_t, a.n_s1, a.n_s2, a.n_s3 = 24, 24, 300, 150, 100
    else:
        a.n_pr, a.n_t, a.n_s1, a.n_s2, a.n_s3 = 400, 400, 5000, 3000, 2000
    os.makedirs(OUT, exist_ok=True)
    import sequence  # noqa: F401  (fail early if SeQUeNCe is missing)
    from importlib.metadata import version
    env = {"python": sys.version.split()[0], "platform": platform.platform(),
           "sequence": version("sequence"), "numpy": np.__version__, "workers": a.workers, "quick": a.quick}
    print("environment:", env, flush=True)
    log = open(os.path.join(OUT, "run_log.txt"), "a", encoding="utf-8")
    for b in [x.strip() for x in a.only.split(",") if x.strip()]:
        t0 = time.time()
        print(f"=== block {b} ===", flush=True)
        res = BLOCKS[b](a)
        res["meta"] = {**env, "block": b, "seconds": time.time() - t0,
                       "n": {"PR": a.n_pr, "T": a.n_t, "S1": a.n_s1, "S2": a.n_s2, "S3": a.n_s3},
                       "finished": time.strftime("%Y-%m-%d %H:%M:%S")}
        with open(os.path.join(OUT, f"{b}.json"), "w", encoding="utf-8") as f:
            json.dump(res, f, default=float)
        msg = f"{time.strftime('%Y-%m-%d %H:%M:%S')}  block {b} done in {time.time() - t0:.0f}s  quick={a.quick}"
        print(msg, flush=True)
        log.write(msg + "\n")
        log.flush()
    log.close()


if __name__ == "__main__":
    main()
