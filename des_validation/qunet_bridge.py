"""
qunet_bridge.py -- the only place where the DES package touches the paper's code.

Used from the paper's code (../code):
  * qunet_tqe.Topology, Realisation, Router : network sampling and the routing
    algorithms under test (algorithms, not physics), and the reference
    closed-form numbers the DES is compared against;
  * qunet_tqe.Physics.clean_F               : the generated-pair fidelity
    (a model INPUT shared by both simulators, not a result);
  * design_law.wcd / fstar / slack / min_T2 : closed-form references.

Nothing here is used inside des_sim.py.
"""
from __future__ import annotations

import math
import os
import sys
from typing import Dict, List, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.normpath(os.path.join(HERE, "..", "code"))
RESULTS_TQE = os.path.normpath(os.path.join(HERE, "..", "results_tqe"))
if CODE not in sys.path:
    sys.path.insert(0, CODE)

import numpy as np  # noqa: E402

import qunet_tqe as Q  # noqa: E402
import design_law as DL  # noqa: E402

SQ2 = math.sqrt(2.0)
SEED0 = 20260922          # = run_study_tqe.SEED0 (common random numbers with block E1)


# ---------------------------------------------------------------------------
# grid paths
# ---------------------------------------------------------------------------
def grid_paths_by_composition(N: int, extra: int) -> Dict[Tuple[int, int], List[int]]:
    """One representative simple corner-to-corner path for every (h, b) that
    exists on the N x N king's-move grid, h = h*..h*+extra, b = number of
    cardinal (1 km) links.  On the all-clean network every path with the same
    (h, b) has the same link multiset, so one representative per class covers
    the best path of each length.  Exhaustive DFS with the admissible prune
    depth + hops_remaining <= h* + extra."""
    T = Q.Topology.grid(N)
    src, dst = 0, N * N - 1
    hop_rem = T.hops_to(dst)
    h0 = hop_rem[src]
    hmax = h0 + extra
    reps: Dict[Tuple[int, int], List[int]] = {}
    vis = [False] * T.n
    vis[src] = True
    path = [src]

    def is_card(u, v):
        return abs(T.elength(u, v) - 1.0) < 1e-9

    def dfs(u, h, b):
        if u == dst:
            reps.setdefault((h, b), list(path))
            return
        for v in T.adj[u]:
            if vis[v] or h + 1 + hop_rem[v] > hmax:
                continue
            vis[v] = True
            path.append(v)
            dfs(v, h + 1, b + (1 if is_card(u, v) else 0))
            path.pop()
            vis[v] = False

    dfs(src, 0, 0)
    return reps


def path_links(T, path: List[int]) -> Tuple[List[float], List[int]]:
    """(link lengths in km, edge ids) along a node path."""
    L, E = [], []
    for a, b in zip(path[:-1], path[1:]):
        e = T.eid[(a, b) if a < b else (b, a)]
        L.append(T.elen[e])
        E.append(e)
    return L, E


# ---------------------------------------------------------------------------
# closed-form references (the paper's model)
# ---------------------------------------------------------------------------
def physics(rounds: int, memory: str, T2: float, kappa: float = 0.002, gen_rate: float = math.inf):
    return Q.Physics(rounds=rounds, memory=memory, T2=T2, kappa=kappa, gen_rate=gen_rate)


def closed_form_path(ph, link_len: List[float]) -> float:
    """Paper's end-to-end fidelity of a clean path with the given link lengths."""
    S = sum(math.log((4 * ph.purify(ph.clean_F(l)) - 1) / 3) for l in link_len)
    return ph.e2e(S, len(link_len), sum(link_len))


def closed_form_slack(N: int, rounds: int, memory: str, T2: float, kappa: float = 0.002,
                      gen_rate: float = math.inf) -> int:
    return DL.slack(Q.Physics(rounds=rounds, memory=memory, T2=T2, kappa=kappa, gen_rate=gen_rate), N)


def closed_form_fstar(N: int, k: int, rounds: int, memory: str, T2: float, kappa: float = 0.002,
                      gen_rate: float = math.inf) -> float:
    ph = Q.Physics(rounds=rounds, memory=memory, T2=T2, kappa=kappa, gen_rate=gen_rate)
    wc, wd = DL.wcd(ph)
    return DL.fstar(ph, N, k, wc, wd)


def closed_form_min_T2(N: int, rounds: int, target: int, memory: str, gen_rate: float,
                       kappa: float = 0.002) -> float:
    """design_law.min_T2 (optimistic storage accounting t = kappa D + C_r / R)."""
    return DL.min_T2(N, rounds, target, memory=memory, gen_rate=gen_rate, kappa=kappa)


def pessimistic_min_T2_depol2(N: int, rounds: int, gen_rate: float, target: int = 1,
                              kappa: float = 0.002) -> float:
    """design_law.py's pessimistic accounting (every link pair waits C_r / R in
    two depolarizing memories, then the delivered pair waits kappa D).  Copied
    from design_law.main() (where it is a nested function) without change of
    logic; test_des_physics.py checks it reproduces design_law.json."""
    def slack_links(r, T2, R_, N, kap=kappa):
        ph = Q.Physics(rounds=r, memory="depol2", T2=T2, kappa=kap)
        t = Q.Physics(rounds=r, gen_rate=R_).t_purify()
        wc, wd = DL.wcd(ph)
        f = math.exp(-2.0 * t / T2)
        wc, wd = wc * f, wd * f
        if DL.fstar(ph, N, 0, wc, wd) < Q.F_TH:
            return -1
        s_ = 0
        while s_ < 40 and DL.fstar(ph, N, s_ + 1, wc, wd) >= Q.F_TH:
            s_ += 1
        return s_

    lo, hi = 1e-4, 1e4
    if slack_links(rounds, hi, gen_rate, N) < target:
        return math.inf
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if slack_links(rounds, mid, gen_rate, N) >= target:
            hi = mid
        else:
            lo = mid
    return hi
