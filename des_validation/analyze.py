"""
analyze.py -- turn <results dir>/*.json (default results_des/) into the summary, LaTeX tables, macros,
figures and a draft paragraph.  Reads whatever blocks exist.

Outputs (written to the same results dir, default results_des/):
  summary.md                  every number and every automatic check (PASS / FAIL / NOTE)
  tab_des_regime.tex          regime calls: closed form vs DES (Mode F, DEJMPS, Mode R fast generation)
  tab_des_gate.tex            BBPSSW vs DEJMPS with gate noise (Mode F)
  tab_des_designchart.tex     5x5 design chart: minimum T2 for s >= 1, closed form vs DES (Mode R)
  tab_des_routing.tex         routing: Mode F (the paper's E1 cell) and Mode R (success, latency, rate)
  des_numbers.tex             \\newcommand macros
  fig_des_profiles.pdf/png    F*(k) on 5x5, r = 3
  fig_des_threshold.pdf/png   DES F*(1) vs T2 against the closed form (5x5, r = 3)
  des_paragraph_draft.md      draft text, every qualitative statement conditional on the numbers

Definitions used throughout
  F*(k)        best mean delivered fidelity over clean paths with h* + k hops (one
               representative path per link composition; all paths of a composition
               are identical in Mode F and differ by < 1 SE in Mode R).
  s (mean)     largest k such that F*(0..k) >= F_th (the fidelity of the average
               delivered state; the paper's definition when fidelity is deterministic).
  s (median)   same with "P(F >= F_th) >= 1/2" in place of "mean >= F_th" (reported in T).
  headroom     paired difference oracle - BFS in success probability; significance by the
               exact McNemar test (two-sided); 95% interval from the exact conditional
               (Clopper-Pearson) interval on the discordant trials, consistent with the test.

Usage:  python analyze.py                      # reads and writes results_des/
        python analyze.py --results DIR        # any other results directory
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results_des")
F_TH = 0.8
Z = 1.959964
H_STAR_5 = 4

CHECKS = []          # (status, name, detail)
LINES = []           # summary.md body
MACROS = {}


def check(ok, name, detail, note_only=False):
    CHECKS.append(("NOTE" if note_only else ("PASS" if ok else "FAIL"), name, detail))


def load(b):
    p = os.path.join(OUT, f"{b}.json")
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None


def wilson(k, n):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    m = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return max(0.0, c - m), min(1.0, c + m)


def mcnemar_exact(b01, b10):
    """Two-sided exact McNemar p-value (binomial test on the discordant pairs)."""
    n = b01 + b10
    if n == 0:
        return 1.0
    k = min(b01, b10)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / 2.0 ** n
    return min(1.0, 2.0 * tail)


def clopper_pearson(k, n, alpha=0.05):
    from scipy.stats import beta
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def paired(a, b):
    """a - b on the same trials.  diff (pp); exact McNemar p (two-sided); and a 95% interval
    consistent with that test: conditional on the n_d discordant trials, theta = P(a-only |
    discordant) has a Clopper-Pearson interval [tL, tU], and diff = (2 theta - 1) n_d / n.
    (The CP interval excludes 1/2 exactly when the doubled-tail exact test has p < 0.05.)"""
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    n = len(a)
    b01, b10 = int(np.sum(a & ~b)), int(np.sum(~a & b))
    nd = b01 + b10
    diff = (b01 - b10) / n
    if nd:
        tL, tU = clopper_pearson(b01, nd)
        lo, hi = (2 * tL - 1) * nd / n, (2 * tU - 1) * nd / n
    else:
        lo = hi = 0.0
    return {"diff": 100 * diff, "lo": 100 * lo, "hi": 100 * hi,
            "a_only": b01, "b_only": b10, "p": mcnemar_exact(b01, b10)}


def macro(name, value):
    MACROS[name] = value


def fmt_T2(x):
    if x is None or not math.isfinite(x):
        return "--"
    return f"{x:#.3g}"


MEMNAME = {"half": "baseline", "deph2": "dephasing", "depol2": "depolarizing"}


def pexp(p):
    """P value for LaTeX, e.g. $1.1\\times10^{-19}$."""
    if p >= 1e-3:
        return f"{p:.2g}"
    if p < 1e-30:                      # beyond this the exact value carries no information
        return "$<10^{-30}$"
    m, e = f"{p:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def fs(s, capped=False):
    """slack for LaTeX: '$-1$', '$\\ge$8', '2'."""
    if capped:
        return f"$\\ge${s}"
    return f"${s}$" if s < 0 else str(s)


def slack_from(prof: dict, ok) -> int:
    """largest k such that ok(prof[0..k]) holds."""
    s = -1
    for k in sorted(prof):
        if k == s + 1 and ok(prof[k]):
            s = k
        else:
            break
    return s


def slack_from_means(prof: dict) -> int:
    return slack_from(prof, lambda F: F >= F_TH)


# ---------------------------------------------------------------------------
# PF: Mode F profiles (+ protocols, + gate noise)
# ---------------------------------------------------------------------------
def analyse_PF(d):
    LINES.append("## PF -- Mode F clean-path profiles (the paper's assumptions; deterministic given success)\n")
    LINES.append("Mode F is an implementation cross-check: it can disagree with the closed form only through an error "
                 "in the closed-form algebra or in the simulator. Mode R (PR, T, S2) is the physics test.\n")
    P = d["paths"]
    best = defaultdict(dict)          # (cfg, mem, scheme, eps) -> k -> F
    depth = {}
    sd_max, n_cf = 0.0, 0
    diff_max = {True: 0.0, False: 0.0}      # keyed by "eps == 0"
    for v in P.values():
        cfg = f"g{v['N']}_r{v['r']}"
        k = v["h"] - v["h_star"]
        key = (cfg, v["memory"], v["scheme"], float(v.get("eps", 0.0)))
        best[key][k] = max(best[key].get(k, -1.0), v["mean"])
        depth[key] = max(depth.get(key, 0), k)
        sd_max = max(sd_max, v["sd"])
        if v["closed_form"] is not None:
            z = float(v.get("eps", 0.0)) == 0.0
            diff_max[z] = max(diff_max[z], abs(v["mean"] - v["closed_form"]))
            n_cf += 1
    check(diff_max[True] < 1e-12, "PF: DES Mode F equals the closed form on every clean path (BBPSSW + twirl, ideal gates)",
          f"max |F_DES - F_closed| = {diff_max[True]:.2e}")
    # with gate noise the paper's code rounds F to 10 decimals before each cached noisy round
    check(diff_max[False] < 1e-9, "PF: same with gate noise (tolerance 1e-9: the paper's code rounds F to 1e-10 per round)",
          f"max |F_DES - F_closed| = {diff_max[False]:.2e}; {n_cf} paths in total")
    check(sd_max < 1e-9, "PF: Mode F is deterministic given success", f"max sd over seeds = {sd_max:.2e}")

    def sl(cfg, mem, sch, eps):
        key = (cfg, mem, sch, eps)
        s = slack_from_means(best[key])
        return s, s == depth[key]

    out, gate = {}, {}
    mism = []
    for ckey, ref in sorted(d["closed_form"].items()):
        cfg, mem, e = ckey.split("|")
        eps = float(e[1:])
        s_cf, Fs = ref["slack"], ref["Fstar"]
        margin = min(abs(F - F_TH) for F in Fs[max(s_cf, 0): max(s_cf, 0) + 2]) if s_cf >= 0 else F_TH - Fs[0]
        s_tw, _ = sl(cfg, mem, "bbpssw_twirl", eps)
        s_dj, dj_cap = sl(cfg, mem, "dejmps", eps)
        if s_tw != s_cf:
            mism.append(ckey)
        row = {"s_cf": s_cf, "margin": margin, "s_F": s_tw, "s_dejmps": s_dj, "dejmps_capped": dj_cap,
               "F0": best[(cfg, mem, "bbpssw_twirl", eps)].get(0), "F0_dejmps": best[(cfg, mem, "dejmps", eps)].get(0)}
        if eps == 0.0:
            s_nt, _ = sl(cfg, mem, "bbpssw", eps)
            row.update({"s_cnot": s_nt, "F0_cnot": best[(cfg, mem, "bbpssw", eps)].get(0)})
            out[(cfg, mem)] = row
        else:
            gate[(cfg, eps)] = row
    check(not mism, "PF: DES regime call s equals the closed-form s in every (configuration, channel, gate error)",
          f"{len(d['closed_form'])} cells" + (f"; MISMATCH in {mism}" if mism else ""))
    LINES.append("| config | channel | s closed | margin delta | s DES Mode F | s DES DEJMPS | s DES CNOT, no twirl | "
                 "F*(0) BBPSSW | F*(0) DEJMPS | F*(0) CNOT, no twirl |")
    LINES.append("|---|---|---|---|---|---|---|---|---|---|")
    for (cfg, mem), v in sorted(out.items()):
        LINES.append(f"| {cfg} | {mem} | {v['s_cf']} | {v['margin']:.4f} | {v['s_F']} | "
                     f"{'>=' if v['dejmps_capped'] else ''}{v['s_dejmps']} | {v['s_cnot']} | {v['F0']:.4f} | "
                     f"{v['F0_dejmps']:.4f} | {v['F0_cnot']:.4f} |")
    LINES.append("\n'>=k': the slack reached the enumeration depth (h* + 8). CNOT without twirl = bilateral CNOT with "
                 "neither twirl nor DEJMPS rotation; it detects bit flips only (an ablation of the twirl).\n")
    if gate:
        LINES.append("Gate noise (two-qubit depolarizing eps after each bilateral CNOT; half memory, T2 = 0.1 s):\n")
        LINES.append("| config | eps | s closed (BBPSSW) | s DES BBPSSW | s DES DEJMPS | F*(0) BBPSSW | F*(0) DEJMPS |")
        LINES.append("|---|---|---|---|---|---|---|")
        for (cfg, eps), v in sorted(gate.items()):
            LINES.append(f"| {cfg} | {eps:g} | {v['s_cf']} | {v['s_F']} | {'>=' if v['dejmps_capped'] else ''}"
                         f"{v['s_dejmps']} | {v['F0']:.4f} | {v['F0_dejmps']:.4f} |")
        LINES.append("")
    n_flip = sum(1 for v in out.values() if v["s_dejmps"] != v["s_cf"])
    check(True, "PF: regime calls under DEJMPS (ideal gates)",
          f"{n_flip}/{len(out)} cells change; min DEJMPS slack {min(v['s_dejmps'] for v in out.values())}", note_only=True)
    o = out.get(("g5_r3", "half"))
    if o:
        if diff_max[True] > 0:
            m, e = f"{diff_max[True]:.0e}".split("e")
            macro("desPFmaxdiff", f"${m}\\times10^{{{int(e)}}}$")
        else:
            macro("desPFmaxdiff", "$0$")
        macro("desDejmpsSfiveThree", fs(o["s_dejmps"], o["dejmps_capped"]))
        macro("desFzeroDejmpsFiveThree", f"{o['F0_dejmps']:.3f}")
    return out, gate, best


# ---------------------------------------------------------------------------
# PR: Mode R at the paper's T2, fast generation
# ---------------------------------------------------------------------------
def analyse_PR(d):
    LINES.append("\n## PR -- Mode R at the paper's T2, fast generation (per-memory storage resolved, kappa kept)\n")
    P = d["paths"]
    groups = defaultdict(lambda: defaultdict(list))
    for v in P.values():
        key = (f"g{v['N']}_r{v['r']}", v["memory"], int(v["R"]))
        groups[key][v["h"] - v["h_star"]].append(v)
    rows, out = [], {}
    for key in sorted(groups):
        info = {k: max(lst, key=lambda x: x["mean"]) for k, lst in groups[key].items()}
        prof = {k: b["mean"] for k, b in info.items()}
        s_des = slack_from_means(prof)
        s_cf = next(iter(groups[key].values()))[0]["closed_form_slack"]
        capped = s_des == max(info)
        decisive = [k for k in (s_des, s_des + 1) if k in info and k >= 0]
        unresolved = [k for k in decisive if abs(info[k]["mean"] - F_TH) < 3 * info[k]["se"]]
        dF = [info[k]["mean"] - info[k]["closed_form"] for k in sorted(info)]
        out[key] = {"s_des": s_des, "s_cf": s_cf, "capped": capped, "unresolved": unresolved, "prof": prof,
                    "info": {k: {kk: info[k][kk] for kk in ("mean", "se", "p_pass", "closed_form", "t_mean")} for k in info},
                    "dF_mean": float(np.mean(dF)), "dF_min": float(min(dF)), "dF_max": float(max(dF))}
        cells = "; ".join(f"k={k}: {info[k]['mean']:.4f}+-{info[k]['se']:.4f} (cf {info[k]['closed_form']:.4f}, "
                          f"P>=.8 {info[k]['p_pass']:.2f})" for k in sorted(info))
        rows.append(f"| {key[0]} | {key[1]} | {key[2]:g} | {s_cf} | {'>=' if capped else ''}{s_des}"
                    f"{' (unresolved k=' + ','.join(map(str, unresolved)) + ')' if unresolved else ''} | "
                    f"{np.mean(dF):+.4f} | {cells} |")
    LINES.append("| config | channel | R (Hz) | s closed (same R) | s DES (mean F) | mean of F_DES - F_closed | "
                 "best-class F*(k) DES +- SE (closed form, P(F>=0.8)) |")
    LINES.append("|---|---|---|---|---|---|---|")
    LINES.extend(rows)
    for R in sorted({k[2] for k in out}, reverse=True):
        cells = [v for k, v in out.items() if k[2] == R]
        agree = sum(1 for v in cells if v["s_des"] == v["s_cf"])
        check(True, f"PR: Mode R calls at R = {R:g} Hz vs closed form at the same R",
              f"{agree}/{len(cells)} agree; F_DES - F_closed between {min(v['dF_min'] for v in cells):+.4f} and "
              f"{max(v['dF_max'] for v in cells):+.4f}", note_only=True)
    return out


# ---------------------------------------------------------------------------
# T: design chart
# ---------------------------------------------------------------------------
def best_points(curve):
    best = {}
    for pt in curve["points"].values():
        g = pt["gi"]
        if g not in best or pt["mean"] > best[g]["mean"]:
            best[g] = pt
    return [best[g] for g in sorted(best)]


def crossing(T2s, ys):
    """first upward zero crossing of ys (already minus the level) on a log-T2 grid."""
    xs = np.log(T2s)
    ys = np.asarray(ys, float)
    if ys[0] >= 0:
        return None, "below_grid"
    for i in range(len(xs) - 1):
        if ys[i] < 0 <= ys[i + 1]:
            return math.exp(xs[i] + (xs[i + 1] - xs[i]) * (-ys[i]) / (ys[i + 1] - ys[i])), "ok"
    return None, "above_grid"


def analyse_T(d, nboot=1000, label="T"):
    LINES.append(f"\n## {label} -- design chart check ({ {'T6': '6x6', 'T7': '7x7'}.get(label, '5x5') }): smallest T2 with s >= 0 / s >= 1, Mode R"
                 + (" (DEJMPS)" if label == "TD" else "") + "\n")
    LINES.append("Mean rule: mean F*(k) >= F_th. Median rule: P(F*(k) >= F_th) >= 1/2. "
                 "CI: 95% bootstrap over deliveries (replicates that do not cross inside the grid are counted as 'dropped').\n")
    rng = np.random.default_rng(7)
    out, rows = {}, []
    for ckey, c in sorted(d["curves"].items()):
        pts = best_points(c)
        T2s = np.array([p["T2"] for p in pts])
        arrs = [np.asarray(p["F"]) for p in pts]
        x, st = crossing(T2s, [a.mean() - F_TH for a in arrs])
        xm, stm = crossing(T2s, [np.mean(a >= F_TH) - 0.5 for a in arrs])
        boots, dropped = [], 0
        if st == "ok":
            for _ in range(nboot):
                m = [a[rng.integers(0, len(a), len(a))].mean() - F_TH for a in arrs]
                xb, sb = crossing(T2s, m)
                if sb == "ok":
                    boots.append(xb)
                else:
                    dropped += 1
        lo, hi = (np.percentile(boots, [2.5, 97.5]) if len(boots) >= 0.9 * nboot else (None, None))
        cf, pe = c["closed_form_T2"], c["pessimistic_T2"]
        ratio = (x / cf) if (x and cf and math.isfinite(cf)) else None
        out[ckey] = {**{k: c[k] for k in ("r", "R", "memory", "kappa", "target")}, "N": c.get("N", 5), "closed_form_T2": cf,
                     "pessimistic_T2": pe, "des_T2": x, "des_lo": lo, "des_hi": hi, "status": st,
                     "des_T2_median": xm, "status_median": stm, "ratio": ratio, "boot_dropped": dropped}
        rows.append(f"| {c['r']} | {c['R']:g} | {c['memory']} | {c['kappa']:g} | s>={c['target']} | {fmt_T2(cf)} | "
                    f"{fmt_T2(pe) if pe else '--'} | {fmt_T2(x) if x else st} | "
                    f"{('[' + fmt_T2(lo) + ', ' + fmt_T2(hi) + ']') if lo is not None else '--'} | "
                    f"{fmt_T2(xm) if xm else stm} | {f'{ratio:.2f}' if ratio else '--'} | {dropped} |")
    LINES.append("| r | R (Hz) | channel | kappa (s/km) | target | T2 closed form, optimistic (s) | T2 pessimistic, depol2 (s) | "
                 "T2 DES, mean rule (s) | 95% CI | T2 DES, median rule (s) | DES / closed | boot dropped |")
    LINES.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    LINES.extend(rows)
    bad = [k for k, v in out.items() if v["status"] != "ok"]
    check(not bad, "T: every DES threshold bracketed by the T2 grid (0.5x-16x of the closed form)",
          f"{len(out) - len(bad)}/{len(out)} bracketed" + (f"; not bracketed: {bad}" if bad else ""))
    noci = [k for k, v in out.items() if v["status"] == "ok" and v["des_lo"] is None]
    check(not noci, "T: bootstrap CI available for every bracketed threshold (>= 90% of replicates cross)",
          f"{len(noci)} without CI" + (f": {noci}" if noci else ""))
    for k, v in out.items():
        if v["pessimistic_T2"] and v["des_T2"]:
            inside = v["closed_form_T2"] <= v["des_T2"] <= v["pessimistic_T2"]
            check(True, f"T: {k}: DES between optimistic and pessimistic accounting?",
                  f"optimistic {fmt_T2(v['closed_form_T2'])} s, DES {fmt_T2(v['des_T2'])} s, "
                  f"pessimistic {fmt_T2(v['pessimistic_T2'])} s -> {'inside' if inside else 'OUTSIDE'}", note_only=True)
    return out


# ---------------------------------------------------------------------------
# S1, S2: routing
# ---------------------------------------------------------------------------
def analyse_S1(d):
    LINES.append("\n## S1 -- Mode F routing, 5x5, r = 3, T2 = 0.1 s, half (the paper's E1 cell, same seeds)\n")
    tr, n = d["trials"], d["n"]
    e1 = None
    p_e1 = os.path.normpath(os.path.join(HERE, "..", "results_tqe", "trials", "E1.npz"))
    if os.path.exists(p_e1):
        z = np.load(p_e1)
        pols = list(z["E1_g5_r3__policies"])
        e1 = {p: z["E1_g5_r3__succ"][:n, pols.index(p)] for p in d["policies"] if p in pols}
    rows, out, total_mis, worst = [], {}, 0, 0.0
    for p in d["policies"]:
        ok_d = np.asarray(tr[p]["ok_des"], bool)
        ok_m = np.asarray(tr[p]["ok_model"], bool)
        Fd, Fm = np.asarray(tr[p]["F_des"], float), np.asarray(tr[p]["F_model"], float)
        sw = np.asarray(tr[p]["swap_failed"], bool)
        valid = (~sw) & (np.asarray(tr[p]["hops"]) > 0)
        mis = int(np.sum(ok_d != ok_m))
        total_mis += mis
        if valid.any():
            worst = max(worst, float(np.max(np.abs(Fd[valid] - Fm[valid]))))
        e1m = int(np.sum(ok_d != e1[p])) if e1 is not None and p in e1 else None
        out[p] = {"sr_des": float(ok_d.mean()), "sr_model": float(ok_m.mean()), "mismatch_model": mis, "mismatch_E1": e1m}
        rows.append(f"| {p} | {100 * ok_d.mean():.2f} | {100 * ok_m.mean():.2f} | "
                    + (f"{100 * e1[p].mean():.2f} | {e1m} |" if e1m is not None else "-- | -- |"))
    LINES.append(f"n = {n} trials (seeds 20260922 + t, identical to the paper's E1).\n")
    LINES.append("| policy | success DES (%) | success qunet_tqe (%) | success paper E1 (%) | DES vs E1 mismatching trials |")
    LINES.append("|---|---|---|---|---|")
    LINES.extend(rows)
    check(total_mis == 0, "S1: DES success equals qunet_tqe success on every trial and policy",
          f"{total_mis} mismatches over {n * len(d['policies'])} (trial, policy) pairs")
    check(worst < 1e-9, "S1: DES fidelity equals qunet_tqe on every path without a swap failure", f"max |dF| = {worst:.2e}")
    if e1 is not None:
        m = sum(v["mismatch_E1"] or 0 for v in out.values())
        check(m == 0, "S1: DES success equals the paper's stored E1 trials", f"{m} mismatches")
    else:
        check(True, "S1: comparison with stored E1 trials", "skipped (../results_tqe/trials/E1.npz not found)", note_only=True)
    for p in ("shortest", "full", "oracle"):
        macro(f"desSone{p.capitalize()}", f"{100 * out[p]['sr_des']:.2f}")
    return out


def analyse_S2(d):
    LINES.append("\n## S2 -- Mode R routing, 5x5, r = 3, R = 1 kHz: headroom, latency, delivery rate\n")
    LINES.append("'oracle' = the paper model's exact oracle path (model at the same R), evaluated in the DES; it need "
                 "not be the DES-optimal path. Latency = t_used = establishment + kappa*D storage. Delivery rate = "
                 "successes / total attempt time (attempts back to back).\n")
    out = {}
    for vkey, v in d["variants"].items():
        LINES.append(f"\n### {vkey}: DES thresholds T2(s>=0) = {v['T2_des_s0']:.4g} s, T2(s>=1) = {v['T2_des_s1']:.4g} s\n")
        for label, run in v["runs"].items():
            tr, pols = run["trials"], d["policies"]
            n = len(tr[pols[0]]["ok_des"])
            want = 0 if label == "s0" else 1
            prof = run.get("profile", {})
            verified = bool(run.get("label_verified", False))
            check(verified, f"S2 {vkey} {label}: DES clean-path profile confirms s = {want} at T2 = {run['T2']:.4g} s",
                  "; ".join(f"F*({k}) = {p['mean']:.4f}+-{p['se']:.4f} (P>=.8 {p['p_pass']:.2f})" for k, p in sorted(prof.items())))
            LINES.append(f"**{label}** (T2 = {run['T2']:.4g} s; intended DES regime s = {want}; verified: {verified}), n = {n}\n")
            LINES.append("DES clean-path profile at this T2: " + "; ".join(
                f"F*({k}) = {p['mean']:.4f} +- {p['se']:.4f}, P(F>=0.8) = {p['p_pass']:.2f}" for k, p in sorted(prof.items())) + "\n")
            LINES.append("| policy | success % [95% CI] | mean hops | latency mean / median (ms) | delivery rate (pairs/s) | "
                         "diagnostic: rate over clean-path attempts only | raw pairs per attempt |")
            LINES.append("|---|---|---|---|---|---|---|")
            res = {}
            for p in pols:
                ok = np.asarray(tr[p]["ok_des"], bool)
                t = np.asarray(tr[p]["t_used"], float)
                h = np.asarray(tr[p]["hops"], float)
                raw = np.asarray(tr[p]["raw"], float)
                nn = np.asarray(tr[p]["n_noisy"])
                valid = np.isfinite(t) & (h > 0)
                clean = valid & (nn == 0)
                lo, hi = wilson(int(ok.sum()), n)
                # attempts for which a policy returns no path (the scalar oracle, when it knows that
                # every path fails) are charged the policy's mean attempt time, like failed attempts
                rate = ok.sum() / (n * t[valid].mean()) if valid.any() else 0.0
                rate_c = ok[clean].sum() / t[clean].sum() if clean.any() else 0.0
                res[p] = {"sr": float(ok.mean()), "lo": lo, "hi": hi, "lat_mean": float(np.mean(t[valid])),
                          "lat_med": float(np.median(t[valid])), "rate": float(rate), "rate_clean": float(rate_c),
                          "hops": float(h[valid].mean()), "raw": float(raw[valid].mean())}
                LINES.append(f"| {p} | {100 * ok.mean():.2f} [{100 * lo:.1f}, {100 * hi:.1f}] | {h[valid].mean():.2f} | "
                             f"{1e3 * res[p]['lat_mean']:.2f} / {1e3 * res[p]['lat_med']:.2f} | {rate:.2f} | {rate_c:.2f} | "
                             f"{res[p]['raw']:.1f} |")
            hd = paired(tr["oracle"]["ok_des"], tr["shortest"]["ok_des"])
            ad = paired(tr["full"]["ok_des"], tr["shortest"]["ok_des"])
            # where do oracle-only successes come from?  detour (h > h*) vs minimum-hop
            ok_o = np.asarray(tr["oracle"]["ok_des"], bool)
            ok_b = np.asarray(tr["shortest"]["ok_des"], bool)
            h_o = np.asarray(tr["oracle"]["hops"])
            only = ok_o & ~ok_b
            n_only_detour = int(np.sum(only & (h_o > H_STAR_5)))
            # P(success) of a clean one-hop detour chosen by the oracle (no swap failure)
            cl_det = (np.asarray(tr["oracle"]["n_noisy"]) == 0) & (h_o == H_STAR_5 + 1) & \
                     ~np.asarray(tr["oracle"]["swap_failed"], bool)
            p_det = float(ok_o[cl_det].mean()) if cl_det.any() else float("nan")
            # latency of CLEAN paths by hop count (distinct paths per trial, all policies pooled)
            # (t_used includes kappa*D, which grows with path length by construction; t_est does not)
            lat_h, est_h, lat_noisy = defaultdict(list), defaultdict(list), []
            has_est = "t_est" in tr[pols[0]]
            for t_i in range(n):
                seen = set()
                for p in pols:
                    h_ = int(tr[p]["hops"][t_i])
                    if h_ <= 0:
                        continue
                    sig = (h_, tr[p]["t_used"][t_i], tr[p]["n_noisy"][t_i])
                    if sig in seen:
                        continue
                    seen.add(sig)
                    if tr[p]["n_noisy"][t_i] == 0:
                        lat_h[h_].append(tr[p]["t_used"][t_i])
                        if has_est:
                            est_h[h_].append(tr[p]["t_est"][t_i])
                    else:
                        lat_noisy.append(tr[p]["t_used"][t_i])

            def msn(d_):
                return {h_: (float(np.mean(v_)), float(np.std(v_, ddof=1) / math.sqrt(len(v_))) if len(v_) > 1 else 0.0,
                             len(v_)) for h_, v_ in sorted(d_.items())}
            lat_clean, est_clean = msn(lat_h), msn(est_h)
            lat_noisy_m = float(np.mean(lat_noisy)) if lat_noisy else float("nan")
            # BFS's own attempts over paths with a noisy link (for the sentence about BFS's latency)
            nb = [tr["shortest"]["t_used"][t_i] for t_i in range(n)
                  if int(tr["shortest"]["hops"][t_i]) > 0 and tr["shortest"]["n_noisy"][t_i] > 0]
            lat_noisy_bfs = float(np.mean(nb)) if nb else float("nan")
            LINES.append(f"\nHeadroom oracle - BFS = {hd['diff']:+.2f} pp, exact 95% conditional on the discordant trials "
                         f"[{hd['lo']:+.2f}, {hd['hi']:+.2f}] (its upper limit equals the estimate when no trial favours BFS), "
                         f"exact McNemar p = {hd['p']:.3g} (oracle-only {hd['a_only']}, BFS-only {hd['b_only']}; "
                         f"{n_only_detour} of the oracle-only successes used a detour). "
                         f"PG - BFS = {ad['diff']:+.2f} pp, p = {ad['p']:.3g}. "
                         f"P(success | oracle path is a clean one-hop detour, swaps ok) = {p_det:.3f}.\n")
            if est_clean:
                LINES.append("Establishment time (without kappa*D) of clean paths by hop count (ms, mean +- SE, n): "
                             + "; ".join(f"h={h_}: {1e3 * m:.2f} +- {1e3 * s_:.2f} ({c})" for h_, (m, s_, c) in est_clean.items()) + "\n")
            LINES.append("Latency incl. kappa*D of clean paths by hop count (ms, mean +- SE, n): "
                         + "; ".join(f"h={h_}: {1e3 * m:.2f} +- {1e3 * s:.2f} ({c})" for h_, (m, s, c) in lat_clean.items())
                         + (f"; paths with a noisy link: {1e3 * lat_noisy_m:.2f} ms ({len(lat_noisy)})" if lat_noisy else "")
                         + (f"; BFS paths with a noisy link: {1e3 * lat_noisy_bfs:.2f} ms ({len(nb)})." if nb else ".") + "\n")
            out[(vkey, label)] = {"T2": run["T2"], "res": res, "headroom": hd, "adaptive": ad, "verified": verified,
                                  "profile": prof, "n_only_detour": n_only_detour, "p_det": p_det,
                                  "n_det": int(cl_det.sum()), "lat_clean": lat_clean, "est_clean": est_clean,
                                  "lat_noisy": lat_noisy_m, "lat_noisy_bfs": lat_noisy_bfs, "n": n}
            check(True, f"S2 {vkey} {label}: headroom oracle - BFS",
                  f"{hd['diff']:+.2f} pp, exact McNemar p = {hd['p']:.3g}", note_only=True)
    return out


# ---------------------------------------------------------------------------
# LaTeX tables
# ---------------------------------------------------------------------------
def dejmps_exact(cfg, mem):
    """Exact DEJMPS slack from the closed form of the paper's code (code/design_law.py, stored in
    ../results_tqe/design_law.json), at the coherence time used for this channel in block PF."""
    p = os.path.normpath(os.path.join(HERE, "..", "results_tqe", "design_law.json"))
    return json.load(open(p))["dejmps"]["des_cells"][f"{cfg}|{mem}"]


def tex_regime(pf, pr):
    L = [r"\begin{table*}[t]", r"\centering",
         r"\caption{Discrete-event check (\textnormal{SeQUeNCe} kernel, density matrices, gate-level purification and swapping) of "
         r"the regime calls at $T_2=0.1$\,s (baseline map), 0.2\,s (per-qubit dephasing) and 0.3\,s (per-qubit "
         r"depolarizing). The last two lie just below the coherence times at which these channels reproduce the baseline regime map "
         r"(Section~\ref{sec:latency}), so that the three channels give different calls. $s$: detour slack; $\delta$: closed-form margin at $R=\infty$. Mode~F uses the model's "
         r"assumptions (an implementation cross-check); DEJMPS replaces twirled BBPSSW (ideal gates, Mode~F); Mode~R "
         r"resolves the storage time of every stored pair at raw-pair rate $R$ (calls on the mean delivered fidelity; closed form at the same $R$ in parentheses).}",
         r"\label{tab:des}", r"\footnotesize",
         r"\begin{tabular}{llcccccc}", r"\toprule",
         r"Grid & Memory & $s$ (closed) & $\delta$ & $s$ Mode F & $s$ DEJMPS & $s$ Mode R, $10^5$\,Hz & $s$ Mode R, $10^4$\,Hz\\",
         r"\midrule"]
    for (cfg, mem), v in sorted(pf.items()):
        N, r = cfg[1], cfg.split("_r")[1]
        ex = dejmps_exact(cfg, mem)
        assert (v["dejmps_capped"] and v["s_dejmps"] <= ex) or (not v["dejmps_capped"] and v["s_dejmps"] == ex), (cfg, mem)
        cells = []
        for R in (100000, 10000):
            x = pr.get((cfg, mem, R)) if pr else None
            if x is None:
                cells.append("--")
            else:
                u = r"$^\dagger$" if x["unresolved"] else ""
                cells.append(f"{fs(x['s_des'], x['capped'])}{u} ({fs(x['s_cf'])})")
        L.append(f"${N}\\times{N}$, $r={r}$ & {MEMNAME[mem]} & {fs(v['s_cf'])} & {v['margin']:.4f} & {fs(v['s_F'])} & "
                 f"{fs(v['s_dejmps'], v['dejmps_capped'])} ({ex}) & " + " & ".join(cells) + r"\\")
    any_u = bool(pr) and any(x["unresolved"] for x in pr.values())
    L += [r"\bottomrule", r"\end{tabular}",
          r"\\[2pt]{\scriptsize DEJMPS: simulation (Bell-diagonal states, no twirl), with the closed form of Supplementary Section~\ref{sec:supp-dejmps} in parentheses; "
          r"$\ge$: slack reached the enumeration depth of the simulation ($h^*+8$ for DEJMPS; BBPSSW profiles are enumerated to $h^*+4$)."
          + (r" $^\dagger$ a deciding mean $F^*(k)$ lies within 3 standard errors of $F_{\rm th}$." if any_u else "") + "}",
          r"\end{table*}"]
    return "\n".join(L)


def tex_gate(gate):
    L = [r"\begin{table}[t]", r"\centering",
         r"\caption{Slack $s$ with gate error $\epsilon$ (two-qubit depolarizing after each of the two local CNOTs of a bilateral CNOT), Mode~F, "
         r"baseline memory map, $T_2=0.1$\,s. B: twirled BBPSSW (closed form = simulation); D: DEJMPS (simulation).}",
         r"\label{tab:desgate}", r"\footnotesize", r"\setlength{\tabcolsep}{4pt}", r"\begin{tabular}{lcccc}",
         r"\toprule", r"Grid & $\epsilon$ (\%) & $s$ (B) & $s$ (D) & $F^*(0)$ (B / D)\\", r"\midrule"]
    for (cfg, eps), v in sorted(gate.items()):
        N, r = cfg[1], cfg.split("_r")[1]
        L.append(f"${N}\\times{N}$, $r={r}$ & {100 * eps:g} & {fs(v['s_F'])} & {fs(v['s_dejmps'], v['dejmps_capped'])} & "
                 f"{v['F0']:.3f} / {v['F0_dejmps']:.3f}\\\\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(L)


def tex_designchart_dejmps(t):
    L = [r"\begin{table}[!t]", r"\centering",
         r"\caption{DEJMPS: smallest coherence time $T_2$ (s) for which the $5\times5$ request is routing-sensitive ($s\ge1$), "
         r"closed form (Section~\ref{sec:supp-dejmps}, optimistic accounting) and discrete-event simulation (Mode~R, mean rule, 95\% bootstrap interval; median rule); $\kappa=2$\,ms/km.}",
         r"\label{tab:desdejmps}", r"\scriptsize", r"\setlength{\tabcolsep}{2.5pt}",
         r"\begin{tabular}{@{}cccccc@{}}", r"\toprule",
         r"$r$ & $R$ (Hz) & Memory & Closed form & Simulation, mean rule & Median rule \\", r"\midrule"]
    for k, v in sorted(t.items(), key=lambda kv: (kv[1]["r"], -kv[1]["R"], kv[1]["memory"])):
        des = fmt_T2(v["des_T2"]) if v["des_T2"] else v["status"].replace("_", " ")
        ci = f" [{fmt_T2(v['des_lo'])}, {fmt_T2(v['des_hi'])}]" if v["des_lo"] is not None else ""
        L.append(f"{v['r']} & {v['R']:g} & {MEMNAME[v['memory']]} & {fmt_T2(v['closed_form_T2'])} & {des}{ci} & "
                 f"{fmt_T2(v['des_T2_median']) if v['des_T2_median'] else '--'}\\\\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(L)


def tex_designchart_six(t):
    L = [r"\begin{table*}[!t]", r"\centering",
         r"\caption{Design-chart check on larger grids: smallest coherence time $T_2$ (s) for which the corner-to-corner "
         r"request is routing-sensitive ($s\ge1$) with BBPSSW: closed form (optimistic accounting), pessimistic all-links-wait "
         r"accounting and discrete-event simulation (Mode~R, mean rule, 95\% bootstrap interval; median rule); $\kappa=2$\,ms/km (Mode~R also counts the light delay of the classical messages).}",
         r"\label{tab:dessix}", r"\footnotesize", r"\setlength{\tabcolsep}{4pt}",
         r"\begin{tabular}{ccccccccc}", r"\toprule",
         r"Grid & $r$ & $R$ (Hz) & Memory & Closed form & Pessimistic & Simulation, mean rule & Median rule & Mean rule/closed \\", r"\midrule"]
    for k, v in sorted(t.items(), key=lambda kv: (kv[1].get("N", 6), kv[1]["r"], -kv[1]["R"], kv[1]["memory"])):
        N = v.get("N", 6)
        des = fmt_T2(v["des_T2"]) if v["des_T2"] else v["status"].replace("_", " ")
        ci = f" [{fmt_T2(v['des_lo'])}, {fmt_T2(v['des_hi'])}]" if v["des_lo"] is not None else ""
        ratio_s = f"{v['ratio']:.2f}" if v["ratio"] else "--"
        pe = _pess_all(v["memory"], v["R"], N, v["r"])
        L.append(f"${N}\\times{N}$ & {v['r']} & {v['R']:g} & {MEMNAME[v['memory']]} & {fmt_T2(v['closed_form_T2'])} & "
                 f"{fmt_T2(pe) if pe else '--'} & {des}{ci} & "
                 f"{fmt_T2(v['des_T2_median']) if v['des_T2_median'] else '--'} & {ratio_s}\\\\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
    return "\n".join(L)


def _pess_all(mem, R, N, r):
    """Pessimistic all-links-wait accounting for every memory channel (../results_tqe/design_law.json,
    written by code/design_law.py); None if unavailable."""
    p = os.path.normpath(os.path.join(HERE, "..", "results_tqe", "design_law.json"))
    if not os.path.exists(p):
        return None
    d = json.load(open(p)).get("latency_all_links", {}).get(f"{mem}|{int(R)}|{N}")
    return d["by_r"][r - 1] if d else None


def tex_designchart(t):
    L = [r"\begin{table*}[t]", r"\centering",
         r"\caption{Smallest coherence time $T_2$ (s) for which the $5\times5$ request with BBPSSW is routing-sensitive ($s\ge1$): "
         r"closed form (optimistic storage accounting $t=\kappa D+C_r/R$), the pessimistic all-links-wait accounting "
         r"(Section~\ref{sec:latency}) and the discrete-event simulation (Mode~R): mean rule (mean $F^*(1)\ge F_{\rm th}$, "
         r"95\% bootstrap interval) and median rule ($F^*(1)\ge F_{\rm th}$ in at least half of the deliveries); $\kappa=2$\,ms/km (Mode~R also counts the light delay of the classical messages, about 0.25\% of $\kappa$).}",
         r"\label{tab:desdesign}", r"\footnotesize",
         r"\begin{tabular}{cccccccc}", r"\toprule",
         r"$r$ & $R$ (Hz) & Memory & Closed form & Pessimistic & Simulation, mean rule & Median rule & Mean rule/closed\\", r"\midrule"]
    for k, v in sorted(t.items(), key=lambda kv: (kv[1]["kappa"] == 0, kv[1]["r"], -kv[1]["R"], kv[1]["memory"])):
        if v["target"] != 1 or v["kappa"] == 0:
            continue
        des = fmt_T2(v["des_T2"]) if v["des_T2"] else v["status"].replace("_", " ")
        ci = f" [{fmt_T2(v['des_lo'])}, {fmt_T2(v['des_hi'])}]" if v["des_lo"] is not None else ""
        ratio_s = f"{v['ratio']:.2f}" if v["ratio"] else "--"
        pe = _pess_all(v["memory"], v["R"], v.get("N", 5), v["r"])
        L.append(f"{v['r']} & {v['R']:g} & {MEMNAME[v['memory']]} & {fmt_T2(v['closed_form_T2'])} & "
                 f"{fmt_T2(pe) if pe else '--'} & {des}{ci} & "
                 f"{fmt_T2(v['des_T2_median']) if v['des_T2_median'] else '--'} & {ratio_s}\\\\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
    return "\n".join(L)


def tex_routing(s1, s2):
    L = [r"\begin{table*}[t]", r"\centering",
         r"\caption{Routing on $5\times5$, $r=3$ with every path evaluated in the discrete-event simulation. Mode~F: the "
         r"$5\times5$, $r=3$ row of Table~\ref{tab:headline}, same trials. Mode~R: raw-pair rate 1\,kHz, $\kappa=2$\,ms/km, at a $T_2$ where the simulation "
         r"gives $s=0$ and one where it gives $s=1$ (both verified on clean paths). Success rates in \%; headroom with "
         r"exact \textnormal{McNemar} $P$. Latency and delivery rate: Supplementary Table~\ref{S-tab:deslatency}.}",
         r"\label{tab:desrouting}", r"\footnotesize", r"\begin{tabular}{lcccccc}", r"\toprule",
         r"Setting & BFS & PG & Dijkstra & Oracle$^\ast$ & Oracle$^\ast-$BFS (pp) & $P$\\", r"\midrule"]
    if s1:
        L.append("Mode F & " + " & ".join(f"{100 * s1[p]['sr_des']:.2f}" for p in ("shortest", "full", "dijkstra", "oracle"))
                 + f" & {100 * (s1['oracle']['sr_des'] - s1['shortest']['sr_des']):+.2f} & --" + r"\\")
    if s2:
        for (vkey, label), v in sorted(s2.items()):
            mem = vkey.split("|")[0]
            hd = v["headroom"]
            L.append(f"Mode R, {MEMNAME[mem]}, $T_2={v['T2']:#.3g}$\\,s, {'$s=0$' if label == 's0' else '$s=1$'} & "
                     + " & ".join(f"{100 * v['res'][p]['sr']:.1f}" for p in ("shortest", "full", "dijkstra", "oracle"))
                     + f" & {hd['diff']:+.1f} & " + pexp(hd["p"]) + r"\\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\\[2pt]{\scriptsize $^\ast$ exact oracle of the scalar model; in Mode~R its path is evaluated in the "
          r"simulation and need not be the best path there.}", r"\end{table*}"]
    return "\n".join(L)


def tex_latency(s2):
    """Mode R mean latency and delivery rate of the routers of Table IV (Supplementary Material)."""
    L = [r"\begin{table}[!t]", r"\centering",
         r"\caption{Mode~R latency and delivery rate of the routers of main-text Table~\ref{tab:desrouting} ($5\times5$, $r=3$, 1\,kHz, "
         r"$\kappa=2$\,ms/km): mean latency (ms, from the start of an attempt until the delivered pair has been stored "
         r"for $\kappa D$, over attempts with a path) / delivery rate (pairs per second of attempt time; attempts without a "
         r"path are charged the mean attempt time).}", r"\label{tab:deslatency}", r"\scriptsize", r"\setlength{\tabcolsep}{2.5pt}", r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
         r"Setting & BFS & PG & Dijkstra & Oracle\\", r"\midrule"]
    for (vkey, label), v in sorted(s2.items()):
        mem = vkey.split("|")[0]
        L.append(f"{MEMNAME[mem]}, {'$s=0$' if label == 's0' else '$s=1$'} & "
                 + " & ".join(f"{1e3 * v['res'][p]['lat_mean']:.1f} / {v['res'][p]['rate']:.1f}"
                              for p in ("shortest", "full", "dijkstra", "oracle")) + r"\\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(L)


# ---------------------------------------------------------------------------
# figures (same style and validated palette as ../code/make_figures.py)
# ---------------------------------------------------------------------------
def figures(pf_best, prd, td):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
    BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#8a8985"
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
        "axes.edgecolor": INK2, "axes.linewidth": 0.6, "xtick.color": INK2, "ytick.color": INK2,
        "axes.labelcolor": INK, "text.color": INK, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.5, "axes.axisbelow": True, "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.4, "lines.markersize": 4.5, "legend.frameon": False,
        "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.02})
    W1 = 3.5
    made = []
    if pf_best is not None:
        fig, ax = plt.subplots(figsize=(W1, 2.3))
        cfg, mem = "g5_r3", "half"
        series = [("bbpssw_twirl", "simulation, BBPSSW + twirl (= closed form)", BLUE, "o"),
                  ("dejmps", "simulation, DEJMPS", AQUA, "^"),
                  ("bbpssw", "simulation, CNOT without twirl", GRAY, "s")]
        for sch, lab, col, mk in series:
            prof = pf_best.get((cfg, mem, sch, 0.0), {})
            ks = sorted(prof)
            ax.plot(ks, [prof[k] for k in ks], color=col, marker=mk, ls="-", label=lab)
        if prd is not None:
            x = prd.get((cfg, mem, 10000))
            if x:
                ks = sorted(x["info"])
                ax.errorbar(ks, [x["info"][k]["mean"] for k in ks], yerr=[Z * x["info"][k]["se"] for k in ks],
                            color=ORANGE, marker="D", ls=(0, (3, 2)), label="simulation, Mode R, $R=10^4$ Hz")
        ax.axhline(F_TH, color=INK2, lw=0.8, ls=(0, (2, 2)))
        ax.text(0.99, F_TH + 0.004, r"$F_{\rm th}$", color=INK2, fontsize=7, ha="right", va="bottom",
                transform=ax.get_yaxis_transform())
        ax.set_xlabel("Extra hops $k$ beyond the minimum ($5\\times5$, $r=3$, $T_2=0.1$ s)")
        ax.set_ylabel("Best clean-path $F^*(k)$")
        ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.12), fontsize=6.5)
        fig.savefig(os.path.join(OUT, "fig_des_profiles.pdf"))
        fig.savefig(os.path.join(OUT, "fig_des_profiles.png"))
        plt.close(fig)
        made.append("fig_des_profiles")
    if td is not None:
        import qunet_bridge as B
        fig, axs = plt.subplots(1, 2, figsize=(2 * W1, 2.3), sharey=True)
        L = [1.0, math.sqrt(2), 1.0, math.sqrt(2), math.sqrt(2)]   # one-hop detour on 5x5: 2 axial + 3 diagonal links
        for ax, R in zip(axs, (1000.0, 100.0)):
            for mem, col, mk in (("half", BLUE, "o"), ("deph2", ORANGE, "s"), ("depol2", AQUA, "^")):
                c = td["curves"].get(f"T|r3|R{int(R)}|{mem}|k0.002|s1")
                if not c:
                    continue
                pts = best_points(c)
                T2s = np.array([p["T2"] for p in pts])
                ax.errorbar(T2s, [p["mean"] for p in pts], yerr=[Z * p["se"] for p in pts], color=col, marker=mk,
                            ls="none", label=f"simulation, {MEMNAME[mem]}")
                fine = np.geomspace(T2s.min(), T2s.max(), 60)
                cf = [B.closed_form_path(B.physics(3, mem, x, kappa=0.002, gen_rate=R), L) for x in fine]
                ax.plot(fine, cf, color=col, lw=1.0, label=f"closed form, {MEMNAME[mem]}")
            ax.axhline(F_TH, color=INK2, lw=0.8, ls=(0, (2, 2)))
            ax.set_xscale("log")
            ax.set_xlabel("Memory coherence time $T_2$ (s)")
            ax.set_title(f"$R={R:g}$ Hz, $5\\times5$, $r=3$, one-hop detour")
        axs[0].set_ylabel("$F^*(1)$")
        axs[1].legend(loc="lower right", fontsize=6, ncol=1)
        fig.savefig(os.path.join(OUT, "fig_des_threshold.pdf"))
        fig.savefig(os.path.join(OUT, "fig_des_threshold.png"))
        plt.close(fig)
        made.append("fig_des_threshold")
    return made


# ---------------------------------------------------------------------------
# draft paragraph: every qualitative statement is conditional on the numbers
# ---------------------------------------------------------------------------
def _cfgname(cfg):
    N, r = cfg[1], cfg.split("_r")[1]
    return f"{N}$\\times${N}, $r={r}$"


def _j(vals, f):
    return ", ".join(f.format(x) for x in vals)


def paragraph(pf, gate, pr, t, s1, s2):
    P = []
    if pf:
        mism = [(k, v) for k, v in pf.items() if v["s_F"] != v["s_cf"]]
        txt = ("**Discrete-event check.** We re-evaluated the regime calls in a discrete-event simulation built on the "
               "SeQUeNCe kernel [wu2021sequence] with full density matrices: purification and swapping are gate circuits "
               "(bilateral CNOTs, Bell-state measurement, Pauli corrections) with outcomes sampled by the Born rule, memory "
               "noise is a Kraus channel applied for the time each pair is actually stored, and fidelity is read from the "
               "density matrix. ")
        if not mism:
            k53 = pf.get(("g5_r3", "half"))
            txt += (f"With the paper's assumptions (Mode F: twirled BBPSSW, instantaneous purification, storage charged to the "
                    f"delivered pair) the simulation reproduces every clean-path fidelity to numerical precision and all "
                    f"{len(pf)} regime calls of Table~\\ref{{tab:des}}"
                    + (f", including the 5$\\times$5, $r=3$ call with margin $\\delta={k53['margin']:.4f}$" if k53 else "")
                    + "; this checks the closed-form algebra and its implementation, not the model. The physical test is "
                      "Mode R, which resolves the storage time of every stored pair.")
        else:
            txt += (f"With the paper's assumptions (Mode F) the simulation DISAGREES with the closed form in {len(mism)} of "
                    f"{len(pf)} cells. INVESTIGATE BEFORE USING ANY RESULT.")
        P.append(txt)
        dj = pf.get(("g5_r3", "half"))
        dj6 = pf.get(("g6_r3", "half"))
        if dj and dj["s_dejmps"] != dj["s_cf"]:
            txt = (f"The regime map is specific to the purification protocol. With DEJMPS in place of twirled BBPSSW (ideal "
                   f"gates), $F^*(0)$ on 5$\\times$5, $r=3$ is {dj['F0_dejmps']:.3f} instead of {dj['F0']:.3f} and the slack is "
                   f"{'at least ' if dj['dejmps_capped'] else ''}{dj['s_dejmps']}")
            if dj6 and dj6["s_cf"] < 0 <= dj6["s_dejmps"]:
                txt += (f"; 6$\\times$6, $r=3$, infeasible under BBPSSW, has slack "
                        f"{'at least ' if dj6['dejmps_capped'] else ''}{dj6['s_dejmps']}")
            txt += "."
            if gate:
                ge = max(e for (_, e) in gate)
                g = gate.get(("g5_r3", ge))
                if g:
                    txt += (f" With gate error $\\epsilon={ge:g}$ the 5$\\times$5, $r=3$ slack is {fs(g['s_F'])} under BBPSSW and "
                            f"{'at least ' if g['dejmps_capped'] else ''}{fs(g['s_dejmps'])} under DEJMPS "
                            f"(Table~\\ref{{tab:desgate}}).")
            txt += (" Proposition~1 needs only monotonicity and is unaffected, but the closed form (\\ref{eq:profile}) and the "
                    "regime map must be recomputed for each protocol")
            if any(v["s_dejmps"] > 5 for v in pf.values()):
                txt += ", and DEJMPS slacks exceed the range $s\\le5$ in which Dijkstra was observed to match the oracle"
            txt += "."
            if all(v["s_cnot"] == -1 for v in pf.values()):
                txt += " Omitting the twirl (bilateral CNOT only) leaves phase errors uncorrected; no configuration is feasible."
            P.append(txt)
    if pr:
        parts = []
        for R in sorted({k[2] for k in pr}, reverse=True):
            cells = {k: v for k, v in pr.items() if k[2] == R}
            agree = sum(1 for v in cells.values() if v["s_des"] == v["s_cf"])
            lo = min(v["dF_min"] for v in cells.values())
            hi = max(v["dF_max"] for v in cells.values())
            word = (f"lower by {abs(hi):.4f}--{abs(lo):.4f}" if hi < 0 else
                    f"higher by {lo:.4f}--{hi:.4f}" if lo > 0 else f"different by {lo:+.4f} to {hi:+.4f}")
            flips = [f"{_cfgname(k[0])} ({k[1]}): $s={v['s_cf']}\\to{v['s_des']}$" for k, v in sorted(cells.items())
                     if v["s_des"] != v["s_cf"]]
            s = (f"at $R=10^{{{int(round(math.log10(R)))}}}$\\,Hz, {agree} of {len(cells)} calls agree with the closed form at "
                 f"the same $R$ and the delivered fidelity is {word}")
            if flips:
                s += " (changed: " + "; ".join(flips) + ")"
            nun = sum(1 for v in cells.values() if v["unresolved"])
            if nun:
                s += f"; {nun} call(s) lie within 3 standard errors of the threshold"
            parts.append(s)
        P.append("In Mode R (Poisson raw-pair generation, nested recurrence on each link, swapping as soon as both links are "
                 "ready, classical messages at $2\\times10^5$\\,km/s, $\\kappa$ kept) at the paper's coherence times: "
                 + "; ".join(parts) + ". The closed form charges storage only to the delivered pair ($t=\\kappa D+C_r/R$), "
                 "whereas the simulation charges every stored pair for its own waiting time.")
    if t:
        rr = [v for v in t.values() if v["target"] == 1 and v["ratio"]]
        if rr:
            lo, hi = min(v["ratio"] for v in rr), max(v["ratio"] for v in rr)
            if lo > 1:
                txt = (f"At slow generation ($R=100$\\,Hz and 1\\,kHz) the design chart is optimistic: the coherence time at "
                       f"which the simulation makes the 5$\\times$5 request routing-sensitive is {lo:.1f}--{hi:.1f} times the "
                       f"closed-form value (Table~\\ref{{tab:desdesign}}).")
            elif hi < 1:
                txt = (f"At slow generation the design chart is conservative: the simulation needs {lo:.1f}--{hi:.1f} times "
                       f"the closed-form coherence time (Table~\\ref{{tab:desdesign}}).")
            else:
                txt = (f"At slow generation the simulated coherence-time requirement is {lo:.1f}--{hi:.1f} times the "
                       f"closed-form value (Table~\\ref{{tab:desdesign}}).")
            dp = [v for v in rr if v["pessimistic_T2"]]
            if dp:
                inside = [v for v in dp if v["closed_form_T2"] <= v["des_T2"] <= v["pessimistic_T2"]]
                txt += (" For depolarizing memories it lies between the optimistic and the pessimistic accounting in "
                        + (f"all {len(dp)} cases, so the two closed forms bracket it." if len(inside) == len(dp)
                           else f"{len(inside)} of {len(dp)} cases."))
            P.append(txt)
    if s1:
        mis = sum(v["mismatch_model"] for v in s1.values())
        if mis == 0:
            P.append(f"Evaluating every router's path in the simulation on the paper's 5$\\times$5, $r=3$ trials (same seeds, "
                     f"Mode F) reproduces the success of every router on every trial (BFS {100 * s1['shortest']['sr_des']:.2f}\\%, "
                     f"PG {100 * s1['full']['sr_des']:.2f}\\%, oracle {100 * s1['oracle']['sr_des']:.2f}\\%).")
        else:
            P.append(f"Mode F routing: {mis} (trial, router) outcomes DIFFER from the paper's simulator. INVESTIGATE.")
    if s2:
        x0 = {k: v for (k, l), v in s2.items() if l == "s0"}
        x1 = {k: v for (k, l), v in s2.items() if l == "s1"}
        unver = [f"{k} {l}" for (k, l), v in s2.items() if not v["verified"]]
        if unver:
            P.append("Mode R routing: the intended regime could not be verified for " + ", ".join(unver)
                     + "; do not use these rows.")
        if x0 and x1:
            mems = list(x1)
            txt = ("With timing resolved (1\\,kHz, " + ", ".join(m.split("|")[0] for m in mems) + " memories), the routing "
                   "headroom (paper-model oracle minus BFS) is " + _j([x0[m]["headroom"]["diff"] for m in mems if m in x0], "{:+.1f}")
                   + "\\,pp at a $T_2$ where the simulation gives $s=0$ and "
                   + _j([x1[m]["headroom"]["diff"] for m in mems], "{:+.1f}") + "\\,pp where it gives $s=1$. ")
            sig0 = [m for m in x0 if x0[m]["headroom"]["p"] < 0.05]
            pos0 = [m for m in sig0 if x0[m]["headroom"]["diff"] > 0]
            pdet = [x0[m]["p_det"] for m in x0 if math.isfinite(x0[m]["p_det"])]
            if not sig0:
                n0 = _j(sorted({x0[m]["n"] for m in x0}), "{}")
                if pdet and max(pdet) > 0:
                    txt += (f"At $s=0$ the headroom is not significant at $n={n0}$ (exact McNemar), but the oracle's clean "
                            f"one-hop detours still succeed in a fraction {_j(pdet, '{:.3f}')} of attempts, so a small "
                            f"headroom is expected and a larger $n$ may resolve it. ")
                else:
                    txt += (f"At $s=0$ the headroom is not significant at $n={n0}$ and no clean detour succeeded, as "
                            f"Proposition~1 predicts. ")
            else:
                txt += (f"At $s=0$ the headroom is significant (exact McNemar $P<0.05$) in {len(sig0)} of {len(x0)} settings"
                        + ("" if len(pos0) == len(sig0) else " (not all positive)") + ". ")
                # explanation only if verified: all oracle-only successes are detours, the detour's mean
                # fidelity is below F_th, and the oracle's clean detours pass in a fraction strictly in (0, 1/2)
                expl = [m for m in pos0
                        if x0[m]["n_only_detour"] == x0[m]["headroom"]["a_only"]
                        and "1" in x0[m]["profile"] and x0[m]["profile"]["1"]["mean"] < F_TH
                        and math.isfinite(x0[m]["p_det"]) and 0 < x0[m]["p_det"] < 0.5]
                if expl and len(expl) == len(pos0):
                    txt += ("In every such case the extra successes are detours whose mean fidelity is below threshold "
                            f"but which succeed in a fraction {_j([x0[m]['p_det'] for m in expl], '{:.3f}')} of attempts, "
                            "because storage times vary from delivery to delivery: near the boundary the regime call is "
                            "soft, and Proposition~1 applies to each delivery's fidelity rather than to the mean. ")
            # latency: one-hop detour on clean links, establishment time (without kappa*D) and with it
            for lab, key in (("establishment time (without $\\kappa D$)", "est_clean"), ("latency including $\\kappa D$", "lat_clean")):
                parts = []
                for m in mems:
                    lc = x1[m].get(key) or {}
                    if H_STAR_5 in lc and H_STAR_5 + 1 in lc:
                        (m0, e0, _), (m1, e1, _) = lc[H_STAR_5], lc[H_STAR_5 + 1]
                        parts.append((m.split("|")[0], m0, m1, m1 - m0, math.sqrt(e0 ** 2 + e1 ** 2)))
                if not parts:
                    continue
                sig = [x for x in parts if abs(x[3]) > 2 * x[4]]
                if sig:
                    txt += (f"On clean links a one-hop detour changes the mean {lab} by "
                            + "; ".join(f"{1e3 * d:+.1f}\\,ms ({n_})" for n_, _, _, d, _ in sig)
                            + (f" (not significant for {', '.join(x[0] for x in parts if x not in sig)})" if len(sig) < len(parts) else "")
                            + ". ")
                else:
                    txt += f"On clean links a one-hop detour does not measurably change the mean {lab}. "
            r = x1[mems[0]]["res"]
            ro, rb = r["oracle"]["rate"], r["shortest"]["rate"]
            txt += (f"Counting every attempt back to back, the delivery rate at $s=1$ is {ro:.1f} pairs/s with the oracle's paths "
                    f"and {rb:.1f} with BFS ({mems[0].split('|')[0]} memory): routing "
                    + ("raises" if ro > rb else "does not raise") + " the delivery rate"
                    + (", by about the ratio of the success probabilities." if ro > rb and
                       abs(ro / rb - r["oracle"]["sr"] / max(r["shortest"]["sr"], 1e-9)) < 0.15 * ro / rb else "."))
            P.append(txt)
    P.append("> DRAFT generated from summary.md. Every number above is in summary.md; check each sentence against it before "
             "use. Labels referenced: tab:des, tab:desgate, tab:desdesign, eq:profile; bib key wu2021sequence.")
    return "\n\n".join(P)


def paper_macros(pf, gate, pr, pr_d, t, s1, s1_d, s2):
    """Every DES number quoted in the paper's text, as a macro (all names letters only)."""
    f2, f3, f4 = (lambda x: f"{x:.2f}"), (lambda x: f"{x:.3f}"), (lambda x: f"{x:.4f}")
    if pf:
        macro("desNcells", str(len(pf)))
        macro("desNgateCells", str(len(gate)) if gate else "0")
        o = pf[("g5_r3", "half")]
        macro("desFzeroBBFiveThree", f3(o["F0"]))
        macro("desFzeroDJFiveThree", f3(o["F0_dejmps"]))
        caps = [v for v in pf.values() if v["dejmps_capped"]]
        macro("desDJminSlack", str(min(v["s_dejmps"] for v in pf.values())))
        macro("desDJallCapped", "yes" if len(caps) == len(pf) else "no")
        o6 = pf[("g6_r3", "half")]
        macro("desFzeroBBSixThree", f3(o6["F0"]))
        macro("desFzeroDJSixThree", f3(o6["F0_dejmps"]))
    if gate:
        for eps, name in ((0.005, "Half"), (0.01, "One")):
            cells = [v for (c, e), v in gate.items() if abs(e - eps) < 1e-12]
            dj = [v["s_dejmps"] for v in cells]
            bb = [v["s_F"] for v in cells]
            capped = any(v["dejmps_capped"] and v["s_dejmps"] == max(dj) for v in cells)
            macro(f"desGate{name}DJmin", str(min(dj)))
            macro(f"desGate{name}DJmax", (r"$\ge$" if capped else "") + str(max(dj)))
            macro(f"desGate{name}BBmin", fs(min(bb)))
            macro(f"desGate{name}BBmax", fs(max(bb)))
    if pr:
        for R, name in ((100000, "Fast"), (10000, "Slow")):
            cells = {k: v for k, v in pr.items() if k[2] == R}
            macro(f"desPR{name}Agree", str(sum(1 for v in cells.values() if v["s_des"] == v["s_cf"])))
            macro(f"desPR{name}N", str(len(cells)))
            macro(f"desPR{name}dFlo", f4(min(abs(v["dF_max"]) for v in cells.values())))
            macro(f"desPR{name}dFhi", f4(max(abs(v["dF_min"]) for v in cells.values())))
            sds = [p["sd"] for p in pr_d["paths"].values() if int(p["R"]) == R]
            macro(f"desPR{name}SDlo", f4(min(sds)))
            macro(f"desPR{name}SDhi", f4(max(sds)))
        # calls that change at 10 kHz: closed-form margin at that rate = best closed-form F*(s_cf) - F_th
        best_cf = defaultdict(float)
        for p in pr_d["paths"].values():
            kk = (f"g{p['N']}_r{p['r']}", p["memory"], int(p["R"]), p["h"] - p["h_star"])
            best_cf[kk] = max(best_cf[kk], p["closed_form"])
        flips = [(k, v) for k, v in pr.items() if k[2] == 10000 and v["s_des"] != v["s_cf"]]
        marg = [best_cf[(k[0], k[1], 10000, v["s_cf"])] - F_TH for k, v in flips]
        macro("desPRSlowFlipN", str(len(flips)))
        macro("desPRSlowFlipAllMinusOne", "yes" if all(v["s_des"] == v["s_cf"] - 1 for _, v in flips) else "no")
        macro("desPRSlowFlipMarginLo", f4(min(marg)))
        macro("desPRSlowFlipMarginHi", f4(max(marg)))
        allm = []
        for k, v in pr.items():
            if k[2] == 10000 and v["s_cf"] >= 0:
                allm.append(best_cf[(k[0], k[1], 10000, v["s_cf"])] - F_TH)
        macro("desPRSlowMarginMax", f4(max(allm)))
        # same one-sided distance (deciding F*(s_cf) - F_th) over the calls that do NOT change
        unch = [best_cf[(k[0], k[1], 10000, v["s_cf"])] - F_TH for k, v in pr.items()
                if k[2] == 10000 and v["s_cf"] >= 0 and v["s_des"] == v["s_cf"]]
        if unch:
            macro("desPRSlowUnchMarginLo", f4(min(unch)))
            macro("desPRSlowUnchMarginHi", f4(max(unch)))
        x = pr[("g5_r3", "half", 10000)]
        macro("desPRSlowFiveThreeFone", f4(x["info"][1]["mean"]))
        macro("desPRSlowFiveThreeFoneCF", f4(x["info"][1]["closed_form"]))
        y = pr[("g5_r3", "half", 100000)]
        macro("desPRFastFiveThreeFone", f4(y["info"][1]["mean"]))
    if t:
        rr = [v for v in t.values() if v["target"] == 1 and v["ratio"] and v["kappa"] > 0]
        macro("desTratioLo", f"{min(v['ratio'] for v in rr):.1f}")
        macro("desTratioHi", f"{max(v['ratio'] for v in rr):.1f}")
        rp = [v for v in rr if v["memory"] != "half"]
        macro("desTratioPhysLo", f"{min(v['ratio'] for v in rp):.1f}")
        macro("desTratioPhysHi", f"{max(v['ratio'] for v in rp):.1f}")
        r0 = [v for v in t.values() if v["target"] == 1 and v["ratio"] and v["kappa"] == 0]
        if r0:
            macro("desTratioNoKapLo", f"{min(v['ratio'] for v in r0):.1f}")
            macro("desTratioNoKapHi", f"{max(v['ratio'] for v in r0):.1f}")
        dp = [v for v in t.values() if v["target"] == 1 and v["pessimistic_T2"] and v["des_T2"]]
        macro("desTbracketIn", str(sum(1 for v in dp if v["closed_form_T2"] <= v["des_T2"] <= v["pessimistic_T2"])))
        macro("desTbracketN", str(len(dp)))
        for mem, nm in (("half", "Half"), ("deph2", "Deph"), ("depol2", "Depol")):
            for R, rn in ((100.0, "Hundred"), (1000.0, "Onek")):
                for r in (3, 4):
                    v = t.get(f"T|r{r}|R{int(R)}|{mem}|k0.002|s1")
                    if v and v["des_T2"]:
                        rw = "Three" if r == 3 else "Four"
                        macro(f"desT{rn}{nm}{rw}", f"{v['des_T2']:.2f}")
                        macro(f"desT{rn}{nm}{rw}Lo", f"{v['des_lo']:.2f}")
                        macro(f"desT{rn}{nm}{rw}Hi", f"{v['des_hi']:.2f}")
                        macro(f"desT{rn}{nm}{rw}CF", f"{v['closed_form_T2']:.2f}")
    if s1:
        macro("desSonePairs", f"{s1_d['n'] * len(s1_d['policies']):,}".replace(",", "{,}"))
        macro("desSoneN", f"{s1_d['n']:,}".replace(",", "{,}"))
    if s2:
        x0 = [v for (k, l), v in s2.items() if l == "s0"]
        x1 = [v for (k, l), v in s2.items() if l == "s1"]
        macro("desStwoN", f"{x1[0]['n']:,}".replace(",", "{,}"))
        macro("desHeadOneLo", f"{min(v['headroom']['diff'] for v in x1):.1f}")
        macro("desHeadOneHi", f"{max(v['headroom']['diff'] for v in x1):.1f}")
        macro("desHeadZeroLo", f"{min(v['headroom']['diff'] for v in x0):.1f}")
        macro("desHeadZeroHi", f"{max(v['headroom']['diff'] for v in x0):.1f}")
        pmax = max(v["headroom"]["p"] for v in x0)
        macro("desHeadZeroPexp", str(int(math.floor(math.log10(pmax))) + 1) if pmax > 0 else "-300")
        macro("desDetourPassLo", f"{100 * min(v['p_det'] for v in x0):.0f}")
        macro("desDetourPassHi", f"{100 * max(v['p_det'] for v in x0):.0f}")
        macro("desZeroFoneLo", f4(min(v["profile"]["1"]["mean"] for v in x0)))
        macro("desZeroFoneHi", f4(max(v["profile"]["1"]["mean"] for v in x0)))
        macro("desAllOnlyDetours", "yes" if all(v["n_only_detour"] == v["headroom"]["a_only"] for v in x0) else "no")
        est = [v["est_clean"][H_STAR_5 + 1][0] - v["est_clean"][H_STAR_5][0] for v in x1]
        lat = [v["lat_clean"][H_STAR_5 + 1][0] - v["lat_clean"][H_STAR_5][0] for v in x1]
        est0 = [v["est_clean"][H_STAR_5 + 1][0] - v["est_clean"][H_STAR_5][0] for v in x0]
        lat0 = [v["lat_clean"][H_STAR_5 + 1][0] - v["lat_clean"][H_STAR_5][0] for v in x0]
        macro("desDetourEstAllLo", f"{1e3 * min(est + est0):.1f}")
        macro("desDetourEstAllHi", f"{1e3 * max(est + est0):.1f}")
        macro("desDetourLatAllLo", f"{1e3 * min(lat + lat0):.1f}")
        macro("desDetourLatAllHi", f"{1e3 * max(lat + lat0):.1f}")
        macro("desNoisyLatLo", f"{1e3 * min(v['lat_noisy_bfs'] for v in x0 + x1):.1f}")
        macro("desNoisyLatHi", f"{1e3 * max(v['lat_noisy_bfs'] for v in x0 + x1):.1f}")
        macro("desCleanLatBFS", f"{1e3 * min(v['lat_clean'][H_STAR_5][0] for v in x0 + x1):.1f}")
        macro("desCleanLatBFSHi", f"{1e3 * max(v['lat_clean'][H_STAR_5][0] for v in x0 + x1):.1f}")
        macro("desDetourEstLo", f"{1e3 * min(est):.1f}")
        macro("desDetourEstHi", f"{1e3 * max(est):.1f}")
        macro("desDetourLatLo", f"{1e3 * min(lat):.1f}")
        macro("desDetourLatHi", f"{1e3 * max(lat):.1f}")
        rb = [v["res"]["shortest"]["rate"] for v in x1]
        ro = [v["res"]["oracle"]["rate"] for v in x1]
        macro("desRateBFSLo", f"{min(rb):.1f}")
        macro("desRateBFSHi", f"{max(rb):.1f}")
        macro("desRateOrLo", f"{min(ro):.1f}")
        macro("desRateOrHi", f"{max(ro):.1f}")
        nd = sum(1 for v in x0 + x1 if v["res"]["dijkstra"]["sr"] > v["res"]["oracle"]["sr"])
        words = {0: "none", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}
        macro("desDijBeatsOracleN", words.get(nd, str(nd)))
        dmax = max([100 * (v["res"]["dijkstra"]["sr"] - v["res"]["oracle"]["sr"]) for v in x0 + x1] + [0.0])
        macro("desDijBeatsMax", f"{dmax:.2f}")
        cap = [(v["res"]["full"]["sr"] - v["res"]["shortest"]["sr"]) / (v["res"]["oracle"]["sr"] - v["res"]["shortest"]["sr"])
               for v in x1]
        macro("desPgCapLo", f"{100 * min(cap):.0f}")
        macro("desPgCapHi", f"{100 * max(cap):.0f}")
        macro("desStwoSettings", words.get(len(x0 + x1), str(len(x0 + x1))))
        macro("desSRbfsOne", f"{100 * min(v['res']['shortest']['sr'] for v in x1):.1f}")
        macro("desSRorOneLo", f"{100 * min(v['res']['oracle']['sr'] for v in x1):.1f}")
        macro("desSRorOneHi", f"{100 * max(v['res']['oracle']['sr'] for v in x1):.1f}")
        macro("desSRbfsOneHi", f"{100 * max(v['res']['shortest']['sr'] for v in x1):.1f}")


def main():
    global OUT
    import argparse
    ap = argparse.ArgumentParser(description="Analyse DES results (reads <dir>/*.json, writes outputs to <dir>).")
    ap.add_argument("--results", default=OUT, help="results directory (default: results_des next to this file)")
    OUT = os.path.abspath(ap.parse_args().results)
    if not os.path.isdir(OUT):
        raise SystemExit(f"results directory not found: {OUT}")
    print(f"results directory: {OUT}")
    pf_d, pr_d, t_d, s1_d, s2_d, td_d, t6_d = (load(b) for b in ("PF", "PR", "T", "S1", "S2", "TD", "T6"))
    if not t6_d and load("T6a") and load("T6b"):      # the two halves of T6 (run_des.py --only T6a,T6b)
        t6_d = load("T6a"); t6_d["curves"].update(load("T6b")["curves"])
    t7_d = None
    if load("T7a") and load("T7b"):                   # 7x7, r = 5, in two halves by rate
        t7_d = load("T7a"); t7_d["curves"].update(load("T7b")["curves"])
    LINES.append("# DES validation -- summary\n")
    metas = [x["meta"] for x in (pf_d, pr_d, t_d, s1_d, s2_d, td_d) if x]
    if metas:
        LINES.append(f"Environment: Python {metas[0]['python']}, SeQUeNCe {metas[0]['sequence']}, {metas[0]['platform']}; "
                     f"workers {metas[0]['workers']}.\n")
        if any(m["quick"] for m in metas):
            LINES.append("**WARNING: at least one block was run with --quick (small n). Not for the paper.**\n")
            print("WARNING: at least one block in this results directory was run with --quick (small n). Not for the paper.")
    n_head = len(LINES)
    pf, gate, pf_best = analyse_PF(pf_d) if pf_d else (None, None, None)
    pr = analyse_PR(pr_d) if pr_d else None
    t = analyse_T(t_d) if t_d else None
    td = analyse_T(td_d, label="TD") if td_d else None
    t6 = analyse_T(t6_d, label="T6") if t6_d else None
    t7 = analyse_T(t7_d, label="T7") if t7_d else None
    s1 = analyse_S1(s1_d) if s1_d else None
    s2 = analyse_S2(s2_d) if s2_d else None
    w = lambda name, text: open(os.path.join(OUT, name), "w", encoding="utf-8").write(text + "\n")
    if pf:
        w("tab_des_regime.tex", tex_regime(pf, pr or {}))
    if gate:
        w("tab_des_gate.tex", tex_gate(gate))
    if t:
        w("tab_des_designchart.tex", tex_designchart(t))
        # thresholds for the design-chart figure of the paper (code/make_figures.py)
        w("design_points.json", json.dumps({k: {f: v[f] for f in ("r", "R", "memory", "kappa", "target",
                                                                     "closed_form_T2", "pessimistic_T2", "des_T2",
                                                                     "des_lo", "des_hi", "des_T2_median")}
                                            for k, v in t.items()}, indent=1, default=float))
    if s1 or s2:
        w("tab_des_routing.tex", tex_routing(s1, s2))
        if s2:
            w("tab_des_latency.tex", tex_latency(s2))
    if td:
        w("tab_des_dejmps.tex", tex_designchart_dejmps(td))
        w("design_points_dejmps.json", json.dumps({k: {f: v[f] for f in ("r", "R", "memory", "kappa", "target",
                                                                            "closed_form_T2", "des_T2", "des_lo", "des_hi",
                                                                            "des_T2_median")}
                                                   for k, v in td.items()}, indent=1, default=float))
        rat = [v["ratio"] for v in td.values() if v["ratio"] and v["memory"] != "half"]
        macro("desDJratioLo", f"{min(rat):.1f}"); macro("desDJratioHi", f"{max(rat):.1f}")
        rat_all = [v["ratio"] for v in td.values() if v["ratio"]]
        macro("desDJratioAllLo", f"{min(rat_all):.1f}"); macro("desDJratioAllHi", f"{max(rat_all):.1f}")
    if t6 and t7:
        t6 = {**t6, **t7}                               # one table and one set of ratios for 6x6 and 7x7
    if t6:
        w("tab_des_six.tex", tex_designchart_six(t6))
        w("design_points_six.json", json.dumps({k: {f: v[f] for f in ("N", "r", "R", "memory", "kappa", "target",
                                                                         "closed_form_T2", "des_T2", "des_lo", "des_hi",
                                                                         "des_T2_median")}
                                                for k, v in t6.items()}, indent=1, default=float))
        rat6 = [v["ratio"] for v in t6.values() if v["ratio"]]
        macro("desSixRatioLo", f"{min(rat6):.1f}"); macro("desSixRatioHi", f"{max(rat6):.1f}")
        macro("desSixN", len(t6))
    made = figures(pf_best, pr, t_d)
    paper_macros(pf, gate, pr, pr_d, t, s1, s1_d, s2)
    w("des_numbers.tex", "\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in sorted(MACROS.items())))
    w("des_paragraph_draft.md", paragraph(pf, gate, pr, t, s1, s2))
    head = ["## Automatic checks\n", "| status | check | detail |", "|---|---|---|"]
    head += [f"| {s} | {n} | {d} |" for s, n, d in CHECKS]
    w("summary.md", "\n".join(LINES[:n_head] + head + [""] + LINES[n_head:]))
    print("\n".join(head))
    n_fail = sum(1 for s, _, _ in CHECKS if s == "FAIL")
    print(f"\nwrote summary.md, tables, macros, figures {made}, des_paragraph_draft.md -- {n_fail} FAIL")


if __name__ == "__main__":
    main()
