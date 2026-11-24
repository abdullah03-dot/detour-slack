"""paper_extra.py -- a few further DES macros quoted in the paper, computed from the
design-point files that analyze.py writes. Run after analyze.py:

    python des_validation/paper_extra.py [--results des_validation/results_des]

Writes <results>/des_extra.tex.
"""
from __future__ import annotations

import argparse
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(HERE, "results_des"))
    a = ap.parse_args()
    rd = a.results
    bb = json.load(open(os.path.join(rd, "design_points.json")))
    dj = json.load(open(os.path.join(rd, "design_points_dejmps.json")))
    six = json.load(open(os.path.join(rd, "design_points_six.json")))
    M = {}

    # DES threshold inside [optimistic closed form, pessimistic all-links-wait accounting]?
    # The pessimistic accounting (every channel, kappa = 2 ms/km) comes from code/design_law.py.
    PLA = json.load(open(os.path.join(HERE, "..", "results_tqe", "design_law.json")))["latency_all_links"]
    ins = lambda v, N: v["closed_form_T2"] <= v["des_T2"] <= PLA[f"{v['memory']}|{int(v['R'])}|{N}"]["by_r"][v["r"] - 1]
    five = [v for v in bb.values() if v["target"] == 1 and v["kappa"] > 0 and v["des_T2"]]
    M["desBrFiveN"] = len(five)
    M["desBrFiveIn"] = sum(ins(v, 5) for v in five)
    d6 = [v for v in six.values() if v["target"] == 1 and v["des_T2"]]
    M["desBrSixN"] = len(d6)
    M["desBrSixIn"] = sum(ins(v, v["N"]) for v in d6)
    # the pessimistic accounting's own choice of rounds on 5x5 (for comparison with the simulation)
    for mem in ("half", "deph2", "depol2"):
        for R_ in (1000, 100):
            br = PLA[f"{mem}|{R_}|5"]["by_r"]
            M[f"pessRstar{ {'half': 'Half', 'deph2': 'Deph', 'depol2': 'Depol'}[mem] }{ {1000: 'Onek', 100: 'Hundred'}[R_] }"] = 1 + br.index(min(br))

    # DEJMPS vs BBPSSW in the same simulation: r = 3, same rate and channel, kappa = 2 ms/km, s >= 1
    key = lambda v: (v["r"], v["R"], v["memory"])
    bbk = {key(v): v for v in bb.values() if v["target"] == 1 and v["kappa"] > 0}
    rat = []
    for v in dj.values():
        if v["target"] == 1 and v["kappa"] > 0 and v["r"] == 3 and v["des_T2"]:
            b = bbk[key(v)]
            rat.append((b["des_T2"] / v["des_T2"], v, b))
    rat.sort(key=lambda x: x[0])
    M["desDJBBLo"] = f"{rat[0][0]:.1f}"
    M["desDJBBHi"] = f"{rat[-1][0]:.1f}"
    M["desDJBBN"] = len(rat)
    hd = [x for x in rat if x[1]["R"] == 100 and x[1]["memory"] == "depol2"][0]
    M["desDJHundredDepol"] = f"{hd[1]['des_T2']:.2f}"
    M["desBBHundredDepol"] = f"{hd[2]['des_T2']:.2f}"

    # thresholds without the kappa*D term, per-qubit channels only (5x5, BBPSSW, s >= 1)
    nk = [v["des_T2"] / v["closed_form_T2"] for v in bb.values()
          if v["target"] == 1 and v["kappa"] == 0 and v["memory"] != "half" and v["des_T2"]]
    M["desTratioNoKapPhysLo"] = f"{min(nk):.1f}"
    M["desTratioNoKapPhysHi"] = f"{max(nk):.1f}"
    T = json.load(open(os.path.join(rd, "T.json")))
    M["desNdeliv"] = T["meta"]["n"]["T"]

    # softness sweep (block S3): headroom across the s = 0 -> 1 boundary
    s3p = os.path.join(rd, "S3.json")
    if os.path.exists(s3p):
        from scipy import stats as sst
        S3 = json.load(open(s3p))
        rows = []
        for key, run in sorted(S3["runs"].items(), key=lambda kv: kv[1]["T2"]):
            tr = run["trials"]
            ob = np.asarray(tr["oracle"]["ok_des"], bool)
            bf = np.asarray(tr["shortest"]["ok_des"], bool)
            pg = np.asarray(tr["full"]["ok_des"], bool)
            w, l = int(np.sum(ob & ~bf)), int(np.sum(~ob & bf))
            pval = 1.0 if w + l == 0 else float(min(1.0, 2 * sst.binom.cdf(min(w, l), w + l, 0.5)))
            pr = run["profile"]
            f1 = pr["1"]
            rows.append({"T2": run["T2"], "factor": run["factor"], "m1": f1["mean"] - 0.8,
                         "m0": pr["0"]["mean"] - 0.8, "pass1": f1["p_pass"],
                         "sd1": f1["se"] * math.sqrt(f1["n"]), "bfs": 100 * bf.mean(), "pg": 100 * pg.mean(),
                         "orc": 100 * ob.mean(), "head": 100 * (ob.mean() - bf.mean()), "p": pval, "n": len(ob)})
        below = [r for r in rows if r["m1"] < 0 and r["m0"] >= 0]
        M["desSoftN"] = f"{rows[0]['n']:,}".replace(",", "{,}")
        M["softHeadMin"] = f"{min(r['head'] for r in below):.2f}"
        M["softHeadMax"] = f"{max(r['head'] for r in below):.2f}"
        M["softBand"] = f"{2 * np.mean([r['sd1'] for r in rows]):.3f}"
        far = min(rows, key=lambda r: r["m1"])
        near = max(below, key=lambda r: r["m1"])
        M["softFarMargin"] = f"{-far['m1']:.3f}"
        M["softFarHead"] = f"{far['head']:.2f}"
        M["softFarP"] = f"{far['p']:.2f}"
        M["softNearMargin"] = f"{-near['m1']:.4f}"
        above = [r for r in rows if r["m1"] >= 0]
        M["softAboveHeadHi"] = f"{max(r['head'] for r in above):.2f}"
        json.dump(rows, open(os.path.join(rd, "soft_sweep.json"), "w"), indent=1)

        def pe(p):
            if p >= 1e-3:
                return f"{p:.3f}"
            e = int(math.floor(math.log10(max(p, 1e-300))))
            return f"$<10^{{{e + 1}}}$" if p > 0 else "$<10^{-300}$"
        L = [r"\begin{table}[!t]", r"\centering",
             r"\caption{Softness of the regime boundary in the discrete-event simulation (Mode~R, 5$\times$5, BBPSSW, "
             r"$r=3$, $R=1$\,kHz, $\kappa=2$\,ms/km, per-qubit depolarizing memories, $n=" + M["desSoftN"] + r"$ paired "
             r"trials per row; Or.$^\ast$: oracle of the scalar model): mean clean-path margins $F^*(k)-F_{\rm th}$ from the simulation, the fraction of "
             r"one-hop-detour deliveries that pass, and the routing headroom (oracle$^\ast-$BFS) with exact "
             r"\textnormal{McNemar} $P$.}", r"\label{tab:dessoft}", r"\scriptsize", r"\setlength{\tabcolsep}{1.5pt}",
             r"\begin{tabular}{@{}cccccccc@{}}", r"\toprule",
             r"$T_2$ (s) & $F^*(0)-F_{\rm th}$ & $F^*(1)-F_{\rm th}$ & pass\,\% & BFS & PG & Or.$^\ast$ & Or.$^\ast-$BFS (pp), $P$ \\",
             r"\midrule"]
        for r in rows:
            L.append(f"{r['T2']:.3g} & {r['m0']:+.4f} & {r['m1']:+.4f} & {100 * r['pass1']:.0f} & {r['bfs']:.2f} & "
                     f"{r['pg']:.2f} & {r['orc']:.2f} & {r['head']:+.2f}, {pe(r['p'])} \\\\")
        L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
        open(os.path.join(rd, "tab_des_soft.tex"), "w").write("% generated by paper_extra.py -- do not edit\n" + "\n".join(L) + "\n")
        for r in rows:
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})

    with open(os.path.join(rd, "des_extra.tex"), "w") as f:
        f.write("% generated by paper_extra.py -- do not edit\n")
        for k, v in M.items():
            f.write("\\newcommand{\\%s}{%s}\n" % (k, v))
    for k, v in M.items():
        print(k, v)


if __name__ == "__main__":
    main()
