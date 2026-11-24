"""
smoke_test.py -- 1-minute check that the environment can run the study.

  python smoke_test.py

Checks: Python and package versions, that ../code (the paper's simulator) is
found, a handful of the physics tests, and that a 2-worker 'spawn' process pool
(the Windows default) runs DES jobs and gives the same answer as a serial run.
"""
from __future__ import annotations

import math
import sys
import time


def main():
    t0 = time.time()
    ok = True
    print(f"Python {sys.version.split()[0]} on {sys.platform}")
    if sys.version_info < (3, 12):
        print("FAIL: SeQUeNCe 1.3 needs Python >= 3.12 (use the qrf-des environment)")
        return 1
    try:
        from importlib.metadata import version
        import numpy, scipy, matplotlib  # noqa: F401
        v = version("sequence")
        print(f"sequence {v}, numpy {numpy.__version__}, scipy {scipy.__version__}, matplotlib {matplotlib.__version__}")
        if v != "1.3.0":
            print(f"WARNING: tested with sequence 1.3.0, found {v}. Run test_des_physics.py before trusting results.")
    except Exception as e:  # noqa: BLE001
        print(f"FAIL: missing package: {e}")
        return 1
    try:
        import qunet_bridge as B
        print(f"paper code found: {B.CODE}")
    except Exception as e:  # noqa: BLE001
        print(f"FAIL: cannot import the paper's code from ../code: {e}")
        return 1

    import test_des_physics as T
    for f in (T.test_twirl, T.test_memory_channels_vs_paper_decay, T.test_swap_pauli_corrections_all_outcomes,
              T.test_swap_werner_composition, T.test_modeR_storage_accounting_exact):
        try:
            f()
            print(f"PASS  {f.__name__}")
        except Exception as e:  # noqa: BLE001
            ok = False
            print(f"FAIL  {f.__name__}: {e}")

    # one Mode F delivery against the closed form (5x5, r = 3, one-hop detour)
    import des_sim as D
    Q = B.Q
    ph = B.physics(3, "half", 0.1)
    L = [1.0, math.sqrt(2), 1.0, math.sqrt(2), math.sqrt(2)]
    F = D.run_path([ph.clean_F(x) for x in L], L, D.mode_F(rounds=3), seed=1).F
    ref = B.closed_form_path(ph, L)
    good = abs(F - ref) < 1e-10
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  Mode F 5x5 r=3 detour: DES {F:.12f}  closed form {ref:.12f}")

    # spawn pool == serial
    import run_des as RD
    jobs = [{"key": f"smoke{i}", "i0": 0, "seeds": [11 + i, 12 + i], "link_F": [ph.clean_F(x) for x in L],
             "link_len": L, "P": {"mode": "R", "gen_rate": 1000.0, "kappa": 0.002, "rounds": 3,
                                  "memory": "deph2", "T2": 0.5}} for i in range(2)]
    par = sorted(RD.pool_map(RD.job_paths, jobs, 2, "smoke"))
    ser = sorted(RD.job_paths(j) for j in jobs)
    good = all(abs(a - b) < 1e-15 for p, s in zip(par, ser) for a, b in zip(p[2], s[2]))
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  2-worker spawn pool reproduces the serial run (Mode R, seeded)")

    print(f"\n{'SMOKE TEST PASSED' if ok else 'SMOKE TEST FAILED'} in {time.time() - t0:.0f}s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
