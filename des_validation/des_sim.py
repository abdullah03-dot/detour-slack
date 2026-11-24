"""
des_sim.py -- discrete-event, density-matrix simulation of one entanglement-
delivery attempt along a fixed path, built on the SeQUeNCe 1.3 kernel.

WHAT IS INDEPENDENT OF THE PAPER'S SIMULATOR (qunet_tqe.py)
(Scope of the claim: in Mode F the inputs -- Werner raw pairs from the paper's
clean_F, the twirl, lambda_s, the memory channels and kappa*D storage of the
delivered pair -- are the paper's model, so Mode F can disagree with the closed
form only through an error in the closed-form algebra or in the simulator: it
is an implementation cross-check.  Mode R changes the physics -- every stored
pair decoheres for its own storage time, purification takes time, messages
take time -- and is the independent test of the regime calls.)
  * Event scheduling: SeQUeNCe Timeline / Event / Process (picosecond clock).
  * Quantum states: SeQUeNCe QuantumManagerDensity (full density matrices).
  * Purification and swapping: gate circuits (CNOT, H, S, measurement, Pauli
    corrections) executed by SeQUeNCe's quantum manager, with measurement
    outcomes sampled from the Born rule.  No closed-form recurrence (BBPSSW
    formula, Werner composition, decay law) is used anywhere in this file.
    Optional gate noise (eps_gate): a two-qubit depolarizing channel after each
    bilateral CNOT, as in the paper's qunet_tqe.noisy_purify_round.
  * Fidelity is read from the density matrix, <Phi+|rho|Phi+>.  SeQUeNCe's own
    `Memory.fidelity` bookkeeping is NOT used: its circuit swapping protocol
    sets it to F1*F2*degradation (the product rule rejected in the paper's
    changelog) and its BBPSSW protocol updates it with the Werner recurrence.

WHAT THIS FILE IMPLEMENTS ITSELF (numpy, unit-tested in test_des_physics.py)
  * Noise channels as Kraus maps: per-qubit depolarizing / dephasing memory
    noise, the pairwise "half" channel, swap depolarization (lambda_s).
  * The isotropic twirl, as an average over the 24 single-qubit Cliffords
    C (x) conj(C) (a unitary 2-design, so the average is the exact twirl).
  * The protocol logic: nested BBPSSW/DEJMPS recurrence per link, swap-ASAP,
    classical-message delays, delivery and storage.
  * The paper's rule that a link whose raw fidelity is <= 1/2 is not purified
    ("no gain"; qunet_tqe.purify): such a link uses its first raw pair as is
    (skip_purify_below_half = True, the default).  Without it, a faulty link
    would be purified at a success probability near 1/2 per round, which costs
    time but cannot rescue it.

TIMING MODEL (all times in seconds; SeQUeNCe clock in ps)
  * Link i (length l_i km) heralds raw pairs as a Poisson process of rate R
    (gen_rate) starting at l_i / c (c = signal_speed, km/s).  A heralded pair
    was emitted l_i / c earlier and decoheres from its emission time.
  * Nested recurrence: two pairs of level j are purified as soon as both exist;
    the outcome reaches both ends after l_i / c; on success the kept pair joins
    level j+1.  A link is ready when it holds a level-r pair; it then stops
    generating.  Pairs left in lower levels are discarded.
  * Swap-ASAP: interior node j performs its Bell measurement as soon as both of
    its links are ready (Bell measurements at different nodes act on different
    qubits and commute).  The Pauli correction is applied at that moment
    (Pauli-frame bookkeeping); for the Pauli memory channels (deph2, depol2)
    this is exactly equivalent to correcting on message arrival.
  * Delivery: t_del = max_j [ t_swap_j + max(dist(n_j, n_0), dist(n_j, n_h)) / c ]
    (h = 1: the link-ready time).  The delivered pair is then stored for
    kappa * D seconds (D = path length) -- the paper's storage charge; set
    kappa = 0 to drop it.  The pair is used at t_used = t_del + kappa * D
    (the latency reported by the study).
  * Every stored pair decoheres for exactly the time it is stored (lazy
    update before every operation on it).

MODE F (faithful to the paper's assumptions) is the parameter point
  gen_rate = inf, signal_speed = inf, kappa = 0.002 s/km, scheme = bbpssw_twirl:
all link pairs exist and are purified at t = 0, and only the delivered pair
is stored (kappa * D).  MODE R (resolved) uses a finite gen_rate and the
physical signal speed 2e5 km/s.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Sequence

import numpy as np

from sequence.kernel.timeline import Timeline
from sequence.kernel.event import Event
from sequence.kernel.process import Process
from sequence.components.circuit import Circuit
from sequence.constants import DENSITY_MATRIX_FORMALISM

F_TH = 0.8                    # success threshold (same as the paper)
C_FIBER = 2.0e5               # km/s, light in fibre
PS = 1e12                     # SeQUeNCe clock unit: 1 ps

# ---------------------------------------------------------------------------
# Elementary matrices
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PHI_PLUS = np.array([1, 0, 0, 1], dtype=complex) / math.sqrt(2.0)
P_PHI = np.outer(PHI_PLUS, PHI_PLUS.conj())
SWAP2 = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)


def werner(F: float) -> np.ndarray:
    """Werner state F|Phi+><Phi+| + (1-F)/3 (I - |Phi+><Phi+|) (input specification)."""
    return F * P_PHI + (1.0 - F) / 3.0 * (np.eye(4) - P_PHI)


def fidelity(rho: np.ndarray) -> float:
    return float(np.real(PHI_PLUS.conj() @ rho @ PHI_PLUS))


def _cliffords() -> List[np.ndarray]:
    """The 24 single-qubit Cliffords modulo global phase (generated by H and S)."""
    H = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2.0)
    S = np.diag([1, 1j]).astype(complex)

    def canon(U):
        idx = np.flatnonzero(np.abs(U.ravel()) > 1e-9)[0]
        ph = U.ravel()[idx] / abs(U.ravel()[idx])
        V = U / ph
        # round and map -0.0 to +0.0 so that equal matrices give equal keys
        return tuple(np.round(V.real, 6).ravel() + 0.0) + tuple(np.round(V.imag, 6).ravel() + 0.0)

    found = {}
    frontier = [I2]
    while frontier:
        nxt = []
        for U in frontier:
            key = canon(U)
            if key in found:
                continue
            found[key] = U
            nxt += [H @ U, S @ U]
        frontier = nxt
    out = list(found.values())
    assert len(out) == 24, len(out)
    return out


_TWIRL_OPS = [np.kron(C, C.conj()) for C in _cliffords()]


def twirl(rho: np.ndarray) -> np.ndarray:
    """Isotropic (U x conj U) twirl of a two-qubit state."""
    acc = np.zeros((4, 4), dtype=complex)
    for U in _TWIRL_OPS:
        acc += U @ rho @ U.conj().T
    return acc / len(_TWIRL_OPS)


def depolarize_qubit(rho: np.ndarray, q: int, c: float) -> np.ndarray:
    """Single-qubit depolarizing channel on qubit q (0 or 1) of a 2-qubit state:
    rho -> c rho + (1-c) Tr_q(rho) (x) I/2, written in Kraus (Pauli) form:
    rho -> (1+3c)/4 rho + (1-c)/4 sum_P P rho P."""
    if c >= 1.0:
        return rho
    acc = (1.0 + 3.0 * c) / 4.0 * rho
    for P in (X, Y, Z):
        Op = np.kron(P, I2) if q == 0 else np.kron(I2, P)
        acc = acc + (1.0 - c) / 4.0 * (Op @ rho @ Op.conj().T)
    return acc


def dephase_qubit(rho: np.ndarray, q: int, c: float) -> np.ndarray:
    """Single-qubit dephasing: off-diagonal elements of qubit q multiplied by c.
    Kraus form rho -> (1+c)/2 rho + (1-c)/2 Z rho Z."""
    if c >= 1.0:
        return rho
    Op = np.kron(Z, I2) if q == 0 else np.kron(I2, Z)
    return (1.0 + c) / 2.0 * rho + (1.0 - c) / 2.0 * (Op @ rho @ Op.conj().T)


W_HALF = None  # set below


def half_channel(rho: np.ndarray, e: float) -> np.ndarray:
    """The paper's default 'half' memory model, F -> 1/2 + (F - 1/2) e, as a
    two-qubit replacement channel rho -> e rho + (1-e) sigma with sigma the
    Werner state of fidelity 1/2 (a fixed point, so the channel is a semigroup
    in t).  It is not a per-qubit channel; see README."""
    if e >= 1.0:
        return rho
    return e * rho + (1.0 - e) * W_HALF


W_HALF = werner(0.5)


def memory_channel(rho: np.ndarray, memory: str, dt: float, T2: float) -> np.ndarray:
    """Apply memory decoherence for a storage time dt to a stored pair (both of
    whose qubits sit in memory for dt)."""
    if dt <= 0.0 or not math.isfinite(T2):
        return rho
    c = math.exp(-dt / T2)
    if memory == "half":
        return half_channel(rho, c)
    if memory == "depol2":
        return depolarize_qubit(depolarize_qubit(rho, 0, c), 1, c)
    if memory == "deph2":
        return dephase_qubit(dephase_qubit(rho, 0, c), 1, c)
    raise ValueError(memory)


# ---------------------------------------------------------------------------
# Gate circuits (run by SeQUeNCe's density-matrix quantum manager)
# ---------------------------------------------------------------------------
def _purif_circuit(scheme: str) -> Circuit:
    """4-qubit circuit on [kept_A, meas_A, kept_B, meas_B]: (DEJMPS rotations),
    bilateral CNOT kept -> meas, Z measurement of both meas qubits."""
    c = Circuit(4)
    if scheme == "dejmps":
        # Rx(+pi/2) on A's qubits, Rx(-pi/2) on B's: Rx(+-pi/2) = H S(dg) H up to phase
        for q in (0, 1):
            c.h(q); c.s(q); c.h(q)
        for q in (2, 3):
            c.h(q); c.sdg(q); c.h(q)
    c.cx(0, 1)
    c.cx(2, 3)
    c.measure(1)
    c.measure(3)
    return c


def _purif_gates(scheme: str) -> Circuit:
    """Same as _purif_circuit without the measurement (used when gate noise is on)."""
    c = Circuit(4)
    if scheme == "dejmps":
        for q in (0, 1):
            c.h(q); c.s(q); c.h(q)
        for q in (2, 3):
            c.h(q); c.sdg(q); c.h(q)
    c.cx(0, 1)
    c.cx(2, 3)
    return c


PURIF_CIRCUITS = {s: _purif_circuit(s) for s in ("bbpssw_twirl", "bbpssw", "dejmps")}
PURIF_GATES = {s: _purif_gates(s) for s in ("bbpssw_twirl", "bbpssw", "dejmps")}
MEAS_13 = Circuit(4)
MEAS_13.measure(1)
MEAS_13.measure(3)

BSM = Circuit(2)          # Bell-state measurement on [q_left_link, q_right_link]
BSM.cx(0, 1)
BSM.h(0)
BSM.measure(0)
BSM.measure(1)

CORR = {}
for _mz in (0, 1):
    for _mx in (0, 1):
        _c = Circuit(1)
        if _mx:
            _c.x(0)
        if _mz:
            _c.z(0)
        CORR[(_mz, _mx)] = _c


# ---------------------------------------------------------------------------
# Parameters and result
# ---------------------------------------------------------------------------
@dataclass
class DESParams:
    rounds: int = 3
    scheme: str = "bbpssw_twirl"     # "bbpssw_twirl" | "bbpssw" (no twirl) | "dejmps"
    memory: str = "half"             # "half" | "deph2" | "depol2"
    T2: float = 0.1                  # s
    kappa: float = 0.002             # s/km storage charge on the delivered pair
    gen_rate: float = math.inf       # heralded raw pairs per second per link
    signal_speed: float = math.inf   # km/s for heralds and classical messages
    lam_swap: float = 0.995          # depolarization of the Bell measurement
    eps_gate: float = 0.0            # two-qubit depolarizing after each bilateral CNOT
    skip_purify_below_half: bool = True   # paper's rule: no purification if raw F <= 1/2
    max_raw_per_link: int = 100000   # safety stop (never reached in practice)

    def validate(self):
        assert 0.0 <= self.eps_gate < 1.0
        assert self.scheme in PURIF_CIRCUITS, self.scheme
        assert self.memory in ("half", "deph2", "depol2"), self.memory
        assert self.rounds >= 0
        if not math.isfinite(self.gen_rate):
            # an infinite rate with delayed messages would generate unboundedly
            # many pairs at t = 0 while waiting for purification outcomes
            assert not math.isfinite(self.signal_speed), \
                "gen_rate = inf requires signal_speed = inf (Mode F)"


def _apply(p: DESParams, kw: dict) -> DESParams:
    for k, v in kw.items():
        if not hasattr(p, k):
            raise TypeError(f"unknown DES parameter {k!r}")
        setattr(p, k, v)
    p.validate()
    return p


def mode_F(**kw) -> DESParams:
    return _apply(DESParams(gen_rate=math.inf, signal_speed=math.inf, kappa=0.002), kw)


def mode_R(gen_rate: float, kappa: float, **kw) -> DESParams:
    return _apply(DESParams(gen_rate=gen_rate, signal_speed=C_FIBER, kappa=kappa), kw)


@dataclass
class DeliveryResult:
    F: float                     # fidelity of the delivered pair (0.25 if a swap failed)
    success: bool                # F >= F_TH and no swap failure
    swap_failed: bool
    t_delivery: float            # s, end-to-end pair established and corrections known
    t_used: float                # s, t_delivery + kappa * D (time at which the pair is used)
    t_link_ready: List[float]
    t_swap: List[float]
    raw_pairs: List[int]         # raw pairs heralded per link
    storage: Optional[Dict[int, tuple]] = None   # key -> (t_birth, t_death), if logged


class _Seg:
    __slots__ = ("left", "right", "t_last")

    def __init__(self, left: int, right: int, t_last: float):
        self.left, self.right, self.t_last = left, right, t_last


# ---------------------------------------------------------------------------
# The simulation of one delivery attempt
# ---------------------------------------------------------------------------
class PathDES:
    """One delivery attempt on a path with h = len(link_F) links."""

    def __init__(self, link_F: Sequence[float], link_len: Sequence[float], P: DESParams,
                 seed: int, swap_fail: Optional[Sequence[bool]] = None, log_storage: bool = False):
        P.validate()
        self.P = P
        self.h = len(link_F)
        assert self.h >= 1 and len(link_len) == self.h
        self.F_raw = [float(f) for f in link_F]
        self.L = [float(x) for x in link_len]
        self.swap_fail = list(swap_fail) if swap_fail is not None else [False] * max(self.h - 1, 0)
        assert len(self.swap_fail) == max(self.h - 1, 0)
        self.rng = np.random.default_rng(seed)
        self.tl = Timeline(formalism=DENSITY_MATRIX_FORMALISM)
        self.qm = self.tl.quantum_manager
        self.circ = PURIF_CIRCUITS[P.scheme]
        self.circ_gates = PURIF_GATES[P.scheme]
        # paper's rule: a link with raw F <= 1/2 is not purified
        self.rounds_link = [0 if (P.skip_purify_below_half and f <= 0.5) else P.rounds for f in self.F_raw]
        self.do_twirl = P.scheme == "bbpssw_twirl"
        h = self.h
        self.buf: List[List[Optional[_Seg]]] = [[None] * (P.rounds + 1) for _ in range(h)]
        self.ready: List[Optional[_Seg]] = [None] * h
        self.t_ready = [math.nan] * h
        self.n_raw = [0] * h
        self.swapped = [False] * h          # index j = interior node j (1..h-1)
        self.t_swap = [math.nan] * h
        self.n_swaps = 0
        self.failed = False
        self.done = False
        self.seg_of: Dict[int, _Seg] = {}
        self.log = {} if log_storage else None
        # cumulative distance along the path (km) for message delays
        self.pos = np.concatenate([[0.0], np.cumsum(self.L)])
        self.result: Optional[DeliveryResult] = None
        self._pending: Dict[tuple, _Seg] = {}     # purified pairs awaiting the outcome message
        self._e2e: Optional[_Seg] = None

    # -- time helpers
    def now(self) -> float:
        return self.tl.now() / PS

    def _delay(self, km: float) -> float:
        c = self.P.signal_speed
        return 0.0 if not math.isfinite(c) else km / c

    def _schedule(self, t: float, method: str, args: list, priority: int):
        self.tl.schedule(Event(int(round(t * PS)), Process(self, method, args), priority))

    # -- state helpers (always through SeQUeNCe's quantum manager)
    def _get_pair(self, seg: _Seg) -> np.ndarray:
        st = self.qm.states[seg.left]
        keys = list(st.keys)
        assert sorted(keys) == sorted([seg.left, seg.right]), (keys, seg.left, seg.right)
        rho = np.array(st.state, dtype=complex)
        if keys[0] != seg.left:
            rho = SWAP2 @ rho @ SWAP2
        return rho

    def _set_pair(self, seg: _Seg, rho: np.ndarray):
        tr = float(np.real(np.trace(rho)))
        assert abs(tr - 1.0) < 1e-9, tr
        self.qm.set([seg.left, seg.right], rho)

    def _two_qubit_depol(self, k1: int, k2: int, eps: float):
        """rho -> (1-eps) rho + eps/16 sum_{a,b} (a x b) rho (a x b)^dag on qubits k1, k2
        of whatever joint state holds them (the paper's gate-noise model)."""
        st = self.qm.states[k1]
        keys = list(st.keys)
        assert k2 in keys
        n = len(keys)
        i1, i2 = keys.index(k1), keys.index(k2)
        rho = np.array(st.state, dtype=complex)
        acc = (1.0 - eps) * rho
        for a in (I2, X, Y, Z):
            for b in (I2, X, Y, Z):
                ops = [I2] * n
                ops[i1], ops[i2] = a, b
                V = ops[0]
                for o in ops[1:]:
                    V = np.kron(V, o)
                acc = acc + eps / 16.0 * (V @ rho @ V.conj().T)
        self.qm.set(keys, acc)

    def _age(self, seg: _Seg, t: float):
        dt = t - seg.t_last
        assert dt >= -1e-12, dt
        if dt > 0.0:
            rho = memory_channel(self._get_pair(seg), self.P.memory, dt, self.P.T2)
            self._set_pair(seg, rho)
        seg.t_last = max(seg.t_last, t)

    def _forget(self, *keys: int):
        t = self.now()
        for k in keys:
            if self.log is not None and k in self.log and self.log[k][1] is None:
                self.log[k] = (self.log[k][0], t)
            self.qm.states.pop(k, None)
            self.seg_of.pop(k, None)

    def _discard_seg(self, seg: _Seg):
        self._forget(seg.left, seg.right)

    # -- run
    def run(self) -> DeliveryResult:
        P = self.P
        for i in range(self.h):
            first = self._delay(self.L[i]) + self._exp()
            self._schedule(first, "_gen", [i], priority=1)
        self.tl.init()
        self.tl.run()
        assert self.result is not None, "simulation ended without a delivery"
        return self.result

    def _exp(self) -> float:
        R = self.P.gen_rate
        return 0.0 if not math.isfinite(R) else float(self.rng.exponential(1.0 / R))

    # -- generation
    def _gen(self, i: int):
        if self.ready[i] is not None or self.done:
            return
        self.n_raw[i] += 1
        if self.n_raw[i] > self.P.max_raw_per_link:
            raise RuntimeError(f"link {i}: more than {self.P.max_raw_per_link} raw pairs")
        t = self.now()
        t_emit = t - self._delay(self.L[i])
        ka, kb = self.qm.new(), self.qm.new()
        self.qm.set([ka, kb], werner(self.F_raw[i]))
        seg = _Seg(ka, kb, t_emit)
        if self.log is not None:
            self.log[ka] = (t_emit, None)
            self.log[kb] = (t_emit, None)
        self._arrive(i, seg, 0)
        if self.ready[i] is None:
            self._schedule(t + self._exp(), "_gen", [i], priority=1)

    def _arrive(self, i: int, seg: _Seg, level: int):
        if self.ready[i] is not None or self.done:
            self._discard_seg(seg)
            return
        if level == self.rounds_link[i]:
            self._link_ready(i, seg)
            return
        partner = self.buf[i][level]
        if partner is None:
            self.buf[i][level] = seg
        else:
            self.buf[i][level] = None
            self._purify(i, partner, seg, level)

    def _arrive_ev(self, i: int, seg_left: int, seg_right: int, level: int):
        self._arrive(i, self._pending.pop((seg_left, seg_right)), level)

    # -- purification
    def _purify(self, i: int, kept: _Seg, meas: _Seg, level: int):
        t = self.now()
        self._age(kept, t)
        self._age(meas, t)
        if self.do_twirl:
            self._set_pair(kept, twirl(self._get_pair(kept)))
            self._set_pair(meas, twirl(self._get_pair(meas)))
        keys = [kept.left, meas.left, kept.right, meas.right]
        if self.P.eps_gate > 0.0:
            # rotations + bilateral CNOTs, gate noise on each CNOT's two qubits, then measurement
            self.qm.run_circuit(self.circ_gates, keys)
            self._two_qubit_depol(kept.left, meas.left, self.P.eps_gate)
            self._two_qubit_depol(kept.right, meas.right, self.P.eps_gate)
            out = self.qm.run_circuit(MEAS_13, keys, float(self.rng.random()))
        else:
            out = self.qm.run_circuit(self.circ, keys, float(self.rng.random()))
        ok = out[meas.left] == out[meas.right]
        # measured qubits are now in a product state with the kept pair: trace them out
        self.qm.discard(meas.left)
        self.qm.remove(meas.left)
        self.qm.discard(meas.right)
        self.qm.remove(meas.right)
        if self.log is not None:
            self.log[meas.left] = (self.log[meas.left][0], t)
            self.log[meas.right] = (self.log[meas.right][0], t)
        if not ok:
            self._discard_seg(kept)
            return
        if self.do_twirl:
            self._set_pair(kept, twirl(self._get_pair(kept)))
        d = self._delay(self.L[i])
        self._pending[(kept.left, kept.right)] = kept
        self._schedule(t + d, "_arrive_ev", [i, kept.left, kept.right, level + 1], priority=0)

    # -- link ready, swapping, delivery
    def _link_ready(self, i: int, seg: _Seg):
        t = self.now()
        self.ready[i] = seg
        self.t_ready[i] = t
        self.seg_of[seg.left] = seg
        self.seg_of[seg.right] = seg
        for lv in range(self.P.rounds + 1):
            if self.buf[i][lv] is not None:
                self._discard_seg(self.buf[i][lv])
                self.buf[i][lv] = None
        if self.h == 1:
            self._deliver_at(t)
            return
        for j in (i, i + 1):                      # interior nodes adjacent to link i
            if 1 <= j <= self.h - 1 and not self.swapped[j] \
                    and self.ready[j - 1] is not None and self.ready[j] is not None:
                self._swap(j)

    def _swap(self, j: int):
        t = self.now()
        self.swapped[j] = True
        self.t_swap[j] = t
        q1 = self.ready[j - 1].right            # qubit at node j from link j-1
        q2 = self.ready[j].left                 # qubit at node j from link j
        segA, segB = self.seg_of[q1], self.seg_of[q2]
        assert segA.right == q1 and segB.left == q2 and segA is not segB
        self._age(segA, t)
        self._age(segB, t)
        # Bell-measurement depolarization lambda_s on the measured qubit of segA
        self._set_pair(segA, depolarize_qubit(self._get_pair(segA), 1, self.P.lam_swap))
        out = self.qm.run_circuit(BSM, [q1, q2], float(self.rng.random()))
        mz, mx = out[q1], out[q2]
        self.qm.discard(q1)
        self.qm.remove(q1)
        self.qm.discard(q2)
        self.qm.remove(q2)
        if self.log is not None:
            self.log[q1] = (self.log[q1][0], t)
            self.log[q2] = (self.log[q2][0], t)
        self.seg_of.pop(q1, None)
        self.seg_of.pop(q2, None)
        # Pauli-frame correction X^mx Z^mz on the far qubit of segB
        if mz or mx:
            self.qm.run_circuit(CORR[(mz, mx)], [segB.right])
        new = _Seg(segA.left, segB.right, t)
        self.seg_of[new.left] = new
        self.seg_of[new.right] = new
        if self.swap_fail[j - 1]:
            self.failed = True
        self.n_swaps += 1
        if self.n_swaps == self.h - 1:
            t_del = max(self.t_swap[k] + self._delay(max(self.pos[k] - self.pos[0], self.pos[-1] - self.pos[k]))
                        for k in range(1, self.h))
            self._e2e = new
            self._deliver_at(t_del)

    def _deliver_at(self, t_del: float):
        if self.h == 1:
            self._e2e = self.ready[0]
        self._schedule(t_del, "_deliver", [], priority=0)

    def _deliver(self):
        t = self.now()
        seg = self._e2e
        self._age(seg, t)
        D = float(self.pos[-1])
        extra = self.P.kappa * D
        if extra > 0:
            rho = memory_channel(self._get_pair(seg), self.P.memory, extra, self.P.T2)
            self._set_pair(seg, rho)
        F = fidelity(self._get_pair(seg))
        if self.log is not None:
            for k in (seg.left, seg.right):
                self.log[k] = (self.log[k][0], t + extra)
        self.done = True
        self.result = DeliveryResult(
            F=0.25 if self.failed else F, success=(not self.failed) and F >= F_TH,
            swap_failed=self.failed, t_delivery=t, t_used=t + extra,
            t_link_ready=list(self.t_ready), t_swap=[self.t_swap[k] for k in range(1, self.h)],
            raw_pairs=list(self.n_raw), storage=self.log)
        self.tl.stop()


def run_path(link_F, link_len, P: DESParams, seed: int, swap_fail=None, log_storage=False) -> DeliveryResult:
    return PathDES(link_F, link_len, P, seed, swap_fail=swap_fail, log_storage=log_storage).run()


def params_dict(P: DESParams) -> dict:
    d = asdict(P)
    return {k: (str(v) if isinstance(v, float) and not math.isfinite(v) else v) for k, v in d.items()}
