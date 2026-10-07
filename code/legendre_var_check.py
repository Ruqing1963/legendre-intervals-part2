"""
legendre_var_check.py -- independent, exact verification of the Legendre-interval variance
anomaly in characteristic 2 (Open Problem 8.2 of the Chen paper).

Dependencies: Python 3 + numpy only (numpy is used purely for array arithmetic; every
number produced is an exact integer or an exact rational derived from integer counts).

Mathematical setup (all exact, no floating-point estimates):
  * n = 2d, l = d-1.  Hayes group G = G_l = {1 + c1 u + ... + c_l u^l} in F_q[u]/(u^{l+1}).
  * A monic F of degree n maps to F*(u) = u^n F(1/u) mod u^{l+1}, i.e. to its top l
    coefficients (a_{n-1}, ..., a_{n-l}).  The fibre of a class g has q^{n-l} = q^{d+1} elements.
  * Legendre interval I_f = {f^2 + s : deg s <= d}  ==  fibre over the class of f^2, and these
    classes are exactly the subgroup H = G^2 of squares (a_{n-j} = 0 for odd j <= l).
  * Psi(g) = sum_{F in fibre(g)} Lambda(F) = #{alpha in F_{q^n} : charpoly(alpha) in fibre(g)}.
    (A root alpha of P^k, deg P = n/k, has charpoly P^k; there are n/k = Lambda(P^k) such roots.)
    So Psi is obtained by running over ALL alpha in F_{q^n} and recording the top l coefficients
    of prod_i (t - alpha^{q^i}); the recurrence only needs the top l coefficients at each step.
  * mean = q^n / |G| = q^{d+1}.   phi(g) = Psi(g) - q^{d+1}.
      Var_G = (1/|G|)   sum_{g in G} phi(g)^2      (full Keating-Rudnick variance)
      Var_H = (1/|H|)   sum_{g in H} phi(g)^2      (second moment of the Legendre intervals about
                                                   the SAME global mean; this is what the paper
                                                   tabulates as var_Lambda_over_q^(d+1))
  * Quadratic characters.  For g in G let p_j(g) be the j-th power sum of the "roots" of g
    (Newton's identities, all signs + in char 2).  p_j is additive on G and p_j(g^2) = 0 for
    odd j (p_{2j} = p_j^2).  The map g -> (p_j(g))_{j odd <= l} is a surjection G -> F_q^m with
    kernel exactly G^2, so the q^m quadratic characters are
          eps_lambda(g) = (-1)^{Tr_{F_q/F_2}( sum_{j odd} lambda_j p_j(g) )}.
    On a root alpha this is the Artin-Schreier character  psi_0( sum_j lambda_j alpha^j ) of F_{q^n}.
  * rho(eps) = sum_g phi(g)^2 eps(g) / sum_g phi(g)^2.   Exact identity (Parseval over G/H):
          Var_H / Var_G = 1 + sum_{eps != 1} rho(eps).
    (This equals the paper's definition Re sum_chi psi(chi) conj(psi(chi eps)) / sum |psi(chi)|^2.)

Usage:  python legendre_var_check.py            # default case list
        python legendre_var_check.py 8,4 2,6    # explicit (q,d) list
"""
from __future__ import annotations

import csv
import itertools
import os
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# primitive polynomials over F_2 (bit N set), verified at runtime
PRIMITIVE = {
    4: 0b10011, 6: 0b1000011, 8: 0b100011101, 10: 0b10000001001, 12: 0b1000001010011,
    14: 0b100000000101011, 16: 0b10001000000001011, 18: 0b1000000000010000001,
    20: 0b100000000000000001001, 22: 0b10000000000000000000011,
    24: 0b1000000000000000010000111, 26: 0b100000000000000000000100111,
}


# ----------------------------------------------------------------------------- GF(2^N) tables
def _clmul_reduce(arr: np.ndarray, s: int, N: int, poly: int) -> np.ndarray:
    """(arr * s) in GF(2^N), arr an array of field elements (ints), s a scalar field element."""
    acc = np.zeros(arr.shape, dtype=np.uint64)
    a = arr.astype(np.uint64)
    b = 0
    while s >> b:
        if (s >> b) & 1:
            acc ^= a << np.uint64(b)
        b += 1
    for bit in range(2 * N - 2, N - 1, -1):
        mask = (acc >> np.uint64(bit)) & np.uint64(1)
        acc ^= mask * np.uint64(poly << (bit - N))
    return acc.astype(np.uint32)


class GF:
    """GF(2^N) with exp/log tables; elements are ints (F_2-polynomial basis)."""

    def __init__(self, N: int):
        self.N = N
        self.size = 1 << N
        self.order = self.size - 1
        poly = PRIMITIVE[N]
        exp = np.zeros(self.order, dtype=np.uint32)
        exp[0] = 1
        L = 1
        while L < self.order:
            s = int(_clmul_reduce(np.array([exp[L - 1]], np.uint32), 2, N, poly)[0])  # x^L
            take = min(L, self.order - L)
            exp[L:L + take] = _clmul_reduce(exp[:take], s, N, poly)
            L += take
        # primitivity check: x^(2^N-1) = 1 and no earlier return to 1
        last = int(_clmul_reduce(np.array([exp[-1]], np.uint32), 2, N, poly)[0])
        assert last == 1 and np.count_nonzero(exp == 1) == 1, "polynomial is not primitive"
        self.exp = exp
        log = np.empty(self.size, dtype=np.int64)
        log[exp] = np.arange(self.order, dtype=np.int64)
        log[0] = -1
        self.log = log
        self.poly = poly

    def mul(self, a: int, b: int) -> int:
        if a == 0 or b == 0:
            return 0
        return int(self.exp[(int(self.log[a]) + int(self.log[b])) % self.order])

    def pow(self, a: int, e: int) -> int:
        if a == 0:
            return 0 if e else 1
        return int(self.exp[(int(self.log[a]) * e) % self.order])

    def mul_vec_by_logs(self, lb: np.ndarray, x: np.ndarray) -> np.ndarray:
        """beta * x where beta is given by its log array lb and x is an array of elements."""
        out = np.zeros(x.shape, dtype=np.uint32)
        mask = x != 0
        lx = self.log[x[mask]]
        out[mask] = self.exp[(lb[mask] + lx) % self.order]
        return out

    def subfield(self, q: int) -> list[int]:
        """Elements of F_q inside GF(2^N) (0 first, then powers of a generator)."""
        step = self.order // (q - 1)
        return [0] + [int(self.exp[j * step]) for j in range(q - 1)]


# ----------------------------------------------------------------------------- Psi over G
def psi_classes(q: int, d: int, gf: GF, track_all: bool = False, chunk: int = 1 << 21):
    """Return (Psi array indexed by class key, extra) where key = sum_j idx(a_{n-j}) q^{j-1}.
    If track_all, also return irreducible counts binned by (a7,a6,a5,a3,a1) (used for d=4 slices)."""
    n, l = 2 * d, d - 1
    N = gf.N
    assert gf.size == q ** n
    fq = gf.subfield(q)
    fqidx = np.zeros(gf.size, dtype=np.int64)
    for i, e in enumerate(fq):
        fqidx[e] = i
    ntrack = n if track_all else l
    psi = np.zeros(q ** l, dtype=np.int64)
    slices = np.zeros(q ** 5, dtype=np.int64) if track_all else None
    qpow = [q ** i for i in range(n)]
    for start in range(0, gf.order, chunk):
        m = np.arange(start, min(start + chunk, gf.order), dtype=np.int64)
        p = [None] * (ntrack + 1)          # p[j] = coefficient of t^{deg-j} ("top-indexed")
        for j in range(1, ntrack + 1):
            p[j] = np.zeros(m.shape, dtype=np.uint32)
        for i in range(n):
            lb = (m * qpow[i]) % gf.order    # log of beta_i = alpha^{q^i}
            beta = gf.exp[lb]
            for j in range(min(ntrack, i + 1), 1, -1):
                p[j] ^= gf.mul_vec_by_logs(lb, p[j - 1])
            p[1] ^= beta
        key = np.zeros(m.shape, dtype=np.int64)
        for j in range(1, l + 1):
            key += fqidx[p[j]] * qpow[j - 1]
        psi += np.bincount(key, minlength=q ** l)
        if track_all:
            # elements of degree exactly n: not in any proper subfield
            deg_mask = np.ones(m.shape, dtype=bool)
            for pr in _prime_factors(n):
                sub = (gf.order) // (q ** (n // pr) - 1)
                deg_mask &= (m % sub) != 0
            # (a7,a6,a5,a3,a1) = (p1,p2,p3,p5,p7) for n = 8
            k2 = (fqidx[p[1]] + fqidx[p[2]] * q + fqidx[p[3]] * q ** 2
                  + fqidx[p[5]] * q ** 3 + fqidx[p[7]] * q ** 4)
            slices += np.bincount(k2[deg_mask], minlength=q ** 5)
    psi[0] += 1  # alpha = 0, charpoly t^n, identity class
    assert psi.sum() == q ** n
    if track_all:
        assert slices.sum() % n == 0
        slices //= n
    return psi, slices


def _prime_factors(n: int) -> list[int]:
    out, f = [], 2
    while f * f <= n:
        if n % f == 0:
            out.append(f)
            while n % f == 0:
                n //= f
        f += 1
    if n > 1:
        out.append(n)
    return out


# ----------------------------------------------------------------------------- group structure
class Hayes:
    def __init__(self, q: int, l: int, gf: GF):
        self.q, self.l, self.gf = q, l, gf
        self.fq = gf.subfield(q)
        self.idx = {e: i for i, e in enumerate(self.fq)}
        self.k = q.bit_length() - 1
        self.elements = list(itertools.product(range(q), repeat=l))  # tuples of F_q indices, c1..cl

    def key(self, g) -> int:
        return sum(ci * self.q ** j for j, ci in enumerate(g))

    def mul(self, g, h):
        fq, gf = self.fq, self.gf
        out = []
        for k in range(1, self.l + 1):
            s = fq[g[k - 1]] ^ fq[h[k - 1]]
            for i in range(1, k):
                s ^= gf.mul(fq[g[i - 1]], fq[h[k - i - 1]])
            out.append(self.idx[s])
        return tuple(out)

    def square(self, g):
        out = []
        for k in range(1, self.l + 1):
            out.append(self.idx[self.gf.mul(self.fq[g[k // 2 - 1]], self.fq[g[k // 2 - 1]])] if k % 2 == 0 else 0)
        return tuple(out)

    def power_sums(self, g):
        """p_1..p_l via Newton (char 2: all signs +), returned as F_q elements (ints)."""
        fq, gf = self.fq, self.gf
        e = [None] + [fq[c] for c in g]
        p = [None] * (self.l + 1)
        for k in range(1, self.l + 1):
            s = e[k] if k % 2 == 1 else 0
            for i in range(1, k):
                s ^= gf.mul(e[i], p[k - i])
            p[k] = s
        return p

    def trace2(self, x: int) -> int:
        """Tr_{F_q/F_2}(x) computed inside GF(2^N)."""
        s, y = 0, x
        for _ in range(self.k):
            s ^= y
            y = self.gf.mul(y, y)
        assert s in (0, 1)
        return s

    def quadratic_characters(self):
        """All q^m characters eps_lambda as arrays over G (indexed by key), plus their labels."""
        odd = [j for j in range(1, self.l + 1) if j % 2 == 1]
        P = {g: self.power_sums(g) for g in self.elements}
        chars = []
        for lam in itertools.product(self.fq, repeat=len(odd)):
            vals = np.empty(self.q ** self.l, dtype=np.int64)
            for g in self.elements:
                s = 0
                for lj, j in zip(lam, odd):
                    s ^= self.gf.mul(lj, P[g][j])
                vals[self.key(g)] = 1 - 2 * self.trace2(s)
            chars.append((tuple(self.idx[x] for x in lam), odd, vals))
        return chars


# ----------------------------------------------------------------------------- analysis
def analyse(q: int, d: int, verbose: bool = True):
    n, l = 2 * d, d - 1
    k = q.bit_length() - 1
    N = k * n
    t0 = time.time()
    gf = GF(N)
    track = (d == 4)
    psi, slices = psi_classes(q, d, gf, track_all=track)
    hay = Hayes(q, l, gf)
    G = hay.elements
    keys = np.array([hay.key(g) for g in G])
    squares = sorted({hay.square(g) for g in G})
    Hkeys = np.array([hay.key(g) for g in squares])
    mean = q ** (d + 1)
    phi = psi.astype(object) - mean           # exact Python ints
    phi2 = np.array([int(x) ** 2 for x in phi], dtype=object)
    S_G = sum(phi2)
    S_H = sum(phi2[Hkeys])
    varG = Fraction(S_G, q ** l)
    varH = Fraction(S_H, len(squares))
    meanH = Fraction(int(sum(psi[Hkeys])), len(squares))
    varH_own = Fraction(sum(int(psi[h] - meanH) ** 2 for h in Hkeys), len(squares))  # about own mean
    ratio = varH / varG if varG else None

    # quadratic characters and rho
    chars = hay.quadratic_characters()
    rows = []
    rho_sum = Fraction(0)
    for lam, odd, vals in chars:
        if not any(lam):
            # sanity: trivial character
            continue
        num = sum(int(a) * int(b) for a, b in zip(phi2, vals))
        rho = Fraction(num, S_G)
        rho_sum += rho
        rows.append((lam, odd, rho))
    # sanity checks on characters
    sq_set = set(squares)
    for lam, odd, vals in chars:
        assert all(vals[hay.key(h)] == 1 for h in squares), "eps not trivial on G^2"
    rng = np.random.default_rng(0)
    for _ in range(50):
        g, h = G[rng.integers(len(G))], G[rng.integers(len(G))]
        gh = hay.mul(g, h)
        for lam, odd, vals in chars[:min(8, len(chars))]:
            assert vals[hay.key(gh)] == vals[hay.key(g)] * vals[hay.key(h)], "eps not multiplicative"
    assert len({tuple(v) for _, _, v in chars}) == len(chars) == q ** len(chars[0][1])
    assert 1 + rho_sum == ratio, (1 + rho_sum, ratio)

    # Fourier transform of Psi restricted to H (H = {1 + c u^2 + ...} is an F_2-space; for d = 4 it is F_q)
    out = {
        "q": q, "d": d, "l": l, "N": N, "|G|": q ** l, "|H|": len(squares),
        "Psi_min": int(psi.min()), "Psi_max": int(psi.max()),
        "VarG/q^(d+1)": varG / mean, "VarH/q^(d+1)": varH / mean,
        "VarH_own/q^(d+1)": varH_own / mean, "meanH/q^(d+1)": meanH / mean,
        "ratio": ratio, "KR": d - 2, "rho_rows": rows, "seconds": round(time.time() - t0, 1),
        "psi_H": [(squares[i], int(psi[h])) for i, h in enumerate(Hkeys)],
    }
    if verbose:
        print(f"\n=== q={q}, d={d}  (n={n}, l={l}, field GF(2^{N}))  [{out['seconds']}s]")
        print(f"  |G|={q**l}  |H|={len(squares)}  mean=q^(d+1)={mean}  Psi in [{psi.min()}, {psi.max()}]")
        print(f"  Var_G/q^(d+1) = {varG/mean} = {float(varG/mean):.6f}   (KR limit {d-2})")
        print(f"  Var_H/q^(d+1) = {varH/mean} = {float(varH/mean):.6f}   (about the global mean)")
        print(f"  Var_H own-mean/q^(d+1) = {varH_own/mean} = {float(varH_own/mean):.6f},  mean_H/q^(d+1) = {meanH/mean}")
        print(f"  ratio Var_H/Var_G = {ratio} = {float(ratio):.6f};   1 + sum rho = {1+rho_sum}  [identity OK]")
        if d == 4:
            print(f"  q(q-1) = {q*(q-1)}  ->  Var_H/q^5 == q(q-1): {varH/mean == q*(q-1)}")
            print(f"  q^2/3  = {Fraction(q*q,3)}  ->  ratio == q^2/3: {ratio == Fraction(q*q,3)}")
        if len(squares) <= 16:
            print("  Psi on H (class c1..cl as F_q indices -> Psi, Psi - mean):")
            for g, v in out["psi_H"]:
                print(f"     {g}: {v:>8d}   {v-mean:>+8d}")
        odd = rows[0][1]
        print(f"  quadratic characters eps_lambda, lambda indexed by odd power sums p_j, j in {odd}:")
        by_support = {}
        for lam, _, rho in rows:
            supp = tuple(j for lj, j in zip(lam, odd) if lj)
            by_support.setdefault(supp, []).append(rho)
        for supp, rs in sorted(by_support.items()):
            tot = sum(rs)
            print(f"     support {supp}: {len(rs):>3d} chars, sum rho = {tot} = {float(tot):+.4f}, "
                  f"values {sorted({float(r) for r in rs})[:6]}")
    out["slices"] = slices
    out["psi"] = psi
    out["hayes"] = hay
    return out


def slice_report(res, show=True):
    """d = 4 only: irreducible counts in each b-slice of each Legendre interval."""
    q, slices, hay = res["q"], res["slices"], res["hayes"]
    print(f"\n  --- d=4 slice structure for q={q}: F = h^2 + t b^2, b = b1 t + b0; a3 = b1^2, a1 = b0^2 ---")
    fq = hay.fq
    table = {}
    for c in range(q):
        for b1 in range(q):
            for b0 in range(q):
                key = 0 + c * q + 0 * q ** 2 + b1 * q ** 3 + b0 * q ** 4
                table[(c, b1, b0)] = int(slices[key])
    for c in range(q):
        n0 = table[(c, 0, 0)]
        nz = [table[(c, b1, b0)] for b1 in range(q) for b0 in range(q) if (b1, b0) != (0, 0)]
        vals = sorted(set(nz))
        print(f"   c=idx{c}: b=0 slice irreducibles={n0}; b!=0 slices: total irreducibles={sum(nz)} "
              f"over {len(nz)} slices, distinct per-slice counts {vals}")
    # E(c): prime-power contribution
    for c in range(q):
        Psi_c = int(res["psi"][c * q])          # key of class (0, c, 0) is c*q
        n_irr = table[(c, 0, 0)] + sum(table[(c, b1, b0)] for b1 in range(q) for b0 in range(q) if (b1, b0) != (0, 0))
        E = Psi_c - 8 * n_irr
        print(f"   c=idx{c}: Psi={Psi_c}, 8*N_irr={8*n_irr}, prime-power part E={E} (q^3={q**3})")
    return table


def _abs_trace_parity(gf: GF, x: np.ndarray) -> np.ndarray:
    """Tr_{GF(2^N)/F_2}(x) in {0,1} for an array x (polynomial basis, so Tr is a bit-mask parity)."""
    mask = 0
    for i in range(gf.N):
        e, s = 1 << i, 0
        y = 1 << i
        for _ in range(gf.N):
            s ^= y
            y = gf.mul(y, y)
        assert s in (0, 1)
        if s:
            mask |= e
    v = (x & np.uint32(mask)).astype(np.uint64)
    for sh in (32, 16, 8, 4, 2, 1):
        v ^= v >> np.uint64(sh)
    return (v & np.uint64(1)).astype(np.int64)


def d4_exponential_sums(res):
    """d = 4: the two numbers A = #{alpha : Tr a = 0, Tr a^3 = 0} and Psi(id), their reduction to
    Artin-Schreier / Gold-type sums over F_{q^8}, and the per-character rho structure predicted by
    the theorem 'phi is supported on {c1 = 0}'."""
    q, psi, hay = res["q"], res["psi"], res["hayes"]
    gf = hay.gf
    fq = hay.fq
    N, order = gf.N, gf.order
    print(f"\n  --- d=4 exponential-sum reduction for q={q} (field GF(2^{N})) ---")
    A = sum(int(psi[c * q]) for c in range(q))                     # sum over H of Psi
    Psi0 = int(psi[0])
    print(f"   A = sum_c Psi(0,c,0) = #{{alpha: Tr a=0, Tr a^3=0}} = {A}   (q^6 = {q**6}, A - q^6 = {A - q**6})")
    print(f"   Psi(id) = #{{alpha: e1=e2=e3=0}} = {Psi0}   (q^5 - q^4 + q^3 = {q**5 - q**4 + q**3})")
    for c in range(1, q):
        print(f"   Psi(0,c,0) for c=idx{c}: {int(psi[c*q])}  = q^5 + ({int(psi[c*q]) - q**5})   [q^3 = {q**3}]")
    # c1 != 0 classes
    off = [int(psi[hay.key(g)]) - q ** 5 for g in hay.elements if g[0] != 0]
    print(f"   classes with c1 != 0: {len(off)} classes, max |Psi - q^5| = {max(abs(x) for x in off)}  (theorem T1 <=> 0)")
    # Gold-type sums W(nu) = sum_beta psi_0(nu (beta^{q+2} + beta^{2q+1}))
    m = np.arange(order, dtype=np.int64)
    t1 = gf.exp[((q + 2) * m) % order]
    t2 = gf.exp[((2 * q + 1) * m) % order]
    base = t1 ^ t2
    Ws = {}
    for nu in fq[1:]:
        lnu = int(gf.log[nu])
        x = np.zeros(order, dtype=np.uint32)
        nz = base != 0
        x[nz] = gf.exp[(gf.log[base[nz]] + lnu) % order]
        par = _abs_trace_parity(gf, x)
        W = int(order - 2 * par.sum()) + 1                      # include beta = 0 (term +1)
        Ws[nu] = W
    Wsum = sum(Ws.values())
    print(f"   W(nu) = sum_beta (-1)^Tr(nu(b^(q+2)+b^(2q+1))), nu in F_q^*: {sorted(Ws.values())}")
    print(f"   check A = q^6 + q^-2 sum_nu W(nu): {q**6 + Fraction(Wsum, q*q)} == {A}: {q**6 + Fraction(Wsum, q*q) == A}")
    # cubic Artin-Schreier sums S(l1,l3) = sum_x psi_0(l3 x^3 + l1 x) over F_{q^8}
    x3 = gf.exp[(3 * m) % order]
    x1 = gf.exp[m]
    Svals = {}
    for l3 in fq[1:]:
        l3log = int(gf.log[l3])
        a = gf.exp[(gf.log[x3] + l3log) % order]
        for l1 in fq:
            if l1:
                b = gf.exp[(gf.log[x1] + int(gf.log[l1])) % order]
                v = a ^ b
            else:
                v = a
            par = _abs_trace_parity(gf, v)
            Svals[(l1, l3)] = int(order - 2 * par.sum()) + 1
    byl3 = {}
    for (l1, l3), v in Svals.items():
        byl3.setdefault(l3, []).append(v)
    print(f"   S(l1,l3) multiset for each l3 != 0: {[sorted(v) for v in byl3.values()]}")
    print(f"   check sum_{{l1}} S(l1,l3) == W(l3): {all(sum(byl3[l3]) == Ws[l3] for l3 in fq[1:])}")
    # predicted rho structure
    SG = sum((int(v) - q ** 5) ** 2 for v in psi)
    SH = sum((int(psi[c * q]) - q ** 5) ** 2 for c in range(q))
    print(f"   sum_G phi^2 = {SG} = 3 q^7 (q-1)? {SG == 3*q**7*(q-1)};   sum_H phi^2 = {SH} = q^7 (q-1)? {SH == q**7*(q-1)}")
    if (q.bit_length() - 1) % 2 == 1:
        pred = Fraction(q - 3, 3 * (q - 1))
        print(f"   k odd => every eps with lambda_3 != 0 should have rho = (q-3)/(3(q-1)) = {pred}")
        rows = [r for lam, odd, r in res["rho_rows"] if lam[1] != 0]
        print(f"   observed: {sorted(set(rows))}  ({len(rows)} characters)")


def main():
    args = sys.argv[1:]
    if args:
        cases = [tuple(int(x) for x in a.split(",")) for a in args]
    else:
        cases = [(2, 4), (4, 4), (8, 4)] + [(2, d) for d in range(3, 13) if d != 4] + [(4, 5), (4, 6)]
    summary, rho_rows = [], []
    for q, d in cases:
        res = analyse(q, d)
        if d == 4:
            slice_report(res)
            d4_exponential_sums(res)
        summary.append({k: (str(v) if isinstance(v, Fraction) else v) for k, v in res.items()
                        if k not in ("rho_rows", "psi_H", "slices", "psi", "hayes")})
        for lam, odd, rho in res["rho_rows"]:
            rho_rows.append({"q": q, "d": d, "odd_j": ";".join(map(str, odd)),
                             "lambda_idx": ";".join(map(str, lam)), "rho": str(rho), "rho_float": float(rho)})
    with open(os.path.join(HERE, "legendre_var_results.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)
    with open(os.path.join(HERE, "legendre_rho_quadratic.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rho_rows[0].keys()))
        w.writeheader()
        w.writerows(rho_rows)
    print("\nwrote legendre_var_results.csv and legendre_rho_quadratic.csv")


if __name__ == "__main__":
    main()
