"""
legendre_d6_tower.py -- d = 6 (n = 12, l = 5): the Legendre profile via the tower F_q c F_{q^6} c F_{q^12}.

Part A (tower identities, brute force over F_{q^12}, q = 2 fully, q = 4 sampled):
   alpha -> (a, n) = (alpha + alpha^{q^6}, alpha^{1+q^6}),  Tr = Tr_{q^6/q}, e_j(a) = elem. sym. of the 6 conjugates:
     e1 = Tr a
     e2 = e2(a) + Tr n
     e3 = e3(a) + Tr a Tr n + Tr(a n)
     e4 = e4(a) + e2(a) Tr n + e1(a) Tr(a n) + Tr(a^2 n) + e2(n)
     e5 = e5(a) + e3(a) Tr n + e2(a) Tr(a n) + e1(a) Tr(a^2 n) + Tr(a^3 n) + e1(a) e2(n) + Tr(a n) Tr n + Tr(a n^2)
     p3 = Tr(a^3 + a n),   p5 = Tr(a^5 + a^3 n + a n^2)
   and e2(beta + beta^q) = Tr_{q^12/q}(beta^{q+1}) + (Tr beta)^2.

Part B (structure theorem, exact data from the paper's Hayes FFT, q = 2..32):
   With U = {alpha : Tr alpha = Tr alpha^3 = Tr alpha^5 = 0}, A6 = |U|, X = sum_U chi(e2), Y = sum_U chi(e4):
     Psi(0,c2,0,c4,0) = q^7 + (A6 - q^9 - X)/q^2 + (X/q)[c2 = 0] + (Y/q)[c2 != 0] chi(c4 / c2^2).
   The script extracts A6, X, Y, checks the formula on every Legendre class, and prints the closed forms.
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "legendre-intervals-explicit-formula", "code"))
import legendre_var_check as LV  # noqa: E402  (GF class)


# ----------------------------------------------------------------------------- Part A
def part_a(q, sample=None):
    k = q.bit_length() - 1
    N = 12 * k
    gf = LV.GF(N)
    order = gf.order
    mul, pw = gf.mul, gf.pow

    def tr6(x):   # Tr_{q^6/q}
        s = 0
        for i in range(6):
            s ^= pw(x, q ** i)
        return s

    def tr12(x):
        s = 0
        for i in range(12):
            s ^= pw(x, q ** i)
        return s

    def esym(roots):
        poly = [1]
        for r in roots:
            new = [0] * (len(poly) + 1)
            for i, c in enumerate(poly):
                new[i] ^= c
                new[i + 1] ^= mul(c, r)
            poly = new
        return poly[1:]  # e1, e2, ...

    rng = np.random.default_rng(1)
    if sample is None:
        alphas = range(1, gf.size)
    else:
        alphas = [int(x) for x in rng.integers(1, gf.size, size=sample)]
    bad = 0
    t0 = time.time()
    for alpha in alphas:
        conj12 = [pw(alpha, q ** i) for i in range(12)]
        e = esym(conj12)
        a = alpha ^ pw(alpha, q ** 6)
        n = mul(alpha, pw(alpha, q ** 6))
        conj6 = [pw(a, q ** i) for i in range(6)]
        ea = esym(conj6)
        Tn, Tan, Ta2n, Ta3n, Tan2 = tr6(n), tr6(mul(a, n)), tr6(mul(mul(a, a), n)), tr6(mul(pw(a, 3), n)), tr6(mul(a, mul(n, n)))
        e2n = esym([pw(n, q ** i) for i in range(6)])[1]
        f1 = ea[0]
        f2 = ea[1] ^ Tn
        f3 = ea[2] ^ mul(ea[0], Tn) ^ Tan
        f4 = ea[3] ^ mul(ea[1], Tn) ^ mul(ea[0], Tan) ^ Ta2n ^ e2n
        f5 = ea[4] ^ mul(ea[2], Tn) ^ mul(ea[1], Tan) ^ mul(ea[0], Ta2n) ^ Ta3n ^ mul(ea[0], e2n) ^ mul(Tan, Tn) ^ Tan2
        p3 = tr12(pw(alpha, 3)); p5 = tr12(pw(alpha, 5))
        g3 = tr6(pw(a, 3) ^ mul(a, n)); g5 = tr6(pw(a, 5) ^ mul(pw(a, 3), n) ^ mul(a, mul(n, n)))
        # e2(beta + beta^q) identity
        beta = alpha
        al2 = beta ^ pw(beta, q)
        e2al = esym([pw(al2, q ** i) for i in range(12)])[1]
        t12 = tr12(beta)
        h2 = tr12(mul(beta, pw(beta, q))) ^ mul(t12, t12)
        if (f1, f2, f3, f4, f5) != tuple(e[:5]) or p3 != g3 or p5 != g5 or e2al != h2:
            bad += 1
    print(f"Part A, q={q}: checked {len(alphas)} elements of F_{{q^12}}, mismatches = {bad}  ({time.time()-t0:.1f}s)")
    return bad == 0


# ----------------------------------------------------------------------------- Part B
def part_b(q):
    import fqlib
    fqlib.MODULI.update({32: (2, 5, [1, 0, 1, 0, 0, 1])})
    import hayes_fft_fq as H
    import shutil, tempfile
    d, l = 6, 5
    N = q ** l
    wd = tempfile.mkdtemp() if N > 4_000_000 else None
    try:
        Psi, ctx = H.all_interval_psi(q, d, work_dir=wd)
    finally:
        if wd:
            shutil.rmtree(wd, ignore_errors=True)
    PsiI = np.rint(Psi).astype(np.int64)
    del Psi
    ind = H._indices_degree(l, l, q, ctx["p"], ctx["k"], ctx["A"], ctx["M"], ctx["NEG"],
                            ctx["FROB"], ctx["FINV"], ctx["axis_of"], ctx["orders"], ctx["strides"])
    A_, M_, INV_ = ctx["A"], ctx["M"], ctx["INV"]
    k = ctx["k"]
    TR = np.zeros(q, np.int64); y = np.arange(q)
    for _ in range(k):
        TR = A_[TR, y]; y = M_[y, y]
    chi = 1 - 2 * TR
    SQ = M_[np.arange(q), np.arange(q)]
    SQRT = np.empty(q, np.int64); SQRT[SQ] = np.arange(q)
    P = np.zeros((q, q), dtype=object)        # P[c2][c4] = Psi(0,c2,0,c4,0)
    for c2 in range(q):
        for c4 in range(q):
            idx = c2 * q + c4 * q ** 3
            P[c2][c4] = int(PsiI[ind[idx]])
    mean = q ** 7
    A6 = sum(int(P[c2][c4]) for c2 in range(q) for c4 in range(q))
    X = sum(int(P[c2][c4]) * int(chi[c2]) for c2 in range(q) for c4 in range(q))
    Y = sum(int(P[c2][c4]) * int(chi[c4]) for c2 in range(q) for c4 in range(q))
    # structure theorem check
    base = Fraction(A6 - q ** 9 - X, q * q)
    ok = True
    for c2 in range(q):
        for c4 in range(q):
            if c2 == 0:
                pred = mean + base + Fraction(X, q)
            else:
                u = M_[c4, INV_[SQ[c2]]]
                pred = mean + base + Fraction(Y, q) * int(chi[u])
            if pred != P[c2][c4]:
                ok = False
    # Fourier check of Psi^(mu,nu) = Y chi(mu / sqrt nu)
    fok = True
    for mu in range(q):
        for nu in range(1, q):
            hat = sum(int(P[c2][c4]) * int(chi[A_[M_[mu, c2], M_[nu, c4]]]) for c2 in range(q) for c4 in range(q))
            if hat != Y * int(chi[M_[mu, INV_[SQRT[nu]]]]):
                fok = False
    eps = 1 if k % 2 else -1
    varH = Fraction(sum((int(P[c2][c4]) - mean) ** 2 for c2 in range(q) for c4 in range(q)), q * q) / mean
    m = Fraction(Y, q ** 3)
    if k % 2:
        varH_pred = Fraction((q - 1) ** 2 * (q + 2) ** 2, q * q) + Fraction(q - 1, q ** 4) * (4 * q * q + m * m)
    else:
        varH_pred = (q - 1) ** 2 + Fraction(q - 1, q ** 4) * m * m
    print(f"\nPart B, q={q} (k={k}):  |U| = A6 = {A6} = q^9 + ({A6 - q**9});  predicted (-1)^(k+1)(q-1)q^5 = {eps*(q-1)*q**5}: {A6 - q**9 == eps*(q-1)*q**5}")
    print(f"   X = sum_U chi(e2) = {X};  predicted (-1)^(k+1) q^6 + q^5 = {eps*q**6 + q**5}: {X == eps*q**6 + q**5}")
    print(f"   Y = sum_U chi(e4) = {Y} = q^3 * {m}   (Y/q^3 integer: {m.denominator == 1})")
    print(f"   structure theorem Psi(0,c2,0,c4,0) = q^7 + (A6-q^9-X)/q^2 + (X/q)[c2=0] + (Y/q)[c2!=0]chi(c4/c2^2): "
          f"{'OK on all classes' if ok else 'FAILED'};  Fourier form Psi^(mu,nu)=Y chi(mu/sqrt(nu)) for nu!=0: {'OK' if fok else 'FAILED'}")
    print(f"   phi_0 (c2=0) = {int(P[0][0]) - mean};  phi(Tr u=0) = {base + Fraction(Y, q)};  phi(Tr u=1) = {base - Fraction(Y, q)}")
    print(f"   Var_H/q^7 = {float(varH):.6f};  closed form in (q, k, m) = {float(varH_pred):.6f}: {varH == varH_pred}")
    return dict(q=q, A6=A6, X=X, Y=Y, m=m)


def part_c(q):
    """Tower reduction of A6, X, Y to sums over a in F_{q^6} (Tr a = 0, a != 0), (lambda, mu) in F_q^2, eps in {0,1}:
         beta' = lambda a + sqrt(mu) sqrt(a) + mu a^3 + eps a^{-2},   beta = e2(a) + a^2 + beta'
         A6 = q^6 + q^4 sum_a sum_{beta'=0} (-1)^eps chi(lambda Tr a^3 + mu Tr a^5)
         X  =       q^4 sum_a chi(e2(a)) sum_{1+beta'=0} (-1)^eps chi(lambda Tr a^3 + mu Tr a^5)
         Y  = G(e2) [ 1 + q^-2 sum_a sum_{lambda,mu,eps} (-1)^eps chi(e4(a) + lambda Tr a^3 + mu Tr a^5 + e2(beta)) ]
       with G(e2) = sum_{n in F_{q^6}} chi(e2(n)).  Checked against the direct values (brute force over F_{q^12})."""
    k = q.bit_length() - 1
    gf6 = LV.GF(6 * k)
    mul, pw = gf6.mul, gf6.pow
    fq = gf6.subfield(q)

    def tr6(x):
        s = 0
        for i in range(6):
            s ^= pw(x, q ** i)
        return s

    def chi(x):   # x in F_q (as element of gf6): (-1)^{Tr_{q/2} x}
        s, y = 0, x
        for _ in range(k):
            s ^= y
            y = mul(y, y)
        return 1 - 2 * s

    def esym6(x):
        poly = [1]
        for i in range(6):
            r = pw(x, q ** i)
            new = [0] * (len(poly) + 1)
            for j, c in enumerate(poly):
                new[j] ^= c
                new[j + 1] ^= mul(c, r)
            poly = new
        return poly[1:]

    sqrt_tab = {mul(s, s): s for s in fq}
    G = sum(chi(esym6(n)[1]) for n in range(gf6.size))
    A6 = q ** 6
    X = 0
    Ysum = 0
    for a in range(1, gf6.size):
        if tr6(a) != 0:
            continue
        ea = esym6(a)
        Ta3, Ta5 = tr6(pw(a, 3)), tr6(pw(a, 5))
        sa = pw(a, 1 << (6 * k - 1))      # sqrt(a) in F_{q^6}: a^(2^(6k-1))
        ainv2 = pw(a, gf6.order - 2)
        base = ea[1] ^ mul(a, a)
        for lam in fq:
            for mu in fq:
                smu = sqrt_tab[mu]
                bp0 = mul(lam, a) ^ mul(smu, sa) ^ mul(mu, pw(a, 3))
                w = chi(mul(lam, Ta3) ^ mul(mu, Ta5))
                for eps in (0, 1):
                    bp = bp0 ^ (ainv2 if eps else 0)
                    sgn = -1 if eps else 1
                    if bp == 0:
                        A6 += q ** 4 * sgn * w
                    if bp == 1:
                        X += q ** 4 * sgn * w * chi(ea[1])
                    beta = base ^ bp
                    Ysum += sgn * chi(ea[3] ^ mul(lam, Ta3) ^ mul(mu, Ta5) ^ esym6(beta)[1])
    Y = G * (1 + Fraction(Ysum, q * q))
    print(f"Part C, q={q}: G(e2) on F_(q^6) = {G} (q^3 = {q**3});  tower reductions give A6 = {A6}, X = {X}, Y = {Y}")
    return A6, X, Y


if __name__ == "__main__":
    part_a(2)
    part_a(4, sample=1500)
    part_c(2)
    part_c(4)
    rows = []
    for q in (2, 4, 8, 16, 32):
        rows.append(part_b(q))
    print("\nsummary: q, A6-q^9, X, Y, m=Y/q^3")
    for r in rows:
        print(f"   {r['q']:3d}  {r['A6']-r['q']**9:>12d}  {r['X']:>14d}  {r['Y']:>16d}  {r['m']}")
