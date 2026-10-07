"""
legendre_q16_hayes.py -- exact Psi for larger q via the paper's Hayes-group FFT
(legendre-intervals-explicit-formula/code/hayes_fft_fq.py, hayes_exact.py), extended to q = 32, 64, 128.

For each (q, d):
  * Psi(x) for all classes (float FFT, rounded; exact CRT version for the smaller cases),
  * Var_G, Var_H (second moment of the Legendre classes about the global mean q^{d+1}), ratio,
  * d = 4: class-by-class comparison with the closed formula of legendre_d4_formula.py,
  * d = 6: Legendre profile, Var_H / q^{d+1} against q^2 and q^4, rho(eps) by support in (p1, p3, p5).

Usage: python legendre_q16_hayes.py 16,4 32,4 64,4 128,4 16,6 32,6
"""
from __future__ import annotations

import os
import shutil
import sys
import tempfile
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(os.path.dirname(HERE), "legendre-intervals-explicit-formula", "code")
sys.path.insert(0, CODE)
sys.path.insert(0, HERE)
import fqlib  # noqa: E402

fqlib.MODULI.update({32: (2, 5, [1, 0, 1, 0, 0, 1]), 64: (2, 6, [1, 1, 0, 0, 0, 0, 1]),
                     128: (2, 7, [1, 1, 0, 0, 0, 0, 0, 1]), 256: (2, 8, [1, 0, 1, 1, 1, 0, 0, 0, 1])})
import hayes_fft_fq as H  # noqa: E402
import hayes_exact as HE  # noqa: E402
from explicit_scan import legendre_mask  # noqa: E402
import legendre_d4_formula as F4  # noqa: E402


def exact_sum_sq(arr):
    """sum of squares of an int64 array without overflow (python ints, chunked)."""
    tot = 0
    for s in range(0, arr.size, 1 << 20):
        chunk = arr[s:s + (1 << 20)].astype(object)
        tot += int(np.dot(chunk, chunk))
    return tot


def power_sums(coef, q, A, M, l):
    """coef: list of arrays c1..cl (F_q ints). Newton in char 2: p_k = sum_{i<k} c_i p_{k-i} + [k odd] c_k."""
    p = [None] * (l + 1)
    for kk in range(1, l + 1):
        s = coef[kk - 1].copy() if kk % 2 == 1 else np.zeros_like(coef[0])
        for i in range(1, kk):
            s = A[s, M[coef[i - 1], p[kk - i]]]
        p[kk] = s
    return p


def run(q, d, exact=False, rho=True):
    t0 = time.time()
    l = d - 1
    N = q ** l
    wd = tempfile.mkdtemp() if N > 4_000_000 else None
    try:
        Psi, ctx = H.all_interval_psi(q, d, work_dir=wd)
    finally:
        if wd:
            shutil.rmtree(wd, ignore_errors=True)
    dev = float(np.max(np.abs(Psi - np.rint(Psi))))
    PsiI = np.rint(Psi).astype(np.int64)
    del Psi
    assert int(PsiI.sum()) == q ** (2 * d)
    if exact:
        PsiE, _, primes = HE.psi_exact(q, d)
        assert np.array_equal(PsiE, PsiI), "FFT result differs from exact CRT result"
        cert = f"exact CRT (primes {primes}) agrees"
    else:
        cert = f"float only, max rounding deviation {dev:.1e}"
    mean = q ** (d + 1)
    phi = PsiI - mean
    mask = legendre_mask(ctx)
    SG = exact_sum_sq(phi)
    SH = exact_sum_sq(phi[mask])
    nH = int(mask.sum())
    varG = Fraction(SG, N) / mean
    varH = Fraction(SH, nH) / mean
    ratio = varH / varG
    print(f"\n=== q={q}, d={d}: |G|={N}, |H|={nH}  [{cert}]  ({time.time()-t0:.1f}s)")
    print(f"  Var_G/q^(d+1) = {float(varG):.6f} (KR {d-2});  Var_H/q^(d+1) = {float(varH):.6f};  ratio = {float(ratio):.6f}")
    print(f"  Psi on H: min {int(PsiI[mask].min())}, max {int(PsiI[mask].max())}, mean q^(d+1) = {mean}")
    if d == 4:
        print(f"  q(q-1) = {q*(q-1)}: Var_H/q^5 == q(q-1): {varH == q*(q-1)};  3-3/q: {varG == Fraction(3*(q-1), q)};  q^2/3: {ratio == Fraction(q*q,3)}")
        # class-by-class comparison with the closed formula (same field encoding)
        A_, M_, NEG_, INV_ = fqlib.field_tables(q)
        psiF, K, I, _ = F4.psi_d4_formula(q, A_, M_, INV_)
        ind = H._indices_degree(l, l, q, ctx["p"], ctx["k"], ctx["A"], ctx["M"], ctx["NEG"],
                                ctx["FROB"], ctx["FINV"], ctx["axis_of"], ctx["orders"], ctx["strides"])
        byc = PsiI[ind]                      # Psi indexed by c1 + c2 q + c3 q^2
        print(f"  closed formula vs FFT, all {N} classes: {'IDENTICAL' if np.array_equal(byc, psiF) else 'MISMATCH'}"
              f"  (max |diff| = {int(np.max(np.abs(byc - psiF)))})")
    if d == 6:
        print(f"  Var_H/q^7 / q^2 = {float(varH)/q**2:.4f};   Var_H/q^7 / q^4 = {float(varH)/q**4:.6f}")
        vals = np.sort(phi[mask])
        print(f"  Legendre deviations (phi = Psi - q^7), {nH} classes: 5 most negative {vals[:5].tolist()}, "
              f"5 most positive {vals[-5:].tolist()}")
        print(f"  share of sum_H phi^2 carried by the identity class: {Fraction(int(phi[0])**2, SH) if SH else 0} "
              f"= {float(Fraction(int(phi[0])**2, SH)) if SH else 0:.4f}")
        uniq, cnt = np.unique(phi[mask], return_counts=True)
        order = np.argsort(-np.abs(uniq))
        print(f"  distinct Legendre deviations: {len(uniq)}; by |phi| (value, count, phi/q^5, share of sum_H phi^2):")
        for i in order[:12]:
            share = Fraction(int(uniq[i]) ** 2 * int(cnt[i]), SH)
            print(f"     {int(uniq[i]):>14d}  x{int(cnt[i]):<5d} phi/q^5 = {int(uniq[i]) / q**5:+.4f}   share {float(share):.4f}")
    if rho and N <= 2_000_000:
        # rho(eps_lambda) via power sums: Phi(p) = sum_{g: P(g) = p} phi(g)^2
        A_ = ctx["A"]; M_ = ctx["M"]
        ind = H._indices_degree(l, l, q, ctx["p"], ctx["k"], ctx["A"], ctx["M"], ctx["NEG"],
                                ctx["FROB"], ctx["FINV"], ctx["axis_of"], ctx["orders"], ctx["strides"])
        idx = np.arange(N, dtype=np.int64)
        coef = []
        r = idx.copy()
        for _ in range(l):
            coef.append(r % q)
            r //= q
        p = power_sums(coef, q, A_, M_, l)
        odd = [j for j in range(1, l + 1) if j % 2 == 1]
        key = np.zeros(N, np.int64)
        for t, j in enumerate(odd):
            key += p[j] * q ** t
        phic = phi[ind].astype(np.float64)
        Phi = np.bincount(key, weights=phic * phic, minlength=q ** len(odd))
        # trace table
        TR = np.zeros(q, np.int64); y = np.arange(q)
        for _ in range(ctx["k"]):
            TR = A_[TR, y]; y = M_[y, y]
        chi = 1 - 2 * TR
        m = len(odd)
        pk = np.arange(q ** m)
        pdig = [(pk // q ** t) % q for t in range(m)]
        by_support = {}
        tot = 0.0
        for lam in range(1, q ** m):
            ld = [(lam // q ** t) % q for t in range(m)]
            s = np.zeros(q ** m, np.int64)
            for t in range(m):
                s = A_[s, M_[ld[t], pdig[t]]]
            rho_val = float(np.dot(Phi, chi[s])) / SG
            tot += rho_val
            supp = tuple(j for t, j in enumerate(odd) if ld[t])
            by_support.setdefault(supp, []).append(rho_val)
        print(f"  quadratic characters (odd power sums {odd}): 1 + sum rho = {1+tot:.6f} (ratio {float(ratio):.6f})")
        for supp, rs in sorted(by_support.items()):
            print(f"     support {supp}: {len(rs):4d} chars, sum rho = {sum(rs):+.4f}, "
                  f"min {min(rs):+.4f}, max {max(rs):+.4f}")
    sys.stdout.flush()
    return dict(q=q, d=d, varG=varG, varH=varH, ratio=ratio)


if __name__ == "__main__":
    args = sys.argv[1:] or ["16,4", "32,4", "64,4", "128,4", "16,6", "32,6"]
    rows = []
    for a in args:
        q, d = map(int, a.split(","))
        rows.append(run(q, d, exact=(q ** (d - 1) <= 1_100_000), rho=(d == 6 and q ** (d - 1) <= 2_000_000)))
    print("\nsummary: q d Var_G/q^(d+1) Var_H/q^(d+1) ratio")
    for r in rows:
        print(f"  {r['q']:4d} {r['d']} {float(r['varG']):.6f} {float(r['varH']):.4f} {float(r['ratio']):.4f}")
