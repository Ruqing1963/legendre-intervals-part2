"""
general_d_scan.py -- the Legendre variance for general even d (Section 9): exact Psi on every class of the Hayes
group G_{d-1} by the Hayes FFT of paper I, and

  * Var_G, Var_{G^2} (second moments about q^{d+1}), their ratio and Var_{G^2}/q^{d+1} against q^2;
  * the identity-class deviation phi(1) = Psi(1) - q^{d+1}, i.e. #{alpha in F_{q^{2d}} : e_1 = ... = e_{d-1} = 0} - q^{d+1};
  * the translation orbit K_d of the identity class, K_d = {(1 + c u)^{2d} mod u^d : c in F_q}: by Lucas' theorem
    |K_d| = 1 if d is a power of 2 and |K_d| = q otherwise; Psi is constant on translation orbits, so phi = phi(1) on K_d;
  * the share of sum_{G^2} phi^2 carried by K_d and by the largest deviations.
Usage: python general_d_scan.py 2,8 4,8 8,8 2,10 4,10
"""
import os
import shutil
import sys
import tempfile
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "legendre-intervals-explicit-formula", "code"))
import fqlib  # noqa: E402

fqlib.MODULI.update({32: (2, 5, [1, 0, 1, 0, 0, 1]), 64: (2, 6, [1, 1, 0, 0, 0, 0, 1]),
                     128: (2, 7, [1, 1, 0, 0, 0, 0, 0, 1])})
import hayes_fft_fq as H  # noqa: E402
from explicit_scan import legendre_mask  # noqa: E402


def exact_sum_sq(arr):
    tot = 0
    for s in range(0, arr.size, 1 << 20):
        chunk = arr[s:s + (1 << 20)].astype(object)
        tot += int(np.dot(chunk, chunk))
    return tot


def binom_mod2(n, j):
    return 1 if (j & n) == j else 0            # Lucas


def run(q, d):
    t0 = time.time()
    l = d - 1
    N = q ** l
    wd = tempfile.mkdtemp() if N > 4_000_000 else None
    try:
        Psi, ctx = H.all_interval_psi(q, d, work_dir=wd)
    finally:
        if wd:
            shutil.rmtree(wd, ignore_errors=True)
    PsiI = np.rint(Psi).astype(np.int64)
    del Psi
    mean = q ** (d + 1)
    phi = PsiI - mean
    mask = legendre_mask(ctx)
    nH = int(mask.sum())
    SG, SH = exact_sum_sq(phi), exact_sum_sq(phi[mask])
    varG = Fraction(SG, N) / mean
    varH = Fraction(SH, nH) / mean
    A_, M_ = ctx["A"], ctx["M"]
    # identity class and the orbit K_d = {(1 + c u)^{2d} mod u^{l+1}}
    id_idx = H.class_index([1] + [0] * l, ctx)
    odd_j = [j for j in range(1, l + 1) if binom_mod2(2 * d, j)]
    Kidx = set()
    for c in range(q):
        series = [1] + [0] * l
        cp = 1
        for j in range(1, l + 1):
            cp = int(M_[cp, c])
            series[j] = cp if binom_mod2(2 * d, j) else 0
        Kidx.add(H.class_index(series, ctx))
    Kidx = sorted(Kidx)
    phiK = [int(phi[i]) for i in Kidx]
    phi_id = int(phi[id_idx])
    print(f"\n=== q={q}, d={d}: |G|={N}, |G^2|={nH}, [G:G^2]=q^{d//2}={q**(d//2)}  ({time.time()-t0:.1f}s)")
    print(f"  Var_G/q^(d+1) = {float(varG):.6f} (KR {d-2}),  Var_G2/q^(d+1) = {float(varH):.6f},  ratio = {float(ratio := varH/varG):.4f}")
    print(f"  Var_G2/q^(d+1) / q^2 = {float(varH)/q**2:.4f},   / q^(d/2) = {float(varH)/q**(d//2):.6f},   / q^(d/2-1) = {float(varH)/q**(d//2-1):.6f}")
    print(f"  phi(1) = Psi(1) - q^(d+1) = {phi_id} = {phi_id/q**d:+.4f} q^d = {phi_id/q**(d-1):+.4f} q^(d-1);"
          f"  Psi(1) = {phi_id + mean}")
    print(f"  Lucas: C(2d,j) odd for j in {odd_j} (j <= d-1)  ->  |K_d| = {len(Kidx)};  phi on K_d all equal to phi(1): {all(v == phi_id for v in phiK)}")
    shareK = Fraction(len(Kidx) * phi_id ** 2, SH) if SH else 0
    print(f"  share of sum_{{G^2}} phi^2 carried by K_d: {float(shareK):.4f};   |K_d| phi(1)^2 / |G^2| / q^(d+1) = {float(Fraction(len(Kidx)*phi_id**2, nH))/mean:.4f}")
    vals, cnt = np.unique(phi[mask], return_counts=True)
    order = np.argsort(-np.abs(vals))
    print(f"  distinct Legendre deviations: {len(vals)}; largest |phi| (value, count, phi/q^(d-1), share of sum_H phi^2):")
    for i in order[:8]:
        print(f"     {int(vals[i]):>16d}  x{int(cnt[i]):<6d}  {int(vals[i])/q**(d-1):+9.4f}   {float(Fraction(int(vals[i])**2*int(cnt[i]), SH)):.4f}")
    sys.stdout.flush()
    return dict(q=q, d=d, varG=varG, varH=varH, ratio=ratio, phi_id=phi_id, K=len(Kidx))


if __name__ == "__main__":
    args = sys.argv[1:] or ["2,4", "4,4", "2,6", "4,6", "2,8", "4,8", "8,8", "2,10", "4,10"]
    rows = []
    for a in args:
        q, d = map(int, a.split(","))
        rows.append(run(q, d))
    print("\nsummary:  q  d  Var_G/q^(d+1)  Var_G2/q^(d+1)  ratio  ratio/q^2  phi(1)/q^(d-1)  |K_d|")
    for r in rows:
        print(f"  {r['q']:3d} {r['d']:3d}  {float(r['varG']):10.4f}  {float(r['varH']):12.4f}  {float(r['ratio']):9.4f}  {float(r['ratio'])/r['q']**2:8.4f}"
              f"  {r['phi_id']/r['q']**(r['d']-1):+10.4f}  {r['K']}")
