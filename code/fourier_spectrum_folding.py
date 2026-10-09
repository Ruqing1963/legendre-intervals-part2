"""
fourier_spectrum_folding.py -- spectral folding of the Legendre variance (Section 8 of the paper).

For H = G^2 < G (Hayes group, |G| = q^l, index |G[2]| = q^{d/2}) and phi = Psi - q^{d+1} with Fourier expansion
phi = sum_chi c_chi chi, restriction Res: G^ -> H^ is onto with fibres of size |G[2]|, and

   Var_G = sum_chi |c_chi|^2 = sum_{psi in H^} E_psi,     E_psi = sum_{chi in fibre(psi)} |c_chi|^2,
   Var_H = sum_{psi in H^} |f_psi|^2,                     f_psi = sum_{chi in fibre(psi)} c_chi,

so Var_H / Var_G <= |G[2]| with equality iff phi is supported on H.  Here f_psi is the psi-Fourier coefficient of
phi restricted to the coset H itself, and E_psi is the average over all cosets gH of the squared psi-coefficient of
phi|_{gH}.  The coherence of a fibre is kappa_psi = |f_psi|^2 / (|G[2]| E_psi) in [0, 1] (Cauchy-Schwarz).

d = 4: all classes from the closed formula (legendre_d4_formula.py), H = G^2 = {(0,c,0)} = F_q, characters (-1)^{Tr(mu c)}.
d = 6: all classes from the Hayes FFT (paper I), H = G^2 = {(c2, c4)} = W_2(F_q) with the Witt law
       (c2,c4)(c2',c4') = (c2+c2', c4+c4'+c2c2'); characters psi_{(a,b)}(x) = i^{Tr_{W_2(F_q)/Z4}((a,b).x)}, of order
       1 (a=b=0), 2 (a=0) or 4 (a != 0).
Usage: python fourier_spectrum_folding.py            (d=4: q = 2..32; d=6: q = 2,4,8,16)
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
sys.path.insert(0, HERE)
import fqlib  # noqa: E402

fqlib.MODULI.update({32: (2, 5, [1, 0, 1, 0, 0, 1]), 64: (2, 6, [1, 1, 0, 0, 0, 0, 1]),
                     128: (2, 7, [1, 1, 0, 0, 0, 0, 0, 1])})
import hayes_fft_fq as H  # noqa: E402
import legendre_d4_formula as F4  # noqa: E402


def field(q):
    A, M, NEG, INV = fqlib.field_tables(q)
    k = q.bit_length() - 1
    x = np.arange(q)
    SQ = M[x, x]
    TR = np.zeros(q, np.int64); y = x.copy()
    for _ in range(k):
        TR = A[TR, y]; y = M[y, y]
    return A, M, INV, SQ, TR, k


def witt_chars(q, A, M, SQ, k):
    """psi[(a,b) index, (x0,x1) index] = i^{Tr_{W2(F_q)/Z4}((a,b)(x0,x1))}."""
    # Witt trace table T[y0, y1] in Z/4
    T = np.zeros((q, q), np.int64)
    for y0 in range(q):
        for y1 in range(q):
            s0, s1 = 0, 0
            c0, c1 = y0, y1
            for _ in range(k):
                s0, s1 = int(A[s0, c0]), int(A[A[s1, c1], M[s0, c0]])
                c0, c1 = int(SQ[c0]), int(SQ[c1])
            assert s0 in (0, 1) and s1 in (0, 1)
            T[y0, y1] = (s0 + 2 * s1) % 4
    psi = np.zeros((q * q, q * q), np.complex128)
    i4 = np.array([1, 1j, -1, -1j])
    for a in range(q):
        for b in range(q):
            for x0 in range(q):
                for x1 in range(q):
                    y0 = int(M[a, x0]); y1 = int(A[M[SQ[a], x1], M[b, SQ[x0]]])
                    psi[a * q + b, x0 * q + x1] = i4[T[y0, y1]]
    return psi


def report(name, q, d, fib_sums, fib_energy, orders, VarH_direct, VarG_direct, extra=""):
    mean = q ** (d + 1)
    idx2 = q ** (d // 2)
    VarH = float(np.sum(np.abs(fib_sums) ** 2)); VarG = float(np.sum(fib_energy))
    print(f"\n=== d={d}, q={q}: |G[2]| = q^{d//2} = {idx2};  Var_G/q^{d+1} = {VarG/mean:.6f} (direct {VarG_direct/mean:.6f}),"
          f"  Var_H/q^{d+1} = {VarH/mean:.6f} (direct {VarH_direct/mean:.6f}),  ratio = {VarH/VarG:.4f} = q^2 x {VarH/VarG/q/q:.4f}")
    for o in sorted(set(orders)):
        sel = orders == o
        sH = float(np.sum(np.abs(fib_sums[sel]) ** 2)); sG = float(np.sum(fib_energy[sel]))
        kap = np.abs(fib_sums[sel]) ** 2 / (idx2 * fib_energy[sel])
        print(f"   fibres over characters of H of order {o}: {int(sel.sum()):5d} fibres;  share of Var_H {sH/VarH:.4f};"
              f"  share of Var_G {sG/VarG:.4f};  coherence kappa: min {kap.min():.4f}, mean {kap.mean():.4f}, max {kap.max():.4f}")
    if extra:
        print("   " + extra)


def run_d4(q):
    A, M, INV, SQ, TR, k = field(q)
    psi, _, _, _ = F4.psi_d4_formula(q, A, M, INV)
    phi = psi.astype(np.float64) - q ** 5
    chi = np.array([[1 - 2 * TR[M[mu, c]] for c in range(q)] for mu in range(q)], np.float64)   # chi[mu, c]
    # cosets of H = {(0,c,0)}: representative (c1,0,c3), members (c1, c, c3 + c1 c)
    Phi = np.zeros((q * q, q))
    for c1 in range(q):
        for c3 in range(q):
            for c in range(q):
                Phi[c1 * q + c3, c] = phi[c1 + c * q + int(A[c3, M[c1, c]]) * q * q]
    Fc = Phi @ chi.T / q                        # coset Fourier coefficients, shape (cosets, mu)
    fib_sums = Fc[0].astype(np.complex128)      # identity coset = H
    fib_energy = np.mean(np.abs(Fc) ** 2, axis=0)
    orders = np.array([1] + [2] * (q - 1))
    VarH_direct = float(np.mean(phi[np.arange(q) * q] ** 2))
    VarG_direct = float(np.mean(phi ** 2))
    pred = [(-1) ** k * q ** 3 - q ** 3 + q * q - (-1) ** k * q * q] + [-q ** 3 + q * q - (-1) ** k * q * q] * (q - 1)
    ok = np.allclose(fib_sums.real, pred)
    report("d4", q, 4, fib_sums, fib_energy, orders, VarH_direct, VarG_direct,
           extra=f"fibre sums: f_0 = {fib_sums[0].real:.0f}, f_mu = {fib_sums[1].real:.0f} (mu != 0);  Theorem 8.2 prediction: {'OK' if ok else 'MISMATCH'};"
                 f"  identity-class share of Var_G: {phi[0]**2/q**3/VarG_direct:.4f}")


def run_d6(q):
    A, M, INV, SQ, TR, k = field(q)
    l = 5
    N = q ** l
    wd = tempfile.mkdtemp() if N > 4_000_000 else None
    try:
        Psi, ctx = H.all_interval_psi(q, 6, work_dir=wd)
    finally:
        if wd:
            shutil.rmtree(wd, ignore_errors=True)
    PsiI = np.rint(Psi).astype(np.int64)
    ind = H._indices_degree(l, l, q, ctx["p"], ctx["k"], ctx["A"], ctx["M"], ctx["NEG"],
                            ctx["FROB"], ctx["FINV"], ctx["axis_of"], ctx["orders"], ctx["strides"])
    byc = PsiI[ind]                             # indexed by c1 + c2 q + c3 q^2 + c4 q^3 + c5 q^4
    A_, M_ = ctx["A"], ctx["M"]
    SQ_ = M_[np.arange(q), np.arange(q)]
    phi = byc.astype(np.float64) - q ** 7
    psi = witt_chars(q, A_, M_, SQ_, k)
    # cosets: rep (c1,0,c3,0,c5); members h=(h2,h4): (c1, h2, c3+c1 h2, h4, c5 + c3 h2 + c1 h4)
    Phi = np.zeros((q ** 3, q * q))
    for c1 in range(q):
        for c3 in range(q):
            for c5 in range(q):
                row = (c1 * q + c3) * q + c5
                for h2 in range(q):
                    cc3 = int(A_[c3, M_[c1, h2]])
                    for h4 in range(q):
                        cc5 = int(A_[A_[c5, M_[c3, h2]], M_[c1, h4]])
                        Phi[row, h2 * q + h4] = phi[c1 + h2 * q + cc3 * q ** 2 + h4 * q ** 3 + cc5 * q ** 4]
    Fc = Phi @ np.conj(psi).T / (q * q)
    fib_sums = Fc[0]
    fib_energy = np.mean(np.abs(Fc) ** 2, axis=0)
    orders = np.array([1 if (a == 0 and b == 0) else (2 if a == 0 else 4) for a in range(q) for b in range(q)])
    VarH_direct = float(np.mean(Phi[0] ** 2)); VarG_direct = float(np.mean(phi ** 2))
    # Y = sum_{G^2} Psi(g) (-1)^{Tr c4}, A6 = sum_{G^2} Psi, X = sum Psi (-1)^{Tr c2}
    chi1 = np.array([1 - 2 * int(ctx_tr(ctx, c)) for c in range(q)])
    PsiH = Phi[0] + q ** 7
    A6 = int(round(PsiH.sum()))
    Yv = int(round(sum(PsiH[h2 * q + h4] * chi1[h4] for h2 in range(q) for h4 in range(q))))
    Xv = int(round(sum(PsiH[h2 * q + h4] * chi1[h2] for h2 in range(q) for h4 in range(q))))
    m = Fraction(Yv, q ** 3)
    s4 = float(np.sum(np.abs(fib_sums[orders == 4]) ** 2))
    extra = (f"A6 = {A6}, X = {Xv}, Y = {Yv} = q^3 * {m};  fibre sums: f_1 = (A6 - q^9)/q^2 = {fib_sums[0].real:.1f},"
             f" f_mu (order 2) = X/q^2 = {fib_sums[1].real:.1f};  order-4 part of Var_H = {s4:.1f} vs q^3 (q-1) m_k^2 = {float(q**3*(q-1)*m*m):.1f}")
    report("d6", q, 6, fib_sums, fib_energy, orders, VarH_direct, VarG_direct, extra=extra)


def ctx_tr(ctx, c):
    """absolute trace of c in the FFT context's field encoding"""
    A_, M_ = ctx["A"], ctx["M"]
    s, y = 0, c
    for _ in range(ctx["k"]):
        s = int(A_[s, y]); y = int(M_[y, y])
    return s


if __name__ == "__main__":
    for q in (2, 4, 8, 16, 32):
        run_d4(q)
    for q in (2, 4, 8, 16):
        run_d6(q)
