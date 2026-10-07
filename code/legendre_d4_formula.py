"""
legendre_d4_formula.py -- closed formula for Psi on the Hayes group G_3 (d = 4, n = 8) over any
F_q, q = 2^k, obtained from the quadratic tower F_q c F_{q^2} c F_{q^4} c F_{q^8}:

  alpha in F_{q^8}  <->  (a, n) = (alpha + alpha^{q^4}, alpha^{1+q^4}) in F_{q^4}^2,
  charpoly(alpha) = N_{F_{q^4}/F_q}(t^2 + a t + n),  so
      e1 = Tr(a),  e2 = e2(a) + Tr(n),  e3 = e3(a) + Tr(a)Tr(n) + Tr(a n)      (Tr = Tr_{q^4/q}).

Counting the pairs (a, n) with t^2 + a t + n irreducible over F_{q^4} (2 roots alpha each) or
a = 0 (alpha in F_{q^4}, one root) gives, for every class (c1, c2, c3):

  Psi(c1, c2, c3) = q^5                                             if c1 != 0   (Theorem T1)
  Psi(0,  c2, c3) = q^5 + q^3 [c3 = 0] + 2 q^3 I(c2, c3) - 2 q^2 K(c2, c3)
      I(c2, c3) = [c2 c3 != 0 and Tr_{q/2}(c2^3 / c3^2) = 1]
      K(c2, c3) = sum_{w : Tr_{q/2} w = 1} sum_{v in F_q} chi( c2 v^2 / w + sqrt(c2) v / w + c3 v^3 / w^2 ),
      chi(x) = (-1)^{Tr_{q/2}(x)}.

Consequences (all exact):  Psi(0,0,0) = q^5 - q^4 + q^3,  Psi(0,c,0) = q^5 + q^3 (k even),
q^5 - q^3 (k odd),  Var_{G^2}(Psi)/q^5 = q(q-1).  The script evaluates the formula for k = 1..KMAX,
checks these identities, Var_G/q^5 = 3 - 3/q, the second-moment identity
sum_{c3 != 0, c2} (K - q I)^2 = q^3 (q-1)/2, and prints the rho(eps) values.

Field arithmetic: the paper's small_q_gf.make_field (polynomial basis, integer encoding), so the
class labelling agrees with hayes_fft_fq / fqlib and allows a class-by-class comparison.
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "legendre-intervals-explicit-formula", "code"))
from small_q_gf import make_field  # noqa: E402

# irreducible moduli over F_2 (low degree first), k = 1..10
MODULI2 = {1: [1, 1], 2: [1, 1, 1], 3: [1, 1, 0, 1], 4: [1, 1, 0, 0, 1], 5: [1, 0, 1, 0, 0, 1],
           6: [1, 1, 0, 0, 0, 0, 1], 7: [1, 1, 0, 0, 0, 0, 0, 1], 8: [1, 0, 1, 1, 1, 0, 0, 0, 1],
           9: [1, 0, 0, 0, 1, 0, 0, 0, 0, 1], 10: [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]}


def _poly(*exps):
    out = [0] * (max(exps) + 1)
    for e in exps:
        out[e] = 1
    return out


# primitive polynomials x^k + ... for the O(q) check (k = 11..22)
MODULI2.update({11: _poly(11, 2, 0), 12: _poly(12, 6, 4, 1, 0), 13: _poly(13, 4, 3, 1, 0), 14: _poly(14, 10, 6, 1, 0),
                15: _poly(15, 1, 0), 16: _poly(16, 12, 3, 1, 0), 17: _poly(17, 3, 0), 18: _poly(18, 7, 0),
                19: _poly(19, 5, 2, 1, 0), 20: _poly(20, 3, 0), 21: _poly(21, 2, 0), 22: _poly(22, 1, 0)})


def field(k):
    q, add, mul, neg, inv = make_field(2, k, MODULI2[k])
    A = np.array(add, np.int64)
    M = np.array(mul, np.int64)
    INV = np.array(inv, np.int64)
    return q, A, M, INV


def tables(q, A, M, INV):
    k = q.bit_length() - 1
    x = np.arange(q)
    SQ = M[x, x]
    SQRT = np.empty(q, np.int64)
    SQRT[SQ] = x
    CUBE = M[SQ, x]
    TR = np.zeros(q, np.int64)
    y = x.copy()
    for _ in range(k):
        TR = A[TR, y]
        y = M[y, y]
    assert set(np.unique(TR)) <= {0, 1}
    W = np.where(TR == 1)[0]
    return SQ, SQRT, CUBE, TR, W


def psi_d4_formula(q, A, M, INV, verbose=False):
    """Psi as an array indexed by idx = c1 + c2 q + c3 q^2 (F_q elements in make_field encoding),
    plus the arrays K, I over (c2, c3)."""
    SQ, SQRT, CUBE, TR, W = tables(q, A, M, INV)
    chi = 1 - 2 * TR                                         # (-1)^{Tr x}
    v = np.arange(q)
    wi = INV[W]                                              # 1/w
    vw1 = M[SQ[v][:, None], wi[None, :]]                     # v^2 / w     shape (q, |W|)
    vw2 = M[v[:, None], wi[None, :]]                         # v / w
    vw3 = M[CUBE[v][:, None], SQ[wi][None, :]]               # v^3 / w^2
    K = np.zeros((q, q), np.int64)                           # K[c2, c3]
    c2 = np.arange(q)
    t12 = A[M[c2[:, None, None], vw1[None]], M[SQRT[c2][:, None, None], vw2[None]]]   # (q, q, |W|)
    for c3 in range(q):
        term = A[t12, M[c3, vw3][None]]
        K[:, c3] = chi[term].sum(axis=(1, 2))
    I = np.zeros((q, q), np.int64)
    for c2v in range(1, q):
        for c3v in range(1, q):
            I[c2v, c3v] = TR[M[CUBE[c2v], INV[SQ[c3v]]]]
    psi = np.full(q ** 3, q ** 5, np.int64)
    for c2v in range(q):
        for c3v in range(q):
            val = q ** 5 + (q ** 3 if c3v == 0 else 0) + 2 * q ** 3 * int(I[c2v, c3v]) - 2 * q ** 2 * int(K[c2v, c3v])
            psi[0 + c2v * q + c3v * q * q] = val
    return psi, K, I, (SQ, SQRT, CUBE, TR, W)


def analyse_formula(k, verbose=True):
    t0 = time.time()
    q, A, M, INV = field(k)
    psi, K, I, (SQ, SQRT, CUBE, TR, W) = psi_d4_formula(q, A, M, INV)
    mean = q ** 5
    assert int(psi.sum()) == q ** 8, "sum over G must be q^8"
    phi = psi.astype(object) - mean
    SG = sum(int(x) ** 2 for x in phi)
    Hidx = [c * q for c in range(q)]                       # classes (0, c, 0)
    SH = sum(int(phi[i]) ** 2 for i in Hidx)
    varG = Fraction(SG, q ** 3) / mean
    varH = Fraction(SH, q) / mean
    ratio = varH / varG
    Aval = sum(int(psi[i]) for i in Hidx)
    # second-moment identity on G_0 \ H
    S2 = 0
    for c3 in range(1, q):
        for c2 in range(q):
            S2 += (int(K[c2, c3]) - q * int(I[c2, c3])) ** 2
    # rho(eps_{lambda1, lambda3}) = r(lambda3) = sum_{x,y} phi(0,x,y)^2 chi(lambda3 y) / SG
    chi = 1 - 2 * TR
    Phi = [sum(int(phi[x * q + y * q * q]) ** 2 for x in range(q)) for y in range(q)]   # Phi(y)
    r = {}
    for lam in range(1, q):
        r[lam] = Fraction(sum(Phi[y] * int(chi[M[lam, y]]) for y in range(q)), SG)
    rho_sum = (q - 1) * 1 + q * sum(r.values())
    kodd = (k % 2 == 1)
    out = dict(k=k, q=q, varG=varG, varH=varH, ratio=ratio, A=Aval,
               psi_id=int(psi[0]), psi_c=sorted({int(psi[c * q]) for c in range(1, q)}),
               S2=S2, Phi=Phi, r=r, rho_sum=rho_sum, seconds=time.time() - t0)
    if verbose:
        print(f"\n=== closed formula, k={k}, q={q}  [{out['seconds']:.1f}s]")
        print(f"  Psi(id) = {psi[0]} = q^5-q^4+q^3: {int(psi[0]) == q**5-q**4+q**3}")
        print(f"  Psi(0,c,0), c!=0: {out['psi_c']}  (q^5{'-' if kodd else '+'}q^3 = {q**5 - q**3 if kodd else q**5 + q**3})")
        print(f"  A = {Aval};  q^6-2q^4+2q^3 = {q**6-2*q**4+2*q**3} (k odd), q^6 = {q**6} (k even)")
        print(f"  Var_H/q^5 = {varH} == q(q-1) = {q*(q-1)}: {varH == q*(q-1)}")
        print(f"  Var_G/q^5 = {varG} == 3-3/q = {Fraction(3*(q-1), q)}: {varG == Fraction(3*(q-1), q)}")
        print(f"  ratio = {ratio} == q^2/3: {ratio == Fraction(q*q, 3)}")
        print(f"  sum_(c3!=0) (K-qI)^2 = {S2} == q^3(q-1)/2 = {q**3*(q-1)//2}: {S2 == q**3*(q-1)//2}")
        print(f"  Phi(y)=sum_x phi(0,x,y)^2 for y!=0: distinct values {sorted(set(Phi[1:]))}  (2q^7 = {2*q**7})")
        vals = sorted(set(r.values()))
        print(f"  r(lambda3) distinct values: {[str(v) for v in vals]}  ((q-3)/(3(q-1)) = {Fraction(q-3, 3*(q-1))})")
        print(f"  1 + sum rho = {1 + rho_sum} == ratio: {1 + rho_sum == ratio}")
    return out


def fourier_check(k):
    """Fourier transform of Psi on G_0 = {c1 = 0} ~ F_q^2 and its predicted values:
         Psi^(mu, 0)  = -q^4                      (mu != 0)
         Psi^(0, nu)  = q^4 (1 - R(0,nu) - 2 L(0,nu))   = -2q^4 / +q^4  (Carlitz values)
         Psi^(mu, nu) = q^4 (1 - R_j - 2 L_j),  j = mu^3/nu^2,  where R_j = #roots in F_q of
                        x^3 + (1+j) x^2 + x + 1 and L_j = #roots with Tr x = 1.
       Then sum_{G_0} phi^2 = q^-2 sum_{(mu,nu) != 0} |Psi^|^2, and Var_G/q^5 = 3 - 3/q is equivalent to
         sum_{j in F_q^*} (1 - 3 L_j - L'_j)^2 = 3q - 5 (k odd) / 3q - 3 (k even),  L'_j = R_j - L_j."""
    q, A, M, INV = field(k)
    psi, K, I, (SQ, SQRT, CUBE, TR, W) = psi_d4_formula(q, A, M, INV)
    chi = 1 - 2 * TR
    P0 = np.array([[int(psi[x * q + y * q * q]) for y in range(q)] for x in range(q)], dtype=object)  # P0[x][y]
    # roots data: for x not in {0,1}, j = (x+1)^3 / x^2
    Lj = np.zeros(q, np.int64); Lpj = np.zeros(q, np.int64)
    for x in range(2, q):
        j = M[CUBE[A[x, 1]], INV[SQ[x]]]
        if TR[x] == 1:
            Lj[j] += 1
        else:
            Lpj[j] += 1
    Rj = Lj + Lpj
    ok = True
    for mu in range(q):
        for nu in range(q):
            if mu == 0 and nu == 0:
                continue
            hat = sum(P0[x][y] * int(chi[A[M[mu, x], M[nu, y]]]) for x in range(q) for y in range(q))
            if nu == 0:
                pred = -q ** 4
            elif mu == 0:
                pred = q ** 4 * (1 - (1 if k % 2 else 3 * int(CUBE[np.arange(1, q)].tolist().count(nu) > 0)) - 2 * (1 if k % 2 else 0))
            else:
                j = M[CUBE[mu], INV[SQ[nu]]]
                pred = q ** 4 * (1 - int(Rj[j]) - 2 * int(Lj[j]))
            if hat != pred:
                ok = False
                print(f"   MISMATCH at (mu,nu)=({mu},{nu}): hat={hat}, pred={pred}")
    print(f"  k={k}: Fourier predictions for all (mu,nu) != 0: {'OK' if ok else 'FAILED'}")
    return ok


def varG_identity_fast(k):
    """O(q) verification of  sum_{j != 0} (1 - 3 L_j - L'_j)^2 = 3q - 5 (k odd) / 3q - 3 (k even)
       using exp/log tables of F_q built directly (polynomial basis)."""
    q = 1 << k
    poly = sum(c << i for i, c in enumerate(MODULI2[k]))
    # multiplication via log tables
    exp = np.zeros(q - 1, np.int64)
    x = 1
    for i in range(q - 1):
        exp[i] = x
        x <<= 1
        if x >> k:
            x ^= poly
    assert len(set(exp.tolist())) == q - 1, "modulus not primitive; pick another"
    log = np.zeros(q, np.int64); log[exp] = np.arange(q - 1)
    xs = np.arange(2, q, dtype=np.int64)
    # trace table
    TR = np.zeros(q, np.int64); y = np.arange(q)
    def sq(a):
        out = np.zeros_like(a); nz = a != 0
        out[nz] = exp[(2 * log[a[nz]]) % (q - 1)]
        return out
    for _ in range(k):
        TR ^= y
        y = sq(y)
    x1 = xs ^ 1
    j = exp[(3 * log[x1] - 2 * log[xs]) % (q - 1)]          # (x+1)^3 / x^2
    Lj = np.bincount(j[TR[xs] == 1], minlength=q)
    Lpj = np.bincount(j[TR[xs] == 0], minlength=q)
    S = int(np.sum((1 - 3 * Lj[1:] - Lpj[1:]) ** 2))
    target = 3 * q - 5 if k % 2 else 3 * q - 3
    return S, target


if __name__ == "__main__":
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    rows = []
    for k in range(1, KMAX + 1):
        rows.append(analyse_formula(k))
    print("\nFourier-side checks (k <= 6):")
    for k in range(1, min(KMAX, 6) + 1):
        fourier_check(k)
    print("\nVar_G identity via roots of x^3+(1+j)x^2+x+1 (O(q) check), k = 1..22:")
    for k in range(1, 23):
        try:
            S, target = varG_identity_fast(k)
            print(f"  k={k:2d} q={1<<k:8d}: sum_j (1-3L_j-L'_j)^2 = {S:8d}, target {target:8d}  {'OK' if S == target else 'FAIL'}")
        except (KeyError, AssertionError) as e:
            print(f"  k={k}: skipped ({e})")
    print("\nsummary (k, q, Var_G/q^5, Var_H/q^5, ratio):")
    for o in rows:
        print(f"  {o['k']:2d} {o['q']:5d}  {float(o['varG']):.6f}  {float(o['varH']):.1f}  {float(o['ratio']):.4f}")
