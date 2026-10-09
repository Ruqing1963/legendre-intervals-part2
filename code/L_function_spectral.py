"""
L_function_spectral.py -- what the five exact values  m_k = Y_k / 8^k = 13, -1, -35, 47, 733  (k = 1..5) do and do
not say about the twisted L-function  L_sigma(T)  of Theorem 7.7  (sum_k Y_k T^k = T d/dT log L_sigma(T)).

1. Linear recurrences: Berlekamp-Massey over Q for (Y_k) and (m_k); exact fits of order d = 1, 2 from the first 2d
   terms with the remaining terms as a test; order 3 is under-determined by five terms (two equations, three unknowns),
   so five terms carry no information beyond "order >= 3", which is the generic linear complexity of five numbers.
2. Weil compatibility: the characteristic roots of every fitted recurrence are compared with 2^{w/2}, w = 0..18
   (Y) and 2^{(w-6)/2} (m): none is a Weil number, as expected for an over-fitted recurrence.
3. Integrality (2-adic) constraints that ARE rigorous: if every Frobenius eigenvalue on H_c^*(X, L) had the
   "supersingular shape" 2^{w/2} * (root of unity) -- the shape that holds for the sums behind A_6 and X -- then
   (a) a single weight w cannot carry Y: Y_2 = -64 forces 2^w | 64, i.e. w <= 6;  (b) m_k odd for all k <= 5 forces
   a nonzero part of weight <= 6 (2-adic valuation of Y_k is exactly 3k).  The paper's Proposition 7.8 (pure weight 9,
   rank 5 excluded by e_5 = 2281 != 2^{15/2}) is reproduced.
4. The point-count dictionary (Proposition 7.9): Y = N_0 - N_1 with N_c = #{alpha in U : e_4(alpha) = c};
   N_1 = (A_6 - Y)/q, N_0 = (A_6 + (q-1) Y)/q, and for odd k  N_1 = q^8 + q^5 - q^4 - q^2 m_k.
Usage: python L_function_spectral.py
"""
from fractions import Fraction
from itertools import product
import math

import numpy as np

M = {1: 13, 2: -1, 3: -35, 4: 47, 5: 733}
Y = {k: 8 ** k * m for k, m in M.items()}
A6 = {k: (1 << k) ** 9 + (-1) ** (k + 1) * ((1 << k) - 1) * (1 << k) ** 5 for k in M}   # exact for k <= 7


def berlekamp_massey(seq):
    """minimal LFSR over Q: returns (L, connection polynomial C(x) = 1 + c1 x + ... + cL x^L) with
       s_n + c1 s_{n-1} + ... + cL s_{n-L} = 0 for n >= L."""
    seq = [Fraction(s) for s in seq]
    C = [Fraction(1)]; B = [Fraction(1)]; L = 0; m = 1; b = Fraction(1)
    for n in range(len(seq)):
        d = seq[n] + sum(C[i] * seq[n - i] for i in range(1, L + 1))
        if d == 0:
            m += 1
            continue
        T = C[:]
        coef = d / b
        C = C + [Fraction(0)] * (len(B) + m - len(C))
        for i in range(len(B)):
            C[i + m] -= coef * B[i]
        if 2 * L <= n:
            L, B, b, m = n + 1 - L, T, d, 1
        else:
            m += 1
    return L, C[:L + 1]


def fit_order(seq, d):
    """solve s_n = sum_{i=1}^d c_i s_{n-i} for n = d+1 .. 2d (exactly), return c and the predictions for n > 2d"""
    if 2 * d > len(seq):
        return None
    A = np.array([[Fraction(seq[n - i]) for i in range(1, d + 1)] for n in range(d, 2 * d)], dtype=object)
    bvec = [Fraction(seq[n]) for n in range(d, 2 * d)]
    # Gaussian elimination over Fractions
    n = d
    Mx = [list(A[r]) + [bvec[r]] for r in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if Mx[r][col] != 0), None)
        if piv is None:
            return None
        Mx[col], Mx[piv] = Mx[piv], Mx[col]
        pv = Mx[col][col]
        Mx[col] = [v / pv for v in Mx[col]]
        for r in range(n):
            if r != col and Mx[r][col] != 0:
                f = Mx[r][col]
                Mx[r] = [a - f * b2 for a, b2 in zip(Mx[r], Mx[col])]
    c = [Mx[r][n] for r in range(n)]
    preds = []
    s = [Fraction(v) for v in seq]
    for idx in range(2 * d, len(seq)):
        preds.append(sum(c[i - 1] * s[idx - i] for i in range(1, d + 1)))
    return c, preds


def roots_of(c):
    """roots of x^d - c1 x^{d-1} - ... - cd"""
    d = len(c)
    poly = [1.0] + [-float(ci) for ci in c]
    return np.roots(poly)


def weil_label(r, scale=1.0):
    a = abs(r) * scale
    if a == 0:
        return "0"
    w = 2 * math.log2(a)
    return f"|.|=2^{w:.3f}" + ("  (Weil weight %d)" % round(w) if abs(w - round(w)) < 1e-6 else "  (not 2^{w/2})")


def main():
    ks = sorted(M)
    mseq = [M[k] for k in ks]; yseq = [Y[k] for k in ks]
    print("m_k =", mseq, "   Y_k =", yseq)
    for name, seq in (("m_k", mseq), ("Y_k", yseq)):
        L, C = berlekamp_massey(seq)
        print(f"\n[{name}] Berlekamp-Massey: minimal order L = {L} (five terms: the generic value is 3; 2L > 5 means the "
              f"polynomial is not unique).  One connection polynomial: {[str(c) for c in C]}")
        for d in (1, 2):
            res = fit_order(seq, d)
            c, preds = res
            print(f"   order {d} from terms 1..{2*d}: c = {[str(x) for x in c]};  predicted terms {2*d+1}..5 = "
                  f"{[str(p) for p in preds]}  vs actual {seq[2*d:]}  -> {'refuted' if any(p != a for p, a in zip(preds, seq[2*d:])) else 'consistent'}")
            for r in roots_of(c):
                print(f"      root {r:.4f}: {weil_label(r, 8.0 if name == 'm_k' else 1.0)}  [as eigenvalue of F on Y]")
        print("   order 3: two equations (n = 4, 5) in three unknowns -> a one-parameter family; no test possible.")

    # ---- rigorous 2-adic constraints under the supersingular shape
    print("\nSupersingular-shape constraints (eigenvalues 2^{w/2} * root of unity):")
    for w in range(5, 12):
        a = Fraction(Y[2], 2 ** w)
        print(f"   pure weight w={w}: Y_2 / 2^w = {a}  -> {'possible' if a.denominator == 1 else 'IMPOSSIBLE (not an algebraic integer)'}")
    v2 = lambda n: (n & -n).bit_length() - 1
    print("   v_2(Y_k) =", [v2(abs(Y[k])) for k in ks], " = 3k exactly, i.e. m_k odd:", [M[k] % 2 for k in ks])
    print("   => under the shape hypothesis a nonzero part of weight <= 6 (the 'q^3' twist itself) is forced.")
    # Newton: e_j of the normalised numbers from p_j = m_j (Proposition 7.8)
    p = [None] + mseq
    e = [1]
    for j in range(1, 6):
        s = sum((-1) ** (i - 1) * e[j - i] * p[i] for i in range(1, j + 1))
        e.append(Fraction(s, j))
    print("   Newton: e_1..e_5 of the eta_i from p_j = m_j:", [str(x) for x in e[1:]], "; |e_5| should be 2^{15/2} =",
          f"{2**7.5:.2f} if b_9 = 5 and all |eta_i| = 2^{3/2}  -> excluded (Prop. 7.8)")

    # ---- the point-count dictionary
    print("\nPoint counts N_c = #{alpha in U : e_4(alpha) = c}:  Y = N_0 - N_1,  A_6 = N_0 + (q-1) N_1")
    for k in ks:
        q = 1 << k
        N1 = Fraction(A6[k] - Y[k], q); N0 = Fraction(A6[k] + (q - 1) * Y[k], q)
        assert N1.denominator == 1 and N0.denominator == 1
        main = q ** 8 + (-1) ** (k + 1) * (q ** 5 - q ** 4)
        print(f"   k={k}: N_1 = {int(N1)} = q^8 +- (q^5 - q^4) - q^2 m_k  [check: {int(N1) == main - q*q*M[k]}],   N_0 = {int(N0)},"
              f"   (N_1 - main)/q^{3.5} = {float(N1 - main) / q**3.5:+.3f},  (N_1 - main)/q^4 = {float(N1 - main) / q**4:+.3f}")
    print("   (epsilon_0 m_k - 1)/12 =", [Fraction((-1) ** (k + 1) * M[k] - 1, 12) for k in ks], " (integers for k <= 5)")


if __name__ == "__main__":
    main()
