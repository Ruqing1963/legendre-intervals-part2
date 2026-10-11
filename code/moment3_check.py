"""
moment3_check.py -- Part IV, Proposition 2.1: the third moment of the sextic family.

For the family F_{a,b}(u) = u^6 + a u^2 + b u + b^2 over F_q x F_q^* (q = 2^k) let N_3 = number of roots in F_{q^3}.
Claim (derived in Part IV, Section 2.2):
    sum_{a,b} N_3(F_{a,b}) = (q-1)^2 + (1/q) * sum_{theta in Theta} G_H(theta),
    Theta = {theta in F_{q^3}^* : Tr_{q^3/q}(1/theta) = 0},
    G_H(theta) = sum_{e in F_{q^3}, Tr_{q^3/q}(e) = 0} (-1)^{Tr_{q^3/2}(theta (e^5 + e^3))}
               = (1/q) sum_{a in F_q} S_3(theta, theta, a),  S_3 = the Part II sums over F_{q^3}.
G_H(theta) is the Gauss sum of the F_2-quadratic form e -> Tr(theta e^5 + theta e^3) restricted to the trace-zero
hyperplane; the pair (lambda_3, lambda_5) = (theta, theta) has j = theta^2.
The script computes both sides exactly for k = 3, 5 (left side by factoring, right side by direct summation) and also
prints the distribution of the values G_H(theta) (which are 0 or +-2^{(2k + r)/2}).
Usage: python moment3_check.py 3 5
"""
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_var_check as LV  # noqa: E402
from quintic_check import Poly  # noqa: E402
from legendre_d6_Wp_pattern import factor_degrees  # noqa: E402

LV.PRIMITIVE.update({3: 0b1011, 5: 0b100101, 7: 0b10000011, 9: 0b1000010001, 15: 0b1000000000000011,
                     21: 0b1000000000000000000101})

FIX3 = {(1, 1, 1, 1, 1, 1): 6, (1, 1, 1, 1, 2): 4, (1, 1, 2, 2): 2, (2, 2, 2): 0, (1, 1, 1, 3): 6, (3, 3): 6,
        (1, 2, 3): 4, (1, 1, 4): 2, (2, 4): 0, (1, 5): 1, (6,): 0}   # fixed points of pi^3 = roots in F_{q^3}


def left_side(k):
    q = 1 << k
    P = Poly(q)
    tot = 0
    for l5 in range(1, q):
        for l3 in range(q):
            B = [1, l5, int(P.M[l3, l3]), 0, 0, 0, int(P.M[l5, l5])]
            tot += FIX3[factor_degrees(P, B, 6)]
    return tot


def right_side(k):
    q = 1 << k
    n = 3 * k
    gf = LV.GF(n)
    order = gf.order
    size = gf.size
    exp, log = gf.exp, gf.log
    x = np.arange(1, size, dtype=np.int64)
    logx = log[x]

    def pw(v, e):           # v array of nonzero elements
        return exp[(e * log[v]) % order]

    def mul(u, v):
        out = np.zeros_like(u)
        m = (u != 0) & (v != 0)
        out[m] = exp[(log[u[m]] + log[v[m]]) % order]
        return out

    def frob_q(v):          # v -> v^q
        return exp[(q * log[v]) % order] if np.ndim(v) else (0 if v == 0 else int(exp[(q * log[v]) % order]))

    def tr_q3_q(v):         # Tr_{q^3/q}
        v1 = np.where(v == 0, 0, exp[(q * log[np.maximum(v, 1)]) % order]); v1[v == 0] = 0
        v2 = np.where(v1 == 0, 0, exp[(q * log[np.maximum(v1, 1)]) % order]); v2[v1 == 0] = 0
        return v ^ v1 ^ v2

    all_el = np.arange(size, dtype=np.int64)
    H = all_el[tr_q3_q(all_el) == 0]                       # trace-zero hyperplane, q^2 elements (incl. 0)
    Hnz = H[H != 0]
    e5e3 = np.zeros_like(H)
    e5e3[H != 0] = pw(Hnz, 5) ^ pw(Hnz, 3)                  # e^5 + e^3 on H
    inv = exp[(-logx) % order]                               # 1/x for x != 0
    Theta = x[tr_q3_q(inv) == 0]
    vals = Counter()
    total = 0
    for th in Theta:
        prod = mul(np.full_like(e5e3, th), e5e3)
        par = LV._abs_trace_parity(gf, prod.astype(np.uint32)).astype(np.int64)
        G = int(np.sum(1 - 2 * par))
        vals[G] += 1
        total += G
    assert total % q == 0, total
    return (q - 1) ** 2 + total // q, vals, len(Theta)


if __name__ == "__main__":
    for k in [int(a) for a in sys.argv[1:]] or [3, 5]:
        q = 1 << k
        L = left_side(k)
        R, vals, nth = right_side(k)
        print(f"k={k}, q={q}: sum N_3 by factoring = {L};  (q-1)^2 + q^-1 sum_Theta G_H = {R};  equal: {L == R}")
        print(f"   |Theta| = {nth} (= q^2 - 1: {nth == q * q - 1});  values of G_H(theta) -> multiplicity: "
              f"{dict(sorted(vals.items()))}")
        nonzero = sum(n for g, n in vals.items() if g != 0)
        print(f"   G_H(theta) in {{0, 2^(k+1)}}: {set(vals) <= {0, 2 ** (k + 1)}};  #theta with G_H != 0 = {nonzero}"
              f" (= (q+1)(q-2)/2 = {(q + 1) * (q - 2) // 2}: {nonzero == (q + 1) * (q - 2) // 2});"
              f"  closed form (q-1)^2 + (q+1)(q-2) = {(q - 1) ** 2 + (q + 1) * (q - 2)}: {L == (q - 1) ** 2 + (q + 1) * (q - 2)}")
