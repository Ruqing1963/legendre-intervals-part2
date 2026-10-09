"""
legendre_d6_Wp_pattern.py -- fine structure of the sums W'(l3,l5) that give X = q^-3 sum W'(l3,l5).

W'(l3,l5) = sum_beta (-1)^{Tr(l3 beta^2 theta beta + l5 beta^4 theta beta + beta^{q+1} + beta)} is the beta-form of
   sum_{l1} S'(l1,l3,l5),  S'(l) = sum_{x in F_{q^12}} (-1)^{Tr(l5 x^5 + l3 x^3 + l1 x) + Tr_{q/2} e_2(x)}.
The x-form polar operator is T'(x) = l5 x^4 + (l5 x)^{1/4} + l3 x^2 + (l3 x)^{1/2} + x + Tr_{q^12/q}(x), and the
radical dimension of W' is 2k + d with d = dim ker(P | K'), K' = ker L' cap F_{q^12},
   L'(x) = T'(x)^4 - (trace part) = l5^4 x^16 + l3^4 x^8 + x^4 + l3^2 x^2 + l5 x,   P = 1 + s^2 + ... + s^10.
For every (l3,l5) we compute the factorisation type over F_q of g(u) = L'(u)/u (degrees of the irreducible factors,
= Frobenius orbit structure on the 15 nonzero elements of ker L') and tabulate it against (d, W'/q^7) using the exact
values of legendre_d6_qf.py.  Usage: python legendre_d6_Wp_pattern.py 3 5 7
"""
import sys
from collections import Counter
from fractions import Fraction

import legendre_d6_qf as QF
from quintic_check import Poly


def factor_degrees(P, g, maxdeg):
    """multiset of degrees of the irreducible factors of the squarefree polynomial g over F_q"""
    q = P.q
    x = [0, 1]
    xq = x
    counts = {}
    rem = g[:]
    degs = []
    for m in range(1, maxdeg + 1):
        xq = P.powmod(xq, q, g)                     # x^{q^m} mod g
        diff = P.trim([int(P.A[a, b]) for a, b in zip(xq + [0] * (len(x) - len(xq)), x + [0] * (len(xq) - len(x)))])
        h = P.gcd(g, diff) if diff else g
        Nm = len(h) - 1 if h else 0                 # number of roots in F_{q^m}
        counts[m] = Nm
        a_m = Nm - sum(e * counts_a for e, counts_a in degs_e(degs) if m % e == 0 and e < m)
        a_m //= m
        degs += [m] * a_m
        if sum(degs) == len(g) - 1:
            break
    return tuple(sorted(degs))


def degs_e(degs):
    c = Counter(degs)
    return list(c.items())


def run(k):
    q = 1 << k
    A6, X, results = QF.run(k, verbose=False)
    P = Poly(q)
    pw = lambda a, e: int(QF.Fq(q).M[a, a]) if False else None  # noqa (unused)
    M, A = P.M, P.A

    def sq(a):
        return int(M[a, a])

    tab = Counter()
    for (l3, l5), (val, rdim, info, n) in results["Wp"].items():
        if l3 == 0 and l5 == 0:
            typ = "0,0"
            g = None
        else:
            # g(u) = l5^4 u^15 + l3^4 u^7 + u^3 + l3^2 u + l5
            l34, l54, l32 = sq(sq(l3)), sq(sq(l5)), sq(l3)
            g = [0] * 16
            g[15] = l54; g[7] = l34; g[3] = 1; g[1] = l32; g[0] = l5
            g = P.trim(g)
            if l5 == 0:                               # g = u (l3^4 u^6 + u^2 + l3^2) -> u * (l3^2 u^3 + u + l3)^2
                g = P.trim([l3, 1, 0, l32])
            typ = factor_degrees(P, g, 15)
            typ = ("l3=0" if l3 == 0 else "l5=0" if l5 == 0 else "mixed", typ)
        d = rdim - 2 * k
        tab[(typ, d, Fraction(val, q ** 7))] += n
    print(f"\nk={k}, q={q}:  (type of g over F_q, d = r - 2k, W'/q^7) -> #(l3,l5)   [X = {X} = q^6+q^5: {X == q**6 + q**5}]")
    for key in sorted(tab, key=str):
        print(f"   {key}  x{tab[key]}")
    return tab


if __name__ == "__main__":
    for k in [int(a) for a in sys.argv[1:]] or [3, 5]:
        run(k)
