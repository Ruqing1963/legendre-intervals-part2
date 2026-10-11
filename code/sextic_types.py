"""
sextic_types.py -- factorisation types over F_q of the sextic family  B_{l3,l5}(u) = l5^2 u^6 + l3^2 u^2 + l5 u + 1
(the hexad of ker L'), counted over (l3, l5) in F_q x F_q^*, compared with the S_6-class counts of Part II, Table 3
(now predictions for the sextic):  1^6: (q-2)(q-8)/720, 2 1^4: (q-2)(q-8)/48, 2^2 1^2: q(q-2)/16, 2^3: (q-2)(q-8)/48,
3 1^3: (q+1)(q-2)/18, 3^2: (q+1)(q-2)/18, 3 2 1: (q+1)(q-2)/6, 4 1^2: (q-2)(q-4)/8, 4 2: q(q-2)/8, 5 1: (q^2+1)/5,
6: (q+1)(q-2)/6   -- for l3 != 0;  and for l3 = 0: (4 1^2): (q-2)/2, (2^3): (q-2)/6, (6): (q+1)/3.
Usage: python sextic_types.py 3 5 7 9
"""
import sys
import time
from collections import Counter
from fractions import Fraction

from quintic_check import Poly
from legendre_d6_Wp_pattern import factor_degrees

PRED = {(1, 1, 1, 1, 1, 1): lambda q: Fraction((q - 2) * (q - 8), 720), (1, 1, 1, 1, 2): lambda q: Fraction((q - 2) * (q - 8), 48),
        (1, 1, 2, 2): lambda q: Fraction(q * (q - 2), 16), (2, 2, 2): lambda q: Fraction((q - 2) * (q - 8), 48),
        (1, 1, 1, 3): lambda q: Fraction((q + 1) * (q - 2), 18), (3, 3): lambda q: Fraction((q + 1) * (q - 2), 18),
        (1, 2, 3): lambda q: Fraction((q + 1) * (q - 2), 6), (1, 1, 4): lambda q: Fraction((q - 2) * (q - 4), 8),
        (2, 4): lambda q: Fraction(q * (q - 2), 8), (1, 5): lambda q: Fraction(q * q + 1, 5), (6,): lambda q: Fraction((q + 1) * (q - 2), 6)}
PRED0 = {(1, 1, 4): lambda q: Fraction(q - 2, 2), (2, 2, 2): lambda q: Fraction(q - 2, 6), (6,): lambda q: Fraction(q + 1, 3)}


def run(k):
    t0 = time.time()
    q = 1 << k
    P = Poly(q)
    M = P.M
    mixed, pure = Counter(), Counter()
    for l5 in range(1, q):
        l52 = int(M[l5, l5])
        for l3 in range(q):
            l32 = int(M[l3, l3])
            B = [1, l5, l32, 0, 0, 0, l52]
            typ = factor_degrees(P, B, 6)
            (mixed if l3 else pure)[typ] += 1
    print(f"\nk={k}, q={q}  [{time.time()-t0:.1f}s]   type of B over F_q -> count (prediction)")
    allok = True
    for typ, n in sorted(mixed.items()):
        pred = PRED.get(typ, lambda q: None)(q)
        ok = (pred is not None and pred == n)
        allok &= ok
        print(f"   l3 != 0  {str(typ):22s} {n:8d}   {str(pred):>12s}  {'OK' if ok else 'MISMATCH'}")
    for typ, n in sorted(pure.items()):
        pred = PRED0.get(typ, lambda q: None)(q)
        ok = (pred is not None and pred == n)
        allok &= ok
        print(f"   l3 = 0   {str(typ):22s} {n:8d}   {str(pred):>12s}  {'OK' if ok else 'MISMATCH'}")
    print(f"   all counts match the Table-3 polynomials: {allok}")


if __name__ == "__main__":
    for k in [int(a) for a in sys.argv[1:]] or [3, 5, 7]:
        run(k)
