"""Tabulate W(l3,l5)/q^7 (and W') against algebraic invariants of (l3, l5):
   j = l3^5 / l5^3 (invariant under (l3,l5) -> (s^3 l3, s^5 l5)),  Tr(j), Tr(1/j), j cube?, j fifth power?, l5 fifth power?"""
import sys
from collections import Counter, defaultdict

import legendre_d6_qf as QF

for k in [int(a) for a in sys.argv[1:]] or [1, 2, 3, 4, 5, 6]:
    q = 1 << k
    A6, X, res = QF.run(k, verbose=False)
    F = QF.Fq(q)
    A, M, INV, TR = F.A, F.M, F.INV, F.TR

    def pw(x, e):
        r = 1
        for _ in range(e):
            r = int(M[r, x])
        return r

    cubes = {pw(x, 3) for x in range(1, q)}
    fifths = {pw(x, 5) for x in range(1, q)}
    print(f"\n### k={k}, q={q}:  |cubes|={len(cubes)}, |fifth powers|={len(fifths)}")
    # pure quintic: W(0,l5) vs l5 fifth power / Tr(l5)
    tab = Counter()
    for (l3, l5), (v, r, info, n) in res["W"].items():
        if l3 == 0 and l5 != 0:
            tab[(v // q ** 7 if v % q ** 7 == 0 else v / q ** 7, l5 in fifths, l5 in cubes)] += n
    print("  W(0,l5)/q^7 by (value, l5 fifth power?, l5 cube?):", dict(tab))
    tab = Counter()
    for (l3, l5), (v, r, info, n) in res["W"].items():
        if l5 == 0 and l3 != 0:
            tab[(v // q ** 7 if v % q ** 7 == 0 else v / q ** 7, l3 in cubes)] += n
    print("  W(l3,0)/q^7 by (value, l3 cube?):", dict(tab))
    # mixed: expand orbits to all pairs and test invariants
    full = {}
    for (l3, l5), (v, r, info, n) in res["W"].items():
        if l3 and l5:
            for s in range(1, q):
                s3 = pw(s, 3); s5 = pw(s, 5)
                a, b = int(M[s3, l3]), int(M[s5, l5])
                for _ in range(k):
                    full[(a, b)] = v
                    a, b = int(F.SQ[a]), int(F.SQ[b])
    assert len(full) == (q - 1) ** 2
    byinv = defaultdict(Counter)
    for (l3, l5), v in full.items():
        j = int(M[pw(l3, 5), INV[pw(l5, 3)]])
        val = v // q ** 7 if v % q ** 7 == 0 else v / q ** 7
        byinv["Tr(j)"][(val, int(TR[j]))] += 1
        byinv["Tr(1/j)"][(val, int(TR[INV[j]]))] += 1
        byinv["j cube"][(val, j in cubes)] += 1
        byinv["j fifth"][(val, j in fifths)] += 1
        byinv["j cube & Tr(j)"][(val, j in cubes, int(TR[j]))] += 1
        byinv["j cube & Tr(1/j)"][(val, j in cubes, int(TR[INV[j]]))] += 1
    for name, c in byinv.items():
        print(f"  mixed W/q^7 by (value, {name}):", dict(sorted(c.items(), key=str)))
