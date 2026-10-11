"""
hexad_traces.py -- trace identities on the hexad used in the sign theorems of Part III.

For every (l3,l5) with the six roots of B(u) = l5^2 u^6 + l3^2 u^2 + l5 u + 1 in F_{q^12} we check:
  (i)  in the classes of order 4 (cycle types 4 1^2 and 4 2) the two roots r5, r6 outside the 4-cycle satisfy
       Tr_{q/2}(r5 + r6) = 0  (pair lemma: two roots in F_q, or a conjugate pair in F_{q^2}, have a sum of trace 0);
  (ii) every class has a root r with Tr_{q^12/2}(r) = 0 (reality of the Witt sums S''(l1));
  (iii) Tr_{q^N/2} vanishes on K = ker L' for N = N_0 (the degree of F_q(K)) when N_0 | 4, and Tr_{q^{2N_0}/2} and
        Tr_{q^{4N_0}/2} vanish on K in all classes (used for the vanishing of q^{(12)} on K).
Usage: python hexad_traces.py 3 5
"""
import itertools
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_d6_qf as QF  # noqa: E402
from hexad_search import Field  # noqa: E402


def run(k):
    Fd = Field(k)
    q = Fd.q
    SQ, SQRT = Fd.SQ, Fd.SQRT
    SQ2 = SQ.dot(SQ) % 2
    SQRT2 = SQRT.dot(SQRT) % 2
    N = Fd.N
    I = np.eye(N, dtype=np.uint8)
    cache = {}

    def S(c):
        if c not in cache:
            cache[c] = Fd.K12.linmap_matrix(lambda a: Fd.scal(c, a))
        return cache[c]

    def trq(a):
        return Fd.apply(Fd.TR12q, a)[0]

    def rel_trace(x, m):
        """Tr_{q^12/q^m}(x) for m | 12, as an element of F_{q^12}."""
        t = Fd.zero
        y = x
        for _ in range(12 // m):
            t = Fd.add(t, y)
            for _ in range(m):
                y = Fd.frob(y)
        return t

    stats = Counter()
    for l5 in range(1, q):
        for l3 in range(q):
            Tp = (S(l5).dot(SQ2) + S(l3).dot(SQ) + I + S(int(Fd.F.SQRT[l3])).dot(SQRT)
                  + S(int(Fd.F.SQRT[int(Fd.F.SQRT[l5])])).dot(SQRT2)) % 2
            K = QF.f2_nullspace(SQ2.dot(Tp) % 2)
            if len(K) != 4:
                continue
            one = Fd.K12.tobits(Fd.one)
            sol = QF.f2_nullspace(np.concatenate([Tp, one[:, None]], axis=1))
            x0 = next(v[:-1] for v in sol if v[-1] == 1)
            coset, kern = [], []
            for cf in itertools.product((0, 1), repeat=4):
                v = x0.copy()
                w = np.zeros_like(x0)
                for c, b in zip(cf, K):
                    if c:
                        v ^= b
                        w ^= b
                coset.append(Fd.K12.frombits(v % 2))
                kern.append(Fd.K12.frombits(w % 2))
            sq3, sq5 = int(Fd.M[l3, l3]), int(Fd.M[l5, l5])

            def Bval(u):
                u2 = Fd.mul(u, u)
                u6 = Fd.mul(Fd.mul(u2, u2), u2)
                return Fd.add(Fd.add(Fd.scal(sq5, u6), Fd.scal(sq3, u2)), Fd.add(Fd.scal(l5, u), Fd.one))
            roots = [u for u in coset if all(v == 0 for v in Bval(u))]
            assert len(roots) == 6
            seen, orb = set(), {}
            for r in roots:
                t = tuple(r)
                if t in seen:
                    continue
                o, x = [], t
                while x not in seen:
                    seen.add(x)
                    o.append(x)
                    x = tuple(Fd.frob(list(x)))
                for y in o:
                    orb[y] = len(o)
            cyc = tuple(sorted(orb[tuple(r)] for r in roots))
            N0 = 1
            for l in set(orb.values()):
                N0 = N0 * l // np.gcd(N0, l)
            if cyc == (2, 2, 2):
                N0 = 2
            # (i)
            if 4 in orb.values():
                others = [r for r in roots if orb[tuple(r)] != 4]
                s = Fd.add(others[0], others[1])
                assert Fd.in_Fq(s)
                stats[(cyc, "Tr_{q/2}(r5+r6)", Fd.F.tr2(int(s[0])))] += 1
            # (ii)
            stats[(cyc, "min over roots of Tr_{q^12/2}(r)", min(Fd.abs_trace(r) for r in roots))] += 1
            # (iii): Tr_{q^m/2}(x) for x in K and m in {N0, 2N0, 4N0} intersect divisors of 12 (K lies in F_{q^m})
            for m in sorted({N0, 2 * N0, 4 * N0} & {1, 2, 3, 4, 6, 12}):
                vals = set()
                for x in kern:
                    tr_m = rel_trace(x, m)                      # Tr_{q^12/q^m} x = (12/m) x since x in F_{q^m}
                    # Tr_{q^m/2}(x) = Tr_{q^12/2}(y) for any y with Tr_{q^12/q^m}(y) = x; use that Tr_{q^12/q^m} is
                    # multiplication by 12/m on F_{q^m}: if 12/m is odd this gives Tr_{q^m/2} x = Tr_{q^12/2} x.
                    if (12 // m) % 2 == 1:
                        vals.add(Fd.abs_trace(x))
                    else:
                        # Tr_{q^m/2}(x) = Tr_{q/2}(Tr_{q^m/q} x); compute Tr_{q^m/q} x by summing m Frobenius conjugates
                        t, y = Fd.zero, x
                        for _ in range(m):
                            t = Fd.add(t, y)
                            y = Fd.frob(y)
                        assert Fd.in_Fq(t)
                        vals.add(Fd.F.tr2(int(t[0])))
                stats[(cyc, f"Tr_(q^{m}/2) on K (N0={N0})", tuple(sorted(vals)))] += 1
    print(f"\nk={k}, q={q}")
    for key, n in sorted(stats.items(), key=str):
        print(f"   {key}  x{n}")


if __name__ == "__main__":
    for k in [int(a) for a in sys.argv[1:]] or [3, 5]:
        run(k)
