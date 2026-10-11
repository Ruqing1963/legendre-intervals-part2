"""
hexad_signs.py -- the radical and the values of Q' through the hexad.

For (l3,l5) with all six roots of B(u) = l5^2 u^6 + l3^2 u^2 + l5 u + 1 in F_{q^12}:
  * the roots satisfy T'(r) = 1 (since L'(r) = r^4 B(r)^2 + B(r) + 1 = 1), so they lie in the coset T'^{-1}(1) = r1 + K;
  * we compute, for every root r: Tr_{q^12/q}(r), Q'(r) = Tr(l5 r^5 + l3 r^3) + Tr_{q/2} e2(r), and the orbit of r;
  * the radical R'' = {x in F_{q^12} : T'(x) = Tr_{q^12/q}(x)} of the x-form of S'(l1), its dimension r, and the set of
    l1 for which Q' + Tr(l1 .) vanishes on R'' (these l1 give |S'(l1)| = 2^{(12k+r)/2}, the others S'(l1) = 0);
  * the exact W'(l3,l5) from legendre_d6_qf, hence the sign pattern: W' = (#nonzero S') * sign * 2^{(12k+r)/2} if all
    nonzero S'(l1) have the same sign.
Usage: python hexad_signs.py 3 5
"""
import itertools
import os
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_d6_qf as QF  # noqa: E402
from hexad_search import Field  # noqa: E402


def wp_table(k):
    """W'(l3,l5) for every pair, from legendre_d6_qf (orbit representatives under (l3,l5) -> (l3^2,l5^2))."""
    _, _, results = QF.run(k, verbose=False)
    q = 1 << k
    F = QF.Fq(q)
    table = {}
    for (l3, l5), (val, rdim, info, n) in results["Wp"].items():
        a, b = l3, l5
        for _ in range(k):
            table[(a, b)] = (val, rdim)
            a, b = int(F.SQ[a]), int(F.SQ[b])
    return table


def run(k):
    Fd = Field(k)
    q = Fd.q
    A, M = Fd.A, Fd.M
    Wp = wp_table(k)
    # T' matrix and Tr_{q^12/q} as an endomorphism (F_q embedded as constants)
    SQRT, SQ = Fd.SQRT, Fd.SQ
    SQRT2 = SQRT.dot(SQRT) % 2
    SQ2 = SQ.dot(SQ) % 2
    N = Fd.N
    I = np.eye(N, dtype=np.uint8)

    def scal_mat(c):
        return Fd.K12.linmap_matrix(lambda a: Fd.scal(c, a))

    sc_cache = {}

    def S(c):
        if c not in sc_cache:
            sc_cache[c] = scal_mat(c)
        return sc_cache[c]

    def frob_pow(a, e):
        for _ in range(e):
            a = Fd.frob(a)
        return a

    def trq(a):                      # Tr_{q^12/q}(a) in F_q
        return Fd.apply(Fd.TR12q, a)[0]

    def Qprime(x, l3, l5):
        x2 = Fd.mul(x, x); x3 = Fd.mul(x2, x); x5 = Fd.mul(x3, x2)
        t = Fd.F.tr2(trq(Fd.add(Fd.scal(l5, x5), Fd.scal(l3, x3))))
        conj = [x]
        for _ in range(11):
            conj.append(Fd.frob(conj[-1]))
        e2 = Fd.zero
        for i in range(12):
            for j in range(i + 1, 12):
                e2 = Fd.add(e2, Fd.mul(conj[i], conj[j]))
        assert Fd.in_Fq(e2)
        return (t + Fd.F.tr2(e2[0])) % 2

    stats = Counter()
    for l5 in range(1, q):
        for l3 in range(q):
            sq3, sq5 = int(M[l3, l3]), int(M[l5, l5])
            Tp = (S(l5).dot(SQ2) + S(l3).dot(SQ) + I + S(int(Fd.F.SQRT[l3])).dot(SQRT) + S(int(Fd.F.SQRT[int(Fd.F.SQRT[l5])])).dot(SQRT2)) % 2
            Lp = SQ2.dot(Tp) % 2                 # L'(x) = T'(x)^4 (fourth power of the value, not the iterate)
            K = QF.f2_nullspace(Lp)
            if len(K) != 4:
                continue
            # particular solution of T'(x) = 1
            one = Fd.K12.tobits(Fd.one)
            aug = np.concatenate([Tp, one[:, None]], axis=1)
            sol = QF.f2_nullspace(np.concatenate([Tp, one[:, None]], axis=1))  # kernel of [T' | 1] gives x with T'x = 1 when last coord = 1
            x0 = None
            for v in sol:
                if v[-1] == 1:
                    x0 = v[:-1]
                    break
            if x0 is None:
                stats[("no T'=1 solution",)] += 1
                continue
            coset = []
            for coeffs in itertools.product((0, 1), repeat=4):
                v = x0.copy()
                for c, b in zip(coeffs, K):
                    if c:
                        v ^= b
                coset.append(Fd.K12.frombits(v % 2))

            def Bval(u):
                u2 = Fd.mul(u, u); u6 = Fd.mul(Fd.mul(u2, u2), u2)
                return Fd.add(Fd.add(Fd.scal(sq5, u6), Fd.scal(sq3, u2)), Fd.add(Fd.scal(l5, u), Fd.one))
            roots = [u for u in coset if all(v == 0 for v in Bval(u))]
            if len(roots) != 6:
                stats[("roots in coset", len(roots))] += 1
                continue
            # orbits of Frobenius on the roots
            rt = [tuple(r) for r in roots]
            seen, cyc, orbit_of = set(), [], {}
            for r in rt:
                if r in seen:
                    continue
                o, x = [], r
                while x not in seen:
                    seen.add(x); o.append(x); x = tuple(Fd.frob(list(x)))
                cyc.append(len(o))
                for y in o:
                    orbit_of[y] = len(o)
            cls = tuple(sorted(cyc))
            # per-root data: (orbit length, Tr_{q^12/q} r == 0 ?, Q'(r))
            rdata = tuple(sorted((orbit_of[tuple(r)], int(trq(r) == 0), Qprime(r, l3, l5)) for r in roots))
            # radical R'' = ker(T' + Tr)
            Rmat = (Tp + Fd.TR12q) % 2
            R = QF.f2_nullspace(Rmat)
            rdim = len(R)
            # Q' + Tr(l1 .) on R'': value of Q' on a basis and Tr_{q^12/q} of basis elements
            Rel = [Fd.K12.frombits(v) for v in R]
            qvals = [Qprime(x, l3, l5) for x in Rel]
            trs = [trq(x) for x in Rel]
            # count l1 with Q'(x) + Tr_{q/2}(l1 * Tr x) = 0 for all basis x -- but Q' is only additive on R'' if R'' is
            # the radical of the bilinear form, which it is; so check on basis vectors
            n_l1 = 0
            for l1 in range(q):
                if all((qv + Fd.F.tr2(int(M[l1, t]))) % 2 == 0 for qv, t in zip(qvals, trs)):
                    n_l1 += 1
            val = Wp[(l3, l5)][0]
            unit = 2 ** ((12 * k + rdim) // 2)
            implied_sign = None
            if n_l1 and val % unit == 0 and abs(val) == n_l1 * unit:
                implied_sign = 1 if val > 0 else -1
            stats[(cls, rdim, n_l1, val // q ** 7 if val % q ** 7 == 0 else val, implied_sign, rdata)] += 1
    print(f"\nk={k}: (class, dim R'', #l1 with S'!=0, W'/q^7, common sign if |W'| = #l1 * 2^((12k+r)/2), "
          f"roots: sorted (orbit length, Tr r = 0?, Q'(r))) -> #pairs")
    for key, n in sorted(stats.items(), key=str):
        print(f"   {key}  x{n}")


if __name__ == "__main__":
    for k in [int(a) for a in sys.argv[1:]] or [3]:
        run(k)
