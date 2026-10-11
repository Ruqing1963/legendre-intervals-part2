"""
part4_brown_invariants.py -- Part IV, Open Problem 2 for the classes of order 4.

Part III, Theorem 6.12: for the classes 4 1^2 and 4 2 (ker L' in F_{q^4}),  S''(l1) = q^4 s_4(l1)  with
    s_4(l1) = sum_{x in F_{q^4}} i^{q^(4)(x)},   q^(4)(x) = Tr_{W_2(F_{q^4})/Z4}[x] + 2 Tr_{q^4/2}(l5 x^5 + l3 x^3 + l1 x),
and s_4 in {0, +-4q^2}; it is nonzero on the hyperplane Tr_{q/2}(l1 tau_0) = kappa_0, tau_0 = r5 + r6.
Conjecture 6.17 says sign s_4 = +1 in class 4 1^2 (r5, r6 in F_q) and -1 in class 4 2 (r5, r6 conjugate in F_{q^2}).
This script, for every pair of these two classes at q = 8 (and a sample at q = 32):
  * finds the hexad in F_{q^12} (roots of B in the coset T'^{-1}(1)) and the 4-cycle / the pair {r5, r6};
  * computes s_4(l1) for all l1 by vectorised Witt sums over F_{2^{4k}} (the same routine as witt_lpoly.py);
  * checks the shift calculus: s_4(l1 + a) = i^{-q^(4)(c)} s_4(l1) for c in F_{q^4} with T'(c) = a, so that one
    representative c per coset of T'(F_q) determines everything;
  * tabulates candidate invariants against the sign:  I_1 = Tr_{q/2}(r5 r6 / tau_0^2) (= 0 iff u^2 + tau_0 u + r5 r6
    splits, i.e. iff class 4 1^2),  I_2 = q^(4)(r5) (Z/4 value of the form at a root outside the 4-cycle, l1 = 0),
    I_3 = q^(4)(r1) at a root of the 4-cycle,  kappa_0 = q^(4)_{l1=0}(r1 + r5)/2.
Usage: python part4_brown_invariants.py 3        (k = 5 with a sample: python part4_brown_invariants.py 5 6)
"""
import itertools
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_d6_qf as QF  # noqa: E402
import legendre_var_check as LV  # noqa: E402
from hexad_search import Field  # noqa: E402
from quintic_check import Poly  # noqa: E402
from legendre_d6_Wp_pattern import factor_degrees  # noqa: E402
from witt_lpoly import witt_sums  # noqa: E402  (also installs the odd-degree primitive polynomials)


def embedding(k, gf):
    """map F_q (fqlib encoding) -> the subfield of order q inside gf (dict int -> int), as in witt_lpoly.run."""
    q = 1 << k
    F = QF.Fq(q)
    step = gf.order // (q - 1)
    for t in range(q - 1):
        cand = int(gf.exp[(t * step) % gf.order])
        tab = {0: 0}
        a_pow, c_pow = 1, 1
        for j in range(q - 1):
            tab[a_pow] = c_pow
            a_pow = int(F.M[a_pow, 2])
            c_pow = int(gf.exp[(gf.log[c_pow] + gf.log[cand]) % gf.order])
        ok = all(tab[int(F.A[a, b])] == (tab[a] ^ tab[b]) for a in range(1, min(q, 64)) for b in range(1, min(q, 64)))
        if ok:
            return tab
    raise RuntimeError("no embedding found")


def run(k, max_pairs=None):
    Fd = Field(k)
    q = Fd.q
    F = Fd.F
    P = Poly(q)
    gf4 = LV.GF(4 * k)
    emb = embedding(k, gf4)
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

    def in_Fq_val(a):
        assert Fd.in_Fq(a), a
        return int(a[0])

    def witt_Z4(x, l1, l3, l5):
        """q^(4)(x) for x in F_{q^4} given as an element of F_{q^12} (Fd encoding): returns value in Z/4.
        Uses the identity q^(12)(x) = q^(4)(x) + 2 Tr_{q^4/2}(x) (Part III, Lemma 6.5) and computes q^(12) in F_{q^12}."""
        # q^(12)(x) = T_{12k}(x) + 2 Tr_{q^12/2}(f(x)) with T_{12k}(x) = Tr_{q^12/2}(x) + 2 e2^abs(x)
        conj = [x]
        for _ in range(12 * k - 1):
            conj.append(Fd.K12.sqr(conj[-1]) if hasattr(Fd.K12, "sqr") else Fd.mul(conj[-1], conj[-1]))
        s0, s1 = Fd.zero, Fd.zero
        for c in conj:                          # Witt addition of Teichmueller lifts: (s0,s1) += (c,0)
            s1 = Fd.add(s1, Fd.mul(s0, c))
            s0 = Fd.add(s0, c)
        T = (in_Fq_val(s0) + 2 * in_Fq_val(s1)) % 4     # s0, s1 in F_2
        x2 = Fd.mul(x, x); x3 = Fd.mul(x2, x); x5 = Fd.mul(x3, x2)
        f = Fd.add(Fd.add(Fd.scal(l5, x5), Fd.scal(l3, x3)), Fd.scal(l1, x))
        q12 = (T + 2 * Fd.abs_trace(f)) % 4
        # Tr_{q^4/2}(x) = Tr_{q/2}(Tr_{q^4/q} x)
        t, y = Fd.zero, x
        for _ in range(4):
            t = Fd.add(t, y)
            y = Fd.frob(y)
        tr4 = F.tr2(in_Fq_val(t))
        return (q12 - 2 * tr4) % 4

    stats = Counter()
    done = Counter()
    for l5 in range(1, q):
        for l3 in range(q):
            B6 = [1, l5, int(P.M[l3, l3]), 0, 0, 0, int(P.M[l5, l5])]
            cls = factor_degrees(P, B6, 6)
            if cls not in ((1, 1, 4), (2, 4)):
                continue
            if max_pairs and done[cls] >= max_pairs:
                continue
            done[cls] += 1
            Tp = (S(l5).dot(SQ2) + S(l3).dot(SQ) + I + S(int(F.SQRT[l3])).dot(SQRT)
                  + S(int(F.SQRT[int(F.SQRT[l5])])).dot(SQRT2)) % 2
            K = QF.f2_nullspace(SQ2.dot(Tp) % 2)
            assert len(K) == 4
            one = Fd.K12.tobits(Fd.one)
            sol = QF.f2_nullspace(np.concatenate([Tp, one[:, None]], axis=1))
            x0 = next(v[:-1] for v in sol if v[-1] == 1)
            coset = []
            for cf in itertools.product((0, 1), repeat=4):
                v = x0.copy()
                for c, b in zip(cf, K):
                    if c:
                        v ^= b
                coset.append(Fd.K12.frombits(v % 2))
            sq3, sq5 = int(Fd.M[l3, l3]), int(Fd.M[l5, l5])

            def Bval(u):
                u2 = Fd.mul(u, u); u6 = Fd.mul(Fd.mul(u2, u2), u2)
                return Fd.add(Fd.add(Fd.scal(sq5, u6), Fd.scal(sq3, u2)), Fd.add(Fd.scal(l5, u), Fd.one))
            roots = [u for u in coset if all(v == 0 for v in Bval(u))]
            assert len(roots) == 6
            orb = {}
            seen = set()
            for r in roots:
                t = tuple(r)
                if t in seen:
                    continue
                o, y = [], t
                while y not in seen:
                    seen.add(y); o.append(y); y = tuple(Fd.frob(list(y)))
                for z in o:
                    orb[z] = len(o)
            four = [list(r) for r in roots if orb[tuple(r)] == 4]
            other = [list(r) for r in roots if orb[tuple(r)] != 4]
            r1, r5, r6 = four[0], other[0], other[1]
            tau0 = Fd.add(r5, r6)
            p56 = Fd.mul(r5, r6)
            t0 = in_Fq_val(tau0)
            # I_1 = Tr_{q/2}(r5 r6 / tau0^2): both r5 r6 and tau0 are in F_q
            p = in_Fq_val(p56)
            inv_t2 = int(F.INV[int(F.M[t0, t0])])
            I1 = F.tr2(int(F.M[p, inv_t2]))
            I2 = witt_Z4(r5, 0, l3, l5)
            I3 = witt_Z4(r1, 0, l3, l5)
            kappa0 = witt_Z4(Fd.add(r1, r5), 0, l3, l5)
            assert kappa0 % 2 == 0
            kappa0 //= 2
            # s_4(l1) for all l1
            s4 = {l1: witt_sums(k, 4, l1, l3, l5, gf4, emb) for l1 in range(q)}
            nz = {l1: v for l1, v in s4.items() if abs(v) > 1e-9}
            signs = {int(round(v.real / (4 * q * q))) for v in nz.values()}
            assert all(abs(v.imag) < 1e-9 and abs(abs(v.real) - 4 * q * q) < 1e-6 for v in nz.values())
            # hyperplane check: nonzero iff Tr(l1 tau0) = kappa0
            hyper_ok = all((F.tr2(int(F.M[l1, t0])) == kappa0) == (l1 in nz) for l1 in range(q))
            # shift calculus with c = r1 (T'(r1) = 1): s_4(l1 + 1) = i^{-q^(4)_{l1}(r1)} s_4(l1)
            shift_ok = True
            for l1 in range(q):
                ph = witt_Z4(r1, l1, l3, l5)
                lhs = s4[(l1 ^ 1)]
                rhs = (1j ** (-ph)) * s4[l1]
                if abs(lhs - rhs) > 1e-6:
                    shift_ok = False
            key = (cls, "sign s_4", tuple(sorted(signs)), "#nonzero", len(nz), "hyperplane ok", hyper_ok, "shift ok", shift_ok,
                   "I1=Tr(r5r6/tau0^2)", I1, "I2=q4(r5)", I2, "I3=q4(r1)", I3, "kappa0", kappa0)
            stats[key] += 1
    print(f"\nk={k}, q={q}: classes 4 1^2 = (1,1,4) and 4 2 = (2,4)")
    for key, n in sorted(stats.items(), key=str):
        print(f"   {key}  x{n}")


if __name__ == "__main__":
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    mp = int(sys.argv[2]) if len(sys.argv) > 2 else None
    run(k, mp)
