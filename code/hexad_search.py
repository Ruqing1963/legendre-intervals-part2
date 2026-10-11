"""
hexad_search.py -- the sextic behind the mixed sums W' (Open Problem 9.1 of Part II).

For (l3, l5) in F_q^2 (q = 2^k, k odd) let L'(x) = l5^4 x^16 + l3^4 x^8 + x^4 + l3^2 x^2 + l5 x and K = ker L' (4-dim over F_2).
The data of Part II (Table 3) say that the Frobenius acts on K through Sp_4(F_2) = S_6, i.e. it preserves one of the
28 symplectic forms on K.  A symplectic form on F_2^4 is the same as a 'hexad up to translation': a 6-set
{r_1..r_6} in the affine space over K with r_i + r_j running over the 15 nonzero elements of K and sum(r_i) = 0
(the backslash-free docstring avoids escape warnings)
(duads = points, synthemes = totally isotropic lines).  Hexads containing 0 are {0} u frame, where a frame is a 5-set of
nonzero elements summing to 0 with any four independent (168 frames, 28 translation classes).

This script, for every (l3, l5) with K inside F_{q^12}:
  * computes K, the Frobenius permutation of K minus 0, the 28 classes, and the Frobenius-stable class(es);
  * builds the canonical sextic B(u) = prod (u + r_i) of the stable class, normalised by the translation that kills
    the u^4 coefficient (then B(u) = u^6 + b3 u^3 + b2 u^2 + b1 u + b0 with b_i in F_q, automatically rational);
  * records the factorisation type of B over F_q (= the S_6 class) and tests which trace form realises the symplectic form.
Usage: python hexad_search.py 3 [5]
"""
import itertools
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_d6_qf as QF  # noqa: E402


class Field:
    """F_{q^12} = F_q[z]/(m) with a few extra helpers."""

    def __init__(self, k):
        self.q = 1 << k
        self.k = k
        self.F = QF.Fq(self.q)
        self.K12 = QF.Fq12(self.F)
        self.A, self.M = self.F.A, self.F.M
        self.N = self.K12.N
        self.zero = [0] * QF.DEG
        self.one = [1] + [0] * (QF.DEG - 1)
        sq = self.K12.linmap_matrix(lambda a: self.K12.mul(a, a))
        self.SQ = sq
        self.SQRT = QF.f2_inverse(sq)
        # q-Frobenius matrix
        s = np.eye(self.N, dtype=np.uint8)
        for _ in range(k):
            s = s.dot(sq) % 2
        self.SIG = s
        tr = np.zeros((self.N, self.N), np.uint8); p = np.eye(self.N, dtype=np.uint8)
        for _ in range(QF.DEG):
            tr = (tr + p) % 2; p = p.dot(s) % 2
        self.TR12q = tr                      # Tr_{q^12/q}

    def add(self, a, b):
        return [int(self.A[x, y]) for x, y in zip(a, b)]

    def mul(self, a, b):
        return self.K12.mul(a, b)

    def scal(self, c, a):
        return [int(self.M[c, x]) for x in a]

    def apply(self, Mx, a):
        return self.K12.frombits(Mx.dot(self.K12.tobits(a)) % 2)

    def sqrt(self, a):
        return self.apply(self.SQRT, a)

    def frob(self, a):
        return self.apply(self.SIG, a)

    def pw(self, a, e):
        r = self.one
        b = a
        while e:
            if e & 1:
                r = self.mul(r, b)
            b = self.mul(b, b)
            e >>= 1
        return r

    def in_Fq(self, a):
        return all(x == 0 for x in a[1:])

    def abs_trace(self, a):
        """Tr_{q^12/2}(a) in F_2"""
        t = self.apply(self.TR12q, a)[0]
        return self.F.tr2(t)

    def Lprime_matrix(self, l3, l5):
        l34 = int(self.M[self.M[l3, l3], self.M[l3, l3]]); l54 = int(self.M[self.M[l5, l5], self.M[l5, l5]])
        l32 = int(self.M[l3, l3])

        def f(a):
            a2 = self.mul(a, a); a4 = self.mul(a2, a2); a8 = self.mul(a4, a4); a16 = self.mul(a8, a8)
            r = self.scal(l54, a16)
            r = self.add(r, self.scal(l34, a8))
            r = self.add(r, a4)
            r = self.add(r, self.scal(l32, a2))
            r = self.add(r, self.scal(l5, a))
            return r
        return self.K12.linmap_matrix(f)

    def poly_from_roots(self, roots):
        """coefficients (low degree first) of prod (u + r)"""
        p = [self.one]
        for r in roots:
            new = [self.zero] * (len(p) + 1)
            for i, c in enumerate(p):
                new[i + 1] = self.add(new[i + 1], c)
                new[i] = self.add(new[i], self.mul(c, r))
            p = new
        return p


def analyse(Fd: Field, l3, l5):
    Mx = Fd.Lprime_matrix(l3, l5)
    ns = QF.f2_nullspace(Mx)
    if len(ns) != 4:
        return dict(dim=len(ns))
    elems = []
    for coeffs in itertools.product((0, 1), repeat=4):
        v = np.zeros(Fd.N, np.uint8)
        for c, b in zip(coeffs, ns):
            if c:
                v ^= b
        elems.append(tuple(Fd.K12.frombits(v)))
    nonzero = [e for e in elems if any(e)]
    idx = {e: i for i, e in enumerate(nonzero)}
    bits = {e: Fd.K12.tobits(list(e)) for e in nonzero}
    # Frobenius permutation and orbit type
    perm = {e: tuple(Fd.frob(list(e))) for e in nonzero}
    seen, orbits = set(), []
    for e in nonzero:
        if e in seen:
            continue
        o, x = [], e
        while x not in seen:
            seen.add(x); o.append(x); x = perm[x]
        orbits.append(len(o))
    orbit_type = tuple(sorted(orbits))

    def add(e, f):
        return tuple(Fd.add(list(e), list(f)))
    zero = tuple(Fd.zero)
    # frames: 4-subsets that are bases, plus their sum
    frames = set()
    for quad in itertools.combinations(nonzero, 4):
        if QF.f2_rank([bits[e] for e in quad]) == 4:
            s = zero
            for e in quad:
                s = add(s, e)
            frames.add(frozenset(quad + (s,)))
    assert len(frames) == 168, len(frames)
    hexads = [frozenset(fr | {zero}) for fr in frames]

    def canon(H):
        return min(tuple(sorted(frozenset(add(h, kap) for h in H))) for kap in elems)
    classes = defaultdict(list)
    for H in hexads:
        classes[canon(H)].append(H)
    assert len(classes) == 28, len(classes)
    stable = []
    for key, Hs in classes.items():
        H = Hs[0]
        sH = frozenset(perm.get(h, zero) if h != zero else zero for h in H)
        if canon(sH) == key:
            stable.append(H)
    out = dict(dim=4, orbit_type=orbit_type, n_stable=len(stable), stable=[])
    # smallest field F_{q^N0} containing K and its absolute trace
    N0 = 1
    for o in orbits:
        N0 = N0 * o // np.gcd(N0, o)
    SIGpow = [np.eye(Fd.N, dtype=np.uint8)]
    for _ in range(N0 - 1):
        SIGpow.append(SIGpow[-1].dot(Fd.SIG) % 2)
    TRN0 = sum(SIGpow) % 2

    def trN0(a):
        return Fd.F.tr2(Fd.apply(TRN0, a)[0])
    out["N0"] = N0
    for H in stable:
        roots = [list(h) for h in H]
        P = Fd.poly_from_roots(roots)                   # degree 6, low first
        assert all(x == 0 for x in P[5]), "e1 != 0"
        c = Fd.sqrt(P[4])
        rts = [Fd.add(r, c) for r in roots]
        Pn = Fd.poly_from_roots(rts)
        assert all(x == 0 for x in Pn[4]) and all(x == 0 for x in Pn[5])
        rational = all(Fd.in_Fq(coef) for coef in Pn[:4])
        sextic = [coef[0] for coef in Pn[:4]] if rational else None      # b0, b1, b2, b3 in F_q
        # symplectic form of the class: omega(duad, duad') = 0 iff disjoint; compare with trace forms over F_{q^N0}
        rlist = sorted(H, key=lambda e: idx.get(e, -1))
        duad = {}
        for i, j in itertools.combinations(range(6), 2):
            duad[add(rlist[i], rlist[j])] = (i, j)
        matches = []
        for a_exp in range(0, 4):
            ok = True
            for x in nonzero:
                xa = list(x)
                for _ in range(a_exp):
                    xa = Fd.mul(xa, xa)
                for y in nonzero:
                    disjoint = len(set(duad[x]) & set(duad[y])) == 0
                    om = 0 if (disjoint or x == y) else 1
                    if trN0(Fd.mul(xa, list(y))) != om:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                matches.append(f"Tr_N0(x^{2**a_exp} y)")
        # factorisation type of the normalised sextic over F_q: orbits of Frobenius on its roots
        rt_t = [tuple(r) for r in rts]
        seen2, cyc = set(), []
        for r in rt_t:
            if r in seen2:
                continue
            o, x = 0, r
            while x not in seen2:
                seen2.add(x); o += 1; x = tuple(Fd.frob(list(x)))
            cyc.append(o)
        out["stable"].append(dict(sextic=sextic, rational=rational, sextic_type=tuple(sorted(cyc)), omega=matches))
    return out


def run(k, pairs=None):
    t0 = time.time()
    Fd = Field(k)
    q = Fd.q
    print(f"\n=== k={k}, q={q}, F_(q^12) with N={Fd.N} bits  (setup {time.time()-t0:.1f}s)")
    stats = Counter()
    rows = []
    pairs = pairs or [(l3, l5) for l3 in range(q) for l5 in range(1, q)]
    for (l3, l5) in pairs:
        r = analyse(Fd, l3, l5)
        r["l3"], r["l5"] = l3, l5
        rows.append(r)
        if r["dim"] != 4:
            stats[("dim", r["dim"])] += 1
            continue
        key = (r["orbit_type"], r["N0"], r["n_stable"],
               tuple(sorted(Counter((s["rational"], s["sextic_type"], tuple(s["omega"])) for s in r["stable"]).items())))
        stats[key] += 1
    print("  (orbit type on K\\0, N0 = [F_q(K):F_q], #stable classes, {(rational, sextic type, omega matches): count}) -> #pairs")
    for key, n in sorted(stats.items(), key=str):
        print(f"   {key}  x{n}")
    print(f"  [{time.time()-t0:.1f}s]")
    return Fd, rows


if __name__ == "__main__":
    import json
    ks = [int(a) for a in sys.argv[1:]] or [3]
    for k in ks:
        Fd, rows = run(k)
        path = os.path.join(os.path.dirname(HERE), "data", f"hexads_k{k}.json")
        with open(path, "w") as fh:
            json.dump(rows, fh, default=lambda o: list(o) if isinstance(o, (tuple, set, frozenset)) else str(o))
        print(f"  wrote {path}")

