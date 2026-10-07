"""
legendre_d6_qf.py -- exact evaluation of A6 and X (d = 6) through F_2-quadratic forms on F_{q^12}.

Hilbert 90: alpha = beta + beta^q maps F_{q^12} q-to-1 onto {Tr alpha = 0}.  With theta = sigma + sigma^{-1}
(sigma = q-Frobenius):  Tr alpha^3 = Tr(beta^2 theta beta),  Tr alpha^5 = Tr(beta^4 theta beta),
e2(alpha) = Tr(beta^{q+1}) + (Tr beta)^2.   Hence (chi = (-1)^{Tr_{q/2}}, Tr = Tr_{q^12/q}, tr = Tr_{q^12/2})

   A6 = q^-3 sum_{l3,l5 in F_q} W(l3,l5),          W(l3,l5) = sum_beta (-1)^{tr(l3 beta^2 theta beta + l5 beta^4 theta beta)}
   X  = q^-3 sum_{l3,l5 in F_q} W'(l3,l5),         W'(l3,l5)= sum_beta (-1)^{tr(l3 beta^2 theta beta + l5 beta^4 theta beta + beta^{q+1} + beta)}

Each W, W' is the exponential sum of an F_2-quadratic form Q (plus a linear form) on F_2^{12k}:
   sum_x (-1)^{Q(x)+L(x)} = 2^{dim R} [Q+L = 0 on R] (-1)^{Arf(Q|V') + Q(x0)} 2^{dim V'/2},
R = radical of the polar form B, V' a complement, B(x0, .) = L on V'.  Everything is F_2-linear algebra
on N = 12k bits, so q = 64, 128 are reachable (far beyond the Hayes FFT at d = 6).

The polar forms are B(beta,gamma) = tr(gamma . T(beta)) with
   T(beta) = l3 (theta beta)^2 + sqrt(l3 theta beta) + l5 (theta beta)^4 + (l5 theta beta)^{1/4}   (+ theta beta for W'),
so the radical is theta^{-1}(ker M), M(y) = l3 y^2 + sqrt(l3 y) + l5 y^4 + (l5 y)^{1/4} (+ y).

Field model: F_{q^12} = F_q[z]/(m(z)), m irreducible of degree 12 found by Rabin's test; F_q from fqlib tables.
Usage: python legendre_d6_qf.py 1 2 3 4 5 6
"""
from __future__ import annotations

import os
import random
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "legendre-intervals-explicit-formula", "code"))
import fqlib  # noqa: E402

fqlib.MODULI.update({32: (2, 5, [1, 0, 1, 0, 0, 1]), 64: (2, 6, [1, 1, 0, 0, 0, 0, 1]),
                     128: (2, 7, [1, 1, 0, 0, 0, 0, 0, 1]), 256: (2, 8, [1, 0, 1, 1, 1, 0, 0, 0, 1])})

DEG = 12


class Fq:
    def __init__(self, q):
        self.q = q
        self.k = q.bit_length() - 1
        A, M, NEG, INV = fqlib.field_tables(q)
        self.A, self.M, self.INV = A, M, INV
        x = np.arange(q)
        self.SQ = M[x, x]
        self.SQRT = np.empty(q, np.int64); self.SQRT[self.SQ] = x
        TR = np.zeros(q, np.int64); y = x.copy()
        for _ in range(self.k):
            TR = A[TR, y]; y = M[y, y]
        self.TR = TR

    def tr2(self, c):
        return int(self.TR[c])


class Fq12:
    """F_q[z]/(m), elements = lists of DEG ints."""

    def __init__(self, F: Fq, seed=0):
        self.F = F
        q, A, M = F.q, F.A, F.M
        rng = random.Random(seed)
        while True:
            m = [rng.randrange(q) for _ in range(DEG)] + [1]
            if self._irreducible(m):
                break
        self.m = m
        self.N = DEG * F.k

    # polynomial helpers over F_q (low degree first)
    def _pmulmod(self, a, b, m):
        A, M = self.F.A, self.F.M
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                Mx = M[x]
                for j, y in enumerate(b):
                    if y:
                        r[i + j] = A[r[i + j], Mx[y]]
        return self._pmod(r, m)

    def _pmod(self, a, m):
        A, M = self.F.A, self.F.M
        a = a[:]
        dm = len(m) - 1
        for top in range(len(a) - 1, dm - 1, -1):
            c = a[top]
            if c:
                for i in range(dm + 1):
                    a[top - dm + i] = A[a[top - dm + i], M[c, m[i]]]
        a = a[:dm]
        while len(a) < dm:
            a.append(0)
        return a

    def _ppow_q(self, a, m, e):
        """a^(q^e) mod m by repeated q-th powering (q = 2^k: k squarings each)."""
        for _ in range(e * self.F.k):
            a = self._pmulmod(a, a, m)
        return a

    def _pgcd_deg(self, a, b):
        def trim(p):
            p = p[:]
            while p and p[-1] == 0:
                p.pop()
            return p
        a, b = trim(a), trim(b)
        while b:
            # a mod b
            a = trim(self._pmod_general(a, b))
            a, b = b, a
        return len(a) - 1

    def _pmod_general(self, a, m):
        A, M, INV = self.F.A, self.F.M, self.F.INV
        a = a[:]
        dm = len(m) - 1
        li = INV[m[-1]]
        for top in range(len(a) - 1, dm - 1, -1):
            c = M[a[top], li]
            if c:
                for i in range(dm + 1):
                    a[top - dm + i] = A[a[top - dm + i], M[c, m[i]]]
        return a[:dm] if dm > 0 else []

    def _irreducible(self, m):
        z = [0, 1] + [0] * (DEG - 2)
        zq = z[:]
        powers = {}
        for e in range(1, DEG + 1):
            zq = self._ppow_q(zq, m, 1)
            powers[e] = zq
        if powers[DEG] != z:
            return False
        for r in (2, 3):
            diff = [int(self.F.A[x, y]) for x, y in zip(powers[DEG // r], z)]
            if self._pgcd_deg(m, diff) > 0:
                return False
        return True

    # field element arithmetic
    def mul(self, a, b):
        return self._pmulmod(a, b, self.m)

    def tobits(self, a):
        k = self.F.k
        v = np.zeros(self.N, np.uint8)
        for i, c in enumerate(a):
            for t in range(k):
                v[i * k + t] = (c >> t) & 1
        return v

    def frombits(self, v):
        k = self.F.k
        return [int(sum(int(v[i * k + t]) << t for t in range(k))) for i in range(DEG)]

    def basis(self):
        k = self.F.k
        out = []
        for i in range(DEG):
            for t in range(k):
                e = [0] * DEG
                e[i] = 1 << t
                out.append(e)
        return out

    def linmap_matrix(self, f):
        """N x N F_2 matrix (columns = images of basis vectors) of an additive map f on elements."""
        B = self.basis()
        Mx = np.zeros((self.N, self.N), np.uint8)
        for j, e in enumerate(B):
            Mx[:, j] = self.tobits(f(e))
        return Mx


# ----------------------------------------------------------------------------- F_2 linear algebra
def f2_inverse(Mx):
    n = Mx.shape[0]
    aug = np.concatenate([Mx.copy() % 2, np.eye(n, dtype=np.uint8)], axis=1)
    r = 0
    for c in range(n):
        piv = np.where(aug[r:, c] == 1)[0]
        if len(piv) == 0:
            raise ValueError("singular")
        p = r + piv[0]
        aug[[r, p]] = aug[[p, r]]
        rows = np.where(aug[:, c] == 1)[0]
        for rr in rows:
            if rr != r:
                aug[rr] ^= aug[r]
        r += 1
    return aug[:, n:]


def f2_nullspace(Mx):
    """basis of {v : Mx v = 0} over F_2, as rows."""
    n, m = Mx.shape
    Amat = Mx.copy() % 2
    pivcols = []
    r = 0
    for c in range(m):
        if r >= n:
            break
        piv = np.where(Amat[r:, c] == 1)[0]
        if len(piv) == 0:
            continue
        p = r + piv[0]
        Amat[[r, p]] = Amat[[p, r]]
        rows = np.where(Amat[:, c] == 1)[0]
        for rr in rows:
            if rr != r:
                Amat[rr] ^= Amat[r]
        pivcols.append(c)
        r += 1
    free = [c for c in range(m) if c not in pivcols]
    basis = []
    for fcol in free:
        v = np.zeros(m, np.uint8)
        v[fcol] = 1
        for i, pc in enumerate(pivcols):
            if Amat[i, fcol]:
                v[pc] = 1
        basis.append(v)
    return basis


def f2_rank(rows):
    if not rows:
        return 0
    Amat = np.array(rows, dtype=np.uint8) % 2
    n, m = Amat.shape
    r = 0
    for c in range(m):
        if r >= n:
            break
        piv = np.where(Amat[r:, c] == 1)[0]
        if len(piv) == 0:
            continue
        p = r + piv[0]
        Amat[[r, p]] = Amat[[p, r]]
        for rr in range(n):
            if rr != r and Amat[rr, c]:
                Amat[rr] ^= Amat[r]
        r += 1
    return r


def quad_exp_sum(B, Qv, Lv):
    """sum_{x in F_2^N} (-1)^{Q(x) + L(x)} for Q with polar matrix B (symmetric, zero diagonal),
    Q(e_i) = Qv[i], L(e_i) = Lv[i].  Returns (value, dim radical, info)."""
    N = B.shape[0]
    U = np.triu(B, 1)

    def Q(v):
        return (int(np.dot(v, Qv)) + int(np.dot(v, U.dot(v) % 2))) % 2

    R = f2_nullspace(B)
    for r in R:
        if (Q(r) + int(np.dot(r, Lv))) % 2:
            return 0, len(R), "Q+L nonzero on radical"
    # complement V'
    comp = []
    cur = list(R)
    rank = len(R)
    for i in range(N):
        e = np.zeros(N, np.uint8); e[i] = 1
        if f2_rank(cur + [e]) > rank:
            cur.append(e); comp.append(e); rank += 1
    g2 = len(comp)
    if g2 == 0:
        return 2 ** len(R), len(R), "radical only"
    C = np.array(comp, dtype=np.uint8)            # rows = complement basis vectors
    Bc = (C.dot(B).dot(C.T)) % 2                  # Gram matrix on V'
    Lc = C.dot(Lv) % 2
    # x0 in V' with B(x0, v) = L(v): Bc^T coords = Lc  (Bc symmetric)
    coords = (f2_inverse(Bc).dot(Lc)) % 2
    x0 = (coords.dot(C)) % 2
    qx0 = Q(x0.astype(np.uint8))
    # Arf invariant of Q on V' via symplectic reduction
    vecs = [C[i].copy() for i in range(g2)]
    arf = 0
    while vecs:
        e = vecs.pop(0)
        # find f with B(e,f)=1
        idx = None
        for j, f in enumerate(vecs):
            if int(e.dot(B).dot(f)) % 2 == 1:
                idx = j; break
        if idx is None:
            raise RuntimeError("degenerate complement")
        f = vecs.pop(idx)
        arf ^= (Q(e) & Q(f))
        new = []
        for u in vecs:
            be = int(u.dot(B).dot(e)) % 2
            bf = int(u.dot(B).dot(f)) % 2
            w = (u + bf * e + be * f) % 2
            new.append(w.astype(np.uint8))
        vecs = new
    sign = -1 if (arf + qx0) % 2 else 1
    return sign * 2 ** (len(R) + g2 // 2), len(R), f"arf={arf}, Q(x0)={qx0}"


# ----------------------------------------------------------------------------- main computation
def run(k, verbose=True):
    t0 = time.time()
    q = 1 << k
    F = Fq(q)
    K = Fq12(F)
    N = K.N
    A, M = F.A, F.M
    I = np.eye(N, dtype=np.uint8)
    SQ = K.linmap_matrix(lambda a: K.mul(a, a))
    SQRT = f2_inverse(SQ)
    SQ2 = SQ.dot(SQ) % 2
    SQRT2 = SQRT.dot(SQRT) % 2
    # sigma = q-Frobenius = SQ^k ; theta = sigma + sigma^11
    sigma = I.copy()
    for _ in range(k):
        sigma = sigma.dot(SQ) % 2
    sig_pows = [I]
    for _ in range(DEG - 1):
        sig_pows.append(sig_pows[-1].dot(sigma) % 2)
    assert np.array_equal(sig_pows[-1].dot(sigma) % 2, I), "sigma^12 != 1"
    theta = (sigma + sig_pows[DEG - 1]) % 2
    TRmat = sum(sig_pows) % 2                      # Tr_{q^12/q} as F_2-matrix (image lies in F_q = first block)

    def scal(c):                                   # scalar multiplication by c in F_q, as N x N F_2 matrix
        blk = np.zeros((k, k), np.uint8)
        for t in range(k):
            blk[:, t] = [(int(M[c, 1 << t]) >> s) & 1 for s in range(k)]
        return np.kron(np.eye(DEG, dtype=np.uint8), blk)

    basis = K.basis()
    bits = [K.tobits(e) for e in basis]
    # trace pairing G[i,j] = tr(e_i e_j)
    G = np.zeros((N, N), np.uint8)
    prods = {}
    for i in range(N):
        for j in range(i, N):
            p = K.mul(basis[i], basis[j])
            trq = K.frombits(TRmat.dot(K.tobits(p)) % 2)[0]     # Tr_{q^12/q}(p) in F_q (coefficient of z^0)
            G[i, j] = G[j, i] = F.tr2(trq)
    # Tr_{q^12/q} of e^2 theta(e), e^4 theta(e), e^{q+1}, e  (F_q-valued vectors)
    def trq_of(a):
        return K.frombits(TRmat.dot(K.tobits(a)) % 2)[0]

    def apply(Mx, a):
        return K.frombits(Mx.dot(K.tobits(a)) % 2)

    u3 = np.zeros(N, np.int64); u5 = np.zeros(N, np.int64); u1 = np.zeros(N, np.int64); u0 = np.zeros(N, np.int64)
    for i, e in enumerate(basis):
        te = apply(theta, e)
        e2 = K.mul(e, e); e4 = K.mul(e2, e2)
        u3[i] = trq_of(K.mul(e2, te))
        u5[i] = trq_of(K.mul(e4, te))
        u1[i] = trq_of(K.mul(e, apply(sigma, e)))
        u0[i] = trq_of(e)
    TRq = F.TR
    SQRTq = F.SQRT
    cache_scal = {}

    def S(c):
        if c not in cache_scal:
            cache_scal[c] = scal(c)
        return cache_scal[c]

    def T_matrix(l3, l5, with_theta_term):
        inner = (S(l3).dot(SQ) + S(SQRTq[l3]).dot(SQRT) + S(l5).dot(SQ2) + S(SQRTq[SQRTq[l5]]).dot(SQRT2)) % 2
        if with_theta_term:
            inner = (inner + I) % 2
        return inner.dot(theta) % 2

    # orbits under Frobenius (l3,l5) -> (l3^2, l5^2) [valid for both W and W'] and, for W only, scaling (s^3 l3, s^5 l5)
    def frob_orbit(l3, l5):
        orb = set()
        a, b = l3, l5
        while (a, b) not in orb:
            orb.add((a, b)); a, b = int(F.SQ[a]), int(F.SQ[b])
        return orb

    def full_orbit(l3, l5):
        orb = set()
        for s in range(1, q):
            s3 = M[M[s, s], s]; s5 = M[M[s3, s], s]
            orb |= frob_orbit(int(M[s3, l3]), int(M[s5, l5]))
        return orb

    results = {"W": {}, "Wp": {}}
    seen = set()
    for l3 in range(q):
        for l5 in range(q):
            if (l3, l5) in seen:
                continue
            orb = full_orbit(l3, l5)
            seen |= orb
            T = T_matrix(l3, l5, False)
            B = (T.T.dot(G)) % 2
            assert np.array_equal(B, B.T)
            Qv = TRq[A[M[l3, u3], M[l5, u5]]].astype(np.uint8)
            val, rdim, info = quad_exp_sum(B, Qv, np.zeros(N, np.uint8))
            results["W"][(l3, l5)] = (val, rdim, info, len(orb))
    seen = set()
    for l3 in range(q):
        for l5 in range(q):
            if (l3, l5) in seen:
                continue
            orb = frob_orbit(l3, l5)
            seen |= orb
            T = T_matrix(l3, l5, True)
            B = (T.T.dot(G)) % 2
            Qv = TRq[A[A[M[l3, u3], M[l5, u5]], u1]].astype(np.uint8)
            Lv = TRq[u0].astype(np.uint8)
            val, rdim, info = quad_exp_sum(B, Qv, Lv)
            results["Wp"][(l3, l5)] = (val, rdim, info, len(orb))
    A6 = sum(v * n for (v, _, _, n) in results["W"].values()) // q ** 3
    X = sum(v * n for (v, _, _, n) in results["Wp"].values()) // q ** 3
    eps = 1 if k % 2 else -1
    print(f"\n=== k={k}, q={q}, N={N} bits  [{time.time()-t0:.1f}s]")
    print(f"  A6 = {A6};  q^9 + (-1)^(k+1)(q-1)q^5 = {q**9 + eps*(q-1)*q**5}: {A6 == q**9 + eps*(q-1)*q**5}")
    print(f"  X  = {X};  (-1)^(k+1) q^6 + q^5 = {eps*q**6 + q**5}: {X == eps*q**6 + q**5}")
    if verbose:
        # distribution of W values (Z = W/q) by type
        from collections import Counter
        cnt = Counter()
        for (l3, l5), (v, rdim, info, n) in results["W"].items():
            typ = "l3=l5=0" if (l3, l5) == (0, 0) else ("l5=0" if l5 == 0 else ("l3=0" if l3 == 0 else "mixed"))
            cnt[(typ, v // q ** 6 if v % q ** 6 == 0 else v / q ** 6, rdim)] += n
        print("  W(l3,l5)/q^6 (dim radical) -> multiplicity:")
        for key in sorted(cnt, key=lambda t: (t[0], -abs(t[1]) if isinstance(t[1], int) else 0)):
            print(f"     {key[0]:8s}  W/q^6 = {key[1]:>6}  r = {key[2]:2d}   x{cnt[key]}")
        cnt2 = Counter()
        for (l3, l5), (v, rdim, info, n) in results["Wp"].items():
            typ = "l3=l5=0" if (l3, l5) == (0, 0) else ("l5=0" if l5 == 0 else ("l3=0" if l3 == 0 else "mixed"))
            cnt2[(typ, v // q ** 6 if v % q ** 6 == 0 else v / q ** 6, rdim)] += n
        print("  W'(l3,l5)/q^6 (dim radical) -> multiplicity:")
        for key in sorted(cnt2, key=lambda t: (t[0], -abs(t[1]) if isinstance(t[1], int) else 0)):
            print(f"     {key[0]:8s}  W'/q^6 = {key[1]:>6}  r = {key[2]:2d}   x{cnt2[key]}")
    return A6, X, results


if __name__ == "__main__":
    ks = [int(a) for a in sys.argv[1:]] or [1, 2, 3, 4, 5]
    for k in ks:
        run(k)
