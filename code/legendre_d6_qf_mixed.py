"""
legendre_d6_qf_mixed.py -- fine structure of W(l3,l5) = sum_{l1} S(l1,l3,l5) for d = 6.

For every orbit representative (l3,l5) this prints
  * r_beta : radical dimension of the beta-form (W),  W/q^7
  * r_x    : dim_F2 ( ker L  cap  F_{q^12} ),  L(x) = l5^4 x^16 + l3^4 x^8 + l3^2 x^2 + l5 x   (x-form radical)
  * Qx0    : whether Q_x = Tr(l5 x^5 + l3 x^3) vanishes on that radical
  * the multiset of S(l1,l3,l5)/q^6 over l1 in F_q (how many l1 contribute, with which signs)
and checks, for the pure quintic terms, Stickelberger's evaluation
  S(0,0,l5) = eps0 q^6 (5 [l5 is a 5th power in F_{q^12}] - 1),  eps0 = (-1)^{k+1}.
For odd k it also extracts J+ = {j = l3^5/l5^3 : W != 0} and scans two-term trace predicates.
Usage: python legendre_d6_qf_mixed.py 1 3 5      (k = 7: pass --no-S to skip the per-l1 sums)
"""
from __future__ import annotations

import sys
import time
from collections import Counter

import numpy as np

import legendre_d6_qf as QF

DEG = 12


def setup(k):
    q = 1 << k
    F = QF.Fq(q)
    K = QF.Fq12(F)
    N = K.N
    A, M = F.A, F.M
    I = np.eye(N, dtype=np.uint8)
    SQ = K.linmap_matrix(lambda a: K.mul(a, a))
    SQRT = QF.f2_inverse(SQ)
    SQ2 = SQ.dot(SQ) % 2
    SQRT2 = SQRT.dot(SQRT) % 2
    sigma = I.copy()
    for _ in range(k):
        sigma = sigma.dot(SQ) % 2
    sig_pows = [I]
    for _ in range(DEG - 1):
        sig_pows.append(sig_pows[-1].dot(sigma) % 2)
    theta = (sigma + sig_pows[DEG - 1]) % 2
    TRmat = sum(sig_pows) % 2

    def scal(c):
        blk = np.zeros((k, k), np.uint8)
        for t in range(k):
            blk[:, t] = [(int(M[c, 1 << t]) >> s) & 1 for s in range(k)]
        return np.kron(np.eye(DEG, dtype=np.uint8), blk)

    basis = K.basis()
    G = np.zeros((N, N), np.uint8)
    for i in range(N):
        for j in range(i, N):
            p = K.mul(basis[i], basis[j])
            G[i, j] = G[j, i] = F.tr2(K.frombits(TRmat.dot(K.tobits(p)) % 2)[0])

    def trq_of(a):
        return K.frombits(TRmat.dot(K.tobits(a)) % 2)[0]

    def apply(Mx, a):
        return K.frombits(Mx.dot(K.tobits(a)) % 2)

    u3 = np.zeros(N, np.int64); u5 = np.zeros(N, np.int64)      # beta-form: Tr(e^2 theta e), Tr(e^4 theta e)
    t1 = np.zeros(N, np.int64); t3 = np.zeros(N, np.int64); t5 = np.zeros(N, np.int64)   # x-form: Tr(e), Tr(e^3), Tr(e^5)
    for i, e in enumerate(basis):
        te = apply(theta, e)
        e2 = K.mul(e, e); e3 = K.mul(e2, e); e4 = K.mul(e2, e2); e5 = K.mul(e4, e)
        u3[i] = trq_of(K.mul(e2, te)); u5[i] = trq_of(K.mul(e4, te))
        t1[i] = trq_of(e); t3[i] = trq_of(e3); t5[i] = trq_of(e5)
    cache = {}

    def S_(c):
        if c not in cache:
            cache[c] = scal(c)
        return cache[c]

    def inner(l3, l5):
        return (S_(l3).dot(SQ) + S_(F.SQRT[l3]).dot(SQRT) + S_(l5).dot(SQ2) + S_(F.SQRT[F.SQRT[l5]]).dot(SQRT2)) % 2

    return dict(q=q, k=k, F=F, K=K, N=N, A=A, M=M, G=G, theta=theta, inner=inner, u3=u3, u5=u5, t1=t1, t3=t3, t5=t5)


def fifth_power_in_q12(F, c):
    """c in F_q^* is a 5th power in F_{q^12} iff k % 4 != 0 or c is a 5th power in F_q."""
    q, k = F.q, F.k
    if k % 4:
        return True
    fifths = {int(QF.np.int64(0))}
    M = F.M
    s = set()
    for x in range(1, q):
        y = 1
        for _ in range(5):
            y = int(M[y, x])
        s.add(y)
    return c in s


def run(k, do_S=True):
    t0 = time.time()
    ctx = setup(k)
    q, F, N, A, M, G = ctx["q"], ctx["F"], ctx["N"], ctx["A"], ctx["M"], ctx["G"]
    TRq, SQ = F.TR, F.SQ
    eps0 = 1 if k % 2 else -1

    def pw(x, e):
        r = 1
        for _ in range(e):
            r = int(M[r, x])
        return r

    def full_orbit(l3, l5):
        orb = set()
        for s in range(1, q):
            s3, s5 = pw(s, 3), pw(s, 5)
            a, b = int(M[s3, l3]), int(M[s5, l5])
            for _ in range(k):
                orb.add((a, b)); a, b = int(SQ[a]), int(SQ[b])
        return orb

    print(f"\n### k={k}, q={q}, N={N}   (setup {time.time()-t0:.1f}s)")
    seen = set()
    rows = []
    for l3 in range(q):
        for l5 in range(q):
            if (l3, l5) in seen or (l3, l5) == (0, 0):
                continue
            orb = full_orbit(l3, l5); seen |= orb
            Tin = ctx["inner"](l3, l5)
            # beta-form
            Tb = Tin.dot(ctx["theta"]) % 2
            Bb = (Tb.T.dot(G)) % 2
            Qb = TRq[A[M[l3, ctx["u3"]], M[l5, ctx["u5"]]]].astype(np.uint8)
            W, rb, _ = QF.quad_exp_sum(Bb, Qb, np.zeros(N, np.uint8))
            # x-form
            Bx = (Tin.T.dot(G)) % 2
            rx = len(QF.f2_nullspace(Bx))
            Qx = TRq[A[M[l3, ctx["t3"]], M[l5, ctx["t5"]]]].astype(np.uint8)
            Svals = None
            if do_S:
                Svals = []
                for l1 in range(q):
                    Lv = TRq[M[l1, ctx["t1"]]].astype(np.uint8)
                    S, _, _ = QF.quad_exp_sum(Bx, Qx, Lv)
                    Svals.append(S)
                assert sum(Svals) == W, "sum over l1 of S must equal W"
            typ = "l5=0" if l5 == 0 else ("l3=0" if l3 == 0 else "mixed")
            rows.append((typ, l3, l5, len(orb), W, rb, rx, Svals))
    # report
    def fmt(v, p):
        return f"{v // q**p}" if v % q**p == 0 else f"{v / q**p:.3f}"
    for typ in ("l5=0", "l3=0", "mixed"):
        cnt = Counter()
        for (t, l3, l5, n, W, rb, rx, Sv) in rows:
            if t != typ:
                continue
            key = (fmt(W, 7), rb, rx)
            if Sv is not None:
                c = Counter(fmt(s, 6) for s in Sv)
                key = key + (tuple(sorted(c.items())),)
            cnt[key] += n
        print(f"  {typ}: (W/q^7, r_beta, r_x[, multiset of S(l1)/q^6 over l1]) -> #pairs")
        for key, n in sorted(cnt.items(), key=str):
            print(f"     {key}  x{n}")
    # Stickelberger check for pure quintic: S(0,0,l5)
    if do_S:
        ok = True
        for (t, l3, l5, n, W, rb, rx, Sv) in rows:
            if t == "l3=0":
                pred = eps0 * q**6 * (5 * int(fifth_power_in_q12(F, l5)) - 1)
                if Sv[0] != pred:
                    ok = False
        print(f"  Stickelberger S(0,0,l5) = eps0 q^6 (5[l5 fifth power in F_q^12] - 1): {'OK' if ok else 'FAILED'}")
    # J+ for odd k
    Wj = None
    if k % 2:
        Wj = {}
        for (t, l3, l5, n, W, rb, rx, Sv) in rows:
            if t == "mixed":
                for s in range(1, q):
                    a, b = int(M[pw(s, 3), l3]), int(M[pw(s, 5), l5])
                    for _ in range(k):
                        j = int(M[pw(a, 5), F.INV[pw(b, 3)]])
                        assert Wj.get(j, W) == W
                        Wj[j] = W
                        a, b = int(SQ[a]), int(SQ[b])
        Jp = {j for j, W in Wj.items() if W}
        print(f"  odd k: W(mixed) is a function of j; |J+| = {len(Jp)} = (q-2)/2 = {(q-2)//2}")
        # two-term scan: {j : Tr(j^a) + Tr(j^b) = c} and {j : Tr(j^a + j^b (j+1)^-1 ...)} limited to a,b < q-1
        hits = []
        Tj = {m: {j for j in range(1, q) if int(TRq[pw(j, m)]) == 1} for m in range(0, q - 1)}
        allj = set(range(1, q))
        for a in range(1, q - 1):
            for b in range(a + 1, q - 1):
                S1 = (Tj[a] ^ Tj[b])            # Tr(j^a)+Tr(j^b) = 1
                if S1 == Jp or (allj - S1) == Jp:
                    hits.append(("Tr(j^a)+Tr(j^b)", a, b))
        # rational: Tr(j^a / (j+1)^b)
        for a in range(0, 16):
            for b in range(1, 16):
                S1 = {j for j in range(1, q) if j != 1 and int(TRq[M[pw(j, a), F.INV[pw(A[j, 1], b)]]]) == 1}
                if S1 == Jp or (allj - {1} - S1) == Jp:
                    hits.append(("Tr(j^a/(j+1)^b)", a, b))
        print(f"  two-term / rational trace predicates matching J+: {hits[:10]}{' ...' if len(hits) > 10 else ''}")
    print(f"  [{time.time()-t0:.1f}s]")
    return rows, Wj


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    do_S = "--no-S" not in sys.argv
    for k in [int(a) for a in args] or [1, 3, 5]:
        run(k, do_S=do_S)
