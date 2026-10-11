"""
moment3_structure.py -- Part IV, Section 2.2: structure of the Gauss sums G_H(theta) behind the third moment.

For theta in Theta = {theta in F_{q^3}^* : Tr_{q^3/q}(1/theta) = 0} let Q_theta(e) = Tr_{q^3/2}(theta e^5 + theta e^3) on
F_{q^3}, H = ker Tr_{q^3/q}, and G_H(theta) = sum_{e in H} (-1)^{Q_theta(e)}.  The polar operator is
T_theta(e) = theta e^4 + theta e^2 + (theta e)^{1/2} + (theta e)^{1/4}; the radical of Q_theta on F_{q^3} is
R = ker T_theta ∩ F_{q^3} (the quintic frame of v^5 + v + c, c = 1/theta, over F_{q^3}), and the radical of the
restriction to H is R_H = {e in H : T_theta(e) in F_q}.
For each theta the script records: dim R, dim R_H, G_H(theta), whether Q_theta vanishes on R_H, the Arf invariant of the
induced form on H/R_H when G_H != 0, and the factorisation type of v^5 + v + c over F_{q^3} (by the number of roots in
F_{q^3}: 1 -> (1,4), 0 -> (2,3), 3 -> (1,1,1,2)), then tabulates the joint distribution.
Usage: python moment3_structure.py 3      (k = 5 needs F_{2^15} tables: about a minute)
"""
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_var_check as LV  # noqa: E402
import legendre_d6_qf as QF  # noqa: E402

LV.PRIMITIVE.update({3: 0b1011, 5: 0b100101, 7: 0b10000011, 9: 0b1000010001, 15: 0b1000000000000011,
                     21: 0b1000000000000000000101})


def run(k):
    q = 1 << k
    n = 3 * k
    gf = LV.GF(n)
    exp, log, order, size = gf.exp, gf.log, gf.order, gf.size

    def mul(a, b):
        return 0 if a == 0 or b == 0 else int(exp[(log[a] + log[b]) % order])

    def pw(a, e):
        if a == 0:
            return 0 if e > 0 else 1
        return int(exp[(e * log[a]) % order])

    def root(a, e):             # a^(1/e) for e a power of 2
        return 0 if a == 0 else int(exp[(log[a] * pow(e, -1, order)) % order])

    def frq(a):
        return pw(a, q)

    def tau(a):
        return a ^ frq(a) ^ frq(frq(a))

    def tr2(a):
        return int(LV._abs_trace_parity(gf, np.array([a], dtype=np.uint32))[0])

    def bits(a):
        return np.array([(a >> i) & 1 for i in range(n)], dtype=np.uint8)

    def frombits(v):
        return int(sum(int(b) << i for i, b in enumerate(v)))

    def linmat(f):              # matrix (columns = images of basis 2^i) of an F_2-linear map
        return np.array([bits(f(1 << i)) for i in range(n)], dtype=np.uint8).T

    Fq = {x for x in range(size) if frq(x) == x}
    Mtau = linmat(tau)
    Hbasis = QF.f2_nullspace(Mtau)                           # basis of H (dim 2k)
    Hel = []
    import itertools
    for cf in itertools.product((0, 1), repeat=len(Hbasis)):
        v = np.zeros(n, dtype=np.uint8)
        for c, b in zip(cf, Hbasis):
            if c:
                v ^= b
        Hel.append(frombits(v))
    Hel = np.array(Hel, dtype=np.int64)
    stats = Counter()
    for th in range(1, size):
        inv = int(exp[(-log[th]) % order])
        if tau(inv) != 0:
            continue
        c = inv                                                # Gold parameter

        def T(e):
            te = mul(th, e)
            return mul(th, pw(e, 4)) ^ mul(th, mul(e, e)) ^ root(te, 2) ^ root(te, 4)

        MT = linmat(T)
        R = QF.f2_nullspace(MT)
        dimR = len(R)
        # R_H = {e in H : T(e) in F_q}: solve on H-coordinates
        img = [T(frombits(b)) for b in Hbasis]
        # T(e) in F_q  <=>  T(e)^q = T(e)  <=>  (frq - id)(T(e)) = 0: linear in e
        Mfq = linmat(lambda a: frq(a) ^ a)
        cols = np.array([bits(x) for x in img], dtype=np.uint8).T           # n x 2k
        Mcond = (Mfq.astype(int).dot(cols.astype(int)) % 2).astype(np.uint8)
        RH = QF.f2_nullspace(Mcond)                                            # coordinates w.r.t. Hbasis
        dimRH = len(RH)
        # Q on H
        e5e3 = np.array([pw(int(e), 5) ^ pw(int(e), 3) for e in Hel], dtype=np.int64)
        vals = np.array([mul(th, int(w)) for w in e5e3], dtype=np.uint32)
        par = LV._abs_trace_parity(gf, vals).astype(np.int64)
        G = int(np.sum(1 - 2 * par))
        # Q on R_H
        qRH = []
        for v in RH:
            e = np.zeros(n, dtype=np.uint8)
            for cc, b in zip(v, Hbasis):
                if cc:
                    e ^= b
            ee = frombits(e)
            qRH.append(tr2(mul(th, pw(ee, 5) ^ pw(ee, 3))))
        vanish = all(x == 0 for x in qRH)
        # type of v^5+v+c over F_{q^3}: number of roots
        nroots = sum(1 for v in range(size) if (pw(v, 5) ^ v ^ c) == 0)
        typ = {1: "(1,4)", 0: "(2,3)", 3: "(1,1,1,2)"}.get(nroots, f"?{nroots}")
        arf = None
        if G != 0:
            arf = 0 if G > 0 else 1
        stats[(typ, "dim R", dimR, "dim R_H", dimRH, "Q|R_H=0", vanish, "G_H", G, "Arf", arf)] += 1
    print(f"\nk={k}, q={q}: theta in Theta (|Theta| = q^2-1), Gold parameter c = 1/theta with Tr_(q^3/q) c = 0")
    for key, m in sorted(stats.items(), key=str):
        print(f"   {key}  x{m}")


if __name__ == "__main__":
    for k in [int(a) for a in sys.argv[1:]] or [3]:
        run(k)
