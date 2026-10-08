"""
kasami_sign.py -- the mixed term W(l3,l5) for odd k through the L-polynomials of the supersingular genus-2 curves
    C_{t,l1} : y^2 + y = x^5 + t x^3 + l1 x,   t = j^{1/5},  l1 in F_q    (representative (l3,l5) = (t,1) of the class j).

For a genus-2 curve over F_q with Frobenius eigenvalues w_i,  S_m(l1) := sum_{x in F_{q^m}} psi(x^5 + t x^3 + l1 x)
= -sum_i w_i^m.  The sums over F_q and F_{q^2} (s1, s2) determine L(T) = prod(1 - w_i T) = 1 - e1 T + e2 T^2 - q e1 T^3
+ q^2 T^4 (e1 = -s1, e2 = (s1^2 + s2)/2) and hence S_12(l1); and W(t,1) = sum_{l1 in F_q} S_12(l1).  This is an
independent check of the quadratic-form computation (legendre_d6_qf_mixed.py) and of Theorem 6.9:
   * type (1,4) (4-cycle):  S_12(l1) = +4 q^6 for the q/2 values l1 with Tr_{q/2}(l1 u1) = kappa,  0 otherwise;
                            W = 2 q^7.  Here u1 is the F_q-root of A(u) = u^5 + t^2 u + 1 and
                            kappa = 1 + Tr_{q/2}(u1^5 + t u1^3);  the nonzero S_12 come from s1 = s2 = 0.
   * types (1,1,1,2), (2,3): every S_12(l1) = +-4 q^6, exactly half of each sign, W = 0.
For every j in F_q^* and every l1 in F_q we tabulate (type, kappa) -> (W/q^7, #S_12 = +4q^6, #S_12 = -4q^6, #S_12 = 0)
and check the l1-predicate of the 4-cycle case.
"""
import sys
from collections import Counter
from fractions import Fraction

import numpy as np

import legendre_var_check as LV


def run(k):
    q = 1 << k
    gf = LV.GF(2 * k)                 # F_{q^2}
    order = gf.order
    step = order // (q - 1)
    Fq = [0] + [int(gf.exp[m * step]) for m in range(q - 1)]
    m_all = np.arange(order, dtype=np.int64)
    x_all = gf.exp[m_all]
    x5 = gf.exp[(5 * m_all) % order]

    def inv(x):
        return gf.pow(x, order - 1)               # x^{-1}  (NOT order-2, which is x^{-2})

    def trq2(x):                                  # Tr_{q/2} on F_q subset F_{q^2}
        s, y = 0, x
        for _ in range(k):
            s ^= y
            y = gf.mul(y, y)
        return s

    # tables: absolute trace parity on F_{q^2}; Tr_{q/2} parity on F_q (indexed by element value)
    def tr_parity(arr):
        return LV._abs_trace_parity(gf, arr.astype(np.uint32))
    trq_tab = np.zeros(gf.size, dtype=np.int8)
    for x in Fq:
        trq_tab[x] = trq2(x)
    Fq_arr = np.array(Fq, dtype=np.uint32)
    Fq_log = gf.log[Fq_arr[1:]]
    Fq5 = np.concatenate(([0], gf.exp[(5 * Fq_log) % order])).astype(np.uint32)
    Fq3 = np.concatenate(([0], gf.exp[(3 * Fq_log) % order])).astype(np.uint32)

    def mul_vec(a, arr):                          # a * arr (a scalar, arr array of elements)
        if a == 0:
            return np.zeros_like(arr)
        return gf.mul_vec_by_logs(np.full(arr.shape, int(gf.log[a]), dtype=np.int64), arr)

    pw = gf.pow
    inv5 = pow(5, -1, q - 1)
    inv2 = pow(2, -1, q - 1)
    tab = Counter()
    pred_ok = True
    for j in Fq[1:]:
        t = pw(j, inv5)                          # j^{1/5} in F_q
        t2 = gf.mul(t, t)
        # roots of A(u) = u^5 + t^2 u + 1 over F_q and F_{q^2}
        r1 = [u for u in Fq[1:] if pw(u, 5) ^ gf.mul(t2, u) ^ 1 == 0]
        vals = x5 ^ mul_vec(t2, x_all) ^ np.uint32(1)
        r2 = int(np.count_nonzero(vals == 0))
        typ = {(1, 1): "(1,4)", (3, 5): "(1,1,1,2)", (0, 2): "(2,3)", (1, 5): "(1,2,2)", (0, 0): "(5)",
               (2, 2): "(1,1,3)"}.get((len(r1), r2), f"?{len(r1)},{r2}")
        kappa = None
        if typ == "(1,4)":
            u1 = r1[0]
            kappa = 1 ^ trq2(pw(u1, 5) ^ gf.mul(t, pw(u1, 3)))
        f1_base = Fq5 ^ mul_vec(t, Fq3)                      # x^5 + t x^3 on F_q
        f2_base = x5 ^ mul_vec(t, gf.exp[(3 * m_all) % order])   # on F_{q^2}^*
        npos = nneg = nzero = 0
        W = 0
        for lam1 in Fq:
            f1 = f1_base ^ mul_vec(lam1, Fq_arr)
            s1 = int(q - 2 * trq_tab[f1].sum())
            f2 = f2_base ^ mul_vec(lam1, x_all)
            s2 = int(order - 2 * tr_parity(f2).sum()) + 1   # include x = 0
            p1, p2 = -s1, -s2
            assert (p1 * p1 - p2) % 2 == 0
            e = [p1, (p1 * p1 - p2) // 2, q * p1, q * q]
            p = [4, p1, p2]
            for m in range(3, 13):
                s = sum((-1) ** (i - 1) * e[i - 1] * p[m - i] for i in range(1, 5) if m - i >= 1)
                if m <= 4:
                    s += (-1) ** (m - 1) * m * e[m - 1]
                p.append(s)
            S12 = -p[12]
            W += S12
            if S12 == 4 * q ** 6:
                npos += 1
            elif S12 == -4 * q ** 6:
                nneg += 1
            elif S12 == 0:
                nzero += 1
            else:
                raise AssertionError(f"unexpected S12={S12}")
            if typ == "(1,4)":
                predicted = trq2(gf.mul(lam1, u1)) == kappa
                if predicted != (S12 != 0) or (predicted and (s1, s2) != (0, 0)):
                    pred_ok = False
        tab[(typ, kappa, Fraction(W, q ** 7), npos, nneg, nzero)] += 1
    print(f"\nk={k}, q={q}:  (type, kappa, W/q^7, #S12=+4q^6, #S12=-4q^6, #S12=0) -> #j")
    for key, n in sorted(tab.items(), key=str):
        print(f"   {key}  x{n}")
    print(f"   4-cycle predicate 'S12 != 0 iff Tr(l1 u1) = kappa, with s1 = s2 = 0': {'OK' if pred_ok else 'FAILS'}")
    return tab, pred_ok


if __name__ == "__main__":
    for k in [int(a) for a in sys.argv[1:]] or [3, 5, 7]:
        run(k)
