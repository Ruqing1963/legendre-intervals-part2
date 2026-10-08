"""
verify_Y_direct.py -- direct point-set evaluation of
    Y_k = sum_{alpha in U} (-1)^{Tr_{q/2} e_4(alpha)},   U = {alpha in F_{q^12}: Tr alpha = Tr alpha^3 = Tr alpha^5 = 0},
by two independent routes, with assertions against the paper's values Y_1 = 104, Y_2 = -64.

Route A (q = 2): enumerate every alpha in F_{2^12}, form the 12 conjugates z_j = alpha^(2^j), compute
          e_1, e_3, e_4, e_5 of (z_0..z_11) by expanding prod (t - z_j), test e_1 = e_3 = e_5 = 0.
          (Equivalent to the twisted F_2-points of X = V(e_1, e_3, e_5) in A^12.)
Route B (q = 2, 4): enumerate monic F = t^12 + c_2 t^10 + c_4 t^8 + e_6 t^6 + ... (e_1 = e_3 = e_5 = 0), weight
          Lambda(F) = #{alpha in F_{q^12}: charpoly(alpha) = F}, computed with the vectorised recurrence of
          legendre_var_check.py (all alpha at once), and sum Lambda(F) chi(c_4).  Route B is the pushforward of
          Route A along the finite map X -> A^9.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import legendre_var_check as LV

EXPECTED = {2: 104, 4: -64, 8: -17920}


def route_a(q):
    k = q.bit_length() - 1
    gf = LV.GF(12 * k)
    mul, pw = gf.mul, gf.pow

    def tr_q_2(x):            # Tr_{q/2} of an element of F_q subset GF(2^{12k})
        s, y = 0, x
        for _ in range(k):
            s ^= y
            y = mul(y, y)
        assert s in (0, 1)
        return s

    Y, size = 0, 0
    for alpha in range(gf.size):
        z = [pw(alpha, q ** j) for j in range(12)] if alpha else [0] * 12
        poly = [1]
        for r in z:
            new = [0] * (len(poly) + 1)
            for i, c in enumerate(poly):
                new[i] ^= c
                new[i + 1] ^= mul(c, r)
            poly = new
        e = poly[1:]               # e_1 .. e_12
        if e[0] == 0 and e[2] == 0 and e[4] == 0:
            size += 1
            Y += 1 - 2 * tr_q_2(e[3])
    return Y, size


def route_b(q):
    res = LV.analyse(q, 6, verbose=False)
    psi, hay = res["psi"], res["hayes"]
    Y = 0
    for c2 in range(q):
        for c4 in range(q):
            g = (0, c2, 0, c4, 0)
            x = hay.fq[c4]
            Y += int(psi[hay.key(g)]) * (1 - 2 * hay.trace2(x))
    return Y, int(sum(int(psi[hay.key((0, c2, 0, c4, 0))]) for c2 in range(q) for c4 in range(q)))


if __name__ == "__main__":
    for q in (2, 4):
        t0 = time.time()
        if q == 2:
            Ya, Ua = route_a(q)
            print(f"q={q}: route A (points of X^sigma(F_q)): |U| = {Ua}, Y = {Ya}   [{time.time()-t0:.1f}s]")
            assert Ya == EXPECTED[q], (Ya, EXPECTED[q])
        t0 = time.time()
        Yb, Ub = route_b(q)
        print(f"q={q}: route B (pushforward to A^9, Lambda-weighted): |U| = {Ub}, Y = {Yb}   [{time.time()-t0:.1f}s]")
        assert Yb == EXPECTED[q], (Yb, EXPECTED[q])
        print(f"   m_{q.bit_length()-1} = Y / q^3 = {Yb // q**3}")
    print("all assertions passed")
