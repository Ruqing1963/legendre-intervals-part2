"""
hexad_fit.py -- fit the coefficients of the canonical sextic B(u) = u^6 + b3 u^3 + b2 u^2 + b1 u + b0 found by
hexad_search.py as Laurent polynomials in (l3, l5) over F_q (Gaussian elimination over F_q with the fqlib tables).
Only pairs with a unique Frobenius-stable hexad class are used for the fit; the formula is then tested on all pairs
(for pairs with several stable classes it must single out exactly one).  Usage: python hexad_fit.py 5 [A B]
"""
import json
import os
import sys
from itertools import product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_d6_qf as QF  # noqa: E402


class FqLin:
    def __init__(self, q):
        self.q = q
        F = QF.Fq(q)
        self.A, self.M, self.INV = F.A, F.M, F.INV

    def pw(self, a, e):
        e %= (self.q - 1)
        r = 1
        for _ in range(e):
            r = int(self.M[r, a])
        return r

    def solve(self, rows, rhs):
        """solve rows . x = rhs over F_q (rows: list of int lists). Returns (solution or None, rank, consistent)."""
        A, M, INV = self.A, self.M, self.INV
        Mx = np.array([r + [b] for r, b in zip(rows, rhs)], dtype=np.int64)
        n, m1 = Mx.shape
        m = m1 - 1
        piv = []
        r = 0
        for c in range(m):
            if r >= n:
                break
            nz = np.nonzero(Mx[r:, c])[0]
            if len(nz) == 0:
                continue
            p = r + nz[0]
            Mx[[r, p]] = Mx[[p, r]]
            inv = int(INV[Mx[r, c]])
            Mx[r] = M[inv, Mx[r]]
            others = np.nonzero(Mx[:, c])[0]
            for o in others:
                if o != r:
                    f = int(Mx[o, c])
                    Mx[o] = A[Mx[o], M[f, Mx[r]]]
            piv.append(c)
            r += 1
        consistent = not np.any(Mx[r:, m])
        x = [0] * m
        if consistent:
            for i, c in enumerate(piv):
                x[c] = int(Mx[i, m])
        return x, r, consistent

    def eval_formula(self, coeffs, monos, l3, l5):
        s = 0
        for cf, (a, b) in zip(coeffs, monos):
            if cf:
                s = int(self.A[s, self.M[cf, self.M[self.pw(l3, a), self.pw(l5, b)]]])
        return s


def main():
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    Amax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    Bmax = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    q = 1 << k
    L = FqLin(q)
    rows = json.load(open(os.path.join(os.path.dirname(HERE), "data", f"hexads_k{k}.json")))
    uniq = [r for r in rows if r.get("dim") == 4 and r["n_stable"] == 1 and r["stable"][0]["sextic"]]
    multi = [r for r in rows if r.get("dim") == 4 and r["n_stable"] > 1]
    print(f"k={k}: {len(uniq)} pairs with a unique stable class, {len(multi)} with several")
    monos = [(a, b) for a in range(0, Amax + 1) for b in range(-Bmax, Bmax + 1)]
    design = [[L.M[L.pw(r['l3'], a), L.pw(r['l5'], b)] for (a, b) in monos] for r in uniq]
    design = [[int(v) for v in row] for row in design]
    formulas = {}
    for i in range(4):
        rhs = [r["stable"][0]["sextic"][i] for r in uniq]
        x, rank, ok = L.solve(design, rhs)
        terms = [(cf, mono) for cf, mono in zip(x, monos) if cf]
        print(f"  b{i}: consistent={ok}, rank {rank}/{len(monos)} unknowns, {len(uniq)} equations; "
              f"{len(terms)} nonzero terms" + ("" if ok else "  (NO Laurent polynomial of this degree)"))
        if ok:
            formulas[i] = x
            print("      " + " + ".join(f"{cf}*l3^{a}*l5^{b}" for cf, (a, b) in terms[:30]))
    if len(formulas) < 4:
        return
    # test on all pairs: which stable classes match the formula?
    from collections import Counter
    stat = Counter()
    for r in rows:
        if r.get("dim") != 4:
            continue
        pred = tuple(L.eval_formula(formulas[i], monos, r['l3'], r['l5']) for i in range(4))
        hits = [s for s in r["stable"] if s["sextic"] and tuple(s["sextic"]) == pred]
        stat[(tuple(r["orbit_type"]), r["n_stable"], len(hits), tuple(tuple(h["sextic_type"]) for h in hits))] += 1
    print("  formula applied to all pairs: (orbit type, #stable, #matching classes, types of the matches) -> count")
    for key, n in sorted(stat.items(), key=str):
        print(f"   {key}  x{n}")


if __name__ == "__main__":
    main()
