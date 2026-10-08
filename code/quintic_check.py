"""
quintic_check.py -- the mixed term for odd k through the Gold-type quintic  v^5 + v = c,  c = j^{-1/2}.

Facts tested here (k odd):
  (1) the kernel equation  u^15 + t^4 u^7 + t^2 u + 1 = 0  (t = j^{1/5}) factors as
        (u^5 + t^2 u + 1)(u^10 + t^2 u^6 + u^5 + 1),   and  u = t^{1/2} v  turns the quintic into  v^5 + v = j^{-1/2};
  (2) v^5 + v + c is never irreducible over F_q (so the whole kernel is F_{q^12}-rational);
  (3) W(l3,l5) != 0  <=>  v^5 + v + c has factorisation type (1,4) over F_q  (Frobenius = 4-cycle);
  (4) #{c in F_q^* : type (1,4)} = (q-2)/2.
Usage: python quintic_check.py          (types for k = 3,5,...,13; correspondence with W for k = 3,5,7)
"""
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "legendre-intervals-explicit-formula", "code"))
import fqlib  # noqa: E402

fqlib.MODULI.update({32: (2, 5, [1, 0, 1, 0, 0, 1]), 64: (2, 6, [1, 1, 0, 0, 0, 0, 1]),
                     128: (2, 7, [1, 1, 0, 0, 0, 0, 0, 1]), 256: (2, 8, [1, 0, 1, 1, 1, 0, 0, 0, 1]),
                     512: (2, 9, [1, 0, 0, 0, 1, 0, 0, 0, 0, 1]), 1024: (2, 10, [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]),
                     2048: (2, 11, [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]),
                     4096: (2, 12, [1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1]),
                     8192: (2, 13, [1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1])})


class Poly:
    """polynomials over F_q as lists (low degree first) using fqlib tables"""

    def __init__(self, q):
        self.q = q
        self.A, self.M, self.NEG, self.INV = fqlib.field_tables(q)

    def trim(self, a):
        a = list(a)
        while a and a[-1] == 0:
            a.pop()
        return a

    def mul(self, a, b):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y:
                        r[i + j] = int(self.A[r[i + j], self.M[x, y]])
        return self.trim(r)

    def mod(self, a, m):
        a = list(a)
        dm = len(m) - 1
        li = int(self.INV[m[-1]])
        for top in range(len(a) - 1, dm - 1, -1):
            c = int(self.M[a[top], li])
            if c:
                for i in range(dm + 1):
                    a[top - dm + i] = int(self.A[a[top - dm + i], self.M[c, m[i]]])
        return self.trim(a[:dm])

    def powmod(self, a, e, m):
        res, base = [1], self.mod(a, m)
        while e:
            if e & 1:
                res = self.mod(self.mul(res, base), m)
            base = self.mod(self.mul(base, base), m)
            e >>= 1
        return res

    def gcd(self, a, b):
        a, b = self.trim(a), self.trim(b)
        while b:
            a, b = b, self.mod(a, b)
        return a

    def divexact_linear(self, a, r):
        """a / (x + r) by synthetic division (assumes a(r) = 0)"""
        n = len(a) - 1
        out = [0] * n
        carry = a[n]
        for i in range(n - 1, -1, -1):
            out[i] = carry
            carry = int(self.A[a[i], self.M[carry, r]])
        assert carry == 0
        return out


_ROOTS = {}


def roots_table(P):
    """roots[c] = list of v in F_q with v^5 + v = c (O(q) histogram instead of O(q^2) search)"""
    q = P.q
    if q in _ROOTS:
        return _ROOTS[q]
    tab = [[] for _ in range(q)]
    for v in range(q):
        v2 = int(P.M[v, v]); v4 = int(P.M[v2, v2]); v5 = int(P.M[v4, v])
        tab[int(P.A[v5, v])].append(v)
    _ROOTS[q] = tab
    return tab


def factor_type(P, c):
    """factorisation type of v^5 + v + c over F_q as a sorted tuple of factor degrees"""
    q = P.q
    f = [c, 1, 0, 0, 0, 1]
    roots = roots_table(P)[c]
    g = f[:]
    for r in roots:
        g = P.divexact_linear(g, r)
    degs = [1] * len(roots)
    d = len(g) - 1
    if d == 0:
        return tuple(sorted(degs))
    # g has no F_q-roots; detect quadratic factors via gcd(g, x^{q^2} - x)
    x = [0, 1]
    xq = P.powmod(x, q, g)
    xq2 = P.powmod(xq, q, g)
    diff = P.trim([int(P.A[a, b]) for a, b in zip(xq2 + [0] * (len(x) - len(xq2)), x + [0] * (len(xq2) - len(x)))])
    h = P.gcd(g, diff)
    dq = len(h) - 1 if h else 0
    nquad = dq // 2
    degs += [2] * nquad
    rest = d - 2 * nquad
    if rest:
        degs.append(rest)        # irreducible of degree 3, 4 or 5
    return tuple(sorted(degs))


def main():
    ks = [3, 5, 7, 9, 11]
    for k in ks:
        q = 1 << k
        P = Poly(q)
        types = Counter()
        for c in range(1, q):
            types[factor_type(P, c)] += 1
        n14 = types[(1, 4)]
        print(f"k={k:2d} q={q:5d}: types of v^5+v+c, c in F_q^*: {dict(sorted(types.items()))};  "
              f"#(1,4) = {n14} vs (q-2)/2 = {(q-2)//2}: {n14 == (q-2)//2};  irreducible quintic: {types[(5,)]}")
    # correspondence with W for k = 3, 5, 7
    import legendre_d6_qf_mixed as MX
    for k in (3, 5, 7):
        q = 1 << k
        rows, Wj = MX.run(k, do_S=False)
        P = Poly(q)
        ok = True
        for j, W in Wj.items():
            c = int(P.M[1, 1])  # placeholder
            inv = int(P.INV[j])
            # sqrt in F_q: x^(q/2)
            c = inv
            for _ in range(k - 1):
                c = int(P.M[c, c])
            typ = factor_type(P, c)
            if (W != 0) != (typ == (1, 4)):
                ok = False
                print(f"   mismatch at j={j}: W={W}, type={typ}")
        print(f"k={k}: [W(j) != 0  <=>  v^5+v+j^(-1/2) has type (1,4)]: {'OK' if ok else 'FAILED'}")


if __name__ == "__main__":
    main()
