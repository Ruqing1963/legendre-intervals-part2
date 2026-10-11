"""
sextic_identity.py -- the hexad of ker L' is the root set of the sextic  B(u) = l5^2 u^6 + l3^2 u^2 + l5 u + 1 = u A(u) + 1.

(a) symbolic: over F_2[l3, l5], reduce L'(r + s) modulo B(r) and B(s) (polynomial division in r, then in s); the
    remainder is zero, so every pairwise sum of roots of B lies in ker L'.  Also checks the identity
    x^4 L'(x) = (x^4 B(x))^2 + x^4 B(x) + x^4  (equivalently L'(x) = x^4 B(x)^2 + B(x) + 1).
(b) numeric: for k = 3, 5 and every (l3, l5) with l5 != 0 such that B splits in F_{q^12}, the 15 pairwise sums of the
    six roots are exactly the nonzero elements of ker L'.
"""
import os
import sys
import itertools

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def symbolic():
    l3, l5, r, s, x = sp.symbols("l3 l5 r s x")
    dom = sp.GF(2)
    B = lambda u: l5 ** 2 * u ** 6 + l3 ** 2 * u ** 2 + l5 * u + 1
    Lp = lambda u: l5 ** 4 * u ** 16 + l3 ** 4 * u ** 8 + u ** 4 + l3 ** 2 * u ** 2 + l5 * u
    # identity L'(x) = x^4 B(x)^2 + B(x) + 1 over F_2
    diff = sp.Poly(sp.expand(Lp(x) - (x ** 4 * B(x) ** 2 + B(x) + 1)), x, l3, l5, domain=dom)
    print("L'(x) == x^4 B(x)^2 + B(x) + 1 over F_2:", diff.is_zero)
    # reduce L'(r+s) modulo B(r), B(s): division by the monic (in r) polynomial B(r)/l5^2 needs the fraction field;
    # instead multiply by l5^(2m) and use pseudo-division over F_2[l3,l5].
    P = sp.Poly(sp.expand(Lp(r + s)), r, s, l3, l5, domain=dom)
    Br = sp.Poly(sp.expand(B(r)), r, s, l3, l5, domain=dom)
    Bs = sp.Poly(sp.expand(B(s)), r, s, l3, l5, domain=dom)
    G = sp.groebner([Br.as_expr(), Bs.as_expr()], r, s, l3, l5, modulus=2, order="lex")
    rem = G.reduce(P.as_expr())[1]
    print("L'(r+s) reduced modulo the ideal (B(r), B(s)) in F_2[r,s,l3,l5]:", sp.expand(rem) == 0, "" if sp.expand(rem) == 0 else rem)
    # sum of the roots and second symmetric function vanish: coefficients of u^5, u^4 are 0 (visible), B' = l5 (separable)
    print("B'(u) =", sp.Poly(sp.diff(B(x), x), x, l3, l5, domain=dom).as_expr())


def numeric(k):
    import hexad_search as HS
    import legendre_d6_qf as QF
    Fd = HS.Field(k)
    q = Fd.q
    M, A, INV = Fd.M, Fd.A, Fd.F.INV
    ok = bad = nosplit = 0
    for l3 in range(q):
        for l5 in range(1, q):
            # roots of B in F_{q^12}: solve via the kernel is circular; find roots by factoring numerically:
            # B(u) = 0 <=> u A(u) = 1; search roots as elements u with L'-related structure is also circular, so scan
            # the (q^12) field is impossible -- instead compute ker L' and test whether its 15 elements are pairwise
            # sums of a 6-set of roots of B: for each kernel element x, the roots r with r + ? ... simpler: take ker L',
            # pick the 168 frames, and check which hexad {0} u frame, translated by sqrt(b4), consists of roots of B.
            Mx = Fd.Lprime_matrix(l3, l5)
            ns = QF.f2_nullspace(Mx)
            if len(ns) != 4:
                nosplit += 1
                continue
            elems = []
            for coeffs in itertools.product((0, 1), repeat=4):
                v = sum((b for c, b in zip(coeffs, ns) if c), start=0 * ns[0])
                elems.append(tuple(Fd.K12.frombits(v % 2)))
            zero = tuple(Fd.zero)
            nonzero = [e for e in elems if e != zero]
            l32, l52 = int(M[l3, l3]), int(M[l5, l5])

            def Bval(u):
                u2 = Fd.mul(u, u); u6 = Fd.mul(Fd.mul(u2, u2), u2)
                return Fd.add(Fd.add(Fd.scal(l52, u6), Fd.scal(l32, u2)), Fd.add(Fd.scal(l5, u), Fd.one))
            found = False
            bits = {e: Fd.K12.tobits(list(e)) for e in nonzero}
            for quad in itertools.combinations(nonzero, 4):
                if QF.f2_rank([bits[e] for e in quad]) != 4:
                    continue
                s5 = zero
                for e in quad:
                    s5 = tuple(Fd.add(list(s5), list(e)))
                H = list(quad) + [s5, zero]
                P = Fd.poly_from_roots([list(h) for h in H])
                c = Fd.sqrt(P[4])
                roots = [Fd.add(list(h), c) for h in H]
                if all(all(v == 0 for v in Bval(rt)) for rt in roots):
                    found = True
                    break
            if found:
                ok += 1
            else:
                bad += 1
    print(f"k={k}: pairs whose kernel is the pairwise-sum set of a hexad of roots of B: {ok}; failures: {bad}; kernel not in F_(q^12): {nosplit}")


if __name__ == "__main__":
    symbolic()
    for k in (3, 5):
        numeric(k)
