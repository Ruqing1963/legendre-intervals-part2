"""
witt_lpoly.py -- the order-4 (Witt) character sums behind X through their L-polynomials.

S''(l1,l3,l5) = sum_{x in F_{q^12}} i^{T(x)} (-1)^{Tr(l5 x^5 + l3 x^3 + l1 x)},  T(x) = Tr_{W_2(F_{q^12})/W_2(F_2)}([x]) in Z/4,
is the trace of Frob^12 on the rank-4 Witt-Artin-Schreier sheaf L_{psi_4}(x, f(x)) over F_q (Swan conductor 5 at infinity).
The same sums over F_{q^m}, m = 1..4,
    s_m = sum_{x in F_{q^m}} i^{T_m(x)} (-1)^{Tr_{q^m/2} f(x)},   T_m(x) = Tr_{W_2(F_{q^m})/Z4}[x],
determine the four Frobenius eigenvalues (Newton over Z[i]) and hence S''_12 = -sum omega_i^12.  We compare with the
exact S'(l1) computed as an F_2-quadratic-form sum on F_{q^12} (x-form of Section 6.7 of Part II), class by class
(class = factorisation type of the sextic B(u) = l5^2 u^6 + l3^2 u^2 + l5 u + 1).
Usage: python witt_lpoly.py 3   [k=5 takes a few minutes]
"""
import itertools
import os
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_d6_qf as QF  # noqa: E402
import legendre_var_check as LV  # noqa: E402
from hexad_search import Field  # noqa: E402
from quintic_check import Poly  # noqa: E402
from legendre_d6_Wp_pattern import factor_degrees  # noqa: E402

# odd degrees needed for F_{q^m}, m = 1, 3 (the GF constructor verifies primitivity)
LV.PRIMITIVE.update({3: 0b1011, 5: 0b100101, 7: 0b10000011, 9: 0b1000010001, 15: 0b1000000000000011,
                     21: 0b1000000000000000000101})
# even k (q = 4, 16): degrees 2, 4, 8, 12, 16 (only added if absent; the GF constructor checks primitivity)
for _d, _p in {2: 0b111, 4: 0b10011, 8: 0b100011101, 12: 0b1000001010011, 16: 0b10001000000001011}.items():
    LV.PRIMITIVE.setdefault(_d, _p)


def witt_sums(k, m, l1, l3, l5, gf, Fq_in_gf):
    """s_m for the field F_{2^{mk}} given as gf = LV.GF(m*k); l's are elements of F_q embedded (Fq_in_gf maps F_q ints)."""
    n = m * k
    order = gf.order
    x = np.arange(1, gf.size, dtype=np.int64)        # all nonzero elements by value
    logx = gf.log[x]
    # f(x) = l5 x^5 + l3 x^3 + l1 x in gf
    def mono(c, e):
        if c == 0:
            return np.zeros_like(x)
        return gf.exp[(gf.log[c] + e * logx) % order]
    fx = mono(Fq_in_gf[l5], 5) ^ mono(Fq_in_gf[l3], 3) ^ mono(Fq_in_gf[l1], 1)
    trf = LV._abs_trace_parity(gf, fx.astype(np.uint32)).astype(np.int64)
    # Witt trace of [x]: iterate over the n conjugates with carry: (s0,s1) += (c,0): s1 ^= s0*c, s0 ^= c
    s0 = np.zeros_like(x); s1 = np.zeros_like(x)
    c = x.copy()
    for _ in range(n):
        prod = np.zeros_like(x)
        mask = (s0 != 0) & (c != 0)
        prod[mask] = gf.exp[(gf.log[s0[mask]] + gf.log[c[mask]]) % order]
        s1 ^= prod
        s0 ^= c
        c = gf.exp[(2 * gf.log[c]) % order]
    # final (s0, s1) in W_2(F_2): s0, s1 in {0,1}
    assert np.all(s0 <= 1) and np.all(s1 <= 1)
    T = (s0 + 2 * s1) % 4
    expo = (T + 2 * trf) % 4
    i4 = np.array([1, 1j, -1, -1j])
    total = i4[expo].sum() + 1          # x = 0 contributes i^0 = 1
    return complex(round(total.real), round(total.imag))


def newton_p12(s):
    """eigenvalues omega_i (4 of them) with p_m = -s_m; returns p_12 (Gaussian integer)."""
    p = {m: -s[m] for m in (1, 2, 3, 4)}
    e1 = p[1]
    e2 = (e1 * p[1] - p[2]) / 2
    e3 = (e2 * p[1] - e1 * p[2] + p[3]) / 3
    e4 = (e3 * p[1] - e2 * p[2] + e1 * p[3] - p[4]) / 4
    e = [1, e1, e2, e3, e4]
    for m in range(5, 13):
        p[m] = e[1] * p[m - 1] - e[2] * p[m - 2] + e[3] * p[m - 3] - e[4] * p[m - 4]
    return p[12], e


def exact_Sprime(Fd, l1, l3, l5, cache):
    """S'(l1,l3,l5) = sum_{x in F_{q^12}} (-1)^{Q'(x) + Tr(l1 x)} by the quadratic-form algorithm (x-form)."""
    key = (l3, l5)
    if key not in cache:
        N = Fd.N
        M = Fd.M
        I = np.eye(N, dtype=np.uint8)
        SQ, SQRT = Fd.SQ, Fd.SQRT
        SQ2 = SQ.dot(SQ) % 2; SQRT2 = SQRT.dot(SQRT) % 2
        scal = lambda c: Fd.K12.linmap_matrix(lambda a: Fd.scal(c, a))
        Tp = (scal(l5).dot(SQ2) + scal(l3).dot(SQ) + I + scal(int(Fd.F.SQRT[l3])).dot(SQRT)
              + scal(int(Fd.F.SQRT[int(Fd.F.SQRT[l5])])).dot(SQRT2)) % 2
        Tpp = (Tp + Fd.TR12q) % 2
        basis = Fd.K12.basis()
        G = np.zeros((N, N), np.uint8)
        for i in range(N):
            for j in range(i, N):
                G[i, j] = G[j, i] = Fd.abs_trace(Fd.mul(basis[i], basis[j]))
        B = (Tpp.T.dot(G)) % 2
        Qv = np.zeros(N, np.uint8); Trv = np.zeros(N, np.int64)
        for i, e in enumerate(basis):
            e2 = Fd.mul(e, e); e3 = Fd.mul(e2, e); e5 = Fd.mul(e3, e2)
            t = Fd.abs_trace(Fd.add(Fd.scal(l5, e5), Fd.scal(l3, e3)))
            conj = [e]
            for _ in range(11):
                conj.append(Fd.frob(conj[-1]))
            s2 = Fd.zero
            for a in range(12):
                for b in range(a + 1, 12):
                    s2 = Fd.add(s2, Fd.mul(conj[a], conj[b]))
            Qv[i] = (t + Fd.F.tr2(s2[0])) % 2
            Trv[i] = Fd.apply(Fd.TR12q, e)[0]
        cache[key] = (B, Qv, Trv)
    B, Qv, Trv = cache[key]
    Lv = np.array([Fd.F.tr2(int(Fd.M[l1, t])) for t in Trv], dtype=np.uint8)
    val, rdim, info = QF.quad_exp_sum(B, Qv, Lv)
    return val, rdim


def run(k, only_classes=None, max_pairs=None):
    q = 1 << k
    Fd = Field(k)
    P = Poly(q)
    gfs = {m: LV.GF(m * k) for m in (1, 2, 3, 4)}
    # embedding of F_q (fqlib encoding) into each gf: both are F_2[x]/(primitive) but with different moduli;
    # build the embedding via a primitive element: map F_q = fqlib tables onto gf's subfield of order q.
    emb = {}
    for m, gf in gfs.items():
        step = gf.order // (q - 1)
        g = int(gf.exp[step])                       # generator of the subfield F_q inside gf
        # find the minimal polynomial relation: fqlib's F_q is Fq(q) with primitive element 2 (x); find an element of
        # gf's subfield with the same minimal polynomial by brute force over the q-1 candidates.
        F = QF.Fq(q)
        # minimal polynomial of x in fqlib encoding: use the field's modulus: compute powers of 2 in fqlib
        def fq_pow(a, e):
            r = 1
            for _ in range(e):
                r = int(F.M[r, a])
            return r
        target = None
        for t in range(q - 1):
            cand = int(gf.exp[(t * step) % gf.order])
            # test: the map 2^j -> cand^j must be additive: check cand^k == sum of lower powers as in fqlib's modulus
            ok = True
            tab = {}
            a_pow = 1; c_pow = 1
            for j in range(q - 1):
                tab[a_pow] = c_pow
                a_pow = int(F.M[a_pow, 2]); c_pow = int(gf.exp[(gf.log[c_pow] + gf.log[cand]) % gf.order]) if c_pow else 0
            tab[0] = 0
            for a in range(1, min(q, 64)):
                for b in range(1, min(q, 64)):
                    if tab[int(F.A[a, b])] != (tab[a] ^ tab[b]):
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                target = tab
                break
        emb[m] = target
    _, _, results = QF.run(k, verbose=False)
    F = QF.Fq(q)
    Wp = {}
    for (a, b), (val, rdim, info, n) in results["Wp"].items():
        for _ in range(k):
            Wp[(a, b)] = val
            a, b = int(F.SQ[a]), int(F.SQ[b])
    stats = Counter()
    eig = Counter()
    sm_stats = Counter()
    pairs = [(l3, l5) for l5 in range(1, q) for l3 in range(q)]
    per_class = Counter()
    unit = 2 ** ((12 * k + 4) // 2)                    # 4 q^6: the modulus of a nonzero Witt sum with 4-dim radical
    for (l3, l5) in pairs:
        B6 = [1, l5, int(P.M[l3, l3]), 0, 0, 0, int(P.M[l5, l5])]
        cls = factor_degrees(P, B6, 6)
        if only_classes and cls not in only_classes:
            continue
        if max_pairs and per_class[cls] >= max_pairs:
            continue
        per_class[cls] += 1
        preds = []
        patt = Counter()
        for l1 in range(q):
            s = {m: witt_sums(k, m, l1, l3, l5, gfs[m], emb[m]) for m in (1, 2, 3, 4)}
            p12, e = newton_p12(s)
            pred = -p12
            preds.append(pred)
            v = pred / unit
            vr = complex(round(v.real * 4) / 4, round(v.imag * 4) / 4)
            patt[(vr, tuple(s[m] for m in (1, 2, 3, 4)))] += 1
            # normalised eigenvalues zeta_i = omega_i / sqrt(q) as roots of unity: exponents a_i/120 (120 = lcm(24, 40))
            roots = np.roots([1, -e[1], e[2], -e[3], e[4]]) / np.sqrt(q)
            raw = [np.angle(z) / (2 * np.pi) * 120 for z in roots]
            # repeated roots are only located to ~1e-4 by np.roots; the exponents are integers, so a loose tolerance is safe
            assert all(abs(a - round(a)) < 0.05 for a in raw) and all(abs(abs(z) - 1) < 1e-3 for z in roots), (raw, roots)
            expo = tuple(sorted(int(round(a)) % 120 for a in raw))
            eig[(cls, vr, expo)] += 1
            # the small-field Gauss sums themselves: s_m / q^{m/2} is 0 or 2^{r_m/2} times an 8th root of unity
            # (r_m = dim of the radical over F_{q^m}); record (r_m, exponent b_m of zeta_8) or '0'
            smn = []
            for m in (1, 2, 3, 4):
                z = s[m] / np.sqrt(q) ** m
                if abs(z) < 1e-9:
                    smn.append("0")
                else:
                    rm = int(round(2 * np.log2(abs(z))))
                    b = np.angle(z) / (2 * np.pi) * 8
                    assert abs(b - round(b)) < 1e-6, (m, s[m], z)
                    smn.append(f"r{rm}:z8^{int(round(b)) % 8}")
            sm_stats[(cls, vr, tuple(smn))] += 1
        tot = sum(preds)
        ok = abs(tot - Wp[(l3, l5)]) < 1e-6
        multiset = tuple(sorted(Counter(v for (v, _) in patt.elements()).items(), key=str))
        stats[(cls, Wp[(l3, l5)] // q ** 7 if Wp[(l3, l5)] % q ** 7 == 0 else Wp[(l3, l5)], ok, multiset)] += 1
    print(f"\nk={k}: (class, W'/q^7, sum over l1 of the L-polynomial predictions == W', multiset of S''(l1)/(4q^6)) -> #(l3,l5)")
    for key, n in sorted(stats.items(), key=str):
        print(f"   {key}  x{n}")
    print(f"  normalised Frobenius eigenvalues zeta_i = omega_i/sqrt(q) of the Witt sheaf, as exponents a_i with zeta_i = exp(2 pi i a_i/120):")
    print(f"  (class, S''(l1)/(4q^6), sorted exponents) -> #(l3,l5,l1)")
    for key, n in sorted(eig.items(), key=str):
        print(f"   {key}  x{n}")
    print(f"  small-field Gauss sums s_m / q^(m/2), m = 1..4, as 'r<radical dim>:z8^<b>' or '0':")
    print(f"  (class, S''(l1)/(4q^6), (s_1, s_2, s_3, s_4)) -> #(l3,l5,l1)")
    for key, n in sorted(sm_stats.items(), key=str):
        print(f"   {key}  x{n}")


if __name__ == "__main__":
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    mp = int(sys.argv[2]) if len(sys.argv) > 2 else None
    run(k, max_pairs=mp)
