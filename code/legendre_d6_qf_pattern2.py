"""For odd k: W(l3,l5) (mixed) is a function of j = l3^5/l5^3 in F_q^*. Find the exponents m with
   {j : W(j) != 0} == {j != 0 : Tr(j^m) = c}  (c = 0 or 1), and test a few other predicates."""
import sys

import legendre_d6_qf as QF

for k in [int(a) for a in sys.argv[1:]] or [3, 5]:
    q = 1 << k
    A6, X, res = QF.run(k, verbose=False)
    F = QF.Fq(q)
    M, INV, TR = F.M, F.INV, F.TR

    def pw(x, e):
        r = 1
        for _ in range(e):
            r = int(M[r, x])
        return r

    Wj = {}
    for (l3, l5), (v, r, info, n) in res["W"].items():
        if l3 and l5:
            for s in range(1, q):
                a, b = int(M[pw(s, 3), l3]), int(M[pw(s, 5), l5])
                for _ in range(k):
                    j = int(M[pw(a, 5), INV[pw(b, 3)]])
                    if j in Wj:
                        assert Wj[j] == v, "W is not a function of j"
                    Wj[j] = v
                    a, b = int(F.SQ[a]), int(F.SQ[b])
    assert len(Wj) == q - 1
    Jplus = {j for j, v in Wj.items() if v != 0}
    print(f"\nk={k}, q={q}: W is a function of j; #{{j: W!=0}} = {len(Jplus)} (= (q-2)/2 = {(q-2)//2})")
    hits = []
    for m in range(1, q - 1):
        for c in (0, 1):
            S = {j for j in range(1, q) if int(TR[pw(j, m)]) == c}
            if S == Jplus:
                hits.append((m, c))
    print(f"   exponents m with {{j: W!=0}} = {{j: Tr(j^m) = c}}: {hits}")
    # other predicates: Tr(j + 1/j), Tr(j^{(q-2)}) = Tr(1/j) etc.
    for name, f in [("Tr(j+1/j)", lambda j: int(TR[F.A[j, INV[j]]])), ("Tr(1/j)", lambda j: int(TR[INV[j]])),
                    ("Tr(j^3)", lambda j: int(TR[pw(j, 3)])), ("Tr(j^5)", lambda j: int(TR[pw(j, 5)])),
                    ("Tr(j^15)", lambda j: int(TR[pw(j, 15)])), ("Tr(j^(1/3))", lambda j: int(TR[pw(j, pow(3, -1, q - 1))])),
                    ("Tr(j^(1/5))", lambda j: int(TR[pw(j, pow(5, -1, q - 1))])),
                    ("Tr(j^(1/15))", lambda j: int(TR[pw(j, pow(15, -1, q - 1))]))]:
        S0 = {j for j in range(1, q) if f(j) == 0}
        print(f"   {name}=0 set equals J+: {S0 == Jplus};  Tr=1 set equals J+: {set(range(1,q)) - S0 == Jplus}")
