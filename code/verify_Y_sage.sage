# verify_Y_sage.sage -- SageMath check of the constant Y_k of the paper (Part II, Sections 4 and 7):
#
#     Y_k = sum_{alpha in U} (-1)^{Tr_{q/2} e_4(alpha)},
#     U   = { alpha in F_{q^12} : Tr alpha = Tr alpha^3 = Tr alpha^5 = 0 },  q = 2^k,  Tr = Tr_{q^12/q},
#
# by direct enumeration of the point set (route A, q = 2) and by the pushforward to the space of monic
# polynomials t^12 + c_2 t^10 + c_4 t^8 + ... with Lambda-weights (route B, q = 2, 4).
# Expected: Y_1 = 104 (m_1 = 13), Y_2 = -64 (m_2 = -1).   Run:  sage verify_Y_sage.sage
#
# Route B for q = 4 factors 4^9 = 262144 polynomials and takes a few minutes.

EXPECTED = {2: 104, 4: -64}

def route_a(k):
    q = 2^k
    K = GF(q^12, 'w')
    k2 = GF(2)
    Y = 0; size = 0
    for a in K:
        conj = [a^(q^j) for j in range(12)]
        R = PolynomialRing(K, 't'); t = R.gen()
        F = prod(t - z for z in conj)
        c = F.list()                      # c[12-m] = e_m (signs irrelevant in characteristic 2)
        e1, e3, e4, e5 = c[11], c[9], c[8], c[7]
        if e1 == 0 and e3 == 0 and e5 == 0:
            size += 1
            # e4 lies in F_q; its absolute trace Tr_{q/2} equals Tr_{F_{q^12}/F_2}(e4)/12 ... use the subfield:
            Fq = GF(q, 'u'); emb = Fq.hom([Fq.gen()])  # placeholder to keep names; compute trace via K
            tr = sum(e4^(2^i) for i in range(k))       # Tr_{q/2}(e4) computed inside K (e4 in F_q)
            assert tr in [K(0), K(1)]
            Y += 1 - 2 * (1 if tr == K(1) else 0)
    return Y, size

def route_b(k):
    q = 2^k
    Fq = GF(q, 'u')
    R = PolynomialRing(Fq, 't'); t = R.gen()
    Y = 0; size = 0
    for coeffs in cartesian_product([Fq]*9):            # (c2, c4, e6, e7, e8, e9, e10, e11, e12)
        c2, c4, e6, e7, e8, e9, e10, e11, e12 = coeffs
        F = t^12 + c2*t^10 + c4*t^8 + e6*t^6 + e7*t^5 + e8*t^4 + e9*t^3 + e10*t^2 + e11*t + e12
        fac = F.factor()
        if len(fac) != 1:
            continue                                   # Lambda(F) = 0 unless F is a prime power
        P, m = fac[0]
        lam = P.degree()                               # Lambda(P^m) = deg P
        size += lam
        tr = c4.trace()                                # Tr_{F_q/F_2}
        Y += lam * (1 - 2 * Integer(tr))
    return Y, size

Ya, Ua = route_a(1)
print("q=2 route A: |U| =", Ua, " Y =", Ya); assert Ya == EXPECTED[2]
for k in (1, 2):
    Yb, Ub = route_b(k)
    print("q=%d route B: |U| = %d  Y = %d  m_k = %d" % (2^k, Ub, Yb, Yb / 2^(3*k))); assert Yb == EXPECTED[2^k]
print("all assertions passed")
