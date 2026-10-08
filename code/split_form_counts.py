"""
split_form_counts.py -- the UNTWISTED F_q-points of X = V(e_1, e_3, e_5) subset A^12 (defined over F_2).

X(F_q) = {z in F_q^12 : e_1(z) = e_3(z) = e_5(z) = 0}.  Grouping z by the multiset of its coordinates
(= a monic polynomial F in F_q[t] of degree 12 splitting completely over F_q), the number of ordered
tuples with that multiset is 12!/prod(m_i!).  We compute
    N_k   = #X(F_q)                                   (q = 2^k)
    T_k   = sum_{z in X(F_q)} (-1)^{Tr_{q/2} e_4(z)}    = Tr(F^k | H_c^*(X, L_psi(e_4)))  (untwisted traces)
for q = 2, 4, 8, to be compared with the twisted traces Y_k = Tr(sigma F^k | ...) of the paper.

CAVEAT (results: #X(F_q) = 2048, 2099200, 1081081856 for q = 2, 4, 8, i.e. about q^11 / 2^{k(k-1)/2}):
for q < 12 no polynomial over F_q has twelve distinct roots, so the untwisted counts are dominated by root
multisets with repetitions and are NOT in the asymptotic regime #X(F_q) ~ q^9 (which needs q >= 12, i.e. q >= 16).
They carry no information about the number of geometric components.  The twisted form (alpha in F_{q^12})
is in the asymptotic regime already for q = 2: |U| = q^9 + (-1)^{k+1}(q-1)q^5.
"""
import os
import sys
from itertools import combinations_with_replacement
from math import factorial

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "legendre-intervals-explicit-formula", "code"))
import fqlib  # noqa: E402

for k in (1, 2, 3):
    q = 1 << k
    A, M, NEG, INV = fqlib.field_tables(q)
    TR = [0] * q
    for x in range(q):
        s, y = 0, x
        for _ in range(k):
            s = int(A[s, y]); y = int(M[y, y])
        TR[x] = s
    N = T = 0
    for ms in combinations_with_replacement(range(q), 12):
        poly = [1]
        for r in ms:
            new = [0] * (len(poly) + 1)
            for i, c in enumerate(poly):
                new[i] = int(A[new[i], c])
                new[i + 1] = int(A[new[i + 1], M[c, r]])
            poly = new
        e = poly[1:]
        if e[0] or e[2] or e[4]:
            continue
        w = factorial(12)
        for x in set(ms):
            w //= factorial(ms.count(x))
        N += w
        T += w * (1 - 2 * TR[e[3]])
    print(f"q={q}: #X(F_q) = {N} = q^9 + ({N - q**9});  untwisted trace T_k = {T};  T_k/q^3 = {T / q**3}")
