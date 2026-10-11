"""
appendix_table_k3.py -- LaTeX rows for the appendix of Part III: all (l3, l5) in F_8 x F_8^* with the factorisation type
of B = l5^2 u^6 + l3^2 u^2 + l5 u + 1 over F_8 and the exact W'(l3, l5)/q^7 of Part II.  F_8 = F_2[a]/(a^3 + a + 1),
integers encode bit vectors in the basis 1, a, a^2.
Usage: python appendix_table_k3.py > ../paper3/table_k3.tex
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import legendre_d6_qf as QF  # noqa: E402
from quintic_check import Poly  # noqa: E402
from legendre_d6_Wp_pattern import factor_degrees  # noqa: E402
from hexad_signs import wp_table  # noqa: E402

NAMES = {0: "0", 1: "1", 2: r"\alpha", 3: r"\alpha+1", 4: r"\alpha^2", 5: r"\alpha^2+1", 6: r"\alpha^2+\alpha", 7: r"\alpha^2+\alpha+1"}
CLS = {(1, 1, 1, 1, 1, 1): "1^6", (1, 1, 1, 1, 2): "2\\,1^4", (1, 1, 2, 2): "2^21^2", (2, 2, 2): "2^3", (1, 1, 1, 3): "3\\,1^3",
       (3, 3): "3^2", (1, 2, 3): "3\\,2\\,1", (1, 1, 4): "4\\,1^2", (2, 4): "4\\,2", (1, 5): "5\\,1", (6,): "6"}


def main():
    k = 3
    q = 8
    P = Poly(q)
    Wp = wp_table(k)
    rows = []
    for l5 in range(1, q):
        for l3 in range(q):
            B = [1, l5, int(P.M[l3, l3]), 0, 0, 0, int(P.M[l5, l5])]
            typ = factor_degrees(P, B, 6)
            val = Wp[(l3, l5)][0]
            w = val // q ** 7
            assert w * q ** 7 == val
            rows.append((l3, l5, CLS[typ], w))
    rows.sort(key=lambda r: (r[2], r[0], r[1]))
    # three columns of rows
    n = len(rows)
    per = (n + 1) // 2
    cols = [rows[i * per:(i + 1) * per] for i in range(2)]
    print(r"\begin{tabular}{@{}cccc@{\qquad\qquad}cccc@{}}")
    print(r"\toprule")
    print(" & ".join([r"$\lambda_3$ & $\lambda_5$ & class & $W'/q^7$"] * 2) + r"\\")
    print(r"\midrule")
    for i in range(per):
        cells = []
        for c in cols:
            if i < len(c):
                l3, l5, cls, w = c[i]
                cells.append(f"${NAMES[l3]}$ & ${NAMES[l5]}$ & ${cls}$ & ${w}$")
            else:
                cells.append(" & & & ")
        print(" & ".join(cells) + r"\\")
    print(r"\bottomrule")
    print(r"\end{tabular}")


if __name__ == "__main__":
    main()
