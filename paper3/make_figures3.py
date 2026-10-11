r"""
make_figures3.py -- figures of Part III (written next to this script):
  fig_hexad_orbits.pdf   : the hexad {r_1..r_6}, its 15 duads = K \ 0 coloured by Frobenius orbit, for four classes;
  fig_eigen_circles.pdf  : normalised Frobenius eigenvalues of the rank-4 Witt sheaf on the unit circle, one panel per
                           class (q = 8 data of witt_lpoly.py, one representative (l3,l5,l1) per class and S'' value);
  fig_class_counts.pdf   : the eleven class counts versus q (data q = 8..512) against the polynomials of Table 2;
  fig_deformation.pdf    : the one-parameter family L_lambda + eps x^4 <-> u A(u) + sqrt(eps): quintic frame at eps = 0,
                           hexad at eps != 0 (schematic).
Run:  python make_figures3.py
"""
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "mathtext.fontset": "cm"})
COL = ["#1b6ca8", "#d1495b", "#edae49", "#00798c", "#66a182", "#8d6a9f", "#2e4057"]


def perm_orbits(perm):
    """orbits of the permutation (dict i->perm(i)) on the 15 duads"""
    duads = list(itertools.combinations(range(6), 2))
    seen, orbits = set(), []
    for d in duads:
        if d in seen:
            continue
        o, x = [], d
        while x not in seen:
            seen.add(x)
            o.append(x)
            x = tuple(sorted((perm[x[0]], perm[x[1]])))
        orbits.append(o)
    return orbits


def fig_hexad_orbits():
    classes = [("$3\\,2\\,1$: orbits (1,2,3,3,6)", {0: 1, 1: 2, 2: 0, 3: 4, 4: 3, 5: 5}),
               ("$4\\,1^2$: orbits (1,2,4,4,4)", {0: 1, 1: 2, 2: 3, 3: 0, 4: 4, 5: 5}),
               ("$5\\,1$: orbits (5,5,5)", {0: 1, 1: 2, 2: 3, 3: 4, 4: 0, 5: 5}),
               ("$2^3$: orbits (1,1,1,2,2,2,2,2,2)", {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4})]
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.0))
    ang = np.pi / 2 + 2 * np.pi * np.arange(6) / 6
    P = np.c_[np.cos(ang), np.sin(ang)]
    for ax, (title, perm) in zip(axes, classes):
        orbits = perm_orbits(perm)
        orbits.sort(key=lambda o: (len(o), o))
        for ci, o in enumerate(orbits):
            for (i, j) in o:
                ax.plot([P[i, 0], P[j, 0]], [P[i, 1], P[j, 1]], color=COL[ci % len(COL)], lw=1.6 if len(o) > 1 else 2.4,
                        ls="-" if len(o) > 1 else "--", alpha=0.9, zorder=1)
        for i in range(6):
            ax.scatter(*P[i], s=120, color="white", edgecolor="black", zorder=3)
            ax.text(P[i, 0], P[i, 1], f"$r_{i+1}$", ha="center", va="center", fontsize=8, zorder=4)
        # Frobenius arrows on the roots
        for i in range(6):
            j = perm[i]
            if j != i:
                mid = 0.82 * (P[i] + 0.5 * (P[j] - P[i]))
                ax.annotate("", xy=1.13 * P[j], xytext=1.13 * P[i],
                            arrowprops=dict(arrowstyle="->", color="gray", lw=0.8, connectionstyle="arc3,rad=0.25"))
        ax.set_aspect("equal")
        ax.set_xlim(-1.35, 1.35)
        ax.set_ylim(-1.35, 1.35)
        ax.axis("off")
        ax.set_title(title, fontsize=9)
    fig.suptitle("$\\mathcal{K}_\\lambda\\setminus 0$ = the 15 duads $r_i+r_j$; colours = Frobenius orbits "
                 "(dashed: fixed duads, in $\\mathbb{F}_q$); grey arrows: Frobenius on the roots", fontsize=9, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(OUT, "fig_hexad_orbits.pdf"))
    plt.close(fig)


def fig_eigen_circles():
    # representative exponent sets a_j (zeta_j = exp(2 pi i a_j / 120)) from data/witt_lpoly_k3_b120.txt
    data = [("$3\\,1^3$,  $S''=+4q^6$", (5, 25, 65, 85), 24), ("$2^3$,  $S''=+4q^6$", (15, 45, 75, 105), 24),
            ("$3^2$,  $S''=+4q^6$", (5, 45, 45, 85), 24), ("$6$,  $S''=+4q^6$", (5, 45, 85, 105), 24),
            ("$2^2 1^2$,  $S''=-4q^6$", (0, 30, 60, 90), 24), ("$3\\,2\\,1$,  $S''=-4q^6$", (10, 50, 70, 110), 24),
            ("$4\\,1^2$,  $S''=+4q^6$", (15, 45, 75, 105), 24), ("$4\\,1^2$,  $S''=0$", (0, 15, 45, 60), 24),
            ("$4\\,2$,  $S''=-4q^6$", (0, 30, 60, 90), 24), ("$4\\,2$,  $S''=0$", (0, 45, 60, 75), 24),
            ("$5\\,1$,  $S''=-q^6$", (9, 33, 57, 81), 40), ("$5\\,1$,  $S''=-q^6$", (3, 27, 51, 99), 40)]
    fig, axes = plt.subplots(3, 4, figsize=(11, 8.2))
    th = np.linspace(0, 2 * np.pi, 400)
    for ax, (title, expo, N) in zip(axes.flat, data):
        ax.plot(np.cos(th), np.sin(th), color="lightgray", lw=1)
        # the N-th roots of unity as ticks
        for a in range(N):
            z = np.exp(2j * np.pi * a / N)
            ax.plot([0.97 * z.real, 1.03 * z.real], [0.97 * z.imag, 1.03 * z.imag], color="lightgray", lw=0.6)
        # 12th powers: the sign of S''
        cnt = {}
        for a in expo:
            cnt[a] = cnt.get(a, 0) + 1
        for a, mult in cnt.items():
            z = np.exp(2j * np.pi * a / 120)
            z12 = z ** 12
            ax.scatter(z.real, z.imag, s=70 + 60 * (mult - 1), color=COL[0] if abs(z12 + 1) < 1e-6 else (COL[1] if abs(z12 - 1) < 1e-6 else COL[2]),
                       edgecolor="black", zorder=3)
            if mult > 1:
                ax.text(1.18 * z.real, 1.18 * z.imag, f"$\\times{mult}$", ha="center", va="center", fontsize=8)
        if N == 40:
            # the missing fifth root: s_1/sqrt(q)
            allf = sorted(set((expo[0] + 24 * t) % 120 for t in range(5)))
            miss = [a for a in allf if a not in expo][0]
            z = np.exp(2j * np.pi * miss / 120)
            ax.scatter(z.real, z.imag, s=70, facecolor="white", edgecolor=COL[3], lw=1.5, zorder=3)
            ax.text(1.2 * z.real, 1.2 * z.imag, "$\\eta$", ha="center", va="center", fontsize=9, color=COL[3])
        ax.set_aspect("equal")
        ax.set_xlim(-1.35, 1.35)
        ax.set_ylim(-1.35, 1.35)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title + f"   ($N={N}$)", fontsize=9)
    fig.suptitle("Normalised Frobenius eigenvalues $\\omega_j/\\sqrt{q}$ of the rank-4 Witt sheaf, $q=8$. "
                 "Blue: $\\zeta^{12}=-1$;  red: $\\zeta^{12}=+1$;  $S''=-q^6\\sum\\zeta_j^{12}$", fontsize=9, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(OUT, "fig_eigen_circles.pdf"))
    plt.close(fig)


def fig_class_counts():
    qs = np.array([8, 32, 128, 512], dtype=float)
    polys = {"$1^6$": lambda q: (q - 2) * (q - 8) / 720, "$2\\,1^4$": lambda q: (q - 2) * (q - 8) / 48,
             "$2^2 1^2$": lambda q: q * (q - 2) / 16, "$2^3$": lambda q: (q - 2) * (q - 8) / 48,
             "$3\\,1^3$": lambda q: (q + 1) * (q - 2) / 18, "$3^2$": lambda q: (q + 1) * (q - 2) / 18,
             "$3\\,2\\,1$": lambda q: (q + 1) * (q - 2) / 6, "$4\\,1^2$": lambda q: (q - 2) * (q - 4) / 8,
             "$4\\,2$": lambda q: q * (q - 2) / 8, "$5\\,1$": lambda q: (q * q + 1) / 5, "$6$": lambda q: (q + 1) * (q - 2) / 6}
    props = {"$1^6$": 1 / 720, "$2\\,1^4$": 1 / 48, "$2^2 1^2$": 1 / 16, "$2^3$": 1 / 48, "$3\\,1^3$": 1 / 18, "$3^2$": 1 / 18,
             "$3\\,2\\,1$": 1 / 6, "$4\\,1^2$": 1 / 8, "$4\\,2$": 1 / 8, "$5\\,1$": 1 / 5, "$6$": 1 / 6}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.0))
    qq = np.linspace(8, 600, 300)
    markers = ["o", "s", "^", "v", "D", "P", "X", "<", ">", "*", "h"]
    for i, (name, f) in enumerate(polys.items()):
        c = COL[i % len(COL)]
        y = np.array([f(q) for q in qq])
        ax1.plot(qq, np.maximum(y, 0.5), color=c, lw=1.2)
        ax1.scatter(qs, np.maximum([f(q) for q in qs], 0.5), marker=markers[i], color=c, edgecolor="black", s=36, zorder=3, label=name)
        dev = np.array([f(q) - props[name] * q * (q - 1) for q in qq])
        ax2.plot(qq, dev, color=c, lw=1.2, label=name)
    ax1.set_xscale("log", base=2)
    ax1.set_yscale("log", base=2)
    ax1.set_xlabel("$q$")
    ax1.set_ylabel("number of $(\\lambda_3,\\lambda_5)$, $\\lambda_3\\neq 0$")
    ax1.set_title("class counts of $B_\\lambda$ (points: exact data $q=8,\\dots,512$; lines: Table 2)")
    ax1.legend(ncol=2, fontsize=7, frameon=False)
    ax2.axhline(0, color="gray", lw=0.6)
    ax2.set_xlabel("$q$")
    ax2.set_ylabel("count $-$ $|C|/720\\cdot q(q-1)$")
    ax2.set_title("deviation from Chebotarev equidistribution (exactly linear in $q$)")
    ax2.legend(ncol=2, fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_class_counts.pdf"))
    plt.close(fig)


def fig_deformation():
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
    # eps = 0: quintic frame {0, u_0..u_4}; eps != 0: hexad {r_1..r_6}; the duads
    ang5 = np.pi / 2 + 2 * np.pi * np.arange(5) / 5
    P5 = np.c_[np.cos(ang5), np.sin(ang5)]
    ax = axes[0]
    pts = np.vstack([[0, 0], P5])
    for i, j in itertools.combinations(range(6), 2):
        ax.plot([pts[i, 0], pts[j, 0]], [pts[i, 1], pts[j, 1]], color=COL[0] if 0 in (i, j) else COL[2], lw=1.2, alpha=0.8)
    ax.scatter(0, 0, s=120, color="white", edgecolor="black", zorder=3)
    ax.text(0, 0, "$0$", ha="center", va="center", fontsize=8, zorder=4)
    for i in range(5):
        ax.scatter(*P5[i], s=120, color="white", edgecolor="black", zorder=3)
        ax.text(P5[i, 0], P5[i, 1], f"$u_{i}$", ha="center", va="center", fontsize=8, zorder=4)
    ax.set_title("$\\varepsilon=0$: $B_{\\lambda,0}=uA_\\lambda(u)$, roots $0,u_0,\\dots,u_4$\n"
                 "$\\mathcal{K}^0_\\lambda\\setminus 0=\\{u_i\\}\\cup\\{u_i+u_j\\}$, group $S_5$", fontsize=8.5)
    ax = axes[1]
    ang6 = np.pi / 2 + 2 * np.pi * np.arange(6) / 6
    P6 = np.c_[np.cos(ang6), np.sin(ang6)]
    for i, j in itertools.combinations(range(6), 2):
        ax.plot([P6[i, 0], P6[j, 0]], [P6[i, 1], P6[j, 1]], color=COL[2], lw=1.2, alpha=0.8)
    for i in range(6):
        ax.scatter(*P6[i], s=120, color="white", edgecolor="black", zorder=3)
        ax.text(P6[i, 0], P6[i, 1], f"$r_{i+1}$", ha="center", va="center", fontsize=8, zorder=4)
    ax.set_title("$\\varepsilon=\\delta^2\\neq 0$: $B_{\\lambda,\\delta}=uA_\\lambda(u)+\\delta$, roots $r_1,\\dots,r_6$\n"
                 "$\\mathcal{K}_\\lambda\\setminus 0=\\{r_i+r_j\\}$, group $S_6\\cong Sp_4(\\mathbb{F}_2)$", fontsize=8.5)
    ax = axes[2]
    # the generalised quadrangle W(2): 15 duads as points on a circle, 15 synthemes as triangles (drawn as 3 chords)
    duads = list(itertools.combinations(range(6), 2))
    ang15 = np.pi / 2 + 2 * np.pi * np.arange(15) / 15
    P15 = np.c_[np.cos(ang15), np.sin(ang15)]
    synthemes = []
    for a, b, c in itertools.combinations(range(15), 3):
        s = set(duads[a]) | set(duads[b]) | set(duads[c])
        if len(s) == 6:
            synthemes.append((a, b, c))
    for si, (a, b, c) in enumerate(synthemes):
        tri = np.array([P15[a], P15[b], P15[c], P15[a]])
        ax.plot(tri[:, 0], tri[:, 1], color=COL[si % len(COL)], lw=0.7, alpha=0.6)
    for d, (i, j) in enumerate(duads):
        ax.scatter(*P15[d], s=70, color="white", edgecolor="black", zorder=3)
        ax.text(1.17 * P15[d, 0], 1.17 * P15[d, 1], f"{i+1}{j+1}", ha="center", va="center", fontsize=6.5)
    ax.set_title("$W(2)$: 15 duads (points), 15 synthemes (isotropic planes)\n"
                 "$\\omega(r_i+r_j,r_k+r_l)=|\\{i,j\\}\\cap\\{k,l\\}|$ mod 2", fontsize=8.5)
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlim(-1.4, 1.4)
        ax.set_ylim(-1.4, 1.4)
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_deformation.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    fig_hexad_orbits()
    fig_eigen_circles()
    fig_class_counts()
    fig_deformation()
    print("figures written to", OUT)
