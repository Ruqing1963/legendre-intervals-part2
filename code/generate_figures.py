"""generate_figures.py -- two vector figures for paper_part2.tex (written to ../paper/):
   fig_gold_cycle.pdf      : factorisation types of v^5 + v + c over F_q (k odd) -- counts (1,4): (q-2)/2,
                             (1,1,1,2): (q-2)/6, (2,3): (q+1)/3, as fractions of q-1, from quintic_check.py
                             (exact for k = 3..11), with the Frobenius cycle type and its effect on W.
   fig_singular_planes.pdf : the singular locus of X = V(e1,e3,e5): 1023 planes P_A stratified by |A| in {2,4,6},
                             all meeting in the diagonal line, with the Lucas-parity reason (schematic)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "paper")
plt.rcParams.update({"font.size": 10, "font.family": "serif", "mathtext.fontset": "cm",
                     "axes.spines.top": False, "axes.spines.right": False})
C1, C2, C3, GRAY, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#898781", "#0b0b0b", "#52514e"

# ------------------------------------------------------------- figure 1: factorisation types (exact counts)
ks = np.array([3, 5, 7, 9, 11])
q = 2.0 ** ks
n14 = (q - 2) / 2
n1112 = (q - 2) / 6
n23 = (q + 1) / 3
tot = q - 1
fig, ax = plt.subplots(figsize=(6.6, 3.6))
x = np.arange(len(ks))
w = 0.26
b1 = ax.bar(x - w, n14 / tot, w, color=C1, label=r"type $(1,4)$, $4$-cycle: $q/2$ sums $S(\lambda_1)=+4q^6$, $W=2q^7$   [$(q-2)/2$ values of $c$]")
b2 = ax.bar(x, n1112 / tot, w, color=C2, label=r"type $(1,1,1,2)$: all $S(\lambda_1)=\pm4q^6$, signs cancel, $W=0$   [$(q-2)/6$]")
b3 = ax.bar(x + w, n23 / tot, w, color=C3, label=r"type $(2,3)$: all $S(\lambda_1)=\pm4q^6$, signs cancel, $W=0$   [$(q+1)/3$]")
ax.axhline(0.5, color=GRAY, lw=0.8, ls="--")
ax.annotate(r"$\to 1/2$", xy=(len(ks) - 0.55, 0.5), xytext=(len(ks) - 0.55, 0.53), fontsize=9, color=MUTED)
ax.axhline(1 / 3, color=GRAY, lw=0.8, ls="--")
ax.annotate(r"$\to 1/3$", xy=(len(ks) - 0.55, 1 / 3), xytext=(len(ks) - 0.55, 0.36), fontsize=9, color=MUTED)
ax.axhline(1 / 6, color=GRAY, lw=0.8, ls="--")
ax.annotate(r"$\to 1/6$", xy=(len(ks) - 0.55, 1 / 6), xytext=(len(ks) - 0.55, 0.195), fontsize=9, color=MUTED)
for bars in (b1, b2, b3):
    for b in bars:
        ax.annotate(f"{b.get_height():.3f}", xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                    xytext=(0, 2), textcoords="offset points", ha="center", fontsize=6.5, color=INK)
ax.set_xticks(x)
ax.set_xticklabels([f"$q=2^{{{k}}}$" for k in ks])
ax.set_ylabel(r"fraction of $c\in\mathbb{F}_q^{\times}$")
ax.set_ylim(0, 0.84)
ax.set_title(r"Factorisation type of $v^5+v+c$ over $\mathbb{F}_q$, $k$ odd, and the mixed term $W$", fontsize=10)
ax.legend(fontsize=7.2, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 1.0), handlelength=1.2)
ax.grid(axis="y", color="#e5e7eb", lw=0.6)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_gold_cycle.pdf")); fig.savefig(os.path.join(OUT, "fig_gold_cycle.png"), dpi=160)
plt.close(fig)

# ------------------------------------------------------------- figure 2: singular planes (schematic)
fig, ax = plt.subplots(figsize=(6.6, 3.9))
ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")
# ambient and X
ax.add_patch(FancyBboxPatch((0.3, 0.4), 9.4, 5.3, boxstyle="round,pad=0.02", fc="none", ec=GRAY, lw=1.0))
ax.text(0.5, 5.45, r"$\mathbb{A}^{12}$,  coordinates $z_0,\dots,z_{11}$", fontsize=9, color=MUTED)
ax.add_patch(FancyBboxPatch((0.9, 0.8), 8.2, 4.3, boxstyle="round,pad=0.02", fc="#f3f4f6", ec=INK, lw=1.0))
ax.text(1.1, 4.85, r"$X=V(e_1,e_3,e_5)$: cone, $\dim 9$, complete intersection, normal", fontsize=9, color=INK)
ax.text(1.1, 4.52, r"smooth where $\#\{z_j\}\geq 3$ (Vandermonde rows $1,z_j^2,z_j^4$ of rank 3)", fontsize=8.5, color=MUTED)
ax.text(1.1, 4.12, r"Lucas: $|A|$ even $\Rightarrow$ $e_1=e_3=e_5\equiv 0$ on $P_A$  (each term has a factor $\binom{\mathrm{even}}{\mathrm{odd}}$)",
        fontsize=8.5, color=INK)
ax.text(1.1, 3.78, r"$|A|$ odd $\Rightarrow$ $e_1=u+v$, so $P_A\cap X=L$;   $\mathrm{Sing}(X)=\bigcup\,P_A$ over even $|A|$: $\dim 2$, codim 7",
        fontsize=8.5, color=INK)
ax.text(1.1, 3.45, r"$1023=66+495+462$ planes ($A$ and its complement give the same plane)", fontsize=8.5, color=MUTED)
# three strata as stacks of parallelograms (y from 1.85 to 2.85, offsets up to +0.3)
strata = [(r"$|A|=2$:  $\binom{12}{2}=66$ planes", C1, 1.3), (r"$|A|=4$:  $\binom{12}{4}=495$ planes", C2, 3.9),
          (r"$|A|=6$:  $\frac{1}{2}\binom{12}{6}=462$ planes", C3, 6.5)]
for label, col, x0 in strata:
    for d in range(3):
        off = 0.15 * d
        poly = Polygon([(x0 + off, 1.85 + off), (x0 + 1.4 + off, 1.85 + off), (x0 + 1.8 + off, 2.75 + off), (x0 + 0.4 + off, 2.75 + off)],
                       closed=True, fc=col, ec="white", lw=1.0, alpha=0.55 if d < 2 else 0.9)
        ax.add_patch(poly)
    ax.text(x0 + 0.05, 1.42, label, fontsize=8, color=INK)
# diagonal line through all planes
ax.plot([1.1, 8.9], [2.3, 2.3], color=INK, lw=1.8)
ax.text(1.1, 1.08, r"thick line: the diagonal $L=\{(u,\dots,u)\}$, contained in every $P_A=\{z_j=u\ (j\in A),\ z_j=v\ (j\notin A)\}$",
        fontsize=7.2, color=INK)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_singular_planes.pdf")); fig.savefig(os.path.join(OUT, "fig_singular_planes.png"), dpi=160)
plt.close(fig)

# ------------------------------------------------------------- figure 3: spectral folding (Section 8.3)
# (a) schematic of the restriction map G^ -> (G^2)^: fibres = cosets of H^perp; the fibres over the characters of G^2 of
#     order <= 2 add their coefficients in phase, the fibres over the order-4 characters cancel.
# (b) d = 6, q = 8 and 16 (fourier_spectrum_folding.py, Table 5): share of Var_{G^2} and of Var_G carried by the
#     fibres over characters of order <= 2 and of order 4.
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.0, 3.3), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axA
ax.set_xlim(0, 10); ax.set_ylim(0, 6.2); ax.axis("off")
ax.text(0.2, 5.85, r"$\hat{\mathcal{G}}$ (coefficients $c_\chi$)", fontsize=8.5, color=INK)
ax.text(7.0, 5.85, r"$\hat{\mathcal{G}}{}^2$ (fibre sums $f_\psi$)", fontsize=8.5, color=INK)
rng = np.random.RandomState(3)
rows = [("coherent", C1, 4.9, [1, 1, 1, 1, 1, 1], r"$\psi$ of order $\leq 2$:  $f_\psi\approx q^4$"),
        ("coherent", C1, 3.7, [1, 1, 1, 1, 1, 1], r"$f_\psi\approx q^4$"),
        ("incoherent", GRAY, 2.3, [1, -1, 1, -1, -1, 1], r"$\psi$ of order $4$:  $f_\psi\approx 0$"),
        ("incoherent", GRAY, 1.1, [-1, 1, -1, 1, 1, -1], r"$f_\psi\approx 0$")]
for kind, col, y, signs, lab in rows:
    for i, s in enumerate(signs):
        x = 0.5 + 0.55 * i
        h = 0.42 * s
        ax.add_patch(plt.Rectangle((x, y), 0.4, h, fc=col, ec="white", lw=0.6, alpha=0.9 if kind == "coherent" else 0.75))
    ax.plot([0.45, 3.85], [y, y], color=GRAY, lw=0.6)
    ax.annotate("", xy=(6.6, y + 0.1), xytext=(4.2, y + 0.1), arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
    ax.text(5.4, y + 0.22, r"$\sum_{\chi\in\mathrm{fibre}}$", fontsize=7.5, ha="center", color=MUTED)
    tot = sum(signs)
    ax.add_patch(plt.Rectangle((6.8, y), 0.5, 0.42 * tot / 6 * 1.6 if tot else 0.02, fc=col, ec="white", lw=0.6))
    ax.text(7.45, y + 0.08, lab, fontsize=7.2, color=INK)
ax.text(0.5, 0.25, r"fibre $=$ coset $\chi H^\perp$, $|H^\perp|=|\mathcal{G}[2]|=q^{3}$", fontsize=7.6, color=MUTED)
ax.text(0.2, 0.2 + 5.3 - 5.3, "", fontsize=1)
ax.set_title(r"(a) folding $\hat{\mathcal{G}}\to\hat{\mathcal{G}}{}^2$", fontsize=9, loc="left")

ax = axB
q8 = {"le2_H": 0.077 + 0.896, "le2_G": 0.036 + 0.228, "o4_H": 0.027, "o4_G": 0.736}
q16 = {"le2_H": 0.062 + 0.935, "le2_G": 0.008 + 0.138, "o4_H": 0.002, "o4_G": 0.854}
labels = [r"order $\leq 2$" + "\n" + r"$q=8$", r"order $\leq 2$" + "\n" + r"$q=16$", r"order $4$" + "\n" + r"$q=8$", r"order $4$" + "\n" + r"$q=16$"]
shareH = [q8["le2_H"], q16["le2_H"], q8["o4_H"], q16["o4_H"]]
shareG = [q8["le2_G"], q16["le2_G"], q8["o4_G"], q16["o4_G"]]
x = np.arange(4); w = 0.36
bH = ax.bar(x - w / 2, shareH, w, color=C1, label=r"share of $\mathrm{Var}_{\mathcal{G}^2}$")
bG = ax.bar(x + w / 2, shareG, w, color=C2, label=r"share of $\mathrm{Var}_{\mathcal{G}}$")
for bars in (bH, bG):
    for b in bars:
        ax.annotate(f"{b.get_height():.3f}", xy=(b.get_x() + b.get_width() / 2, b.get_height()), xytext=(0, 2),
                    textcoords="offset points", ha="center", fontsize=6.3, color=INK)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7.2)
ax.set_ylim(0, 1.15); ax.set_ylabel("share", fontsize=8)
ax.legend(fontsize=7, frameon=False, loc="upper right")
ax.grid(axis="y", color="#e5e7eb", lw=0.6)
ax.set_title(r"(b) $d=6$: energy by fibre type", fontsize=9, loc="left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_spectral_folding.pdf")); fig.savefig(os.path.join(OUT, "fig_spectral_folding.png"), dpi=160)
plt.close(fig)

# ------------------------------------------------------------- figure 4: general even d (Section 9.2)
# exact data: Table 6 / general_d_scan.py, legendre_q16_hayes.py, Theorem 3.x (d = 4) and Table 4 (d = 6).
ratio = {4: [(2, 4 / 3), (4, 16 / 3), (8, 64 / 3), (16, 256 / 3), (32, 1024 / 3)],
         6: [(2, 2.7028), (4, 2.4354), (8, 16.6571), (16, 51.8240), (32, 261.5810)],
         8: [(2, 0.3265), (4, 2.2544), (8, 2.2755)],
         10: [(2, 1.0993), (4, 1.3608)]}
# |phi(1)| / q^{3d/4}: d=4: (q^4-q^3)/q^3 = q-1; d=6: |phi_0|/q^{4.5} with phi_0 = 32, -768, 35840, -983040, 34537472;
# d=8: 40, 6336, 491008 / q^6; d=10: 48, 19848 / q^7.5
thr = {4: [(q, q - 1) for q in (2, 4, 8, 16, 32)],
       6: [(2, 32 / 2 ** 4.5), (4, 768 / 4 ** 4.5), (8, 35840 / 8 ** 4.5), (16, 983040 / 16 ** 4.5), (32, 34537472 / 32 ** 4.5)],
       8: [(2, 40 / 2 ** 6), (4, 6336 / 4 ** 6), (8, 491008 / 8 ** 6)],
       10: [(2, 48 / 2 ** 7.5), (4, 19848 / 4 ** 7.5)]}
cols = {4: C2, 6: C1, 8: C3, 10: "#b26dd6"}
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.0, 3.3))
ax = axA
qq = np.array([2, 4, 8, 16, 32, 64])
ax.plot(qq, qq ** 2 / 3, ls="--", color=C2, lw=0.8, alpha=0.6)
ax.plot(qq, qq ** 2 / 4, ls="--", color=C1, lw=0.8, alpha=0.6)
for d, pts in ratio.items():
    xs, ys = zip(*pts)
    ax.plot(xs, ys, marker="o", ms=4, lw=1.4, color=cols[d], label=rf"$d={d}$")
ax.set_xscale("log", base=2); ax.set_yscale("log", base=2)
ax.set_xlabel(r"$q=2^k$", fontsize=8.5); ax.set_ylabel(r"$\mathrm{Var}_{\mathcal{G}^2}/\mathrm{Var}_{\mathcal{G}}$", fontsize=8.5)
ax.text(36, 36 ** 2 / 3 * 1.15, r"$q^2/3$", fontsize=7.5, color=C2)
ax.text(36, 36 ** 2 / 4 / 1.9, r"$q^2/4$", fontsize=7.5, color=C1)
ax.text(9, 2.28 * 1.35, r"$\approx 2.28$", fontsize=7.5, color=C3)
ax.legend(fontsize=7.5, frameon=False, loc="upper left")
ax.grid(color="#e5e7eb", lw=0.6, which="major")
ax.set_title(r"(a) amplification against $q$", fontsize=9, loc="left")

ax = axB
w = 0.19
for i, d in enumerate((4, 6, 8, 10)):
    for q, v in thr[d]:
        pos = np.log2(q) + (i - 1.5) * w
        ax.bar(pos, v, w, color=cols[d], label=rf"$d={d}$" if q == thr[d][0][0] else None)
ax.axhline(1.0, color=INK, lw=0.9, ls="--")
ax.text(1.05, 1.08, "threshold $1$", fontsize=7.5, color=INK)
ax.set_yscale("log", base=2)
ax.set_xticks([1, 2, 3, 4, 5]); ax.set_xticklabels([r"$2$", r"$4$", r"$8$", r"$16$", r"$32$"])
ax.set_xlabel(r"$q$", fontsize=8.5); ax.set_ylabel(r"$|\varphi(1)|/q^{3d/4}$", fontsize=8.5)
ax.legend(fontsize=7.5, frameon=False, loc="upper left", ncol=2)
ax.grid(axis="y", color="#e5e7eb", lw=0.6)
ax.set_title(r"(b) identity class against the threshold", fontsize=9, loc="left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_even_d_bifurcation.pdf")); fig.savefig(os.path.join(OUT, "fig_even_d_bifurcation.png"), dpi=160)
plt.close(fig)

# ------------------------------------------------------------- figure 5: S_6 classes and the signs of W' (Section 6.7, Table 3)
# asymptotic density of each class = leading coefficient of the exact count / q^2 (Table 3); W'/q^7 by class.
classes = [("$1^6$", 1 / 720, 4), ("$2\\,1^4$", 1 / 48, 0), ("$2^2 1^2$", 1 / 16, 0), ("$2^3$", 1 / 48, 4),
           ("$3\\,1^3$", 1 / 18, 4), ("$3^2$", 1 / 18, 4), ("$3\\,2\\,1$", 1 / 6, 0), ("$4\\,1^2$", 1 / 8, 2),
           ("$4\\,2$", 1 / 8, -2), ("$5\\,1$", 1 / 5, -1), ("$6$", 1 / 6, 4)]
assert abs(sum(c[1] for c in classes) - 1) < 1e-12 and abs(sum(c[1] * c[2] for c in classes) - 1) < 1e-12
valcol = {4: C1, 2: C3, 0: GRAY, -2: C2, -1: "#c0392b"}
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.0, 3.4), gridspec_kw={"width_ratios": [1.05, 1]})
ax = axA
x = np.arange(len(classes))
dens = [c[1] for c in classes]
bars = ax.bar(x, dens, 0.7, color=[valcol[c[2]] for c in classes], edgecolor="white", lw=0.5)
for b, c in zip(bars, classes):
    ax.annotate(f"1/{round(1/c[1])}", xy=(b.get_x() + b.get_width() / 2, b.get_height()), xytext=(0, 2),
                textcoords="offset points", ha="center", fontsize=6.3, color=INK)
ax.set_xticks(x); ax.set_xticklabels([c[0] for c in classes], fontsize=7.5, rotation=0)
ax.set_ylabel(r"density of the class ($\times q^2$ pairs)", fontsize=8)
ax.set_ylim(0, 0.25)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=valcol[v], label=rf"$W'/q^7={v}$") for v in (4, 2, 0, -2, -1)], fontsize=6.8, frameon=False,
          loc="upper left", ncol=2, title=r"value of $W'$", title_fontsize=7)
ax.grid(axis="y", color="#e5e7eb", lw=0.6)
ax.set_title(r"(a) Frobenius classes in $S_6\cong Sp_4(\mathbb{F}_2)$", fontsize=9, loc="left")
ax = axB
contrib = [c[1] * c[2] for c in classes]
cum = np.cumsum(contrib)
ax.bar(x, contrib, 0.7, color=[valcol[c[2]] for c in classes], edgecolor="white", lw=0.5)
ax.plot(x, cum, color=INK, lw=1.2, marker="o", ms=3, label="running total")
ax.axhline(1.0, color=INK, lw=0.8, ls="--")
ax.axhline(0, color=GRAY, lw=0.6)
ax.text(0.1, 1.03, r"$\sum=1$: leading term $q^9$ of $\sum_{\rm mixed}W'=q^9-3q^8+q^7$", fontsize=7.2, color=INK)
ax.annotate(r"$5$-cycles: $-\frac{1}{5}$", xy=(9, -0.2), xytext=(5.6, -0.42), fontsize=7.5, color="#c0392b",
            arrowprops=dict(arrowstyle="->", color="#c0392b", lw=0.8))
ax.annotate(r"$6$-cycles: $+\frac{2}{3}$", xy=(10, 0.667), xytext=(6.4, 0.78), fontsize=7.5, color=C1,
            arrowprops=dict(arrowstyle="->", color=C1, lw=0.8))
ax.set_xticks(x); ax.set_xticklabels([c[0] for c in classes], fontsize=7.5)
ax.set_ylabel(r"density $\times\ W'/q^7$", fontsize=8)
ax.set_ylim(-0.5, 1.15)
ax.legend(fontsize=7, frameon=False, loc="center left")
ax.grid(axis="y", color="#e5e7eb", lw=0.6)
ax.set_title(r"(b) contributions to $\sum_{\rm mixed} W'/q^9$", fontsize=9, loc="left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_S6_classes.pdf")); fig.savefig(os.path.join(OUT, "fig_S6_classes.png"), dpi=160)
plt.close(fig)

# ------------------------------------------------------------- figure 6: spectral bounds on m_k (Section 7.2, Props 7.8 and 7.10)
mk = {1: 13, 2: -1, 3: -35, 4: 47, 5: 733}
ks = np.array(sorted(mk)); q = 2.0 ** ks
norm = np.array([abs(mk[k]) for k in ks]) / q ** 1.5
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.0, 3.3))
ax = axA
ax.plot(ks, norm, marker="o", ms=5, lw=1.2, color=C1, label=r"$|m_k|/q^{3/2}$")
for k, v in zip(ks, norm):
    ax.annotate(f"{v:.2f}", xy=(k, v), xytext=(0, 5), textcoords="offset points", ha="center", fontsize=7, color=INK)
for b, ls in ((4, ":"), (5, "--")):
    ax.axhline(b, color=C2, lw=0.9, ls=ls)
    ax.text(1.05, b + 0.1, rf"pure weight $9$, rank $b_9={b}$: $|m_k|\leq {b}\,q^{{3/2}}$", fontsize=7, color=C2)
ax.axhspan(4.05, 5, color=C2, alpha=0.08)
ax.text(2.55, 2.5, "rank $5$ excluded by Newton:\n" + r"$e_5=2281\neq 2^{15/2}=181.0$", fontsize=7.2, color="#c0392b", ha="center")
ax.set_xticks(ks); ax.set_xlabel(r"$k$  ($q=2^k$)", fontsize=8.5); ax.set_ylabel(r"$|m_k|/q^{3/2}$", fontsize=8.5)
ax.set_ylim(0, 6)
ax.legend(fontsize=7.5, frameon=False, loc="upper center")
ax.grid(color="#e5e7eb", lw=0.6)
ax.set_title(r"(a) the five values against the pure weight-$9$ bound", fontsize=9, loc="left")
ax = axB
v2 = [3 * k for k in ks]
ax.plot(ks, v2, marker="s", ms=5, lw=0, color=C1, label=r"$v_2(Y_k)=3k$ (exact, $m_k$ odd)")
ax.plot(ks, 3 * ks, color=C1, lw=0.9, ls="-")
kk = np.linspace(1, 5, 50)
ax.plot(kk, 3.5 * kk, color="#c0392b", lw=1.0, ls="--", label=r"$v_2\geq 7k/2$: any part of weight $\geq 7$")
ax.fill_between(kk, 3.5 * kk, 3.5 * kk + 6, color="#c0392b", alpha=0.08)
ax.fill_between(kk, 3 * kk - 3, 3 * kk, color=C3, alpha=0.08)
ax.text(1.1, 15.3, "weight $\\geq 7$ alone: impossible\n($Y_2=-64$ is not $2^7\\times$ algebraic integer)", fontsize=7, color="#c0392b")
ax.text(3.05, 5.2, "a part of weight $\\leq 6$\nis forced", fontsize=7, color=C3)
ax.set_xticks(ks); ax.set_xlabel(r"$k$", fontsize=8.5); ax.set_ylabel(r"$2$-adic valuation", fontsize=8.5)
ax.set_ylim(0, 22)
ax.legend(fontsize=6.8, frameon=False, loc="upper left")
ax.grid(color="#e5e7eb", lw=0.6)
ax.set_title(r"(b) $2$-adic constraints (supersingular shape)", fontsize=9, loc="left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_spectral_bounds.pdf")); fig.savefig(os.path.join(OUT, "fig_spectral_bounds.png"), dpi=160)
plt.close(fig)
print("wrote fig_gold_cycle, fig_singular_planes, fig_spectral_folding, fig_even_d_bifurcation, fig_S6_classes, fig_spectral_bounds (pdf/png)")
