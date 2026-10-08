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
print("wrote fig_gold_cycle.pdf/png and fig_singular_planes.pdf/png")
