"""Figures for paper_part2.tex (vector PDF + PNG preview). All numbers are the exact values computed in the paper.
   fig_profile.pdf : Legendre deviations phi = Psi - q^{d+1} on G^2 for q = 16, d = 4 and d = 6 (Theorem 4.2).
   fig_scaling.pdf : Var_{G^2}/q^{d+1} and Var_G/q^{d+1} against q (log-log), d = 4 and d = 6, with q^2 and q^4 references."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.abspath(__file__))     # figures are written next to this script (paper/)

plt.rcParams.update({"font.size": 10, "font.family": "serif", "mathtext.fontset": "cm",
                     "axes.spines.top": False, "axes.spines.right": False})
# validated categorical slots 1-3 of the reference palette (blue, orange, aqua); 4th line is neutral gray
C_D4, C_D6, C_REF, C_REF2 = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
INK, MUTED = "#0b0b0b", "#52514e"

# ------------------------------------------------------------------ Figure 1: profiles, q = 16
q = 16
# d = 4: 16 Legendre classes (0,c,0): identity -(q^4-q^3), others +q^3 (k even)
d4_levels = [(-(q**4 - q**3), 1, r"$c=0$: $-(q^4-q^3)$"), (q**3, 15, r"$c\neq0$: $+q^3$")]
# d = 6: 256 classes: c2=0 -> -q^5+q^4 (16), Tr(c4/c2^2)=0 -> +q^2 m_4 (120), Tr=1 -> -q^2 m_4 (120), m_4 = 47
d6_levels = [(-983040, 16, r"$c_2=0$: $-q^5+q^4$"), (12032, 120, r"$\mathrm{Tr}(c_4/c_2^2)=0$: $+q^2m_4$"),
             (-12032, 120, r"$\mathrm{Tr}(c_4/c_2^2)=1$: $-q^2m_4$")]

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.4))
rng = np.random.default_rng(3)
for ax, levels, title, color, lin in [(axes[0], d4_levels, r"$d=4$: $\varphi$ on the $q$ Legendre classes", C_D4, 1e3),
                                      (axes[1], d6_levels, r"$d=6$: $\varphi$ on the $q^2$ Legendre classes", C_D6, 3e3)]:
    ax.set_yscale("symlog", linthresh=lin, linscale=0.6)
    ax.axhline(0, color=MUTED, lw=0.8, ls="--")
    for i, (lvl, cnt, lab) in enumerate(levels):
        x = rng.normal(0, 0.06, cnt)
        ax.axhline(lvl, color=color, lw=1.2, alpha=0.9)
        ax.scatter(x, [lvl] * cnt, s=10, color=INK, alpha=0.55, zorder=3, linewidths=0)
        ax.annotate(f"{lab}  ({cnt})", xy=(0.27, lvl), fontsize=8.5, color=INK, va="center")
    ax.set_xlim(-0.3, 1.3)
    ax.set_xticks([])
    ax.set_title(title, fontsize=10)
    ax.grid(axis="y", color="#e5e7eb", lw=0.6)
axes[0].set_yticks([-60000, -10000, 0, 4096]); axes[0].set_yticklabels([r"$-6\times10^4$", r"$-10^4$", "0", r"$+4096$"])
axes[1].set_yticks([-983040, -100000, -12032, 0, 12032]); axes[1].set_yticklabels([r"$-9.8\times10^5$", r"$-10^5$", r"$-12032$", "0", r"$+12032$"])
axes[0].set_ylabel(r"$\varphi=\Psi-q^{5}$  (symlog),  $q=16$")
axes[1].set_ylabel(r"$\varphi=\Psi-q^{7}$  (symlog),  $q=16$")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_profile.pdf")); fig.savefig(os.path.join(OUT, "fig_profile.png"), dpi=160)
plt.close(fig)

# ------------------------------------------------------------------ Figure 2: scaling
q4 = np.array([2, 4, 8, 16, 32, 64, 128])
varH4 = q4 * (q4 - 1.0)                       # Theorem 3.4
varG4 = 3 - 3.0 / q4                          # Theorem 3.9
q6 = np.array([2, 4, 8, 16, 32])
varH6 = np.array([15.5625, 9.011719, 79.093506, 225.505600, 1100.884360])
varG6 = np.array([5.757812, 3.700378, 4.748336, 4.351375, 4.208574])
qs = np.logspace(np.log10(2), np.log10(128), 100)

fig, ax = plt.subplots(figsize=(7.0, 4.0))
ax.plot(qs, qs**2, color=C_REF, lw=1.2, ls="--", label=r"reference $q^{2}$")
ax.plot(qs, qs**4 / 8**4 * 79.093506, color=C_REF2, lw=1.2, ls=":", label=r"reference $q^{4}$ (through the $d=6$, $q=8$ point)")
ax.plot(q4, varH4, "o-", color=C_D4, lw=1.4, ms=5, label=r"$d=4$: $\mathrm{Var}_{\mathcal{G}^2}/q^{5}=q(q-1)$ (Thm 3.4)")
ax.plot(q6, varH6, "s-", color=C_D6, lw=1.4, ms=5, label=r"$d=6$: $\mathrm{Var}_{\mathcal{G}^2}/q^{7}$ (exact, Thm 7.1)")
ax.plot(q4, varG4, "o", mfc="white", color=C_D4, ms=5, label=r"$d=4$: $\mathrm{Var}_{\mathcal{G}}/q^{5}=3-3/q$")
ax.plot(q6, varG6, "s", mfc="white", color=C_D6, ms=5, label=r"$d=6$: $\mathrm{Var}_{\mathcal{G}}/q^{7}$")
ax.axhline(2, color=C_D4, lw=0.8, ls="-.", alpha=0.6)
ax.axhline(4, color=C_D6, lw=0.8, ls="-.", alpha=0.6)
ax.annotate("Keating–Rudnick limit $d-2=2$", xy=(45, 2), xytext=(45, 1.35), fontsize=8, color=MUTED)
ax.annotate("Keating–Rudnick limit $d-2=4$", xy=(45, 4), xytext=(45, 5.6), fontsize=8, color=MUTED)
ax.set_xscale("log", base=2); ax.set_yscale("log", base=10)
ax.set_xticks(q4); ax.set_xticklabels([str(v) for v in q4])
ax.set_xlabel(r"$q=2^k$")
ax.set_ylabel(r"normalised variance $\mathrm{Var}/q^{d+1}$")
ax.grid(True, which="major", color="#e5e7eb", lw=0.6)
ax.legend(fontsize=7.8, frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_scaling.pdf")); fig.savefig(os.path.join(OUT, "fig_scaling.png"), dpi=160)
plt.close(fig)
print("wrote fig_profile.pdf/png, fig_scaling.pdf/png")
