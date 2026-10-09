# 第十三部分：终稿——Figure 3（S₆ 类与 W′ 符号）、Figure 6（m_k 的谱约束）、摘要/引言一致性、参考文献独立分页

日期：2026-10-08。代码：`code/generate_figures.py`（现生成 6 幅矢量图：fig_gold_cycle、fig_singular_planes、fig_S6_classes、fig_spectral_bounds、fig_spectral_folding、fig_even_d_bifurcation；加上 `paper/make_figures.py` 的 fig_profile、fig_scaling 共 8 幅）。论文 33 页，0 errors、0 warnings（含 hyperref）、0 overfull/underfull。

## 新图与编号（按论文最终顺序）
| 图 | 文件 | 位置 | 内容 |
|---|---|---|---|
| 1 | fig_profile | §4 | Legendre 偏差剖面（q=16） |
| 2 | fig_gold_cycle | §6.6 | v⁵+v+c 的 S₅ 分解型分布 |
| **3** | **fig_S6_classes** | §6.7 | (a) S₆ ≅ Sp₄(F₂) 的 11 个共轭类的渐近密度（1/720, 1/48, 1/16, 1/48, 1/18, 1/18, 1/6, 1/8, 1/8, 1/5, 1/6，按 W′/q⁷ 着色）；(b) 密度 × W′/q⁷ 与累计和：5-轮换 −1/5，6-轮换 +2/3，两个 3 阶类 +4/9，2³ +1/12，恒等 +1/180，两个 4 阶类 ±1/4 相消；总和恰为 1 = Σ_mixed W′ 的首项系数 q⁹（Prop 6.18(c)：q⁹ − 3q⁸ + q⁷） |
| 4 | fig_scaling | §7.1 | 归一化方差对 q |
| 5 | fig_singular_planes | §7.2 | X 的 1023 张奇异平面 |
| **6** | **fig_spectral_bounds** | §7.2 | (a) \|m_k\|/q^{3/2} = 4.60, 0.13, 1.55, 0.73, 4.05 对纯权重 9、秩 b₉ 的界 b₉；b₉ ≥ 5 被迫、b₉ = 5 被 Newton 恒等式 e₅ = 2281 ≠ 2^{15/2} 排除；(b) 2-adic 图：v₂(Y_k) = 3k（m_k 奇）与"权重 ≥ 7 的分量赋值 ≥ 7k/2"虚线，数据严格在其下 ⟹ 权重 ≤ 6 分量被迫；Y₂ = −64 排除任何单一权重 ≥ 7 |
| 7 | fig_spectral_folding | §8.3 | 谱折叠示意与能量份额 |
| 8 | fig_even_d_bifurcation | §9.2 | d = 4..10 的放大分岔与阈值 |

## 一致性微调
- 摘要："The conjectured closed forms in the remaining cases (A₆ = q⁹ − (q−1)q⁵ for even k, and X = (−1)^{k+1}q⁶ + q⁵ for all k) are verified exactly for q ≤ 128 …"——不再让读者误以为奇数 k 的 A₆ 未证。
- 引言 Main results 第三条：明确 "this *proves* A₆ = q⁹ + (q−1)q⁵ for all odd k (Corollary 6.13)"，并把剩余猜想限定为 A₆（偶 k）与 X（所有 k）。
- 参考文献前 `\clearpage`：Table 7 完整占第 32 页，第 33 页为独立的参考文献与作者信息页。

## 绘图注意事项（供复现）
- matplotlib mathtext：无 `\widehat`、`\tfrac`、`\ge`；`\mathcal{G}` 必须带花括号。
- 全部图用 Computer Modern 数学字体（`mathtext.fontset = cm`），serif 正文字体，矢量 PDF 加 160 dpi 预览 PNG。
