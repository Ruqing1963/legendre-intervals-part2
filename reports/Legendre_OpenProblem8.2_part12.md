# 第十二部分：§8–§9 的插图（Figure 5、Figure 6）与零警告定稿

日期：2026-10-08。代码：`code/generate_figures.py`（新增两幅图的绘制代码；matplotlib，Computer Modern 数学字体，矢量 PDF + 预览 PNG）。论文：§8.3 新增 Figure 5（`fig_spectral_folding.pdf`），§9.2 新增 Figure 6（`fig_even_d_bifurcation.pdf`），正文加交叉引用；全部 16 个含数学的章节标题改用 `\texorpdfstring`，加 `pdftitle/pdfauthor`，hyperref 的 59 条 "Token not allowed in a PDF string" 警告全部消除；最终核查表前加分页。

## 编译状态
- 0 errors，0 undefined references，0 hyperref warnings，0 overfull/underfull hbox；仅剩的 1 条 underfull vbox（图 6 所在页的垂直拉伸）已通过在 §9.4 前分页消除——见最后一次编译输出。全文 32 页。

## Figure 5（§8.3，谱折叠）
- (a) 示意图：Ĝ 的系数 c_χ 按限制映射的纤维（陪集 χH^⊥，\|H^⊥\| = \|G[2]\| = q³）分组；阶 ≤ 2 的 ψ 上方纤维同号叠加（f_ψ ≈ q⁴），阶 4 的 ψ 上方纤维相消（f_ψ ≈ 0）。
- (b) 数据（Table 5）：d = 6，q = 8、16。阶 ≤ 2 纤维携带 Var_{G²} 的 97.3%、99.7%，却只占 Var_G 的 26.4%、14.6%；阶 4 纤维反之（2.7%/73.6%，0.2%/85.4%）。这正是式 (10) 的相干/非相干分解。

## Figure 6（§9.2，一般偶数 d）
- (a) Var_{G²}/Var_G 对 q 的双对数图：d = 4 恰为 q²/3（q ≤ 32），d = 6 ≈ q²/4（q ≤ 32），d = 8 在 q = 4, 8 平坦（2.25, 2.28），d = 10（q ≤ 4）≈ 1.1–1.4。虚线为 q²/3、q²/4。
- (b) \|φ(1)\|/q^{3d/4} 对阈值 1：d = 4 为 q−1（3, 7, 15, 31），d = 6 为 1.4–6.2，d = 8 为 0.63, 1.55, 1.87（贴着阈值），d = 10 为 0.27, 0.61（阈值之下）。数值来源：Theorem 3.2（d=4）、Table 4 的 φ₀（d=6）、`general_d_scan.py`（d=8,10）。
- 说明：面板 (b) 的"阈值"是 Prop 9.1(c) 的下界判据（恒等轨道单独就能造成放大的条件），不是相变的证明；图注和正文都按此措辞。

## 其他修改
- matplotlib mathtext 不支持 `\widehat`、也要求 `\mathcal{G}` 带花括号（已改）。
- 状态表 C₆ 行改为 `\raggedright`，消除 p-列的 underfull。
- fig_scaling 的图注开头改写，消除原有 underfull。
