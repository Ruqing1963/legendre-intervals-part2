# 第五部分：常数 $Y$ 的代数几何定义——坐标、簇、维数与 $\ell$-进上同调

日期：2026-10-08。承接第三、四部分及论文 `paper_part2.tex` §4、§7。代码：`code/verify_Y_direct.py`（已运行，两条独立路线均断言通过）、`code/verify_Y_sage.sage`（SageMath 版本；本机无 Sage，未运行）、`code/split_form_counts.py`（非扭形式的点数，见其中 caveat）。

**本部分的目的**是把论文中"九维几何对象上的 Frobenius 迹"这一说法变成无歧义的定义。核心结论先行：

- $Y_k$ 是一个**固定的、定义在 $\mathbb F_2$ 上的** 9 维仿射簇 $X=V(e_1,e_3,e_5)\subset\mathbb A^{12}$ 上，Artin–Schreier 层 $\mathcal L=e_4^*\mathcal L_{\psi_0}$ 关于**扭 Frobenius** $\sigma\circ F^k$（$\sigma$ 为坐标的 12-循环置换）的 Lefschetz 迹。"九维"指 $X$ 的 Krull 维数（$X$ 是 $\mathbb A^{12}$ 中余维 3 的完全交），不是纤维维数，也不是 Betti 数。
- $X$ 几何整、正规、Cohen–Macaulay，是锥；奇异轨迹是 1023 张过原点的平面之并（维数 2，余维 7），其方程完全显式。
- $Y_k$ 不是任何单个簇的 $\mathbb F_{2^k}$-点计数：定义 $U_k$ 的"相对迹"条件依赖于中间域 $\mathbb F_q$，$U_k$ 是 $X$ 的 $\sigma$-扭形式的 $\mathbb F_q$-点。正因为是扭迹，$(Y_k)_k$ 仍是一个有理函数（扭 $L$-函数）的对数导数系数。
- Deligne 权：$H^j_c$ 的 Frobenius 特征根权 $\le j$；归一化 $2^{3k}=q^3$ 对应 Tate 扭 $(3)$，即权 6；若上同调集中在 $j=9$ 且纯，则 $m_k=Y_k/q^3$ 为权 3，$|m_k|\le b_9\,q^{3/2}$。数据要求 $b_9\ge5$，且"纯权 3 且 $b_9=5$"被 $|e_5|=2^{15/2}\ne2281$ 排除。

---

## 1. 定义与坐标化

### 1.1 算术定义（论文 §4）

**定义 1.1.** 对 $q=2^k$，令 $\operatorname{Tr}=\operatorname{Tr}_{\mathbb F_{q^{12}}/\mathbb F_q}$，$\psi_q(x)=(-1)^{\operatorname{Tr}_{\mathbb F_q/\mathbb F_2}(x)}$，
$$U_k=\{\alpha\in\mathbb F_{q^{12}}:\ \operatorname{Tr}\alpha=\operatorname{Tr}\alpha^3=\operatorname{Tr}\alpha^5=0\},\qquad
Y_k=\sum_{\alpha\in U_k}\psi_q\bigl(e_4(\alpha)\bigr),$$
其中 $e_m(\alpha)$ 是 $\alpha$ 的 12 个共轭 $\alpha,\alpha^q,\dots,\alpha^{q^{11}}$ 的第 $m$ 个初等对称多项式（$\in\mathbb F_q$）。

### 1.2 坐标化：分裂形式

**定义 1.2（分裂形式）.** 在 $\mathbb A^{12}_{\mathbb F_2}=\operatorname{Spec}\mathbb F_2[z_0,\dots,z_{11}]$ 上令
$$I=\langle e_1(z),\,e_3(z),\,e_5(z)\rangle,\qquad X:=V(I)=\{z:\ e_1=e_3=e_5=0\},$$
$e_m(z)=\sum_{j_1<\dots<j_m}z_{j_1}\cdots z_{j_m}$。

**命题 1.3.** $I=\langle p_1,p_3,p_5\rangle$，其中 $p_m=\sum_jz_j^m$。

*证明.* 特征 2 下的 Newton 恒等式：$p_1=e_1$；$p_3=e_1^3+e_1e_2+e_3$；$p_5=e_5+e_2e_3+e_1(\cdots)$（由 $5e_5=e_4p_1-e_3p_2+e_2p_3-e_1p_4+p_5$ 及 $p_2=e_1^2$、$p_4=e_1^4$）。故 $e_3\equiv p_3\pmod{e_1}$，$e_5\equiv p_5\pmod{e_1,e_3}$，两理想相等。$\square$

**定义 1.4（扭 Frobenius 与扭形式）.** 令 $\sigma\in\operatorname{Aut}(\mathbb A^{12})$，$\sigma(z)_j=z_{j-1}$（下标模 12）。$\sigma$ 定义在 $\mathbb F_2$ 上，保持 $I$（$e_m$ 对称），且与绝对 Frobenius $F:z\mapsto(z_j^2)$ 交换。$X$ 关于 $\mathbb F_q$ 的 $\sigma$-扭形式 $X^{(\sigma)}$ 是由 Frobenius $F_q^{(\sigma)}:=\sigma\circ F^k$ 定义的 $\mathbb F_q$-结构（$H^1(\mathbb F_q,S_{12})$ 中由 12-循环代表的上循环；它恰是 Weil 限制 $\operatorname{Res}_{\mathbb F_{q^{12}}/\mathbb F_q}\mathbb A^1\cong\mathbb A^{12}_{\mathbb F_q}$ 的扭）。

**命题 1.5.** 映射 $\alpha\mapsto z(\alpha)=(\alpha,\alpha^q,\dots,\alpha^{q^{11}})$ 是双射 $U_k\xrightarrow{\ \sim\ }X^{(\sigma)}(\mathbb F_q)=\{z\in X(\bar{\mathbb F}_2):\sigma F^k(z)=z\}$，且 $e_4(z(\alpha))=e_4(\alpha)$。因此
$$Y_k=\sum_{z\in X^{(\sigma)}(\mathbb F_q)}\psi_q\bigl(e_4(z)\bigr).$$

*证明.* $F^k(z)_j=z_j^q$，$\sigma F^k(z)_j=z_{j-1}^q$；固定点即 $z_j=z_{j-1}^q$ 对所有 $j$，亦即 $z_j=z_0^{q^j}$ 且 $z_0^{q^{12}}=z_0$，故 $z_0\in\mathbb F_{q^{12}}$。条件 $e_1=e_3=e_5=0$ 即 $p_1=p_3=p_5=0$（命题 1.3），即 $\operatorname{Tr}\alpha=\operatorname{Tr}\alpha^3=\operatorname{Tr}\alpha^5=0$。$\square$

**注 1.6（为何不是固定簇的点计数）.** $U_k$ 中的条件 $\operatorname{Tr}_{q^{12}/q}$ 依赖于 $k$；用 $\mathbb F_q$-基把 $\mathbb F_{q^{12}}$ 写成 $\mathbb F_q^{12}$ 得到的簇 $V(p_1,p_3,p_5)\subset\mathbb A^{12}_{\mathbb F_q}$（论文 §7 的 "$V$"）正是 $X^{(\sigma)}$，它是 $X\otimes\mathbb F_q$ 的非平凡扭形式，不是某个 $\mathbb F_2$-簇的基变换。所以 $(Y_k)_k$ **不是**一个固定簇的 $\mathbb F_{2^k}$-点计数序列，而是固定簇 $X$ 上**扭迹** $\operatorname{Tr}(\sigma F^k)$ 的序列（§3）。

### 1.3 商模型与 Hayes 群

**命题 1.7.** 令 $\pi:\mathbb A^{12}\to\mathbb A^{12}/S_{12}=\operatorname{Spec}\mathbb F_2[e_1,\dots,e_{12}]\cong\mathbb A^{12}$，$\Lambda:=V(e_1,e_3,e_5)\cong\mathbb A^9$，坐标 $(e_2,e_4,e_6,e_7,e_8,e_9,e_{10},e_{11},e_{12})$。则 $X=\pi^{-1}(\Lambda)$，$\pi|_X:X\to\Lambda$ 有限平坦，次数 $12!$，在 $\Lambda$ 的稠密开集上 étale。$\Lambda(\mathbb F_q)$ 与首一多项式 $F=t^{12}+e_2t^{10}+e_4t^8+\dots+e_{12}\in\mathbb F_q[t]$（$a_{11}=a_9=a_7=0$）一一对应，而
$$\#\bigl(\pi^{-1}(F)\cap X^{(\sigma)}(\mathbb F_q)\bigr)=\Lambda(F),\qquad\text{故}\qquad
Y_k=\sum_{F\in\Lambda(\mathbb F_q)}\Lambda(F)\,\psi_q(a_8(F)).$$

*证明.* $\mathbb F_2[z]$ 是 $\mathbb F_2[e]$ 上秩 $12!$ 的自由模，故 $\pi$ 有限平坦，$X=\pi^{-1}(\Lambda)$ 亦然。通有 $F$ 可分：$F'=e_7t^4+e_9t^2+e_{11}=(\sqrt{e_7}t^2+\sqrt{e_9}t+\sqrt{e_{11}})^2\ne0$。$\pi^{-1}(F)$ 中的 $\sigma F^k$-固定点是 $F$ 的根的排序 $(z_j)$ 满足 $z_j=z_0^{q^j}$，即 $\alpha=z_0\in\mathbb F_{q^{12}}$ 且 $\operatorname{charpoly}(\alpha)=F$；这样的 $\alpha$ 恰有 $\Lambda(F)$ 个（$F=P^m$ 时为 $P$ 的 $12/m$ 个根，否则为 0）。$\square$

这正是论文 §4 中 $Y=\sum_{c_2,c_4}\Psi(0,c_2,0,c_4,0)\chi(c_4)$ 的几何来源：复合 $X\to\Lambda\to\mathbb A^2$，$(e_2,e_4)$，纤维维数依次为 $0$ 与 $7$，$9=2+7$。

### 1.4 Artin–Schreier 覆盖与层

**定义 1.8.** $f:=e_4|_X$，$\mathbb F_2$-系数的 4 次齐次多项式。Artin–Schreier 覆盖
$$Z=\{(z,y)\in X\times\mathbb A^1:\ y^2+y=e_4(z)\}\subset\mathbb A^{13},\qquad \varpi:Z\to X,$$
$\mathcal L:=e_4^*\mathcal L_{\psi_0}$，$\mathcal L_{\psi_0}$ 为 $\mathbb A^1_{\mathbb F_2}$ 上对应 $\psi_0=(-1)^{(\cdot)}$ 的 Artin–Schreier 层；对 $\mathbb F_q$，$\psi_q=\psi_0\circ\operatorname{Tr}_{q/2}$，故 $\mathcal L\otimes\mathbb F_q$ 的 Frobenius 迹函数正是 $\psi_q\circ e_4$。

**命题 1.9.** (i) $\varpi$ 是 étale 的 2 次覆盖（$d(y^2+y)=dy\ne0$），$Z$ 在 $\varpi^{-1}(X_{\rm sm})$ 上非奇异。(ii) $f$ 在仿射簇 $X$ 上无极点；在射影闭包 $\overline X\subset\mathbb P^{12}$ 上，$f=e_4(z)/z_\infty^4$ 的极点因子为 $4X_\infty$，$X_\infty=\overline X\cap\{z_\infty=0\}\cong\mathbb P(X)\subset\mathbb P^{11}$（$X$ 是锥，$X_\infty$ 是锥底）。因为极点阶 $4$ 被 $p=2$ 整除，$\mathcal L$ 沿 $X_\infty$ 的 Swan 导子严格小于 4（可用 $y\mapsto y+g$ 消去 $f$ 的 2 次幂部分），这是任何 Katz–Laumon 型纯性论证必须处理的点。(iii) $\mathcal L$ 在 $X_{\bar{\mathbb F}_2}$ 上几何非平凡：限制到对角线 $L=\{(u,\dots,u)\}\subset X$ 得 $e_4=\binom{12}4u^4=u^4$，而 $u^4\notin\{g^2+g+c:g\in\bar{\mathbb F}_2[u]\}$（$g=u^2+au+b$ 给 $u^4+(a^2+1)u^2+au+\cdots$，需 $a=0$ 与 $a^2+1=0$ 矛盾），故 $Z$ 几何连通。

---

## 2. 维数、奇异轨迹、整性

**命题 2.1（维数）.** $X$ 是 $\mathbb A^{12}$ 中由正则序列 $(e_1,e_3,e_5)$（次数 $1,3,5$）定义的完全交，$\dim X=9$（Krull 维数），纯维，Cohen–Macaulay，且是以原点为顶点的仿射锥（方程齐次）。"九维"的确切含义是 $\dim X=9$；等价地，$X\to\Lambda\cong\mathbb A^9$ 有限（纤维维数 0），$X\to\mathbb A^2$（$e_2,e_4$）纤维维数 7。Betti 数 $b_j=\dim H^j_c(X_{\bar{\mathbb F}_2},\mathcal L)$ 是另一回事（§3）。

*证明.* 命题 1.7：$X$ 有限平坦地映到 $\mathbb A^9$，故纯 9 维；$\mathbb F_2[z]/(e_1,e_3,e_5)$ 是 $\mathbb F_2[e]/(e_1,e_3,e_5)\cong\mathbb F_2[e_2,e_4,e_6,\dots,e_{12}]$ 上的自由模，故 $(e_1,e_3,e_5)$ 是正则序列。$\square$

**命题 2.2（奇异轨迹）.** 记 $\Delta_2=\{z:\ \#\{z_0,\dots,z_{11}\}\le2\}$。则
$$\operatorname{Sing}(X)=X\cap\Delta_2=\bigcup_{\substack{A\subset\{0,\dots,11\}\\ |A|\ \text{偶},\ 0<|A|<12}}P_A,\qquad P_A=\{z:\ z_j=u\ (j\in A),\ z_j=v\ (j\notin A)\}\cong\mathbb A^2,$$
共 $\binom{12}2+\binom{12}4+\tfrac12\binom{12}6=66+495+462=1023$ 张过原点的平面（$A$ 与其补集给出同一平面），都含对角线 $L=\{(u,\dots,u)\}$。于是 $\dim\operatorname{Sing}(X)=2$，余维 7。

*证明.* $X$ 是完全交，由 Jacobi 判别法，$X$ 在 $z$ 光滑当且仅当 $de_1,de_3,de_5$ 在 $z$ 线性无关。$\partial e_m/\partial z_j=e_{m-1}(z\setminus z_j)$，在 $X$ 上（$e_1=e_3=e_5=0$）：
$$\partial_je_1=1,\qquad \partial_je_3=e_2+z_j^2,\qquad \partial_je_5=e_4+e_2z_j^2+z_j^4 .$$
行化简后三行为 $(1)_j,(z_j^2)_j,(z_j^4)_j$，即 $w_j=z_j^2$ 的 Vandermonde 行，秩 3 当且仅当 $\{w_j\}$ 至少 3 个不同值，当且仅当 $\{z_j\}$ 至少 3 个不同值（平方在特征 2 下单射）。反之，设 $z$ 只取值 $u,v$，$|A|=r$ 个 $v$、$s=12-r$ 个 $u$。$e_m(z)=\sum_i\binom{s}{m-i}\binom ri u^{m-i}v^i$。若 $r$ 偶（则 $s$ 偶），对奇数 $m$ 每一项含 $\binom{\text{偶}}{\text{奇}}$，由 Lucas 定理为偶，故 $e_1=e_3=e_5=0$，整张平面 $P_A\subset X$。若 $r$ 奇，$e_1=su+rv=u+v$，$e_1=0$ 迫使 $u=v$，落在 $L$ 上。$\square$

**命题 2.3（整性）.** $X$ 几何整（几何不可约且约化）、正规。

*证明.* 约化：通有纤维 étale（命题 1.7）故 $X$ 通有约化；完全交是 Cohen–Macaulay 的，无嵌入分支，故 $X$ 约化。正规：Cohen–Macaulay 给出 Serre 条件 $S_2$，命题 2.2 给出奇异轨迹余维 $7\ge2$，即 $R_1$，由 Serre 判据 $X$ 正规。不可约：$X$ 是锥，故连通；连通的正规簇不可约。全部论证对 $\bar{\mathbb F}_2$ 同样成立，故 $X$ 几何整。$\square$

**推论 2.3′（单值群）.** $X_{\bar{\mathbb F}_2}$ 的不可约分支与 $S_{12}/G$ 一一对应，$G$ 为通有多项式 $F_{\rm gen}=t^{12}+e_2t^{10}+\dots+e_{11}t+e_{12}$ 在 $\bar{\mathbb F}_2(e_2,e_4,e_6,\dots,e_{12})$ 上的 Galois 群（$X$ 是有序根组的簇，$\pi$ 的通有纤维是分裂代数）。由命题 2.3 分支只有一个，故 $G=S_{12}$：族 $\{F:e_1=e_3=e_5=0\}$ 的几何单值群是整个 $S_{12}$，尽管 $e_1=e_3=e_5=0$ 是特征 2 下 Legendre 区间的特殊约束。

**注 2.4.** 由命题 2.3，$H^{18}_c(X_{\bar{\mathbb F}_2},\mathcal L)=H^0(X_{\bar{\mathbb F}_2},\mathcal L^\vee)^\vee(-9)=0$（命题 1.9(iii)，$\mathcal L$ 几何非平凡于几何整簇）。这是 $Y_k=o(q^9)$ 的上同调原因；更强的 $Y_k=O(q^{9/2})$ 需要 $H^j_c=0$（$9<j<18$）或其扭迹相消，本文未证（§3.3）。

**注 2.5（非扭形式的点数并不反映维数）.** `split_form_counts.py` 给出 $\#X(\mathbb F_q)=2048,\ 2\,099\,200,\ 1\,081\,081\,856$（$q=2,4,8$），远大于 $q^9$。这不是分支数的信号：$q<12$ 时 $\mathbb F_q$ 上没有 12 个两两不同的根，$X(\mathbb F_q)$ 全由带重根的多项式的排序贡献（权 $12!/\prod m_i!$），渐近 $\#X(\mathbb F_q)\sim q^9$ 要到 $q\ge16$ 才开始。扭形式 $X^{(\sigma)}(\mathbb F_q)=U_k$ 则从 $q=2$ 起就有 $|U_k|=q^9+(-1)^{k+1}(q-1)q^5$（第四部分猜想，$q\le128$ 精确）。

---

## 3. $\ell$-进上同调与权

### 3.1 扭 Lefschetz 迹公式

**命题 3.1.** 取 $\ell\ne2$。$\sigma$ 作用于 $X$ 并保持 $e_4$，故提升为 $\mathcal L$ 的自同构；$\sigma$ 与 $F$ 交换。对每个 $k\ge1$，
$$Y_k=\sum_{j=0}^{18}(-1)^j\operatorname{Tr}\bigl(\sigma\circ F^k\ \big|\ H^j_c(X_{\bar{\mathbb F}_2},\mathcal L)\bigr),$$
其中 $F$ 是 $\mathbb F_2$ 上的几何 Frobenius 在上同调上的作用。

*证明.* 命题 1.5 把 $Y_k$ 写成扭形式 $X^{(\sigma)}/\mathbb F_q$ 上 $\mathcal L$ 的迹函数之和，$X^{(\sigma)}$ 的几何 Frobenius 在 $X_{\bar{\mathbb F}_2}$ 的上同调上作用为 $\sigma F^k$（$\sigma$ 有限阶且与 $F$ 交换）。对此应用 Grothendieck–Lefschetz 迹公式（SGA 4½, Rapport, Thm 3.2；有限阶自同构与 Frobenius 复合的情形见 Deligne–Lusztig, Ann. Math. 103 (1976), Thm 3.2 的证明）。$\square$

**推论 3.2（$L$-函数）.** $\sigma$ 与 $F$ 可同时三角化，$\sigma$ 的特征根为 12 次单位根 $\zeta$，故
$$Y_k=\sum_i\pm\zeta_i\lambda_i^k,\qquad L_\sigma(T):=\prod_j\det\bigl(1-\sigma FT\,\big|\,H^j_c\bigr)^{(-1)^{j+1}},\qquad T\frac{d}{dT}\log L_\sigma(T)=\sum_{k\ge1}Y_kT^k,$$
$L_\sigma\in\mathbb Q(\zeta_{12})(T)$，而因 $Y_k\in\mathbb Z$，实际 $L_\sigma\in\mathbb Q(T)$。$(Y_k)_k$ 是线性递推序列；$m_k=Y_k/8^k=\sum_i\pm\zeta_i(\lambda_i/8)^k$。$(13,-1,-35,47,733)$ 是这个有理函数的前五个对数导数系数。

### 3.2 权

**命题 3.3（Deligne, Weil II）.** $\mathcal L$ 是纯权 0 的光滑层，故 $F$ 在 $H^j_c(X_{\bar{\mathbb F}_2},\mathcal L)$ 上的特征根 $\lambda$ 是代数整数且对每个复嵌入 $|\lambda|\le2^{w/2}$，$w\le j$（混合，权 $\le j$）。$\sigma$ 不改变绝对值。归一化 $2^{3k}=q^3$ 对应 Tate 扭 $(3)$，即权 6：$m_k$ 是 $H^\bullet_c(X,\mathcal L)(3)$ 上 $\sigma F^k$ 的迹，其中 $H^j_c(3)$ 的权 $\le j-6$。

**推论 3.4.** 若 $H^j_c(X_{\bar{\mathbb F}_2},\mathcal L)=0$ 对 $j\ne9$（集中于中间维数 $=\dim X$）且 $H^9_c$ 纯权 9，则 $m_k=\pm\sum_{i=1}^{b_9}\zeta_i\eta_i^k$，$|\eta_i|=2^{3/2}$，从而 $|m_k|\le b_9\,q^{3/2}$，$C_6=\lim\operatorname{Var}_{\mathcal G^2}/q^9=1$（论文定理 7.1）。

### 3.3 数据对 $b_9$ 与纯性的约束

| $k$ | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| $m_k$ | 13 | $-1$ | $-35$ | 47 | 733 |
| $|m_k|/2^{3k/2}$ | 4.60 | 0.13 | 1.55 | 0.73 | 4.05 |

- $|m_1|/2^{3/2}=4.6$ 与 $|m_5|/2^{15/2}=4.05$ 要求推论 3.4 中 $b_9\ge5$；第三部分写的 "$|m_k|\le4q^{3/2}$" 不成立，特此更正。
- Newton 恒等式由 $m_1..m_5$ 给出 $e_1..e_5=13,85,361,1069,2281$。若 $b_9=5$ 且纯权 3，则 $|e_5|=\prod|\eta_i|=2^{15/2}\approx181\ne2281$，矛盾（此论证只用绝对值，对 $\sigma$-扭 $\zeta_i$ 同样有效）。故：**纯权 3 时 $b_9\ge6$，或存在低权（$j<9$ 经 Tate 扭后权 $<3$）或高次（$9<j<18$）的贡献**。
- 增长率 $|m_k|\asymp2^{3k/2}$（而非 $2^{2k}$ 或更高）在 $k\le5$ 的范围内与"集中于 $j=9$ 且纯"相容，但这不是证明：Katz 关于 $\sum\psi(f)$ 上同调集中于中间维数的定理要求 $\deg f$ 与 $p$ 互素且首项形式非退化，这里 $\deg e_4=4=2^2$，$X$ 又非光滑（命题 2.2），均不满足现成定理的假设。

### 3.4 与 $A_6,X$ 的对照

$A_6=|U_k|$ 与 $X=\sum_U\chi(e_2)$ 对应于 $\mathcal G^2\cong W_2(\mathbb F_q)$ 的**特征**（$\psi_{(0,1)}$），经 Parseval 化为 $q^2$ 条超奇异曲线 $C_\lambda$ 的 $H^1$ 上 Frobenius 迹之和（第四部分定理 7），特征根为 $\sqrt q\cdot$单位根，故为 $q$ 的幂的 $\pm$ 组合。$Y$ 对应的 $(-1)^{\operatorname{Tr}e_4}$ 在 $W_2(\mathbb F_q)$ 上**不是特征**（Witt 加法进位 $x_1x_1'$），所以 $Y$ 不能按特征分解到曲线，而是 9 维簇 $X$ 上层 $\mathcal L$ 的整体扭迹；$m_5=733$ 为素数、$m_1=13$，排除了 $q$ 的幂的少项组合。

---

## 4. 验证

**4.1 直接点集计算（已运行）.** `verify_Y_direct.py`：
- 路线 A（$q=2$）：遍历 $\mathbb F_{2^{12}}$ 的 4096 个元素，构造 12 个共轭，展开 $\prod(t-z_j)$ 得 $e_1,e_3,e_4,e_5$，筛 $e_1=e_3=e_5=0$，求和 $(-1)^{e_4}$——这正是 $X^{(\sigma)}(\mathbb F_2)$ 上 $\mathcal L$ 的迹函数之和。结果 $|U|=544$，$Y_1=104$。
- 路线 B（$q=2,4$）：命题 1.7 的推出形式，$\sum_F\Lambda(F)\psi(a_8)$，$\Lambda$ 由 `legendre_var_check.py` 对全部 $\alpha$ 的向量化特征多项式递推给出。结果 $(544,104)$、$(259072,-64)$，即 $m_1=13$，$m_2=-1$，与论文一致，断言全部通过。

**4.2 SageMath 脚本.** `verify_Y_sage.sage` 实现同样两条路线（路线 B 对 $q=4$ 需分解 $4^9$ 个多项式，数分钟），含断言；本机未安装 Sage，脚本未在此运行，其 $q=2$ 的路线 A 与 Python 路线 A 逐句对应。

---

## 5. 结论与尚未完成的事

| 问题 | 回答 | 状态 |
|---|---|---|
| $Y$ 的坐标与理想 | $z\in\mathbb A^{12}$，$I=\langle e_1,e_3,e_5\rangle=\langle p_1,p_3,p_5\rangle$ | 定理 |
| AS 覆盖 | $y^2+y=e_4(z)$，$\deg e_4=4$，极点 $4X_\infty$ | 定理 |
| "九维" | $\dim X=9$（Krull），完全交余维 3；纤维维数 $0$（$X\to\mathbb A^9$）与 $7$（$X\to\mathbb A^2$） | 定理 |
| 连通/奇异 | 几何整、正规的锥；$\operatorname{Sing}=1023$ 张平面之并，维数 2 | 定理 |
| 迹公式 | $Y_k=\sum(-1)^j\operatorname{Tr}(\sigma F^k\mid H^j_c(X,\mathcal L))$，$\sigma$ 为 12-循环 | 定理 |
| 权 | $w\le j$；$2^{3k}$ 对应 Tate 扭 (3)；若集中于 $j=9$ 且纯则 $m_k$ 权 3 | 定理（条件部分标明） |
| 集中于 $j=9$、$b_9$ 的值 | 未证；数据要求 $b_9\ge5$，排除"纯且 $b_9=5$" | 开放 |

要完成 $Y$ 的"最终闭式"，需要：(a) 证明 $H^j_c(X,\mathcal L)=0$（$j\ne9$）或至少 $\sigma F$-迹在 $j\ne9$ 相消——由于 $p\mid\deg e_4$ 与 $X$ 的奇异性，需要沿 $X_\infty$ 的 Swan 导子分析（命题 1.9(ii)）和奇异轨迹（命题 2.2）处的局部计算；(b) 计算 $b_9$ 与 $\sigma F$ 在 $H^9_c$ 上的特征多项式——由推论 3.2，这等价于算出有限多个 $Y_k$（$k\le2b_9$ 左右）后用 Berlekamp–Massey 识别线性递推；目前只有 $k\le5$，而 $k=6$（$\mathbb F_{2^{72}}$）已超出现有枚举方法。
