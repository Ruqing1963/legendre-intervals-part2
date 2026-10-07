# 第四部分：$A_6$、$X$ 的超奇异曲线/二次型刻画，及 $Y$ 的上同调定位

日期：2026-10-07。承接第三部分。代码：`legendre_d6_qf.py`（$\mathbb F_2$-二次型的精确指数和：根基 + Arf 不变量，$N=12k\le84$ 比特），`legendre_d6_qf_pattern.py`、`legendre_d6_qf_pattern2.py`（取值模式分析）；输出 `qf_k456.txt`、`qf_k7.txt`、`qf_pattern_k7.txt`。

**先说结论（诚实版）。**
- $A_6$、$X$ 的闭式现已由**与 FFT 完全独立的第二种精确方法**验证到 $q=128$（$k=1,\dots,7$）；
- 两者都严格归结为 $\mathbb F_{q^{12}}$ 上有限多个（$q^2-1$ 个）超奇异 Artin–Schreier 曲线族指数和的求和，每一项的绝对值由一个显式线性化多项式的根的个数决定、符号由 Arf 不变量决定（§1–2，**定理**）；
- 其中"纯三次"部分 $W(\lambda_3,0)$ 对所有 $k$ **已完整证明**（§3，引理 E）；"纯五次"与"混合"部分的根基维数与符号分布已精确算出并给出组合恒等式（§3），但**对一般 $k$ 的证明未完成**：混合项的非零集合不是任何简单的迹条件（§3.3 的扫描），需要逐项分析 15 次线性化多项式根的子域归属与 Arf 不变量。因此 $A_6,X$ 的闭式**仍不是无条件定理**；本部分把缺口压缩成一个完全显式、有限、可机械验证的"符号引理"。
- $Y$：给出 Grothendieck–Lefschetz 迹公式下的上同调表述与权重分析（§4）。要更正第三部分的一处粗略说法：$|m_k|\le4q^{3/2}$ 在 $k=1,5$ **不成立**（$13>4\sqrt8$，$733>4\cdot32^{3/2}$），因此若 $Y$ 来自纯权 9 的上同调，其秩至少为 5；并且 $m_1..m_5$ 的 Newton 恒等式给出整数初等对称函数，但 5 个绝对值为 $2^{3/2}$ 的代数整数之积不可能是整数，故"纯权 3、秩 5"也被排除——真实结构是秩 $\ge6$ 或含混合权。

---

## 1. 定理 7：$A_6-q^9$ 是超奇异曲线族的 Frobenius 迹和 **[定理]**

$\mathcal G^2$ 的指示函数展开为 $\mathcal G/\mathcal G^2$ 的特征：$\sum_{g\in\mathcal G^2}\Psi(g)=\dfrac{|\mathcal G^2|}{|\mathcal G|}\sum_{\varepsilon^2=1}\psi_{12}(\varepsilon)$。由第一部分命题 2.1，$\varepsilon=\varepsilon_\lambda$，$\lambda=(\lambda_1,\lambda_3,\lambda_5)\in\mathbb F_q^3$，且 $\psi_{12}(\varepsilon_\lambda)=\sum_{\alpha\in\mathbb F_{q^{12}}}\psi_0(\lambda_5\alpha^5+\lambda_3\alpha^3+\lambda_1\alpha)=:S(\lambda)$。故
$$A_6-q^9=\frac1{q^3}\sum_{\lambda\ne0}S(\lambda).\tag{1.1}$$
$S(\lambda)$ 是 Artin–Schreier 曲线 $C_\lambda:\ y^2+y=\lambda_5x^5+\lambda_3x^3+\lambda_1x$ 的仿射点数偏差：$\#C_\lambda(\mathbb F_{q^{12}})^{\rm aff}=q^{12}+S(\lambda)$，且若 $\omega_1,\dots,\omega_{2g}$ 是 $C_\lambda/\mathbb F_q$ 的 Frobenius 特征值（$g=2,1,0$ 按 $\lambda_5\ne0$、$\lambda_5=0\ne\lambda_3$、仅 $\lambda_1\ne0$），则 $S(\lambda)=-\sum_j\omega_j^{12}$。

**二次型结构。** $\operatorname{Tr}_{q^{12}/2}(\lambda_5x^5+\lambda_3x^3)=\operatorname{Tr}\bigl(x\cdot(\lambda_5x^4+\lambda_3x^2)\bigr)$ 是 $\mathbb F_2$-二次型 $Q_\lambda$（$x\cdot L(x)$ 型，$L$ 线性化），$\operatorname{Tr}(\lambda_1x)$ 是线性型。故对每个 $m\ge1$，$\sum_{x\in\mathbb F_{q^m}}(-1)^{\operatorname{Tr}(\lambda_5x^5+\lambda_3x^3+\lambda_1x)}\in\{0,\pm2^{(mk+r_m)/2}\}$；由此 $\omega_j/\sqrt q$ 全为单位根（van der Geer–van der Vlugt："$y^2+y=xL(x)$ 型曲线超奇异"），即 **$\omega_j=\sqrt q\,\zeta_j$**——这正是用户提示的形式。极化形式与根基：
$$B_\lambda(x,y)=\operatorname{Tr}\bigl(y\,T_\lambda(x)\bigr),\quad T_\lambda(x)=\lambda_5x^4+(\lambda_5x)^{1/4}+\lambda_3x^2+(\lambda_3x)^{1/2},\qquad
R_\lambda=\ker T_\lambda=\{x:\ \lambda_5^4x^{16}+\lambda_3^4x^8+\lambda_3^2x^2+\lambda_5x=0\},$$
$r_\lambda=\dim_{\mathbb F_2}R_\lambda=\log_2\bigl(1+\#\{\text{非零根}\in\mathbb F_{q^{12}}\text{ of }\lambda_5^4x^{15}+\lambda_3^4x^7+\lambda_3^2x+\lambda_5\}\bigr)\in\{0,\dots,4\}$，且
$$S(\lambda)=\begin{cases}0&\text{若 }\operatorname{Tr}(\lambda_1x)\ne Q_\lambda(x)\text{ 于某 }x\in R_\lambda\\ \pm2^{(12k+r_\lambda)/2}&\text{否则（符号 = Arf 不变量）}.\end{cases}$$

**Hilbert 90 形式。** 对 $\lambda_1$ 求和：$\sum_{\lambda_1}S(\lambda_1,\lambda_3,\lambda_5)=W(\lambda_3,\lambda_5):=\sum_{\beta\in\mathbb F_{q^{12}}}(-1)^{\operatorname{Tr}(\lambda_3\beta^2\theta\beta+\lambda_5\beta^4\theta\beta)}$（$\alpha=\beta+\beta^q$，$\theta=\sigma+\sigma^{-1}$），于是
$$A_6=q^9+\frac1{q^3}\sum_{(\lambda_3,\lambda_5)\ne0}W(\lambda_3,\lambda_5),\qquad X=\frac1{q^3}\sum_{\lambda_3,\lambda_5}W'(\lambda_3,\lambda_5),\quad W'=\sum_\beta(-1)^{\operatorname{Tr}(\lambda_3\beta^2\theta\beta+\lambda_5\beta^4\theta\beta+\beta^{q+1}+\beta)}.\tag{1.2}$$
$W,W'$ 的极化形式为 $\operatorname{Tr}(\gamma\cdot T(\beta))$，$T(\beta)=M(\theta\beta)$，$M(y)=\lambda_3y^2+\sqrt{\lambda_3y}+\lambda_5y^4+(\lambda_5y)^{1/4}$（$W'$ 另加 $+y$），根基 $=\theta^{-1}(\ker M\cap\operatorname{im}\theta)$，$\ker\theta=\mathbb F_{q^2}$，$\operatorname{im}\theta=\ker\operatorname{Tr}_{q^{12}/q^2}$，故根基维数 $=2k+\dim(\ker M\cap\operatorname{im}\theta)$。

---

## 2. 精确算法与 $q\le128$ 的验证

`legendre_d6_qf.py` 把 $\mathbb F_{q^{12}}=\mathbb F_q[z]/(m)$（$m$ 由 Rabin 判别随机搜得）视为 $\mathbb F_2^{12k}$，预计算 $\beta\mapsto\beta^2,\beta^{1/2},\sigma,\theta$ 与迹配对矩阵 $G_{ij}=\operatorname{Tr}(e_ie_j)$，对每个 $(\lambda_3,\lambda_5)$ 组装极化矩阵 $B=T^{\mathsf T}G$ 与 $Q(e_i)$，然后：求根基 $R$（$\mathbb F_2$ 零空间）、检验 $Q+L\equiv0$ 于 $R$、取补空间 $V'$、解 $B(x_0,\cdot)=L|_{V'}$、辛 Gram–Schmidt 求 Arf，得 $\sum=2^{\dim R}(-1)^{\operatorname{Arf}+Q(x_0)}2^{\dim V'/2}$。利用对称 $(\lambda_3,\lambda_5)\sim(s^3\lambda_3,s^5\lambda_5)\sim(\lambda_3^2,\lambda_5^2)$（对 $W'$ 只有后者）压缩轨道。

| $k$ | $q$ | $N$ | $A_6=q^9+(-1)^{k+1}(q-1)q^5$ | $X=(-1)^{k+1}q^6+q^5$ | 用时 |
|---|---|---|---|---|---|
| 1–3 | 2–8 | 12–36 | ✓ | ✓ | <1 s |
| 4 | 16 | 48 | ✓ | ✓ | 2 s |
| 5 | 32 | 60 | ✓ | ✓ | 9 s |
| 6 | 64 | 72 | ✓ | ✓ | 48 s |
| 7 | 128 | 84 | ✓ | ✓ | 225 s |

（$k\le5$ 与 Hayes FFT 的值逐一相同；$k=6,7$ 为新数据。）

**$W(\lambda_3,\lambda_5)/q^7$ 的分布**（$q^2-1$ 个非平凡对；"纯三次" $\lambda_5=0$，"纯五次" $\lambda_3=0$，"混合"两者非零；括号内为根基维数）：

| $k$ | 纯三次（$q-1$ 个） | 纯五次（$q-1$ 个） | 混合（$(q-1)^2$ 个） |
|---|---|---|---|
| 1 | $0$ $(4)$ | $2$ $(4)$ | $0$ $(6)$ ×1 |
| 3 | $0$ $(8)$ | $2$ $(8)$ | $2$ $(8)$ ×21, $0$ $(10)$ ×28 |
| 5 | $0$ $(12)$ | $2$ $(12)$ | $2$ $(12)$ ×465, $0$ $(14)$ ×496 |
| 7 | $0$ $(16)$ | $2$ $(16)$ | $2$ $(16)$ ×8001, $0$ $(18)$ ×8128 |
| 2 | $-2$ $(6)$ | $0$ $(8)$ | $-4$ $(8)$ ×3, $1$ $(4)$ ×6 |
| 4 | $-2$ $(10)$ | $-4$ $(12)$ ×3, $1$ $(8)$ ×12 | $-4$ $(12)$ ×75, $1$ $(8)$ ×90, $0$ $(12)$ ×60 |
| 6 | $-2$ $(14)$ | $0$ $(16)$ | $-4$ $(16)$ ×1386, $1$ $(12)$ ×1638, $0$ $(16)$ ×945 |

**可读出的规律（$k\le7$ 精确）：**
- $k$ 奇：$W(\lambda_3,0)=0$，$W(0,\lambda_5)=2q^7$；混合项取 $2q^7$ 恰在 $(q-1)(q-2)/2$ 个对上、否则为 $0$。总和 $(q-1)\cdot2q^7+\tfrac{(q-1)(q-2)}2\cdot2q^7=(q-1)q^8$ ✓。
- $k$ 偶：$W(\lambda_3,0)=-2q^7$；$W(0,\lambda_5)=0$（$k\equiv2\bmod4$），或（$k\equiv0\bmod4$，此时 $5\mid q-1$）$-4q^7$ 于 $(q-1)/5$ 个五次剩余、$+q^7$ 于其余；混合项 $\in\{-4q^7,\,q^7,\,0\}$，每个固定 $\lambda_5$ 下三类个数 $(n_-,n_+,n_0)$ 为 $(1,2,0),(5,6,4),(22,26,15)$（$k=2,4,6$），满足 $n_+-4n_-=2-q$，正是 $\sum W=-(q-1)q^8$ 所需。

---

## 3. 已证部分与剩余缺口

### 3.1 引理 E（纯三次项，所有 $k$）**[定理]**
$$W(\lambda_3,0)=\begin{cases}0&k\text{ 奇}\\-2q^7&k\text{ 偶}\end{cases}\qquad(\lambda_3\in\mathbb F_q^\times).$$
**证明.** $\ker M=\{0\}\cup\{y:y^3=\lambda_3^{-1}\}$，三个立方根 $y$ 都在 $\mathbb F_{q^{12}}$ 中（$3\mid(q^{12}-1)/(q-1)$）。记 $\zeta=y^{q-1}\in\mu_3$：$\zeta=1$ 时 $y\in\mathbb F_q$；$\zeta\ne1$ 时，$k$ 奇（$q\equiv2\bmod3$，$\zeta^q=\zeta^2$）得 $y\in\mathbb F_{q^2}$，$k$ 偶（$\zeta^q=\zeta$）得 $y\in\mathbb F_{q^3}$。三者都满足 $\operatorname{Tr}_{q^{12}/q^2}(y)=0$，故 $\ker M\subset\operatorname{im}\theta$，根基 $R=\theta^{-1}(\ker M)$，$\dim R=2k+2$，$|W|\in\{0,2^{7k+1}=2q^7\}$。$Q(\beta)=\operatorname{Tr}_{q^{12}/2}(\lambda_3\beta^2\theta\beta)$ 在 $R$ 上：$\beta\in\ker\theta=\mathbb F_{q^2}$ 时为 $0$；对 $\theta\beta=y$，不同解相差 $c\in\mathbb F_{q^2}$，$\operatorname{Tr}(\lambda_3c^2y)=0$（$c^2y$ 落在 $\mathbb F_{q^2}$ 或 $\mathbb F_{q^6}$，指数为偶的子域），故只需一个特解 $\beta_y$。
- $y\in\mathbb F_{q^2}$（$k$ 奇，$\zeta\ne1$）：在 $\mathbb F_{q^4}$ 上 $\theta\beta=\sigma(\operatorname{Tr}_{q^4/q^2}\beta)$，取 $\beta_y\in\mathbb F_{q^4}$ 使 $\operatorname{Tr}_{q^4/q^2}\beta_y=y^q$。则 $Q(\beta_y)=\operatorname{Tr}_{q^4/2}(\lambda_3y\beta_y^2)=\operatorname{Tr}_{q^2/2}\bigl(\lambda_3y(\operatorname{Tr}_{q^4/q^2}\beta_y)^2\bigr)=\operatorname{Tr}_{q^2/2}(\lambda_3y^{2q+1})=\operatorname{Tr}_{q^2/2}(\zeta^2)=k\cdot\operatorname{Tr}_{4/2}(\omega)=k\equiv1$。故 $Q|_R\ne0$，$W=0$。
- $y\in\mathbb F_q$：同法 $Q(\beta_y)=\operatorname{Tr}_{q^2/2}(\lambda_3y^3)=\operatorname{Tr}_{q^2/2}(1)=2k\equiv0$。$y\in\mathbb F_{q^3}$（$k$ 偶，$\lambda_3$ 非立方）：$\operatorname{Tr}_{q^3/q}y=y(1+\zeta+\zeta^2)=0$，$\theta$ 在 $\mathbb F_{q^3}$ 的迹零部分可逆，故 $\beta_y\in\mathbb F_{q^3}$，$Q(\beta_y)=\operatorname{Tr}_{q^{12}/2}(\cdot\in\mathbb F_{q^3})=0$。故 $k$ 偶时 $Q|_R=0$，$|W|=2q^7$。
- 符号：$W=\sum_{\lambda_1\in\mathbb F_q}S(\lambda_1,\lambda_3,0)$，每项 $|S|=2q^6$（$r=2$），$|W|=2q^7=q\cdot2q^6$ 迫使 $q$ 项同号，等于 $S(0,\lambda_3,0)$ 的符号；由 Hasse–Davenport（第二部分 §3）$S(0,\lambda_3,0)=(-1)^{6k+1}2^{6k+1}=-2q^6$。故 $W=-2q^7$。$\square$

（顺带得到：$d=6$ 时 Carlitz 和的 $\lambda_1$-无关性在 $k$ 奇时**失效**——$\mathbb F_q\not\subset T(\mathbb F_{q^6})$，因为立方根落在 $\mathbb F_{q^2}$ 而 $[\mathbb F_{q^6}:\mathbb F_{q^2}]=3$ 为奇——这解释了 §2 表中 $k$ 奇、$\lambda_5=0$ 列的 $0$。）

### 3.2 剩余缺口的精确表述 **[符号引理，未证]**
记 $\varepsilon_0=(-1)^{k+1}$。$A_6$ 的闭式等价于
$$\sum_{\lambda_5\ne0}W(0,\lambda_5)+\sum_{\lambda_3\lambda_5\ne0}W(\lambda_3,\lambda_5)=(q-1)q^7(\varepsilon_0q+2).\tag{3.1}$$
每一项由 $\ker M_{\lambda_3,\lambda_5}$（15 次多项式 $\lambda_5^4y^{15}+\lambda_3^4y^7+\lambda_3^2y+\lambda_5$ 在 $\mathbb F_{q^{12}}$ 中的根）与 $\operatorname{im}\theta$ 的交、$Q$ 在根基上的取值、以及 Arf 不变量决定。纯五次项的根基是 $\{y^{15}=\lambda_5^{-3}\}$ 的子集，其子域归属按 $k\bmod4$（$5\mid q^4-1$）分类，可仿引理 E 处理；**混合项**是真正的困难：§2 的数据表明其非零集合在 $k$ 奇时是轨道不变量 $j=\lambda_3^5/\lambda_5^3\in\mathbb F_q^\times$ 的函数（`legendre_d6_qf_pattern2.py` 验证 $W$ 只依赖 $j$），恰有 $(q-2)/2$ 个 $j$ 非零，但该集合**不是** $\{\operatorname{Tr}(j^m)=c\}$ 型（$k=3$ 时碰巧等于 $\{\operatorname{Tr}j=0\}$，$k=5$ 时对所有 $m\le q-2$ 均不成立），也不是 $\{\operatorname{Tr}(j+1/j)=c\}$ 等；它由 15 次多项式的根落入 $\ker\operatorname{Tr}_{q^{12}/q^2}$ 的方式刻画，没有初等闭式。因此我不认为（3.1）能用"超奇异曲线点数公式"一步证出：点数公式给出 $|S(\lambda)|$ 的候选值，但哪些 $\lambda$ 取 $0$、哪些取 $\pm$，是 15 次线性化多项式根的精细算术。

### 3.3 $X$ 的情形
同样的刻画对 $W'$ 成立（§2 表之后的分布见 `qf_k*.txt`），$X$ 的闭式等价于 $\sum_{\lambda_3,\lambda_5}W'=q^3(\varepsilon_0q^6+q^5)$，验证到 $k=7$。$\lambda_3=\lambda_5=0$ 项 $W'(0,0)=\sum_\beta(-1)^{\operatorname{Tr}(\beta^{q+1}+\beta)}$ 的数据为 $W'(0,0)/q^6=2,-4,8,-16,32,-64,128$（$k=1..7$），即 $W'(0,0)=(-1)^{k+1}q^7=\varepsilon_0q^7$：这是 $\mathbb F_{q^{12}}$ 上 Hermite 型二次型 $\operatorname{Tr}(\beta^{q+1})$ 带线性项的 Gauss 和（根基 $=\mathbb F_{q^2}$，维数 $2k$），其值可由 Coulter 关于 $\sum\chi(ax^{2^\alpha+1}+bx)$ 的定理直接给出；其余项与 $A_6$ 同难。

---

## 4. $Y$ 的上同调定位 **[现代表述 + 诚实的限制]**

### 4.1 代数簇与层
把 $\mathbb F_{q^{12}}$ 视为 $\mathbb F_q$ 上 12 维向量空间，$\alpha\mapsto(e_1,\dots,e_{12})(\alpha)$ 与 $p_j(\alpha)$ 都是 $\mathbb F_q$-多项式映射（$p_j$ 为 $j$ 次型，$e_j$ 为 $j$ 次型）。令
$$V=\{p_1=p_3=p_5=0\}\subset\mathbb A^{12}_{\mathbb F_q},\qquad \dim V=9,$$
$V(\mathbb F_q)=U$。则
$$Y=\sum_{P\in V(\mathbb F_q)}\psi_q\bigl(e_4(P)\bigr)=\sum_{i}(-1)^i\operatorname{Tr}\bigl(\mathrm{Frob}_q\,\big|\,H^i_c(V_{\bar{\mathbb F}_q},\ e_4^*\mathcal L_\psi)\bigr)$$
（Grothendieck–Lefschetz 迹公式，$\mathcal L_\psi$ 为 $\mathbb A^1$ 上的 Artin–Schreier 层）。Deligne（Weil II）：$H^i_c$ 的权 $\le i$。$Y$ 的"预期"大小是 $q^{\dim V/2}=q^{9/2}$，对应 $H^9_c$ 的纯权 9 部分；归一化 $m_k=Y/q^3$ 则对应权 $9-6=3$。

对比 $A_6$、$X$：它们是 $\mathcal G^2$ 上**特征**的和（$\chi(e_2)=\psi_{(0,1)}$ 是 $W_2(\mathbb F_q)$ 的特征），经 Parseval 变成 $q^2$ 个**曲线**（$C_\lambda$，亏格 $\le2$）的 $H^1$ 的 Frobenius 迹，而且这些曲线超奇异（§1），所以值是 $q$ 的幂的 $\pm$ 组合。$Y$ 则不是特征和：$(-1)^{\operatorname{Tr}e_4}$ 在 $W_2(\mathbb F_q)$ 上不是特征（Witt 加法的进位 $x_1x_1'$ 破坏可乘性；$\psi_{(1,0)}=i^{\operatorname{Tr}x_1}(-1)^{\operatorname{Tr}x_2+\cdots}$ 才是特征），因此 $Y$ 不能线性分解为单个 $L$-函数的零点幂和，而是 9 维簇 $V$ 上四次型层的"整体"迹。

### 4.2 权重与秩：数据说了什么
$$m_k=Y_k/q^3:\quad 13,\ -1,\ -35,\ 47,\ 733\ (k=1..5),\qquad |m_k|/q^{3/2}=4.60,\ 0.13,\ 1.55,\ 0.73,\ 4.05 .$$
- 若 $Y_k=\pm\sum_{j=1}^{b}\omega_j^k$（$|\omega_j|=2^{9/2}$，即 $H^9_c$ 纯权 9、秩 $b$，且为 $\mathbb F_2$ 上同一对象的 $\mathbb F_{2^k}$-迹），则 $|m_k|\le b\,q^{3/2}$，数据要求 $b\ge5$。第三部分写的 "$|m_k|\lesssim4q^{3/2}$" 在 $k=1,5$ 不成立，特此更正。
- Newton 恒等式：由 $m_1..m_5$ 得 $e_1..e_5=13,85,361,1069,2281$，**全为整数**——与"$m_k$ 是固定代数整数的幂和"相容；但 5 个绝对值 $2^{3/2}$ 的数之积 $=\pm2^{15/2}\notin\mathbb Z$，故**纯权 3 且秩 5 不可能**；最简单的相容情形是秩 $\ge6$（偶秩，$\prod=\pm2^{9}$ 或 $\pm2^{3b/2}$）或权 3 与权 $\le2$ 的混合（例如 $H^8_c$ 的权 8 部分贡献 $q^{-3}\cdot q^4=q$ 量级的项）。
- 一个重要的技术点：$V$ 本身依赖 $q$（它是 $\operatorname{Res}_{\mathbb F_{q^{12}}/\mathbb F_q}\mathbb A^1$ 中由相对迹定义的子簇），所以 $(Y_k)_k$ 并非字面上同一个 $\mathbb F_2$-簇的 $\mathbb F_{2^k}$-点计数。要得到真正的 $L$-函数，应改用 Hayes 群的模空间表述：$\Psi(g)=\#\{(F,\alpha):F\in\text{fibre}(g),\ F(\alpha)=0,\ \alpha\in\mathbb F_{q^{12}}\}$，参数空间 $\{(c_2,c_4,c_6,\dots,c_{12})\}=\mathbb A^9$ 定义在 $\mathbb F_2$ 上，$Y_k$ 是 $\mathbb A^9$ 上某个（关于 $(c_2,c_4)$ 的 Fourier 变换后的）$\ell$-进复形的 $\mathbb F_{2^k}$-迹——在这个意义下 $(13,-1,-35,47,733,\dots)$ 是一个有理 $L$-函数的系数，由有限多个数据决定，而 Weil II 保证其零点/极点的绝对值是 $2^{w/2}$。

### 4.3 为什么"纯权迹"就是最终闭式
在 Deligne 之后，一个指数和族被视为"已解出"，当且仅当把它写成了某个 $\ell$-进层 $\mathcal F$（已知秩、已知权、已知单值群）在 Frobenius 下的迹：此时 (i) 其在所有扩域上的值由有限个数（$L$-函数的系数）确定；(ii) 绝对值由权给出（这里就是 $m_k=O(q^{3/2})$，足以推出第三部分推论 $C_6=1$）；(iii) 统计分布由单值群给出（Katz–Sarnak）。$A_6,X$ 是退化情形：对应的层（$C_\lambda$ 的 $H^1$）超奇异，Frobenius 特征值是 $\sqrt q$ 乘单位根，于是迹是 $q$ 的幂的有限 $\pm$ 组合，"闭式"退化为初等公式；$Y$ 是非退化情形：$m_5=733$ 为素数，$m_1=13$，不可能是 $\pm2^{a}$ 的少项组合，所以不存在 $q$ 的初等公式，$L$-函数本身就是答案。

**但必须说清楚本报告没有做到的事**：我没有确定 $\mathcal F$——既没有证明 $H^i_c(V,e_4^*\mathcal L_\psi)$ 集中在 $i=9$（这需要 $e_4|_V$ 在无穷远处的非退化性/Katz–Laumon 型估计，而 $e_4$ 是 $\mathbb F_{q^{12}}$ 的坐标上的四次型，$V$ 由 1、3、5 次型切出，相应的奇点分析未做），也没有算出秩 $b$。能负责任地说的是：$Y$ 是上述簇上上述层的迹（这是迹公式，无条件）；数据与"纯权 9、偶秩 $\ge6$"或"权 9 为主、含低权修正"相容；第三部分 §3 的塔归约给出了另一个可能更有用的几何模型——$Y/G-1$ 是 $a\in\mathbb F_{q^6}$（$\operatorname{Tr}a=0$，$a\ne0$）上的和，被加项对 $(\lambda,\sqrt\mu)$ 求和后是 $\mathbb F_q$ 上三次 AS 和，即 **$Y$ 可写成 $\operatorname{Res}_{\mathbb F_{q^6}/\mathbb F_q}$ 的 5 维参数空间上一个椭圆型 AS 曲线族的 $H^1$ 的"相对迹的和"**；秩的计算应从这里入手。

---

## 5. 小结

| 项 | 第三部分状态 | 本部分状态 |
|---|---|---|
| $A_6,X$ 闭式 | $k\le5$（FFT） | $k\le7$（独立二次型方法）；严格归结为 (3.1)；纯三次项全证（引理 E） |
| 超奇异结构 | 提及 | 定理 7：$S(\lambda)=-\sum\omega_j^{12}$，$\omega_j=\sqrt q\zeta_j$，根基 = 15 次线性化多项式的核 |
| 混合项 | — | 非零集合不是迹条件；给出精确计数 |
| $Y$ | $q^3m_k$ | 迹公式表述；权 3 分析；更正 $|m_k|\le4q^{3/2}$；秩 $\ge5$，纯权 3 秩 5 排除 |
