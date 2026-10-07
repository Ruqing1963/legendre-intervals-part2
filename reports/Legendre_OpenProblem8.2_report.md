# 函数域 Legendre 区间（特征 2）方差异动：Open Problem 8.2 的结构分析

> **更新（同日，第二部分）**：`Legendre_OpenProblem8.2_part2.md` 用 $\mathbb F_{q^4}\subset\mathbb F_{q^8}$ 的二次塔给出了 $d=4$ 对**所有** $q=2^k$ 的完整证明：(T1)、引理 1（$\Psi(\mathrm{id})=q^5-q^4+q^3$）、(T3)（符号 $(-1)^k$）、$\operatorname{Var}_{\mathcal G^2}/q^5=q(q-1)$、引理 2（Carlitz 和的 $\lambda_1$-无关性）以及 $\operatorname{Var}_{\mathcal G}/q^5=3-3/q$ 均成为定理；§3.6 的全部遗留项已关闭。并用论文的 Hayes 群 FFT 计算了 $q=16,32$ 的 $d=6$，判定偶数 $d$ 的标度为 $q^2$。下文 §3.2、§3.6 中标注"精确验证 $k\le3$"的结论请以第二部分为准。
> **第三部分** `Legendre_OpenProblem8.2_part3.md`：用 $\mathbb F_{q^{12}}\supset\mathbb F_{q^6}$ 塔与平移/伸缩对称给出 $d=6$ Legendre 剖面的结构定理（三个值，由 $c_2=0$ 与 $\operatorname{Tr}(c_4/c_2^2)$ 决定）及 $\operatorname{Var}_{\mathcal G^2}/q^7$ 的闭式，归结为三个常数 $A_6,X,Y$；$A_6,X$ 初等（$k\le5$ 验证），$Y=q^3m_k$ 为真正的指数和。
> **第四部分** `Legendre_OpenProblem8.2_part4.md`：$A_6,X$ 归结为 $q^2-1$ 个超奇异 Artin–Schreier 曲线的二次型指数和（根基 = 15 次线性化多项式的核，符号 = Arf），用独立的 $\mathbb F_2$-线性代数方法精确验证到 $q=128$；纯三次项对所有 $k$ 证毕，混合项的符号引理未证（其非零集合不是迹条件）。$Y$ 的 Grothendieck–Lefschetz 表述与权重分析；更正：$|m_k|\le4q^{3/2}$ 在 $k=1,5$ 不成立。

日期：2026-10-07。配套代码：`legendre_var_check.py`（仅依赖 numpy，所有数值均为精确整数/有理数）；输出：`legendre_var_results.csv`、`legendre_rho_quadratic.csv`、`full_run_output.txt`。
论文仓库（`legendre-intervals-explicit-formula/`）中的 `data/char2_variance.csv`、`data/char2_correlations.csv`、`data/exact_scan.csv` 被用作独立交叉验证；本报告的全部数值都由独立代码重新计算并与之逐项核对。

**学术纪律标注**：下文每条结论都标明属于
- **[定理]**：给出完整证明；
- **[精确验证]**：对 $q\in\{2,4,8\}$（或所列 $d$）由精确整数计算验证，对一般 $k$ 未证明；
- **[猜想/启发]**：仅由数据或结构类比支持。

---

## 0. 结论速览

| 量 ($d=4$, $n=8$, $\ell=3$) | $q=2$ | $q=4$ | $q=8$ | 一般 $q=2^k$ |
|---|---|---|---|---|
| $\operatorname{Var}_{\mathcal G}(\Psi)/q^5$ | $3/2$ | $9/4$ | $21/8$ | $3-3/q$ **[精确验证 $k\le3$]** |
| $\operatorname{Var}_{\mathcal G^2}(\Psi)/q^5$ | $2$ | $12$ | $56$ | $q(q-1)$ **[精确验证 $k\le3$]** |
| 比率 | $4/3$ | $16/3$ | $64/3$ | $q^2/3$ |
| $\Psi(\text{单位类})-q^5$ | $-8$ | $-192$ | $-3584$ | $-(q^4-q^3)$ **[精确验证]** |
| $\Psi(0,c,0)-q^5$，$c\ne0$ | $-8$ | $+64$ | $-512$ | $-q^3$（$k$ 奇）, $+q^3$（$k$ 偶）**[精确验证；模 Carlitz 和的 $\lambda_1$-无关性可证]** |
| $A=\#\{\alpha:\operatorname{Tr}\alpha=\operatorname{Tr}\alpha^3=0\}$ | $48$ | $4096$ | $254976$ | $q^6-2q^4+2q^3$（$k$ 奇）, $q^6$（$k$ 偶） |
| $\rho(\varepsilon)$，$\varepsilon$ 仅依赖 $p_1$（$q-1$ 个） | $1$ | $1$ | $1$ | $1$ **[定理，模 (R1)]** |
| $\rho(\varepsilon)$，$\lambda_3\neq0$（$q^2-q$ 个） | $-1/3$ | $-1/9,\ 5/9$ | $5/21$ | $k$ 奇：$\frac{q-3}{3(q-1)}$ |

三个核心结构事实：

1. **[定理]** $\mathcal G/\mathcal G^2\cong\mathbb F_q^{\,m}$，$m=\lfloor d/2\rfloor$，同构由**奇数幂和** $p_1,p_3,\dots$ 给出；非平凡二次特征恰为 Artin–Schreier 特征 $\alpha\mapsto\psi_0(\lambda_1\alpha+\lambda_3\alpha^3+\cdots)$。
2. **[定理，模一个显式特征和 (R1)；$k\le 3$ 已精确验证]** $d=4$ 时偏差 $\varphi=\Psi-q^5$ **完全支撑在子群 $\{c_1=0\}$ 上**（$q^2$ 个类）；Legendre 子群 $\mathcal G^2$ 是其中 $q$ 个类。方差"爆炸"本质上是**单位类 $t^8$ 周围区间的素数亏损 $-(q^4-q^3)$** 这一件事：$q(q-1)=(q-1)^2+(q-1)$，前一项来自单位类，后一项来自其余 $q-1$ 个 Legendre 类的 $\pm q^3$。恒等式本身对一般 $k$ 只差两步：单位类计数 (T2) 与 Carlitz 和的 $\lambda_1$-无关性（§3.6）。
3. **[精确验证]** $\sum_{\chi\ne1}\psi_8(\chi)=-q^7+q^6$，即 $\mathbb E_\chi[\operatorname{Tr}\Theta_\chi^8]\to-1\neq0$：在 $p=2$、$\ell=3$ 处 Frobenius 类并非在 $U(2)$ 中 Haar 等分布，这就是 $\operatorname{Var}_{\mathcal G}/q^5\to3$ 而非 KR 的 $d-2=2$ 的根源。

---

## 1. 设定与精确恒等式 **[定理]**

$n=2d$，$\ell=d-1$，$\mathcal G=\mathcal G_\ell=\{1+c_1u+\dots+c_\ell u^\ell\}\subset(\mathbb F_q[u]/u^{\ell+1})^\times$，$|\mathcal G|=q^\ell$。
首一多项式 $F=t^n+a_{n-1}t^{n-1}+\cdots$ 映到 $F^*(u)=u^nF(1/u)\bmod u^{\ell+1}$，即 $(c_1,\dots,c_\ell)=(a_{n-1},\dots,a_{n-\ell})$。每个类的纤维有 $q^{n-\ell}=q^{d+1}$ 个多项式。

**Legendre 区间 = $\mathcal G^2$ 的纤维。** $f^2+s$（$\deg s\le d$）的奇数位系数 $a_{n-1},a_{n-3},\dots$ 全为 $0$，而 $(1+c_1u+\cdots)^2=1+c_1^2u^2+c_2^2u^4+\cdots$，所以 $\mathcal I_f$ 恰是 $f^2$ 所在类的纤维，且这些类恰为 $\mathcal G^2$。当 $d=4$：$\mathcal G^2=\{1+cu^2\}$，$|\mathcal G^2|=q$，指数 $q^2$。

**$\Psi$ 的根计数表示。** 对 $\alpha\in\mathbb F_{q^n}$ 记 $\mathrm{charpoly}(\alpha)=\prod_{i<n}(t-\alpha^{q^i})$。若 $F=P^k$，$\deg P=n/k$，则恰有 $n/k=\Lambda(F)$ 个 $\alpha$ 以 $F$ 为特征多项式。故
$$\Psi(g)=\sum_{F\in\text{fibre}(g)}\Lambda(F)=\#\{\alpha\in\mathbb F_{q^n}:\ \mathrm{charpoly}(\alpha)\in\text{fibre}(g)\}=\#\{\alpha:\ (e_1,\dots,e_\ell)(\alpha)=g\},$$
其中 $e_j$ 是 $\alpha$ 的 $n$ 个共轭的初等对称函数。这是代码的计算原理（对所有 $\alpha$ 向量化地递推特征多项式的前 $\ell$ 个系数），也是后文所有代数推导的出发点。均值 $q^n/|\mathcal G|=q^{d+1}$。

**方差的定义（与论文一致）。** 记 $\varphi(g)=\Psi(g)-q^{d+1}$。论文 Table 4 的 `var_Lambda_over_q^(d+1)` 是 **关于全局均值 $q^{d+1}$ 的二阶矩**：
$$\operatorname{Var}_{\mathcal G}=\frac1{|\mathcal G|}\sum_{g\in\mathcal G}\varphi(g)^2,\qquad \operatorname{Var}_{\mathcal G^2}=\frac1{|\mathcal G^2|}\sum_{g\in\mathcal G^2}\varphi(g)^2 .$$
（例：$q=2,d=4$ 时两个 Legendre 区间的 $\Psi$ 都等于 $24$，关于自身均值的方差为 $0$，但关于全局均值 $32$ 的二阶矩为 $64=2\cdot q^5$，这正是 Table 4 的 "2.00"。代码同时输出两种定义。）

**Parseval 恒等式。** 设 $\psi_n(\chi)=\sum_{\deg F=n}\Lambda(F)\chi(F)$，则 $\varphi(g)=|\mathcal G|^{-1}\sum_{\chi\ne1}\psi_n(\chi)\bar\chi(g)$，于是
$$\sum_{\chi}\psi_n(\chi)\overline{\psi_n(\chi\varepsilon)}=|\mathcal G|\sum_g\varphi(g)^2\,\bar\varepsilon(g),\qquad
\rho(\varepsilon)=\frac{\sum_g\varphi(g)^2\varepsilon(g)}{\sum_g\varphi(g)^2},\qquad
\frac{\operatorname{Var}_{\mathcal G^2}}{\operatorname{Var}_{\mathcal G}}=1+\sum_{\varepsilon\ne1,\ \varepsilon|_{\mathcal G^2}=1}\rho(\varepsilon).$$
（$\varphi$ 均值为零，所以 $\chi=1$ 与 $\chi\varepsilon=1$ 的项自动消失，论文中的定义与此等价。）代码对每个 $(q,d)$ 都精确核验了此恒等式。

**$L$-函数。** 对 $\chi\ne1$，$L(u,\chi)=\sum_F\chi(F)u^{\deg F}$ 是次数 $\le\ell-1=d-2$ 的多项式（$\deg F\ge\ell$ 时类等分布），$\psi_n(\chi)=-\sum_{j}\gamma_j^{\,n}$，$|\gamma_j|=\sqrt q$。因此
$$\frac{\operatorname{Var}_{\mathcal G}}{q^{d+1}}=\frac1{|\mathcal G|}\sum_{\chi\ne1}\bigl|\operatorname{Tr}\Theta_\chi^{\,n}\bigr|^2,\qquad \Theta_\chi\in U(d-2)\ \text{(}\chi\text{ 本原时)},$$
Keating–Rudnick 的预测 $d-2$ 就是 $\int_{U(d-2)}|\operatorname{Tr}U^n|^2dU=d-2$。

---

## 2. 二次特征的结构 **[定理]**

**命题 2.1.** 设 $p_j(g)$ 为 $g\in\mathcal G_\ell$ 对应"根"的 $j$ 次幂和（由 Newton 恒等式从 $c_1,\dots,c_\ell$ 定义，特征 2 下全部符号为 $+$：$p_k=\sum_{i<k}c_ip_{k-i}+[k\text{ 奇}]\,c_k$）。则
1. $p_j(gh)=p_j(g)+p_j(h)$（根的并集），$p_{2j}=p_j^2$，故 $p_j(g^2)=2p_j(g)=0$ 对奇数 $j$；
2. $\Pi:\ g\mapsto(p_j(g))_{j\text{ 奇},\,j\le\ell}\in\mathbb F_q^{\,m}$，$m=\#\{j\le\ell\text{ 奇}\}=\lfloor d/2\rfloor$，是满同态（$p_j=c_j+(\text{仅含 }c_{<j}\text{ 的项})$，可逐个解出奇数位 $c_j$），核的阶为 $q^{\ell-m}=|\mathcal G^2|$，故 $\ker\Pi=\mathcal G^2$，$\mathcal G/\mathcal G^2\cong\mathbb F_q^{\,m}$；
3. 全部 $q^m$ 个二次特征为 $\varepsilon_\lambda(g)=(-1)^{\operatorname{Tr}_{\mathbb F_q/\mathbb F_2}\left(\sum_{j\text{ 奇}}\lambda_jp_j(g)\right)}$，$\lambda\in\mathbb F_q^{\,m}$。

**推论 2.2（Artin–Schreier 解读）。** 若 $g=\mathrm{class}(\mathrm{charpoly}\,\alpha)$，则 $p_j(g)=\operatorname{Tr}_{\mathbb F_{q^n}/\mathbb F_q}(\alpha^j)$，于是
$$\varepsilon_\lambda(\alpha)=\psi_0\Bigl(\sum_{j\text{ 奇}}\lambda_j\alpha^j\Bigr),\qquad \psi_0=\text{ }\mathbb F_{q^n}\text{ 的典范加法特征}.$$
二次扭转 = 用奇次 Artin–Schreier 多项式 $y^2+y=\sum\lambda_j x^j$ 的覆盖去扭转。$d=4$ 时 $p_1=c_1$，$p_3=c_3+c_1c_2+c_1^3$。

**推论 2.3（奇偶 $d$ 的第一个结构差异）。** $\varepsilon_\lambda$ 本原（不通过 $\mathcal G_{\ell-1}$ 分解）当且仅当 $\lambda_\ell\ne0$，这要求 $\ell$ 为奇数，即 **$d$ 为偶数**。$d$ 为奇数时，所有二次特征都是非本原的（导子 $\le u^{\ell}$）。$d$ 偶数时有 $q^m-q^{m-1}$ 个本原二次特征。

注：提示中给出的 $\varepsilon(g)=(-1)^{\operatorname{Tr}(\lambda_1c_1+\lambda_3c_3)}$ **不是**特征（$c_3$ 不可加：$c_3(gh)=c_3+c_3'+c_1c_2'+c_2c_1'$），必须用 $p_3$ 代替 $c_3$。

---

## 3. $d=4$ 的完整刻画（任务一、任务三）

记 $n=8$，类 $g=(c_1,c_2,c_3)=(a_7,a_6,a_5)$，$\mathcal G^2=\{(0,c,0)\}$，$\mathcal G_0:=\{c_1=0\}\cong(\mathbb F_q^2,+)$（$(1+xu^2+yu^3)(1+x'u^2+y'u^3)=1+(x+x')u^2+(y+y')u^3$），$\mathcal G^2\subset\mathcal G_0\subset\mathcal G$，阶 $q\mid q^2\mid q^3$。

### 3.1 切片结构与素幂贡献 **[定理]**

Prop. 8.1 的分解 $F=h^2+tb^2$，$h=t^4+h_3t^3+h_2t^2+h_1t+h_0$，$b=b_1t+b_0$：
$$F=t^8+h_3^2t^6+h_2^2t^4+b_1^2t^3+h_1^2t^2+b_0^2t+h_0^2 .$$
因此纤维 $(0,c,0)$ 恰是 $\{h_3^2=c\}$，$h_2,h_1,h_0$ 自由（$q^3$ 个 $h$），$b$ 自由（$q^2$ 个），共 $q^5$。**$b=0$ 切片包含 $q^3$ 个完全平方（不是提示所说的 $q$ 个）**，它们对 $\Psi$ 的贡献并非 $0$：$\Lambda(h^2)=\Lambda(h)$，故
$$E(c):=\sum_{b=0}\Lambda(F)=\sum_{\deg h=4,\ h_3=\sqrt c}\Lambda(h)=q^3\quad\text{（对所有 }c\text{ 精确成立）}.$$
证明：$h\mapsto h_3$ 是 $\mathcal G_1\cong\mathbb F_q$ 上的类映射，其非平凡特征 $\chi_\lambda(h)=\psi_q(\lambda h_3)$ 的 $L$-函数为 $1+\bigl(\sum_a\psi_q(\lambda a)\bigr)u=1$，所以 $\sum_{\deg h=4}\Lambda(h)\chi_\lambda(h)=0$，$\Lambda$ 在 $h_3$ 的每个值上精确等分布，每值 $q^4/q=q^3$。（直接验证：不可约四次式中 $h_3=s\ne0$ 的有 $q^3/4$ 个，$h_3=0$ 的有 $(q^3-q^2)/4$ 个，再加 $P^2$（$\deg P=2$）与 $P^4$（$\deg P=1$）全落在 $c=0$，合计 $q^3$。）代码输出 `prime-power part E = q^3` 对 $q=2,4,8$ 均成立。

于是 $\Psi(0,c,0)=q^3+8N_8(c)$，$N_8(c)$ = 纤维中不可约八次式个数，全部来自 $b\ne0$ 的 $q^2-1$ 个切片。切片内部的不可约个数并不均匀（代码 `slice structure` 输出：$q=4$ 时 $c=0$ 的 15 个切片的不可约数 $\in\{0,8\}$，$c\ne0$ 时 $\in\{6,8,12,16\}$），所谓"正交干涉"在切片层面没有干净的结构；干净的结构出现在类（Fourier）层面，见下。

### 3.2 三个精确事实

**(T1) 偏差支撑在 $\mathcal G_0$ 上**：$\Psi(g)=q^5$ 对所有 $c_1\ne0$ 的 $q^3-q^2$ 个类。**[定理，模 (R1)；$q\le8$ 精确验证]**

证明（约化）。$\Psi(s,c_2,c_3)=\#\{\alpha:\ \operatorname{Tr}\alpha=s,\ e_2=c_2,\ e_3=c_3\}$。
(i) 伸缩 $\alpha=s\beta$ 把 $(s,c_2,c_3)$ 变成 $(1,c_2/s^2,c_3/s^3)$，故只需 $s=1$，即 $(e_2,e_3)$ 在 $V_1=\{\operatorname{Tr}\beta=1\}$ 上等分布。
(ii) 平移 $\beta\mapsto\beta+\gamma$（$\gamma\in\mathbb F_q$）保持 $V_1$，且 $e_2\mapsto e_2+\gamma$，$e_3\mapsto e_3+\gamma^2$（$8\gamma=28\gamma^2=56\gamma^3=0$，$7\gamma=\gamma$，$21\gamma^2=\gamma^2$）。在群语言里这是自同构 $\tau_\gamma(g)=g\cdot(1+\gamma su^2+\gamma s(\gamma+s)u^3)$。
(iii) 对特征和 $S(\mu,\nu)=\sum_{\beta\in V_1}\psi_q(\mu e_2+\nu e_3)$，(ii) 给出 $S=\psi_q(\mu\gamma+\nu\gamma^2)S$，故 $S=0$ 除非 $\operatorname{Tr}((\mu+\sqrt\nu)\gamma)=0\ \forall\gamma$，即 $\mu=\sqrt\nu$。
(iv) 剩余的 $q-1$ 个和：用 $e_3=p_3+e_2+1$（$e_1=1$），
$$\textbf{(R1)}\qquad S(\sqrt\nu,\nu)=\psi_q(\nu)\sum_{\beta\in\mathbb F_{q^8},\ \operatorname{Tr}\beta=1}\psi_0(\nu\beta^3)\,\psi_q\bigl((\nu+\sqrt\nu)\,e_2(\beta)\bigr)=0\quad(\nu\in\mathbb F_q^\times).$$
(T1) $\iff$ (R1)。代码对 $q=2,4,8$ 直接验证 `max |Psi - q^5| = 0` 于全部 $c_1\ne0$ 的类。$q=2$ 时 (R1) 退化为 $\#\{\beta\in\mathbb F_{256}:\operatorname{Tr}\beta=1,\operatorname{Tr}\beta^3=0\}=64$，等价于 Carlitz 型三次和 $\sum_x(-1)^{\operatorname{Tr}(x^3)}=\sum_x(-1)^{\operatorname{Tr}(x^3+x)}=-32$（代码 `S(l1,l3)` 输出核实）。对一般 $k$，(R1) 是一个"三次 AS 项 $\times$ 二次型"的混合指数和，尚无现成的闭式定理，这是通向一般 $k$ 证明的**唯一障碍**（见 3.5）。

**(T2) 单位类的亏损**：$\Psi(0,0,0)=\#\{\alpha\in\mathbb F_{q^8}:e_1=e_2=e_3=0\}=q^5-q^4+q^3=q^3\,\Phi_6(q)$。**[精确验证 $q=2,4,8$]**

**(T3) 其余 Legendre 类**：$\Psi(0,c,0)-q^5=-q^3$（$k$ 奇：$q=2,8$），$=+q^3$（$k=2$：$q=4$），对所有 $c\ne0$。**[精确验证 $k\le3$；一般 $k$ 见 3.6 (R2)]**
由伸缩对称 $\alpha\mapsto s\alpha$（$c\mapsto s^2c$，$\mathbb F_q^\times$ 中每个元都是平方）知 $\Psi(0,c,0)$ 在 $c\ne0$ 上为常数 **[定理]**，因此 Legendre 子群上的 $\Psi$ 由两个整数决定：
$$A:=\sum_c\Psi(0,c,0)=\#\{\alpha:\operatorname{Tr}\alpha=0,\ \operatorname{Tr}\alpha^3=0\},\qquad \Psi(0,0,0),$$
且 $\Psi(0,c,0)=\frac{A-\Psi(0,0,0)}{q-1}$（$c\ne0$）。数据：$A=q^6-2q^4+2q^3$（$q=2$：48；$q=8$：254976），$A=q^6$（$q=4$）。给定 (T2)，这两个 $A$ 值分别精确给出 $-q^3$ 与 $+q^3$：
$\frac{q^6-2q^4+2q^3-(q^5-q^4+q^3)}{q-1}=q^3(q^2-1)$，$\frac{q^6-(q^5-q^4+q^3)}{q-1}=q^3(q^2+1)$。

### 3.3 方差闭式的证明（在 T1–T3 之下）**[定理 given T2,T3]**

$$\sum_{g\in\mathcal G^2}\varphi(g)^2=(q^4-q^3)^2+(q-1)\,q^6=q^6\bigl[(q-1)^2+(q-1)\bigr]=q^7(q-1),$$
与 (T3) 的符号无关。故
$$\boxed{\ \frac{\operatorname{Var}_{\mathcal G^2}(\Psi)}{q^5}=\frac{q^7(q-1)}{q\cdot q^5}=q(q-1)\ }$$
注意这个恒等式**只用到 (T2) 与 $A$ 的值**（不需要 (T1)）：$\sum_{\mathcal G^2}\varphi^2=(\Psi(\mathrm{id})-q^5)^2+\frac{(A-\Psi(\mathrm{id})-(q-1)q^5)^2}{q-1}$，代入 (T2) 与 3.6 中 $A$ 的两个闭式都得到 $q^7(q-1)$。
其中 $(q-1)^2q^6$ 来自单位类，$(q-1)q^6$ 来自其余 $q-1$ 个 Legendre 类：**$q(q-1)$ 中的因子 $(q-1)$ 来自单位类亏损 $q^4-q^3=q^3(q-1)$，因子 $q$ 来自 $(q-1)+1$。** 提示中"$b=0$ 与 $b\ne0$ 切片正交干涉产生 $q(q-1)$"的图像不成立；真正的机制是单位类（即 $t^8$ 周围的区间 $\{t^8+s:\deg s\le4\}$）的素数亏损占方差的 $(q-1)/q$。

同样 **[精确验证]** $\sum_{\mathcal G}\varphi^2=\sum_{\mathcal G_0}\varphi^2=3q^7(q-1)$，即 $\operatorname{Var}_{\mathcal G}/q^5=3(q-1)/q$，比率 $=q^2/3$。

### 3.4 $\rho(\varepsilon)$ 的微观结构（任务三）**[定理 given T1；数值与论文 Table 4 完全吻合]**

由 (T1)，$\sum_g\varphi^2\varepsilon(g)=\sum_{g\in\mathcal G_0}\varphi^2\varepsilon(g)$。在 $\mathcal G_0$ 上 $p_1=0$，$p_3=c_3=:y$，故 $\varepsilon_{\lambda_1,\lambda_3}|_{\mathcal G_0}=\psi_q(\lambda_3y)$ **与 $\lambda_1$ 无关**。于是：

1. $\rho(\varepsilon_{\lambda_1,0})=1$ 对全部 $q-1$ 个仅依赖 $p_1=a_{n-1}$ 的特征（导子 $u^2$，$\psi_8\equiv0$）。这解释了论文 `char2_correlations.csv` 中 $q=2,4,8$ 各恰有 $1,3,7$ 个 $\rho=1.0$ 的特征。
2. $\rho(\varepsilon_{\lambda_1,\lambda_3})=r(\lambda_3):=\dfrac{\sum_{x,y}\varphi(0,x,y)^2\psi_q(\lambda_3y)}{\sum\varphi^2}$ 对 $\lambda_3\ne0$，共 $q(q-1)$ 个特征，取值只依赖 $\lambda_3$。
3. 伸缩自同构 $\sigma_s:(c_1,c_2,c_3)\mapsto(sc_1,s^2c_2,s^3c_3)$ 保持 $\Psi$，故 $r(\lambda_3)=r(s^3\lambda_3)$：$r$ 在 $\mathbb F_q^\times/(\mathbb F_q^\times)^3$ 上为常数。$k$ 奇（$3\nmid q-1$）时 $r$ 为常数；$k$ 偶时至多 3 个值（$q=4$：$5/9,-1/9,-1/9$，对应 $\lambda_3\in\{1\},\{\omega\},\{\omega^2\}$）。
4. $k$ 奇时的显式值：$\Phi(y):=\sum_x\varphi(0,x,y)^2$ 在 $y\ne0$ 上常数，$=\frac{3q^7(q-1)-q^7(q-1)}{q-1}=2q^7$，故
$$r=\frac{q^7(q-1)-2q^7}{3q^7(q-1)}=\frac{q-3}{3(q-1)}:\qquad q=2:\ -\tfrac13,\quad q=8:\ \tfrac5{21}=0.2381 .$$
论文 Table 4 中 $q=8,d=4$ 的 56 个 $\rho=0.2381$ 与 7 个 $\rho=1.0$ 正是 $56\cdot\frac5{21}+7=\frac{64}{3}-1$。论文把"阶 2 坐标"归为 7 个特征，其 $\rho$ 之和 $7\cdot\frac5{21}=\frac53=1.6667$，"其余"$=\frac{64}{3}-1-\frac53=18.6667$，与提示引用的 1.67 / 18.67 一致——所以**相关性并非"弥散"：除 $\rho=1$ 的非本原特征外，所有本原二次特征的 $\rho$ 相同**（$k$ 奇），只是论文的坐标系（$(\mathbb Z/4)^k\times(\mathbb Z/2)^k$ 的生成元）与 $p_j$-坐标不对齐，掩盖了这一点。
5. **"哪些二次特征强相关"的判据**：$\rho(\varepsilon)=1\iff\varepsilon$ 通过 $p_1$ 分解（导子 $u^2$）$\iff\varepsilon$ 在 $\mathcal G_0$ 上平凡；其余特征的 $\rho$ 由 $\lambda_3$ 的三次剩余类决定。
6. 平均值 **[定理 given 3.3]**：$\rho_{\rm avg}=\dfrac{q^2/3-1}{q^2-1}\to\dfrac13$。注意这里的 $1/3$ 来自 $\operatorname{Var}_{\mathcal G}/q^5\to3$，是 $d=4$ 的数值巧合，**不是** $1/(d-1)$ 的一般规律：论文数据 $q=8,d=6$ 给 $\sum\rho/(q^3-1)=15.66/511\approx0.03\ne1/5$。

### 3.5 为何 $\operatorname{Var}_{\mathcal G}/q^5\to3$ 而不是 KR 的 $2$ **[精确验证 + 结构解释]**

$\sum_{\chi}\psi_8(\chi)=|\mathcal G|\Psi(\text{id})$ 与 (T2) 给出 $\sum_{\chi\ne1}\psi_8(\chi)=-q^7+q^6$，即对 $q^3$ 个特征平均 $\psi_8\approx-q^4=-q^{n/2}$：$\mathbb E_\chi\operatorname{Tr}\Theta_\chi^8\to-1$。Haar 测度下该期望为 $0$，所以 Frobenius 类在 $U(2)$ 中**不**等分布。$q=2$ 时可完全手算（$\mathcal G\cong\mathbb Z/4\times\mathbb Z/2$，$\chi=(\chi(1+u),\chi(1+u^3))$，$L(u,\chi)=1+c_1u+c_2u^2$，$c_1=\sum_a\chi(1+au)$，$c_2=\sum_{a,b}\chi(1+au+bu^2)$）：

| $\chi$ | $L(u,\chi)$ | 零点 $\gamma$ | $\psi_8=-\sum\gamma^8$ | $|\operatorname{Tr}\Theta^8|^2$ |
|---|---|---|---|---|
| $(\pm i,+1)$ | $1+(1\pm i)u$ | $-(1\pm i)$ | $-16$ | 1 |
| $(\pm i,-1)$ | $1+(1\pm i)u\pm2iu^2$ | $\sqrt2e^{\mp i75^\circ},\sqrt2e^{\pm i165^\circ}$ | $+16$ | 1 |
| $(-1,+1)=\varepsilon_{1,0}$ | $1$ | — | $0$ | 0 |
| $(-1,-1)=\varepsilon_{0,1}$ | $1+2u^2$ | $\pm i\sqrt2$ | $-32$ | 4 |
| $(1,-1)=\varepsilon_{1,1}$ | $1+2u+2u^2$ | $-1\pm i$ | $-32$ | 4 |

校验：$\sum_{\chi\ne1}|\psi_8|^2=3072=|\mathcal G|^2\operatorname{Var}_{\mathcal G}=64\cdot48$ ✓，$\sum_{\chi\ne1}\psi_8=-64=-q^7+q^6$ ✓。两个本原二次特征的零点满足 $\gamma^8=16$（$\Theta^8=\mathbf 1$），贡献 $|\operatorname{Tr}\Theta^8|^2=4$，是 KR 均值 $2$ 的两倍，这就是过剩的来源。结构上，二次特征的 $L$-函数即 Artin–Schreier 曲线 $y^2+y=\lambda_1x+\lambda_3x^3$（亏格 1）的 $L$-多项式，其 Frobenius 特征值是 $\sqrt q\cdot$（低阶单位根），$\Theta^{8}$ 退化为恒等——这是特征 2 下超奇异椭圆曲线 $y^2+y=x^3$ 的经典事实（$q=2^k$ 上其 Frobenius 满足 $\pi^2=\pm 2^k\cdot\zeta$，$\pi^{8}$ 或 $\pi^{24}$ 为实数）。

关于 KR 定理的适用范围：$d=4$ 对应 $N=\ell-1=2$ 个非平凡零点，处在 KR–Katz 等分布定理通常要求的范围（$N\ge3$，即 $h\le n-5$）之外，且 Katz 的 Witt 向量单值群计算以 $p$ 为奇数为标准假设；本数据表明 $p=2,\ \ell=3$ 时等分布确实不成立。**未查阅文献原文核对具体假设条款，此处仅陈述数据事实。**

### 3.6 通向一般 $k$ 的完整证明：剩余障碍与工具 **[明确未完成]**

全部 $d=4$ 结论归结为 $\mathbb F_{q^8}$ 上三个显式指数和：

- **(R1)**（给出 T1）：$\sum_{\operatorname{Tr}\beta=1}\psi_0(\nu\beta^3)\psi_q((\nu+\sqrt\nu)e_2(\beta))=0$；
- **(R2)**（给出 $A$）：用 Hilbert 90 写 $\alpha=\beta+\beta^q$，则 $\operatorname{Tr}(\alpha^3)=\operatorname{Tr}(\beta^{q+2}+\beta^{2q+1})=:Q'(\beta)$ 是 $\beta$ 的 **Gold 型二次型**，
$$A=q^6+q^{-2}\sum_{\nu\in\mathbb F_q^\times}W(\nu),\qquad W(\nu)=\sum_{\beta\in\mathbb F_{q^8}}(-1)^{\operatorname{Tr}(\nu(\beta^{q+2}+\beta^{2q+1}))}=\sum_{\lambda_1\in\mathbb F_q}\sum_{x}\psi_0(\nu x^3+\lambda_1x).$$
  代码核实：$q=2$：$W=-64$；$q=4$：$W\in\{-2048,1024,1024\}$，和为 $0$；$q=8$：七个 $W$ 全为 $-65536=-2q^5$。三次 AS 和 $S(\lambda_1,\lambda_3)=\sum_x\psi_0(\lambda_3x^3+\lambda_1x)$ 在全部 $(q,\lambda_1,\lambda_3)$ 上取值 $-2^{4k+1}$（$\lambda_3$ 为 $\mathbb F_{q^8}$ 中的立方）或 $+2^{4k}$（非立方），与 Carlitz（1979, *Explicit evaluation of certain exponential sums*, Math. Scand.）对 $\mathbb F_{2^{2m}}$ 的闭式 $\sum_x\psi_0(ax^3)=(-1)^{m+1}2^{m+1}$（$a$ 立方）、$(-1)^m2^m$（非立方）在 $m=4k$ 处一致（我已对 $N=2,4$ 手工复核该闭式），且**与 $\lambda_1\in\mathbb F_q$ 无关**。
  由此得到 $A$ 的闭式 **[定理，模 $\lambda_1$-无关性]**：$\lambda_3\in\mathbb F_q^\times$ 在 $\mathbb F_{q^8}$ 中是立方 $\iff$ $k$ 奇，或 $k$ 偶且 $\lambda_3$ 是 $\mathbb F_q$ 中的立方（因 $(q^8-1)/(q-1)\equiv 8\equiv2\pmod 3$ 当 $q\equiv1\pmod3$）。于是
  $k$ 奇：$W(\nu)=q\cdot(-2q^4)$，$A=q^6-q^{-2}(q-1)\cdot2q^5=q^6-2q^4+2q^3$；
  $k$ 偶：$\sum_\nu W(\nu)=q\Bigl[\tfrac{q-1}3(-2q^4)+\tfrac{2(q-1)}3q^4\Bigr]=0$，$A=q^6$。
  两者与数据（$48,\ 254976,\ 4096$）完全一致，并通过 (T2) 给出 (T3) 的符号规律。
  $\lambda_1$-无关性的部分证明：$x\mapsto x+c$（$c\in\mathbb F_{q^8}$）给出 $S(b,\lambda_3)=\psi_0(\lambda_3c^3)\,S(b+\lambda_3c^2+\sqrt{\lambda_3c},\lambda_3)$；$c\in\mathbb F_q$ 时常数项 $\operatorname{Tr}_{q/2}(8\lambda_3c^3)=0$，而 $T(c)=\lambda_3c^2+\sqrt{\lambda_3c}$ 是 $\mathbb F_q$ 上的 $\mathbb F_2$-线性映射，核为 $\{0\}\cup\{c^3=\lambda_3^{-1}\}$，故 $S(\cdot,\lambda_3)$ 在 $\mathbb F_q$ 的一个指数 $\le4$ 的子群的陪集上为常数。补全需要 Carlitz 同文中对一般 $b$ 的 $\sum_x\psi_0(ax^3+bx)$ 公式（现成工具，本报告未逐条核对其分情形条件）；
- **(R3)**（给出 $\Psi(\text{id})$）：$\#\{\beta:\ Q_a(\beta)=0,\ Q'(\beta)=0\}/q$，其中 $Q_a(\beta)=e_2(\beta+\beta^q)$ 也是 $\beta$ 的二次型。故 $\mathcal G_0$ 上**全部** $\Psi(0,x,y)$ 由**二次型束** $\{\mu Q_a+\nu Q'\}_{(\mu,\nu)\in\mathbb F_q^2}$ 的秩（根基维数）与 Arf 不变量分布决定：$\sum_\beta(-1)^{\operatorname{Tr}(\mu Q_a+\nu Q')}\in\{0,\pm2^{(8k+r_{\mu\nu})/2}\}$。

**当前状态小结**：对一般 $k$，方差恒等式 $\operatorname{Var}_{\mathcal G^2}/q^5=q(q-1)$ 已归结为 **(R3)（即 T2）与 (R2) 中的 $\lambda_1$-无关性** 两点；$\rho$ 的微观结构与 $\operatorname{Var}_{\mathcal G}/q^5=3-3/q$ 还额外需要 **(R1)（即 T1）**。三者在 $k\le3$ 均已精确验证。

**工具**：(a) Carlitz 的三次和闭式（已用于 (R2)）与 Gold 函数二次型的秩定理（Lahtonen–McGuire–Ward 2007）；(b) 二次型束的根基是线性化多项式的核，可用算子 $\sigma+\sigma^{-1}$（$\sigma=$ $q$-Frobenius，$\ker(\sigma+\sigma^{-1})=\mathbb F_{q^2}$）在 $\mathbb F_{q^8}/\mathbb F_q$ 上显式分析，预期结果按 $k$ 的奇偶与 $\nu$ 的三次剩余类分情形——这正是数据中 $k=1$ 与 $k\ge2$ 的 $A$ 不同、$k$ 偶时 $r(\lambda_3)$ 取 3 个值的来源；(c) 几何上 (R1) 是 AS 曲线 $y^2+y=\nu x^3$ 与二次型超曲面的纤维积，可用 Deligne Weil II 给出 $O(q^{7/2})$ 的界但不能给出"精确为 $0$"；精确消失需要 (b) 的代数计算或 Katz 式的单值群论证。**本报告未完成 (a)(b) 对一般 $k$ 的推演。**

---

## 4. 大 $q$ 极限与奇偶 $d$ 分裂（任务二）

### 4.1 数据（论文 `char2_variance.csv`，比率 $=\operatorname{Var}_{\mathcal G^2}/\operatorname{Var}_{\mathcal G}$；$q=2,4$ 的 $d\le12$、$d\le6$ 行已由本代码独立复算一致）

| $d$ | $q=2$ | $q=4$ | $q=8$ |
|---|---|---|---|
| 4 | 1.333 | 5.333 | 21.333 |
| 5 | 1.067 | 5.428 | 4.306 |
| 6 | 2.703 | 2.435 | 16.657 |
| 7 | 1.987 | 1.788 | 1.988 |
| 8 | 0.327 | 2.254 | — |
| 9 | 0.813 | 0.936 | — |

$\operatorname{Var}_{\mathcal G^2}/q^{d+1}$：$d=6$ 时 $q=2,4,8\mapsto15.56,\ 9.01,\ 79.09$；$d=8$ 时 $q=2,4\mapsto1.59,\ 13.25$。

### 4.2 可以证明的结构差异 **[定理]**

1. （推论 2.3）$d$ 偶 $\iff$ 存在本原二次特征（$\lambda_\ell\ne0$）。$d=4$ 的分析显示全部"过剩"来自本原二次特征（$\rho=\frac{q-3}{3(q-1)}\to\frac13$，不随 $q$ 衰减），非本原的贡献 $\rho=1$ 只有 $q-1$ 个，占比 $O(1/q)$。
2. Prop. 8.1 的度数读法：$d=2m$ 时 $\deg b\le m-1$，$F=h^2+tb^2$ 的 $t^{d+1}$（奇次）系数只能来自 $tb^2$ 的 $t\cdot b_m^2t^{2m}$，而 $b_m=0$，故 **$a_{d+1}=a_{n-\ell}=0$ 被度数强制**；$d$ 奇时 $d+1$ 偶，$a_{d+1}=h_{(d+1)/2}^2$ 自由。换言之 $d$ 偶时 Legendre 子群把 Hayes 群**最外层**坐标 $c_\ell$（通过 $p_\ell$）钉死，这恰是本原二次特征能看到的坐标；$d$ 奇时最外层坐标不受约束，二次特征全部非本原。这是提示中"$\deg a-\deg b=1$（不对称）vs $\deg a=\deg b$（对称）"的精确含义。
3. $d=4$ 的主项是单位类亏损 $-(q^4-q^3)\approx-q^{d}$，相对亏损 $1/q$。

### 4.3 关于主导幂次的结论 **[基于数据的判断]**

提示问 $d=2m$ 时 $\operatorname{Var}_{\mathcal G^2}/q^{d+1}$ 是否 $\asymp q^{2(m-1)}$。$d=4$：$q(q-1)\asymp q^2$ ✓。$d=6$（$m=3$，预测 $q^4$）：$q=8$ 时观测 $79.1\ll q^4=4096$，仅 $\approx1.2q^2$；$q=4$：$9.0$ vs $q^4=256$。**$q^{2(m-1)}$ 不被数据支持。** 三个数据点与 $\asymp q^2$ 相容但不足以确定幂次；而且 $d=6$ 的单位类亏损不再是 $q^d$ 量级（`exact_scan.csv`：$q=4,d=6$ 的最小 $\Psi=15616=q^7-3q^4$，亏损 $3q^4=3q^{d-2}$），所以 $d=4$ 的"单类主导"机制在 $d=6$ 不再成立，过剩转由更多本原二次特征的中等相关承担（$q=8,d=6$：511 个特征，$\sum\rho=15.66$）。**猜想**：对偶数 $d$，$\lim_{q\to\infty}\operatorname{Var}_{\mathcal G^2}/(q^{d+1}\cdot q^2)$ 存在且为正常数（$d=4$ 时为 $1$）；奇数 $d$ 时 $\operatorname{Var}_{\mathcal G^2}/\operatorname{Var}_{\mathcal G}$ 有界（$q=8$：$4.3,\ 2.0$）。验证需要 $q=16,\ d=4,6$（$\mathbb F_{2^{32}},\mathbb F_{2^{48}}$）的计算，超出当前精确枚举的范围，建议改用论文的 Hayes 群 FFT（`hayes_exact.py`）。

另注意 $q=4$ 时 $d=5$ 的比率 5.43 **高于** $d=4$，$q=2$ 时奇偶交替完全不可见：奇偶分裂是大 $q$ 现象，小 $q$ 被 $O(1)$ 的涨落淹没。

本代码给出的 $q=4$ Legendre 类剖面（`full_run_output.txt`）也显示 $d=4$ 的"单类主导"在相邻 $d$ 已经变形：$d=5$ 时 $\Psi$ 只依赖 $c_4$（$c_4=0$：$-360$，否则 $+184$，各 4 个类）；$d=6$ 时只有 $c_2=0$ 的 4 个类有大亏损 $-768=-3q^4$，其余 12 个类为 $\pm16$。即 $d=6$ 的亏损类是一个 $q$-族而非单点，且量级降为 $q^{d-2}$。

---

## 5. $q=2$ 固定、$d\to\infty$ 的"退相干"（任务四）**[启发式]**

论文数据（本代码独立复算 $d\le12$ 一致）：比率 $1.33,1.07,2.70,1.99,0.33,0.81,1.10,1.79,2.13,1.14,1.54,0.98,1.15$（$d=4..16$），即 $\sum_\varepsilon\rho(\varepsilon)$ 在 $[-0.67,1.7]$ 间振荡，振幅缓慢减小但到 $d=12$ 仍有 $1.13$。

启发式：$\rho(\varepsilon)$ 是 $|\mathcal G|=2^{d-1}$ 项的归一化相关；若 $\psi(\chi),\psi(\chi\varepsilon)$ 相位独立，$\rho(\varepsilon)=O(2^{-(d-1)/2})$，$2^{\lfloor d/2\rfloor}-1$ 个特征的和若符号随机则 $\to0$（$\sim2^{-d/4}$），若符号相干则 $O(1)$ 不消失。数据介于两者之间：单个 $|\rho|$ 的最大值从 $d=6$ 的 $0.65$ 降到 $d=14$ 的 $0.06$（符合 $2^{-(d-1)/2}$ 的衰减），但求和后仍有 $O(1)$ 涨落，说明存在**部分相干**。数论图像：$d$ 小时二次特征与 Frobenius 的代数关系（如 3.5 中 $\Theta^8=\mathbf 1$）是"刚性"的，少数特征控制全部方差（低维相干）；$d$ 增大后零点个数 $d-2$ 增多，二次扭转对应的 AS 曲线 $y^2+y=\sum\lambda_jx^j$ 亏格增大，其 Frobenius 不再落在有限阶子群，相关性被平均掉（退相干）。**是否 $\lim_{d\to\infty}$ 比率 $=1$：数据相容但未证明，且非单调；更稳妥的表述是 Cesàro 意义下趋于 1。** 对 $q=2$ 这一极限等价于"二次扭转渐近独立"，与 Katz–Sarnak 在 $q$ 固定、$d\to\infty$ 时并无现成定理（KR 的所有极限都是 $q\to\infty$），属真正的开放问题。

---

## 6. 验证代码（任务五）

`legendre_var_check.py`：
- `GF(N)`：$\mathbb F_{2^N}$ 的 exp/log 表（倍增法生成，运行时验证本原性），$N=2dk\le26$；
- `psi_classes`：对全部 $\alpha\in\mathbb F_{q^n}$ 向量化递推 $\prod_i(t-\alpha^{q^i})$ 的前 $\ell$ 个系数，`bincount` 得到 $\Psi$（$\alpha=0$ 单独计入单位类）；断言 $\sum\Psi=q^n$；$d=4$ 时额外记录 $(a_7,a_6,a_5,a_3,a_1)$ 以输出 $b$-切片的不可约计数；
- `Hayes`：群乘法、平方、Newton 幂和、二次特征 $\varepsilon_\lambda$（运行时断言：在 $\mathcal G^2$ 上平凡、乘性、两两不同、个数 $q^m$）；
- `analyse`：$\operatorname{Var}_{\mathcal G}$、$\operatorname{Var}_{\mathcal G^2}$（两种均值）、$\rho(\varepsilon)$、恒等式 $1+\sum\rho=$ 比率的精确核验、$d=4$ 的 $q(q-1)$ 与 $q^2/3$ 检验；
- `d4_exponential_sums`：$A$、$\Psi(\text{id})$、Gold 型和 $W(\nu)$、三次 AS 和 $S(\lambda_1,\lambda_3)$ 及其关系式的核验，$k$ 奇时 $\rho=\frac{q-3}{3(q-1)}$ 的核验。

运行：`python legendre_var_check.py`（默认 $d=4,q=2,4,8$；$q=2,d=3..12$；$q=4,d=5,6$），或 `python legendre_var_check.py 8,4`。$q=8,d=4$ 需枚举 $\mathbb F_{2^{24}}$，约 1–2 分钟，内存 <1 GB。全部输出见 `full_run_output.txt`。

---

## 7. 对提示中若干表述的更正

1. "$b=0$ 切片有 $q$ 个完全平方，素数贡献为 $0$"：应为 $q^3$ 个完全平方，$\Lambda$-贡献恰为 $q^3$（§3.1）。
2. "$q(q-1)$ 由 $b=0$ 与 $b\ne0$ 切片正交干涉产生"：切片层面无此结构；$q(q-1)=(q-1)^2+(q-1)$ 来自单位类亏损 $q^3(q-1)$ 与其余 $q-1$ 个类的 $\pm q^3$（§3.3）。
3. "$\lim\rho_{\rm avg}=1/3=1/(d-1)$"：$1/3$ 成立，但来自 $\operatorname{Var}_{\mathcal G}/q^5\to3$；$1/(d-1)$ 不推广（§3.4.6）。
4. "相关性弥散在整个二次特征群上"：恰相反，除 $q-1$ 个 $\rho=1$ 的非本原特征外，$k$ 奇时所有本原二次特征的 $\rho$ 相同（§3.4）。
5. "$\varepsilon(g)=(-1)^{\operatorname{Tr}(\lambda_1c_1+\lambda_3c_3)}$"：不是群特征，应用 $p_3=c_3+c_1c_2+c_1^3$（§2）。
6. "偶数 $d=2m$ 时主导幂次 $q^{2(m-1)}$"：$d=6$ 数据否定（§4.3）。
