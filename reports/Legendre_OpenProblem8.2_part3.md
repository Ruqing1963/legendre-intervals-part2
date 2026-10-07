# 第三部分：$d=6$ 的 Legendre 剖面——$\mathbb F_{q^{12}}\supset\mathbb F_{q^6}$ 塔

日期：2026-10-07。承接第一、二部分。代码：`legendre_d6_tower.py`（塔恒等式的穷举核验、结构定理的逐类核验、三个常数的提取、塔归约公式的核验）；数据来源仍是论文的 Hayes 群 FFT（`hayes_fft_fq.py`，$q=2,\dots,32$，$q\le16$ 经 CRT 精确化）。

**结论概览。** 设 $n=12$，$\ell=5$，Legendre 子群 $\mathcal G^2=\{(0,c_2,0,c_4,0)\}$（$q^2$ 个类），
$$U=\{\alpha\in\mathbb F_{q^{12}}:\ \operatorname{Tr}\alpha=\operatorname{Tr}\alpha^3=\operatorname{Tr}\alpha^5=0\},\qquad
A_6=|U|,\quad X=\sum_{\alpha\in U}\chi(e_2(\alpha)),\quad Y=\sum_{\alpha\in U}\chi(e_4(\alpha)).$$

| 结论 | 状态 |
|---|---|
| **定理 5（结构）**：$\Psi(0,c_2,0,c_4,0)=q^7+\dfrac{A_6-q^9-X}{q^2}+\dfrac Xq[c_2=0]+\dfrac Yq[c_2\ne0]\,\chi\!\left(\dfrac{c_4}{c_2^2}\right)$ | **定理**（§2，平移与伸缩对称） |
| **定理 6（塔归约）**：$A_6,X,Y$ 化为 $\mathbb F_{q^6}$ 上的显式有限和（含 $e_2$ 的 Gauss 和） | **定理**（§3，$q=2,4,8$ 数值核验） |
| $A_6=q^9+(-1)^{k+1}(q-1)q^5$，$X=(-1)^{k+1}q^6+q^5$ | **精确验证 $k\le5$**，证明路线见 §4 |
| $Y=q^3m_k$，$m_k=13,-1,-35,47,733$（$k=1..5$） | 真正的指数和，无初等闭式（§4） |
| $\operatorname{Var}_{\mathcal G^2}/q^7$ 的闭式与 $C_6=\lim\operatorname{Var}_{\mathcal G^2}/q^9=1$ | **定理**（模 $A_6,X$ 的闭式与 $Y=o(q^{11/2})$） |

与 $d=4$ 的本质区别：$d=4$ 时三个相应的常数全部初等（故 $\Psi$ 全部初等）；$d=6$ 时 $A_6,X$ 初等但 $Y$ 不是——它是 $U$ 上一个四次型的特征和，这正是 Table 4 中 $d=6$ 列"不规则"数字（如 $q=32$ 的 $733$）的来源，也是 Legendre 区间在 $d\ge6$ 开始呈现真正算术随机性的标志。

---

## 1. 塔恒等式 **[引理 C]**

对 $\alpha\in\mathbb F_{q^{12}}$ 令 $a=\alpha+\alpha^{q^6}$，$n=\alpha^{1+q^6}$（$\in\mathbb F_{q^6}$），$\operatorname{Tr}=\operatorname{Tr}_{q^6/q}$，$e_j(a)$ 为 $a$ 的 6 个共轭的初等对称函数。$\alpha$ 的 12 个共轭是 $\prod_{i<6}(t^2+a_it+n_i)$ 的根，展开得
$$\begin{aligned}
e_1&=\operatorname{Tr}a,\\
e_2&=e_2(a)+\operatorname{Tr}n,\\
e_3&=e_3(a)+\operatorname{Tr}a\operatorname{Tr}n+\operatorname{Tr}(an),\\
e_4&=e_4(a)+e_2(a)\operatorname{Tr}n+e_1(a)\operatorname{Tr}(an)+\operatorname{Tr}(a^2n)+e_2(n),\\
e_5&=e_5(a)+e_3(a)\operatorname{Tr}n+e_2(a)\operatorname{Tr}(an)+e_1(a)\operatorname{Tr}(a^2n)+\operatorname{Tr}(a^3n)+e_1(a)e_2(n)+\operatorname{Tr}(an)\operatorname{Tr}n+\operatorname{Tr}(an^2),
\end{aligned}$$
（从 $k$ 个因子取 $a_i$、从 $j$ 个因子取 $n_i$ 的组合；$e_2(a\setminus j)=e_2(a)+a_je_1(a)+a_j^2$ 等恒等式在特征 2 下化简）。幂和：由 $s_k=\alpha^k+\alpha^{q^6k}$ 满足 $s_k=as_{k-1}+ns_{k-2}$ 得 $s_3=a^3+an$，$s_5=a^5+a^3n+an^2$，故
$$p_1=\operatorname{Tr}a,\qquad p_3=\operatorname{Tr}(a^3+an),\qquad p_5=\operatorname{Tr}(a^5+a^3n+an^2).$$
另外 $e_2(\beta+\beta^q)=\operatorname{Tr}_{q^{12}/q}(\beta^{q+1})+(\operatorname{Tr}\beta)^2$（与 $d=4$ 同样的循环指标计算）。全部恒等式由 `part_a` 对 $q=2$（全体 4095 个元素）与 $q=4$（1500 个抽样）穷举核验，零不符。

**纤维结构**（同第二部分引理 A）：$t^2+at+n$ 在 $\mathbb F_{q^6}$ 上不可约（$a\ne0$，$\operatorname{Tr}_{q^6/2}(n/a^2)=1$）时 2 个 $\alpha$；$a=0$（$\alpha\in\mathbb F_{q^6}$，$n=\alpha^2$）时 1 个；其余 0 个。特别地 $\mathbb F_{q^6}\subset U$（12 个共轭成对出现，奇数幂和为 0），$e_2|_{\mathbb F_{q^6}}=(\operatorname{Tr}\alpha)^2$，$e_4|_{\mathbb F_{q^6}}=e_2(\alpha)^2$。

---

## 2. 定理 5：Legendre 剖面的结构 **[定理]**

$\mathcal G^2$ 上 $c_1=c_3=c_5=0$ 等价于 $p_1=p_3=p_5=0$（Newton：$p_3=c_3+c_1c_2+c_1^3$，$p_5=c_5+(\text{含 }c_1,c_3\text{ 的项})$），故 $\Psi(0,c_2,0,c_4,0)=\#\{\alpha\in U:e_2=c_2,e_4=c_4\}$，而 $\mathcal G^2\cong(\mathbb F_q^2,+)$（$(1+xu^2+yu^4)(1+x'u^2+y'u^4)\equiv1+(x+x')u^2+(y+y'+xx')u^4$——注意 $u^4$ 位并非直接相加，但 $e_4$-坐标仍可用 Witt 坐标线性化；下面的论证直接在 $\alpha$ 上进行，不依赖此）。

**两种对称。** (i) 平移 $\alpha\mapsto\alpha+c$（$c\in\mathbb F_q$）保持 $U$：$\operatorname{Tr}(\alpha+c)=\operatorname{Tr}\alpha$，$\operatorname{Tr}(\alpha+c)^3=\operatorname{Tr}\alpha^3+c(\operatorname{Tr}\alpha)^2$，$\operatorname{Tr}(\alpha+c)^5=\operatorname{Tr}\alpha^5+c(\operatorname{Tr}\alpha)^4$。用 $e_j(\alpha+c)=\sum_i\binom{12-i}{j-i}c^{j-i}e_i$ 并取模 2（$\binom{12}{2}=66,\binom{11}{1}=11$；$\binom{12}4=495,\binom{11}3=165,\binom{10}2=45,\binom91=9$ 均为奇数），在 $U$ 上（$e_1=e_3=0$）：
$$e_2(\alpha+c)=e_2(\alpha),\qquad e_4(\alpha+c)=e_4(\alpha)+c^2e_2(\alpha)+c^4 .$$
(ii) 伸缩 $\alpha\mapsto s\alpha$（$s\in\mathbb F_q^\times$）保持 $U$，$(e_2,e_4)\mapsto(s^2e_2,s^4e_4)$。

**证明。** 令 $F(\lambda)=\sum_{\alpha\in U}\chi(\lambda e_2+e_4)$。由 (i)，$F(\lambda)=\sum_U\chi(\lambda e_2+e_4+c^2e_2+c^4)=\chi(c^4)F(\lambda+c^2)=\chi(c)F(\lambda+c^2)$，故 $F(\lambda)=\chi(\sqrt\lambda)F(0)=\chi(\lambda)Y$。由 (ii)，对 $\nu\ne0$ 取 $s=\nu^{-1/4}$：$\widehat\Psi_H(\mu,\nu):=\sum_U\chi(\mu e_2+\nu e_4)=F(\mu\nu^{-1/2})=\chi(\mu/\sqrt\nu)\,Y$；对 $\mu\ne0$：$\widehat\Psi_H(\mu,0)=\sum_U\chi(\mu e_2)=X$。Fourier 反演：
$$\Psi(0,c_2,0,c_4,0)=\frac1{q^2}\Bigl[A_6+X\bigl(q[c_2=0]-1\bigr)+Y\sum_{\nu\ne0}\chi(\nu c_4)\sum_\mu\chi\bigl(\mu(\nu^{-1/2}+c_2)\bigr)\Bigr]
=q^7+\frac{A_6-q^9-X}{q^2}+\frac Xq[c_2=0]+\frac Yq[c_2\ne0]\chi\Bigl(\frac{c_4}{c_2^2}\Bigr).\ \square$$

**推论 5.1（三值剖面）。** 记 $\varphi=\Psi-q^7$，$B:=(A_6-q^9-X)/q^2$：
- $c_2=0$ 的 $q$ 个类：$\varphi_0=B+X/q$；
- $c_2\ne0$，$\operatorname{Tr}_{q/2}(c_4/c_2^2)=0$ 的 $q(q-1)/2$ 个类：$\varphi_A=B+Y/q$；
- $\operatorname{Tr}(c_4/c_2^2)=1$ 的 $q(q-1)/2$ 个类：$\varphi_B=B-Y/q$。
$$\frac{\operatorname{Var}_{\mathcal G^2}(\Psi)}{q^7}=\frac{q\varphi_0^2+q(q-1)(B^2+Y^2/q^2)}{q^9}.$$
这解释了第二部分 §6 观察到的"恰三个值、$c_2=0$ 的 $q$ 个类主导"。`part_b` 对 $q=2,4,8,16,32$ 的全部 $q^2$ 个 Legendre 类逐一核验了定理 5 的公式与 $\widehat\Psi_H(\mu,\nu)=Y\chi(\mu/\sqrt\nu)$。

**注（一般 $d$）。** 同样的论证对任意 $d$ 成立：$U_d=\{p_j=0,\ j\text{ 奇}\le\ell\}$ 在 $\mathbb F_q$-平移与伸缩下不变，$e_{2j}$ 的平移规律由 $\binom{2d-i}{2j-i}\bmod2$ 决定（Lucas 定理），Legendre 剖面总可归结为 $\lfloor d/2\rfloor$ 个变量上的少数特征和。$d=4$ 时 $U_4=\{p_1=p_3=0\}$，$e_2(\alpha+c)=e_2$，得到的正是 (T3)（$\Psi(0,c,0)$ 在 $c\ne0$ 上常值）。

---

## 3. 定理 6：三个常数的塔归约 **[定理]**

设 $S(a)=\{n\in\mathbb F_{q^6}:\ \operatorname{Tr}(an)=\operatorname{Tr}a^3,\ \operatorname{Tr}(\sqrt a\,n)^2+\operatorname{Tr}(a^3n)=\operatorname{Tr}a^5,\ \operatorname{Tr}_{q^6/2}(a^{-2}n)=1\}$（由引理 C：$p_3=0\iff\operatorname{Tr}(an)=\operatorname{Tr}a^3$；$p_5=0\iff\operatorname{Tr}(an^2)+\operatorname{Tr}(a^3n)=\operatorname{Tr}a^5$，而 $\operatorname{Tr}(an^2)=\operatorname{Tr}(\sqrt a\,n)^2$）。则对 $a\ne0$，$\operatorname{Tr}a=0$：
$$A_6=q^6+2\sum_a|S(a)|,\qquad X=2\sum_a\chi(e_2(a))\sum_{n\in S(a)}\chi(\operatorname{Tr}n),\qquad
Y=G+2\sum_a\sum_{n\in S(a)}\chi\bigl(e_4(a)+e_2(a)\operatorname{Tr}n+\operatorname{Tr}(a^2n)+e_2(n)\bigr),$$
其中 $G=\sum_{n\in\mathbb F_{q^6}}\chi(e_2(n))$（$\mathbb F_{q^6}$ 上 $e_2$ 的 Gauss 和；$X$ 中 $\mathbb F_{q^6}$ 的贡献 $\sum\chi((\operatorname{Tr}\alpha)^2)=0$）。

用特征展开示性函数：$[\operatorname{Tr}(an)=\operatorname{Tr}a^3]=q^{-1}\sum_\lambda\chi(\lambda\operatorname{Tr}(an)+\lambda\operatorname{Tr}a^3)$；$[\operatorname{Tr}(\sqrt an)^2+\operatorname{Tr}(a^3n)=\operatorname{Tr}a^5]=q^{-1}\sum_\mu\chi(\operatorname{Tr}(\sqrt\mu\sqrt a\,n)+\mu\operatorname{Tr}(a^3n)+\mu\operatorname{Tr}a^5)$（$\chi(\mu t^2)=\chi(\sqrt\mu\,t)$）；$[\text{不可约}]=\frac12\sum_{\varepsilon\in\{0,1\}}(-1)^\varepsilon\chi(\varepsilon\operatorname{Tr}(a^{-2}n))$。令
$$\beta'=\lambda a+\sqrt\mu\sqrt a+\mu a^3+\varepsilon a^{-2},\qquad \beta=e_2(a)+a^2+\beta',\qquad w=\chi(\lambda\operatorname{Tr}a^3+\mu\operatorname{Tr}a^5).$$

**引理 D（$e_2$ 的 Gauss 和）。** $e_2$ 在 $\mathbb F_{q^6}$ 上的极化形式为 $B(x,y)=\operatorname{Tr}x\operatorname{Tr}y+\operatorname{Tr}(xy)$，非退化（根基：$y+\operatorname{Tr}y=0\Rightarrow y\in\mathbb F_q\Rightarrow 6y=0$）。对任意 $\beta$，$\sum_n\chi(e_2(n)+\operatorname{Tr}(\beta n))=\chi(e_2(\beta))\,G$。证：$n_0=\beta+\operatorname{Tr}\beta$ 满足 $B(n_0,n)=\operatorname{Tr}(\beta n)$，而 $e_2(\beta+c)=e_2(\beta)+c\operatorname{Tr}\beta+c^2$（$c\in\mathbb F_q$，用 $\binom51,\binom62$ 为奇），代 $c=\operatorname{Tr}\beta$ 得 $e_2(n_0)=e_2(\beta)$。$\square$

于是（$\sum_n\chi(\operatorname{Tr}(\gamma n))=q^6[\gamma=0]$）：
$$\boxed{\begin{aligned}
A_6&=q^6+q^4\sum_{a}\ \sum_{(\lambda,\mu,\varepsilon):\ \beta'=0}(-1)^\varepsilon w,\\
X&=q^4\sum_a\chi(e_2(a))\sum_{(\lambda,\mu,\varepsilon):\ \beta'=1}(-1)^\varepsilon w,\\
Y&=G\Bigl[1+q^{-2}\sum_a\sum_{\lambda,\mu\in\mathbb F_q}\sum_{\varepsilon}(-1)^\varepsilon\chi\bigl(e_4(a)+\lambda\operatorname{Tr}a^3+\mu\operatorname{Tr}a^5+e_2(\beta)\bigr)\Bigr].
\end{aligned}}$$
`part_c` 用这三式独立算出 $q=2,4,8$ 的 $(A_6,X,Y)=(544,96,104),\ (259072,-3072,-64),\ (134447104,294912,-17920)$，与 `part_b` 从 Hayes FFT 提取的值完全一致；同时得到 $G=8,-64,512$，即 **$G=(-1)^{k+1}q^3$**（$k\le3$）。

**$A_6$ 与 $X$ 的方程只在"薄层"上有解。** $\beta'=0$ 且 $\varepsilon=0$：$\lambda a+\sqrt\mu\sqrt a+\mu a^3=0$，除以 $\sqrt a$ 后平方得 $\mu^2a^5+\lambda^2a+\mu=0$（$\mu\ne0$），故 $a$ 的次数 $\le5$，即 $a\in\mathbb F_{q^2}\cup\mathbb F_{q^3}$——与 $d=4$ 的"情形 II 薄层"完全平行（那里是 $a\in\mathbb F_{q^2}$）。$\varepsilon=1$：乘以 $a^2$ 并平方得 $\mu^2a^{10}+\lambda^2a^6+\mu a^5+1=0$；令 $u=\mu a^5+\lambda a^3+1$，它等价于
$$u^2+u=\lambda a^3+1,$$
即点 $(a,u)$ 落在**超奇异椭圆曲线** $E_\lambda:\ y^2+y=\lambda x^3+1$ 的 $\mathbb F_{q^6}$-点上，并满足有理性约束 $\mu=(u+\lambda a^3+1)/a^5\in\mathbb F_q$。由 $u\in\mathbb F_{q^6}$ 得 $\operatorname{Tr}_{q^6/2}(\lambda a^3+1)=0$，故 $\chi(\lambda\operatorname{Tr}a^3)=1$，权 $w=(-1)^{\operatorname{Tr}_{q^6/2}(u)}$。$X$ 的方程 $\beta'=1$ 同理化为 $E_\lambda$ 型曲线上带约束的点计数。**这就是 $A_6,X$ 初等（$\pm$ 若干 $q$ 的幂）的原因**：超奇异曲线的点数是初等的。

**$Y$ 为何不初等。** $Y$ 的归约式中多出 $e_2(\beta)$——$\beta$ 关于 $(\lambda,\sqrt\mu,\varepsilon)$ 是仿射的，$e_2(\beta)$ 因而是 $(\lambda,s=\sqrt\mu)$ 的二次型加上 $s^3B(\sqrt a,a^3)$、$s^4e_2(a^3)$ 等高次项（$\chi(s^4\cdot)$ 可线性化，$\chi(s^3\cdot)$ 不能）：对 $\lambda$ 求和后留下 $\mathbb F_q$ 上的**三次** AS 和，再对 $a\in\mathbb F_{q^6}$ 求和。$\sum_a$ 不再坍缩到薄层，结果是一个 5 维族上的指数和，大小 $q^{9/2}$ 量级且不规则。

---

## 4. 三个常数的值

由 `part_b`（$q=2..32$ 精确）：

| $k$ | $q$ | $A_6-q^9$ | $X$ | $Y$ | $m_k=Y/q^3$ |
|---|---|---|---|---|---|
| 1 | 2 | $32=(q-1)q^5$ | $96=q^6+q^5$ | $104$ | $13$ |
| 2 | 4 | $-3072=-(q-1)q^5$ | $-3072=-q^6+q^5$ | $-64$ | $-1$ |
| 3 | 8 | $229376=(q-1)q^5$ | $294912=q^6+q^5$ | $-17920$ | $-35$ |
| 4 | 16 | $-15728640$ | $-15728640$ | $192512$ | $47$ |
| 5 | 32 | $1040187392$ | $1107296256$ | $24018944$ | $733$ |

**[精确验证 $k\le5$；一般 $k$ 为猜想]**
$$A_6=q^9+(-1)^{k+1}(q-1)q^5,\qquad X=(-1)^{k+1}q^6+q^5,\qquad G=(-1)^{k+1}q^3 .$$
等价表述：$U$ 上 $e_2$ 的分布为 $\#\{\alpha\in U:e_2=x\}=q^8-2q^4[k\text{ 奇}]$（$x\ne0$），$=q^8+(q-1)q^4(q+2)$ 或 $q^8-(q-1)q^5$（$x=0$，$k$ 奇/偶）。
证明路线（未完成）：(a) 用 §3 的薄层刻画，把 $A_6,X$ 写成 $E_\lambda(\mathbb F_{q^6})$ 上满足 $\mathbb F_q$-有理性约束的点的带权计数，再用 $y^2+y=x^3+c$ 的超奇异性（Frobenius 的 $\pi^{2}$ 为 $\pm q^{?}\cdot$单位根）给出闭式；或 (b) Hilbert 90：$\alpha=\beta+\beta^q$，$\operatorname{Tr}\alpha^3=\operatorname{Tr}(\beta^2\theta\beta)$，$\operatorname{Tr}\alpha^5=\operatorname{Tr}(\beta^4\theta\beta)$（$\theta=\sigma+\sigma^{-1}$），$A_6=q^9+q^{-3}\sum_{(\lambda_3,\lambda_5)\ne0}W(\lambda_3,\lambda_5)$，$W$ 为 $\mathbb F_2^{12k}$ 上二次型的指数和，其根基为 $\theta^{-1}(\ker M\cap V)$，$M(y)=\lambda_3y^2+\sqrt{\lambda_3y}+\lambda_5y^4+(\lambda_5y)^{1/4}$，$V=\ker\operatorname{Tr}_{q^{12}/q^2}$；$\lambda_5=0$ 的部分已由 Carlitz（第二部分 §3）给出 $W(\lambda_3,0)=-2q^7$，剩余 $\sum_{\lambda_5\ne0}W=(q-1)q^7\bigl((-1)^{k+1}q+2\bigr)$ 需五次 Gauss 和（$\mathbb F_{16}\ni\mu_5$，纯 Gauss 和）与 Arf 不变量。两条路线都是机械的但篇幅可观。

**$Y$。** $m_k=13,-1,-35,47,733$：$|m_k|\le4q^{3/2}$，无明显闭式（$733$ 为素数）。$Y$ 是 $U$ 上四次型 $e_4$ 的特征和，可写成定理 6 中的 5 维族指数和；期望它是某个曲面/曲线族 $L$-函数的 Frobenius 迹，但本报告未能识别。$Y$ 的符号规律：$k=1,4,5$ 为正，$k=2,3$ 为负，不随 $k$ 奇偶交替。

---

## 5. $d=6$ 的最终闭式

把 §4 的值代入推论 5.1（$B=(A_6-q^9-X)/q^2=-2q^3[k\text{ 奇}]$）：

$$\varphi_0=\begin{cases}q^5+q^4-2q^3&k\text{ 奇}\\-q^5+q^4&k\text{ 偶}\end{cases},\qquad
\varphi_{A,B}=\begin{cases}-2q^3\pm q^2m_k&k\text{ 奇}\\ \pm q^2m_k&k\text{ 偶}\end{cases}$$
（第二部分 §6 表中 $c_2=0$ 列的闭式由此得到证明——模 $A_6,X$ 的闭式——其余两列即 $\varphi_{A,B}$：$q=8$：$-3264=-1024-2240$，$1216=-1024+2240$，$2240=64\cdot35$ ✓）。

$$\boxed{\frac{\operatorname{Var}_{\mathcal G^2}(\Psi)}{q^7}=\begin{cases}\dfrac{(q-1)^2(q+2)^2}{q^2}+\dfrac{(q-1)(4q^2+m_k^2)}{q^4}&k\text{ 奇}\\[2ex](q-1)^2+\dfrac{(q-1)m_k^2}{q^4}&k\text{ 偶}\end{cases}}$$
`part_b` 核验：$15.5625,\ 9.011719,\ 79.093506,\ 225.505600,\ 1100.884360$ 与 FFT 值完全相等。

**推论（偶数 $d$ 标度律，$d=6$）。** 若 $A_6,X$ 的闭式对所有 $k$ 成立且 $m_k=o(q^{5/2})$（数据：$|m_k|\le4q^{3/2}$），则
$$\lim_{q\to\infty}\frac{\operatorname{Var}_{\mathcal G^2}(\Psi)}{q^{9}}=1=C_6,\qquad \frac{\operatorname{Var}_{\mathcal G^2}}{\operatorname{Var}_{\mathcal G}}\sim\frac{q^2}{4}$$
（$\operatorname{Var}_{\mathcal G}/q^7\to4$ 为 KR 值，$q=16,32$ 时 $4.35,4.21$）。与 $d=4$（$C_4=1$，比率 $q^2/3$）合起来：**偶数 $d$ 的 Legendre 方差异动是 $q^2$ 律，主项由 $f_{d-1}=0$（$c_2=0$）的 $q^{d/2-1}$ 个区间贡献，其偏差为 $\pm q^{d-1}(1+O(1/q))$，符号 $(-1)^{k+1}$（$d=4$ 时 $-(q^4-q^3)$ 恒负，$d=6$ 时为 $(-1)^{k+1}q^5+q^4-\cdots$）。**

---

## 6. 与 $d=4$ 的对照及尚未完成的部分

| | $d=4$ | $d=6$ |
|---|---|---|
| $U_d$ | $\{p_1=p_3=0\}\subset\mathbb F_{q^8}$ | $\{p_1=p_3=p_5=0\}\subset\mathbb F_{q^{12}}$ |
| Legendre 坐标 | $e_2$ | $(e_2,e_4)$ |
| 平移作用 | $e_2$ 不变 | $e_2$ 不变，$e_4\mapsto e_4+c^2e_2+c^4$ |
| 剖面 | 2 个值（$c=0$ / $c\ne0$） | 3 个值（$c_2=0$ / $\operatorname{Tr}(c_4/c_2^2)=0,1$） |
| 常数 | $A=\#U$，$\Psi(\mathrm{id})$：均初等 | $A_6,X$ 初等（猜想），$Y$ 非初等 |
| 薄层 | $a\in\mathbb F_{q^2}\setminus\mathbb F_q$ | $a\in\mathbb F_{q^2}\cup\mathbb F_{q^3}$ 及曲线 $y^2+y=\lambda x^3+1$ |

未完成：(1) $A_6,X$（等价地 $G$）闭式的证明（§4 两条路线）；(2) $Y$ 的算术解释。(2) 是真正的开放点：按 Katz–Sarnak 的图像，$d=6$ 有 4 个零点，$e_4$ 这一"第二 Witt 坐标"对应的特征和不再退化到有限阶 Frobenius，$m_k$ 应是某个 $\ell$-进层的 Frobenius 迹；确定该层（及其秩，由 $|m_k|\lesssim 4q^{3/2}$ 推测在 $\mathbb F_q$ 上为曲面型或秩 $\le 8$ 的曲线族）需要 Deligne 的理论，超出本报告范围。
