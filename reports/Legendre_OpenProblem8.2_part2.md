# 第二部分：$d=4$ 的完整证明（任意 $q=2^k$）与偶数 $d$ 标度律的判定

日期：2026-10-07。承接 `Legendre_OpenProblem8.2_report.md`（下称"第一部分"）。
新代码：`legendre_d4_formula.py`（闭式及其全部校验）、`legendre_q16_hayes.py`（调用论文 `hayes_fft_fq.py` / `hayes_exact.py`，扩展到 $q=32,64,128$）；输出 `hayes_q16_output.txt`、`hayes_d6_profile.txt`。

**本部分把第一部分 §3.6 的全部遗留项关闭**：下列结论对**所有** $q=2^k$ 成立并给出完整证明（不依赖任何未核对的文献定理；Hasse–Davenport 关系仅在一个可选的旁注中用到）。

| 结论（$d=4$，$n=8$） | 状态 |
|---|---|
| (T1) $c_1\ne0\Rightarrow\Psi(g)=q^5$ | **定理**（§2） |
| 引理 1：$\#\{\alpha\in\mathbb F_{q^8}:e_1=e_2=e_3=0\}=q^5-q^4+q^3$ | **定理**（§2） |
| (T3) $\Psi(0,c,0)-q^5=(-1)^k q^3$（$c\ne0$） | **定理**（§2） |
| $\operatorname{Var}_{\mathcal G^2}(\Psi)/q^5=q(q-1)$ | **定理**（§2） |
| 引理 2：$S(\lambda_1,\lambda_3)$ 与 $\lambda_1\in\mathbb F_q$ 无关；Carlitz 值 | **定理**（§3） |
| $\operatorname{Var}_{\mathcal G}(\Psi)/q^5=3-3/q$，比率 $=q^2/3$ | **定理**（§4） |
| $\rho(\varepsilon)$：$q-1$ 个为 $1$，其余 $=\frac{q-3}{3(q-1)}$（$k$ 奇） | **定理**（§4.4） |
| 偶数 $d$ 标度律：$d=6$ 时 $\operatorname{Var}_{\mathcal G^2}/q^{7}\asymp q^2$，**不是** $q^4$ | **精确数据判定**（$q\le32$，§6） |

全部闭式均与独立计算逐类核对：$q=2,4,8$ 对第一部分的根枚举，$q=16,32,64,128$ 对论文的 Hayes 群 FFT（$2^{21}$ 个类逐一相等），$\operatorname{Var}_{\mathcal G}$ 恒等式另以 $O(q)$ 算法核到 $q=2^{22}$。

---

## 1. 二次塔参数化 **[引理 A]**

设 $\sigma_4:\alpha\mapsto\alpha^{q^4}$ 为 $\mathbb F_{q^8}/\mathbb F_{q^4}$ 的非平凡自同构。对 $\alpha\in\mathbb F_{q^8}$ 令
$$a=\alpha+\alpha^{q^4}\in\mathbb F_{q^4},\qquad n=\alpha^{1+q^4}\in\mathbb F_{q^4},$$
则 $\alpha,\alpha^{q^4}$ 是 $t^2+at+n$ 的两根，而 $\alpha$ 的 8 个 $\mathbb F_q$-共轭是 4 个 $\mathbb F_q$-共轭二次式 $t^2+a_it+n_i$（$a_i=a^{q^i},n_i=n^{q^i}$，$i<4$）的根。展开 $\prod_{i<4}(t^2+a_it+n_i)$ 得（记 $\operatorname{Tr}=\operatorname{Tr}_{q^4/q}$，$e_j(a)$ 为 $a$ 的四个共轭的初等对称函数）
$$e_1(\alpha)=\operatorname{Tr}(a),\qquad e_2(\alpha)=e_2(a)+\operatorname{Tr}(n),\qquad e_3(\alpha)=e_3(a)+\operatorname{Tr}(a)\operatorname{Tr}(n)+\operatorname{Tr}(an).$$
（$e_3$：从四个因子中取三个 $a_i$，或取一个 $n_j$ 与另一个因子的 $a_i$：$\sum_{i\ne j}a_in_j=\operatorname{Tr}a\operatorname{Tr}n-\operatorname{Tr}(an)$。）

**纤维结构。** $(a,n)\in\mathbb F_{q^4}^2$ 的原像：
- 若 $t^2+at+n$ 在 $\mathbb F_{q^4}$ 上不可约（$\iff a\ne0$ 且 $\operatorname{Tr}_{q^4/2}(n/a^2)=1$），恰有 2 个 $\alpha$；
- 若 $a=0$（$\iff\alpha\in\mathbb F_{q^4}$），$n=\alpha^2$，恰有 1 个 $\alpha$；
- 若 $a\ne0$ 且二次式可约，原像为空（其根 $\beta\in\mathbb F_{q^4}$ 满足 $\beta+\beta^{q^4}=0\ne a$）。

于是对任意类 $(c_1,c_2,c_3)$：
$$\Psi(c_1,c_2,c_3)=2\,\#\Bigl\{(a,n):\ a\ne0,\ \operatorname{Tr}_{q^4/2}\tfrac n{a^2}=1,\ (e_1,e_2,e_3)=(c_1,c_2,c_3)\Bigr\}+\#\{\alpha\in\mathbb F_{q^4}:(e_1,e_2,e_3)=(c_1,c_2,c_3)\}.\tag{1.1}$$

两个事实反复使用：(i) 对 $x\in\mathbb F_{q^4}$，$\operatorname{Tr}_{q^8/2}(x)=\operatorname{Tr}_{q^4/2}(2x)=0$；对 $x\in\mathbb F_{q^2}$，$\operatorname{Tr}_{q^4/q}(x)=0$；对 $x\in\mathbb F_q$，$\operatorname{Tr}_{q^4/q}(x)=4x=0$。(ii) 对 $a\in\mathbb F_{q^4}$，$\mathbb F_q$-线性泛函 $n\mapsto\operatorname{Tr}(n)$ 与 $n\mapsto\operatorname{Tr}(an)$ 线性无关 $\iff a\notin\mathbb F_q$；而 $\mathbb F_2$-泛函 $n\mapsto\operatorname{Tr}_{q^4/2}(a^{-2}n)$ 属于它们生成的 $\mathbb F_2$-空间 $\{n\mapsto\operatorname{Tr}_{q^4/2}(cn):c\in\mathbb F_q+\mathbb F_qa\}$ $\iff a^{-2}\in\mathbb F_q+\mathbb F_qa\iff a\in\mathbb F_{q^2}$（若 $a\notin\mathbb F_{q^2}$，$a^{-2}=x+ya$ 将给出 $a$ 的次数 $\le3$ 的方程，与 $[\mathbb F_q(a):\mathbb F_q]=4$ 矛盾；若 $a\in\mathbb F_{q^2}\setminus\mathbb F_q$，$\mathbb F_q+\mathbb F_qa=\mathbb F_{q^2}\ni a^{-2}$）。

---

## 2. 主定理：$\Psi$ 在 $\mathcal G_3$ 上的闭式

**定理 1.** 对任意 $q=2^k$，记 $\chi(x)=(-1)^{\operatorname{Tr}_{q/2}(x)}$，$W=\{w\in\mathbb F_q:\operatorname{Tr}_{q/2}w=1\}$。则
$$\Psi(c_1,c_2,c_3)=q^5\quad(c_1\ne0),$$
$$\Psi(0,c_2,c_3)=q^5+q^3[c_3=0]+2q^3\,I(c_2,c_3)-2q^2\,K(c_2,c_3),$$
$$I(c_2,c_3)=\bigl[c_2c_3\ne0,\ \operatorname{Tr}_{q/2}(c_2^3/c_3^2)=1\bigr],\qquad
K(c_2,c_3)=\sum_{w\in W}\sum_{v\in\mathbb F_q}\chi\Bigl(\frac{c_2v^2}{w}+\frac{\sqrt{c_2}\,v}{w}+\frac{c_3v^3}{w^2}\Bigr).$$

**证明.** 按 (1.1) 分四种情形计算贡献。固定目标 $(c_1,c_2,c_3)$。

*情形 I：$a\notin\mathbb F_{q^2}$。* 条件 $e_1=c_1$ 即 $\operatorname{Tr}a=c_1$；给定这样的 $a$，对 $n$ 的条件为两条 $\mathbb F_q$-线性方程 $\operatorname{Tr}(n)=c_2+e_2(a)$，$\operatorname{Tr}(an)=c_3+e_3(a)+c_1(c_2+e_2(a))$，由事实 (ii) 无关，解集为 $q^2$ 元仿射子空间；不可约条件 $\operatorname{Tr}_{q^4/2}(a^{-2}n)=1$ 是其上非常值的仿射 $\mathbb F_2$-泛函（事实 (ii)），恰取一半：$q^2/2$ 个 $n$，每个 2 个根。$a$ 的个数：$\#\{a\in\mathbb F_{q^4}:\operatorname{Tr}a=c_1\}=q^3$，减去 $\mathbb F_{q^2}$ 中的（它们的迹全为 0）：$q^3-q^2[c_1=0]$。贡献 $2\cdot\frac{q^2}2(q^3-q^2[c_1=0])=q^5-q^4[c_1=0]$。

*情形 II：$a\in\mathbb F_{q^2}\setminus\mathbb F_q$。* 此时 $\operatorname{Tr}a=0$，只对 $c_1=0$ 有贡献。$a$ 的四个共轭为 $a,a^q,a,a^q$，令 $T=a+a^q\in\mathbb F_q^\times$（$T=0\iff a\in\mathbb F_q$），$N=a^{1+q}$，则 $e_2(a)=a^2+a^{2q}+4a^{1+q}=T^2$，$e_3(a)=2(a^2a^q+aa^{2q})=0$。对 $n$ 的条件：$\operatorname{Tr}(n)=c_2+T^2$，$\operatorname{Tr}(an)=c_3$（$q^2$ 个解）。不可约条件：由事实 (ii) 写 $a^{-2}=x+ya$，则 $\operatorname{Tr}_{q^4/2}(a^{-2}n)=\operatorname{Tr}_{q/2}\bigl(x\operatorname{Tr}(n)+y\operatorname{Tr}(an)\bigr)=\operatorname{Tr}_{q/2}\bigl(x(c_2+T^2)+yc_3\bigr)$ 在解集上为常数。由 $a^2=Ta+N$ 得 $a^{-1}=(a+T)/N$，$a^{-2}=(N+T^2+Ta)/N^2$，即 $x=(N+T^2)/N^2$，$y=T/N^2$。而 $\operatorname{Tr}_{q/2}(xT^2)=\operatorname{Tr}_{q/2}(T^2/N+T^4/N^2)=\operatorname{Tr}_{q/2}(z+z^2)=0$（$z=T^2/N$）。故该情形贡献 $2q^2$ 乘以满足
$$\tau(a):=\operatorname{Tr}_{q/2}\Bigl(\frac{(N+T^2)c_2+Tc_3}{N^2}\Bigr)=1$$
的 $a\in\mathbb F_{q^2}\setminus\mathbb F_q$ 的个数。这样的 $a$ 与不可约二次式 $t^2+Tt+N$（$T\ne0$，$\operatorname{Tr}_{q/2}(N/T^2)=1$）二对一对应。令 $N=T^2w$（$w\in W$），$v=1/T$：
$$\tau=\operatorname{Tr}_{q/2}\Bigl(\frac{(w+1)c_2}{T^2w^2}+\frac{c_3}{T^3w^2}\Bigr)=\operatorname{Tr}_{q/2}\Bigl(\frac{c_2v^2}{w}+\frac{c_2v^2}{w^2}+\frac{c_3v^3}{w^2}\Bigr),\qquad \operatorname{Tr}_{q/2}\frac{c_2v^2}{w^2}=\operatorname{Tr}_{q/2}\frac{\sqrt{c_2}\,v}{w}.$$
贡献 $=4q^2M$，$M=\#\{(v,w)\in\mathbb F_q^\times\times W:\tau=1\}$。用 $\chi$ 写 $[\tau=1]=\frac{1-\chi(\cdot)}2$ 并补上 $v=0$ 项（每个 $w$ 贡献 $\chi(0)=1$）：$M=\frac{q(q-1)}4-\frac12\bigl(K-\frac q2\bigr)=\frac{q^2}4-\frac K2$，故贡献 $=q^4-2q^2K(c_2,c_3)$。

*情形 III：$a\in\mathbb F_q^\times$。* $\operatorname{Tr}a=0$（需 $c_1=0$），$e_2(a)=6a^2=0$，$e_3(a)=4a^3=0$。条件：$\operatorname{Tr}(n)=c_2$，$a\operatorname{Tr}(n)=c_3$，即 $c_3=ac_2$；不可约条件 $\operatorname{Tr}_{q^4/2}(a^{-2}n)=\operatorname{Tr}_{q/2}(a^{-2}c_2)=1$。若 $c_2=0$ 则 $c_3=0$ 且迹为 0，无解；若 $c_2\ne0$ 则 $a=c_3/c_2$ 需非零，条件化为 $\operatorname{Tr}_{q/2}(c_2^3/c_3^2)=1$，此时 $n$ 有 $q^3$ 个（一条 $\mathbb F_q$-条件），各 2 根。贡献 $2q^3I(c_2,c_3)$。

*情形 IV：$a=0$，$\alpha\in\mathbb F_{q^4}$。* $e_1=0$，$e_2=\operatorname{Tr}(\alpha^2)=(\operatorname{Tr}\alpha)^2$，$e_3=0$。贡献 $[c_1=0][c_3=0]\cdot\#\{\alpha\in\mathbb F_{q^4}:(\operatorname{Tr}\alpha)^2=c_2\}=q^3[c_1=0][c_3=0]$（平方是双射）。

相加：$c_1\ne0$ 时只有情形 I，$\Psi=q^5$；$c_1=0$ 时 $\Psi=q^5-q^4+q^4-2q^2K+2q^3I+q^3[c_3=0]$。$\square$

**推论 1（T1）.** $\varphi=\Psi-q^5$ 支撑在 $\mathcal G_0=\{c_1=0\}$ 上。因此 $\rho(\varepsilon_{\lambda_1,0})=1$ 对全部 $q-1$ 个仅依赖 $p_1$ 的二次特征，且 $\rho(\varepsilon_{\lambda_1,\lambda_3})$ 与 $\lambda_1$ 无关（第一部分 §3.4 的论证现在无条件成立）。

**推论 2（引理 1）.** $K(0,0)=\sum_{w\in W}\sum_v\chi(0)=q^2/2$，$I(0,0)=0$，故
$$\Psi(0,0,0)=\#\{\alpha\in\mathbb F_{q^8}:e_1=e_2=e_3=0\}=q^5+q^3-q^4=q^5-q^4+q^3 .$$
等价地：首项系数 $a_7=a_6=a_5=0$ 的不可约八次首一多项式恰有 $(q^5-q^4)/8=q^4(q-1)/8$ 个（素幂部分 $E=q^3$，第一部分 §3.1）。

**推论 3（T3）.** 对 $c\ne0$：令 $y=\sqrt c\,v$，则 $\chi(cv^2/w+\sqrt c\,v/w)=\chi((y^2+y)/w)$，而 $y\mapsto y^2+y$ 二对一映满 $H_0=\{\operatorname{Tr}z=0\}$，故 $\sum_y\chi((y^2+y)/w)=2\sum_{z\in H_0}\chi(z/w)=q[w=1]$（$w\ne1$ 时 $z\mapsto\operatorname{Tr}(z/w)$ 在 $H_0$ 上非平凡）。于是 $K(c,0)=q[1\in W]=q[k\text{ 奇}]$，
$$\Psi(0,c,0)=q^5+q^3-2q^3[k\text{ 奇}]=\begin{cases}q^5+q^3&k\text{ 偶}\\ q^5-q^3&k\text{ 奇}\end{cases}\qquad(c\ne0).$$
并得 $A=\sum_c\Psi(0,c,0)=\#\{\alpha:\operatorname{Tr}\alpha=\operatorname{Tr}\alpha^3=0\}=q^6$（$k$ 偶）或 $q^6-2q^4+2q^3$（$k$ 奇）。

**定理 2（Legendre 方差闭式，所有 $k$）.**
$$\sum_{g\in\mathcal G^2}\bigl(\Psi(g)-q^5\bigr)^2=(q^4-q^3)^2+(q-1)q^6=q^7(q-1),\qquad
\boxed{\frac{\operatorname{Var}_{\mathcal G^2}(\Psi)}{q^5}=q(q-1)}\quad\text{对所有 }q=2^k .$$
证明即推论 2、3。$\square$

**数值核对。** `legendre_d4_formula.py`：$k=1,\dots,6$ 的闭式与第一部分（$q\le8$，根枚举）完全一致；`legendre_q16_hayes.py`：$q=16,32,64,128$ 的闭式与论文 FFT 算法（$q\le64$ 另经模两素数 CRT 精确化）在全部 $q^3$ 个类上逐一相等（`IDENTICAL, max |diff| = 0`）。

---

## 3. 引理 2：三次 Artin–Schreier 和对 $\lambda_1\in\mathbb F_q$ 的无关性

记 $\psi_0$ 为 $\mathbb F_{q^8}$ 的典范加法特征，$S(b,a)=\sum_{x\in\mathbb F_{q^8}}\psi_0(ax^3+bx)$。

**引理 2.** 对 $a\in\mathbb F_q^\times$，$b\in\mathbb F_q$：$S(b,a)=S(0,a)$。

**证明.** 对任意 $c\in\mathbb F_{q^8}$，代换 $x\mapsto x+c$：
$$\operatorname{Tr}_{q^8/2}\bigl(a(x+c)^3+b(x+c)\bigr)=\operatorname{Tr}\bigl(ax^3+(b+ac^2+\sqrt{ac})x\bigr)+\operatorname{Tr}(ac^3+bc),$$
其中用了 $\operatorname{Tr}(acx^2)=\operatorname{Tr}((\sqrt{ac}\,x)^2)=\operatorname{Tr}(\sqrt{ac}\,x)$。因此
$$S(b,a)=\psi_0(ac^3+bc)\,S\bigl(b+T_a(c),a\bigr),\qquad T_a(c):=ac^2+\sqrt{ac}\ \ (\mathbb F_2\text{-线性}).$$
**断言**：存在 $c\in\mathbb F_{q^4}$ 使 $T_a(c)=b$。若然，则 $ac^3+bc\in\mathbb F_{q^4}$，由事实 (i) 其绝对迹为 0，于是 $S(0,a)=S(T_a(c),a)=S(b,a)$（取上式中 $b\to0$）。

断言的证明：$T_a$ 保持 $\mathbb F_{q^4}$（$a\in\mathbb F_q$）。关于迹型 $\langle u,v\rangle=\operatorname{Tr}_{q^4/2}(uv)$ 计算伴随：$\operatorname{Tr}(T_a(c)u)=\operatorname{Tr}(ac^2u)+\operatorname{Tr}(\sqrt{ac}\,u)=\operatorname{Tr}\bigl(c\,\sqrt{au}\bigr)+\operatorname{Tr}\bigl(c\,au^2\bigr)$（用 $\operatorname{Tr}(x^2)=\operatorname{Tr}(x)$ 两次），故 $T_a^*(u)=\sqrt{au}+au^2=T_a(u)$：$T_a$ 自伴，$\operatorname{im}T_a=(\ker T_a)^\perp$。$\ker T_a|_{\mathbb F_{q^4}}=\{0\}\cup\{u\in\mathbb F_{q^4}:u^3=a^{-1}\}$。对这样的 $u\ne0$，$u^3\in\mathbb F_q$，故 $\zeta:=u^{q-1}$ 满足 $\zeta^3=1$。若 $\zeta=1$，$u\in\mathbb F_q$，$\operatorname{Tr}_{q^4/2}(bu)=\operatorname{Tr}_{q/2}(4bu)=0$。若 $\zeta\ne1$：$k$ 偶时 $\zeta^q=\zeta$，$u^{q^4}=\zeta^4u\ne u$，与 $u\in\mathbb F_{q^4}$ 矛盾，此情形不出现；$k$ 奇时 $\zeta^q=\zeta^2$，$u^{q^2}=\zeta^{1+q}u=\zeta^3u=u$，即 $u\in\mathbb F_{q^2}$，于是 $\operatorname{Tr}_{q^4/2}(bu)=\operatorname{Tr}_{q^2/2}(2bu)=0$。总之 $b\perp\ker T_a$，即 $b\in\operatorname{im}(T_a|_{\mathbb F_{q^4}})$。$\square$

**推论（Carlitz 值的两种推导）.** (a) 由定理 1 的 Fourier 变换（§4，(4.3)）：$S(0,\nu)=\sum_{\lambda_1}S(\lambda_1,\nu)/q=\sum_{\alpha:\operatorname{Tr}\alpha=0}\psi_0(\nu\alpha^3)=\widehat\Psi(0,\nu)$，而 $\widehat\Psi(0,\nu)=-2q^4$（$\nu$ 在 $\mathbb F_{q^8}$ 中为立方：$k$ 奇时恒成立，$k$ 偶时当且仅当 $\nu$ 是 $\mathbb F_q$ 中的立方）或 $+q^4$（否则）。这给出 Carlitz 1979 定理在 $\mathbb F_{2^{8k}}$、$a\in\mathbb F_q^\times$ 情形的一个初等证明。(b) 旁注：经典路线用 $\mathbb F_4$ 上的三次 Gauss 和 $g_4=\sum_{y\in\mathbb F_4^\times}\chi_3(y)(-1)^{\operatorname{Tr}y}=1-\omega-\omega^2=2$ 与 Hasse–Davenport 提升 $g_{4^m}=(-1)^{m-1}2^m$，得 $\sum_{x\in\mathbb F_{2^{2m}}}\psi_0(ax^3)=2\operatorname{Re}(\bar\chi_3(a))\,(-1)^{m-1}2^m$，即立方 $(-1)^{m+1}2^{m+1}$、非立方 $(-1)^m2^m$；$m=4k$ 时与 (a) 一致。第一部分 §3.6 中 (R2) 的所有数值（$W(\nu)=-64;\ -2048,1024,1024;\ 7\times(-65536)$）均由此解释。

---

## 4. 定理 3：$\operatorname{Var}_{\mathcal G}(\Psi)/q^5=3-3/q$（所有 $k$）

### 4.1 $\mathcal G_0\cong\mathbb F_q^2$ 上的 Fourier 变换
$\mathcal G_0=\{(0,x,y)\}$ 在乘法下同构于 $(\mathbb F_q^2,+)$。令 $\widehat\Psi(\mu,\nu)=\sum_{x,y}\Psi(0,x,y)\chi(\mu x+\nu y)$。由推论 1 与 Parseval，
$$\sum_{g\in\mathcal G}\varphi(g)^2=\sum_{x,y}\varphi(0,x,y)^2=\frac1{q^2}\sum_{(\mu,\nu)\ne(0,0)}|\widehat\Psi(\mu,\nu)|^2\tag{4.1}$$
（$\sum_{\mathcal G_0}\varphi=\#\{\operatorname{Tr}\alpha=0\}-q^2\cdot q^5=0$）。对定理 1 的三项分别变换。

*$q^3[y=0]$ 项*：$q^3\sum_x\chi(\mu x)=q^4[\mu=0]$。

*$K$ 项*：用恒等式 $\sum_{x\in\mathbb F_q}\chi(x\alpha+\sqrt x\,\beta)=\sum_s\chi(s(\sqrt\alpha+\beta))=q[\sqrt\alpha=\beta]$（$x=s^2$）得
$$\widehat K(\mu,\nu)=\sum_{w\in W}\sum_v\Bigl[\sum_y\chi\bigl(y(\tfrac{v^3}{w^2}+\nu)\bigr)\Bigr]\Bigl[\sum_x\chi\bigl(x(\tfrac{v^2}{w}+\mu)+\sqrt x\,\tfrac vw\bigr)\Bigr]
=q^2\,\#\Bigl\{(v,w)\in\mathbb F_q\times W:\ v^3=\nu w^2,\ v^2(w+1)=\mu w^2\Bigr\}.$$
（第二个条件来自 $\sqrt{v^2/w+\mu}=v/w\iff\mu=v^2/w^2+v^2/w$。）

*$I$ 项*：$I=[x,y\ne0]\frac{1-\chi(x^3/y^2)}2$，
$$\widehat I(\mu,\nu)=\tfrac12\sum_{x,y\ne0}\chi(\mu x+\nu y)-\tfrac12J(\mu,\nu),\qquad
J(\mu,\nu)=\sum_{x,y\ne0}\chi\Bigl(\frac{x^3}{y^2}+\mu x+\nu y\Bigr)\overset{x=yt}{=}\sum_{t\ne0}\sum_{y\ne0}\chi\bigl(y(t^3+\mu t+\nu)\bigr)=q\,R^*(\mu,\nu)-(q-1),$$
$R^*(\mu,\nu)=\#\{t\in\mathbb F_q^\times:t^3+\mu t+\nu=0\}$。

### 4.2 三个频率区
(a) $\nu=0,\ \mu\ne0$：$R^*=\#\{t\ne0:t^2=\mu\}=1$，$J=1$，$\widehat I=\frac12(-(q-1)-1)=-\frac q2$；$\widehat K=q^2\#\{(v,w):v=0,\ 0=\mu w^2\}=0$。故
$$\widehat\Psi(\mu,0)=2q^3\cdot(-\tfrac q2)=-q^4 .\tag{4.2}$$
(b) $\mu=0,\ \nu\ne0$：$R^*=\#\{t:t^3=\nu\}=:R_0(\nu)$（$k$ 奇为 1；$k$ 偶为 $3[\nu\text{ 立方}]$），$\widehat I=-\frac q2R_0$；$\widehat K=q^2\#\{(v,w):v^3=\nu w^2,\ w=1\}=q^2[1\in W]\#\{v:v^3=\nu\}=q^2[k\text{ 奇}]$。故
$$\widehat\Psi(0,\nu)=q^4-q^4R_0(\nu)-2q^4[k\text{ 奇}]=\begin{cases}-2q^4&k\text{ 奇}\\ -2q^4\ (\nu\text{ 立方}),\ \ +q^4\ (\nu\text{ 非立方})&k\text{ 偶}\end{cases}\tag{4.3}$$
(c) $\mu\nu\ne0$：由 $v^3=\nu w^2$ 与 $v^2(w+1)=\mu w^2$ 相除得 $v=\nu(w+1)/\mu$，代回得 $\nu^2(w+1)^3=\mu^3w^2$。令 $j:=\mu^3/\nu^2\in\mathbb F_q^\times$，则
$$\widehat K(\mu,\nu)=q^2L_j,\qquad L_j:=\#\{w\in W:\ (w+1)^3=jw^2\}.$$
对 $R^*$：代换 $t=(\nu/\mu)s$ 把 $t^3+\mu t+\nu=0$ 化为 $s^3+j(s+1)=0$；再令 $s=1+1/w$（$s\ne1$ 自动成立），$s^3+js+j=0\iff w^3+(1+j)w^2+w+1=0\iff(w+1)^3=jw^2$。故 $R^*(\mu,\nu)=R_j:=\#\{w\in\mathbb F_q:(w+1)^3=jw^2\}$，且 $L_j$ 恰是其中 $\operatorname{Tr}w=1$ 的根数；记 $L'_j=R_j-L_j$。于是 $\widehat I=\frac q2(1-R_j)$，
$$\widehat\Psi(\mu,\nu)=q^4(1-R_j)-2q^4L_j=q^4\bigl(1-3L_j-L'_j\bigr).\tag{4.4}$$
对每个 $j\in\mathbb F_q^\times$，$\#\{(\mu,\nu):\mu^3/\nu^2=j\}=q-1$（对每个 $\nu$ 数 $\mu^3=j\nu^2$ 的解：$k$ 奇恒 1 个；$k$ 偶时恰 $(q-1)/3$ 个 $\nu$ 使 $j\nu^2$ 为立方，各 3 个解）。

代入 (4.1)：
$$\sum_{\mathcal G}\varphi^2=q^6(q-1)\Bigl[1+\kappa_k+\Sigma\Bigr],\qquad \kappa_k=\begin{cases}4&k\text{ 奇}\\ \tfrac13\cdot4+\tfrac23\cdot1=2&k\text{ 偶}\end{cases},\qquad \Sigma:=\sum_{j\in\mathbb F_q^\times}(1-3L_j-L'_j)^2 .\tag{4.5}$$

### 4.3 关键恒等式
**引理 B.** $\Sigma=3q-4+(-1)^k$，即 $3q-5$（$k$ 奇）、$3q-3$（$k$ 偶）。

**证明.** 令 $C_j(x)=x^3+(1+j)x^2+x+1$，$c(x)=2-\chi(x)$（$\operatorname{Tr}x=1\Rightarrow3$，$=0\Rightarrow1$），则 $3L_j+L'_j=\sum_{C_j(x)=0}c(x)$。每个 $x\in\mathbb F_q\setminus\{0,1\}$ 恰是一个 $C_j$（$j=(x+1)^3/x^2\ne0$）的根，$x\in\{0,1\}$ 不是任何 $C_j$（$j\ne0$）的根。于是
$$\Sigma=(q-1)-2\sum_{x\ne0,1}c(x)+\sum_{x\ne0,1}c(x)^2+\mathrm{OFF},\qquad \mathrm{OFF}=\sum_{\substack{x\ne x'\\ \text{同一 }C_j\text{ 的根}}}c(x)c(x').$$
用 $\sum_{x\ne0,1}\chi(x)=-1-(-1)^k$：$\sum c=2(q-2)+1+(-1)^k$，$\sum c^2=\sum(5-4\chi)=5(q-2)+4+4(-1)^k$，故 $\Sigma=2q-1+2(-1)^k+\mathrm{OFF}$。

对角外项：$x\ne x'$ 是同一 $C_j$ 的两根 $\iff$ $x,x'\notin\{0,1\}$，$x\ne x'$，且第三根 $x''=1/(xx')$ 使 $e_2=1$，即
$$x^2x'^2+xx'+x+x'=0 .\tag{4.6}$$
（反之，(4.6) 的解给出 $e_2=e_3=1$ 的三次式 $C_j$，$j=x+x'+x''+1\ne0$，因为 $j=0$ 的 $C_0=(x+1)^3$ 只有根 1；$x''\in\{0,1\}$ 亦被排除：$x''=1\iff xx'=1$，代入 (4.6) 得 $x=x'$。）(4.6) 的点上 $x+x'=p^2+p$（$p=xx'$），故 $\chi(x)\chi(x')=\chi(p^2+p)=1$，且由对称性 $\sum\chi(x)=\sum\chi(x')$，得 $\mathrm{OFF}=5P-4\sum_{\text{点}}\chi(x)$，$P$ 为 (4.6) 在 $x,x'\notin\{0,1\}$ 上的点数。固定 $x\notin\{0,1\}$，(4.6) 是 $x'$ 的二次方程 $x^2x'^2+(x+1)x'+x=0$，有 $1+\chi\bigl(x^3/(x+1)^2\bigr)$ 个根（且根自动 $\notin\{0,1\}$）。因此
$$P=\sum_{x\ne0,1}\Bigl(1+\chi\tfrac{x^3}{(x+1)^2}\Bigr),\qquad \sum_{\text{点}}\chi(x)=\sum_{x\ne0,1}\chi(x)+\sum_{x\ne0,1}\chi\Bigl(x+\frac{x^3}{(x+1)^2}\Bigr).$$
令 $y=x+1$：$\frac{x^3}{(x+1)^2}=y+1+\frac1y+\frac1{y^2}$，$\chi$ 值为 $\chi(y)(-1)^k$（因 $\operatorname{Tr}(\frac1y+\frac1{y^2})=0$）；$x+\frac{x^3}{(x+1)^2}=\frac{x}{(x+1)^2}=\frac1y+\frac1{y^2}$，$\chi$ 值为 1。于是 $P=(q-2)+(-1)^k(-1-(-1)^k)=q-3-(-1)^k$，$\sum_{\text{点}}\chi(x)=-1-(-1)^k+(q-2)=q-3-(-1)^k$，$\mathrm{OFF}=P=q-3-(-1)^k$，$\Sigma=3q-4+(-1)^k$。$\square$

### 4.4 结论
**定理 3.** 对所有 $q=2^k$：$\sum_{\mathcal G}\varphi^2=q^6(q-1)(1+\kappa_k+3q-4+(-1)^k)=3q^7(q-1)$，即
$$\boxed{\frac{\operatorname{Var}_{\mathcal G}(\Psi)}{q^5}=3-\frac3q,\qquad \frac{\operatorname{Var}_{\mathcal G^2}}{\operatorname{Var}_{\mathcal G}}=\frac{q^2}3 .}$$
（$k$ 奇：$1+4+3q-5=3q$；$k$ 偶：$1+2+3q-3=3q$。）

**定理 4（$\rho$ 的完整刻画）.** (i) $\rho(\varepsilon_{\lambda_1,0})=1$（$q-1$ 个）；(ii) $\rho(\varepsilon_{\lambda_1,\lambda_3})=r(\lambda_3)$ 与 $\lambda_1$ 无关，且 $r(s^3\lambda_3)=r(\lambda_3)$；(iii) $k$ 奇时 $r\equiv\frac{q-3}{3(q-1)}$；(iv) $\rho_{\rm avg}=\frac{q^2/3-1}{q^2-1}\to\frac13$。证明：(i)(ii) 由推论 1 与伸缩对称；(iii) $k$ 奇时立方映射满，故 $\Phi(y)=\sum_x\varphi(0,x,y)^2$ 在 $y\ne0$ 上为常数 $=\frac{3q^7(q-1)-q^7(q-1)}{q-1}=2q^7$，于是 $r=\frac{q^7(q-1)-2q^7}{3q^7(q-1)}$；(iv) 由定理 2、3。$k$ 偶时 $r$ 在 $\lambda_3$ 的两个三次剩余类上取两个值（$q=4$：$5/9,-1/9$；$q=16$：$14/45,5/18$；$q=64$：$127/378,239/756$），它们由 (4.2)–(4.4) 的显式 Fourier 数据决定。$\square$

**数值核对。** (4.2)–(4.4) 对全部 $(\mu,\nu)$ 在 $k\le6$ 逐点核对（`Fourier predictions ... OK`）；引理 B 以 $O(q)$ 算法核到 $k=22$（$q=4194304$）；定理 3 的值 $3-3/q$ 由论文 FFT 算法在 $q=16,32,64,128$ 独立得到（$2.8125,\ 2.90625,\ 2.953125,\ 2.9765625$）。

---

## 5. 关于 $d=4$ 的一句话总结
在 $\mathbb F_{q^8}$ 上，$\alpha\mapsto(\alpha+\alpha^{q^4},\alpha^{1+q^4})$ 把"前三个系数"问题化为 $\mathbb F_{q^4}$ 上二次式的计数；所有非平凡性集中在 $a\in\mathbb F_{q^2}\setminus\mathbb F_q$ 这一薄层（情形 II），其贡献由 $\mathbb F_q$ 上的小三次和 $K$ 控制。Legendre 子群 $\{(0,c,0)\}$ 恰好是 $K$ 退化为 $q[k\text{ 奇}]$ 的位置，故 $\Psi$ 取两个值 $q^5-q^4+q^3$ 与 $q^5\pm q^3$，方差 $q(q-1)q^5$；全群方差则由三次族 $w^3+(1+j)w^2+w+1$ 的根的迹分布决定，经引理 B 给出 $3-3/q$。KR 极限 $2$ 失效的根源是 (4.2)–(4.4)：$\widehat\Psi$ 的 $q^2-1$ 个非平凡频率中有 $O(q^2)$ 个取值 $\asymp q^4$（而非 $O(q^{7/2})$），对应 $\sum_\chi\psi_8(\chi)=-q^7+q^6$。

---

## 6. 偶数 $d$ 标度律的判定（$d=6$，$q\le32$）

`legendre_q16_hayes.py` 调用论文的 Hayes 群 FFT（$q=16$：$|\mathcal G|=2^{20}$，并经 CRT 精确化；$q=32$：$|\mathcal G|=2^{25}$，浮点舍入偏差 $0$）。连同第一部分的数据：

| $q$ | $\operatorname{Var}_{\mathcal G}/q^7$ (KR 4) | $\operatorname{Var}_{\mathcal G^2}/q^7$ | 比率 | $\operatorname{Var}_{\mathcal G^2}/q^{7}\,/\,q^2$ | $\operatorname{Var}_{\mathcal G^2}/q^{7}\,/\,q^4$ |
|---|---|---|---|---|---|
| 2 | 5.758 | 15.56 | 2.70 | 3.89 | 0.97 |
| 4 | 3.700 | 9.01 | 2.44 | 0.563 | 0.0352 |
| 8 | 4.748 | 79.09 | 16.66 | 1.236 | 0.0193 |
| 16 | 4.351 | 225.51 | 51.82 | 0.881 | 0.00344 |
| 32 | 4.209 | 1100.88 | 261.58 | 1.075 | 0.00105 |

**判定**：$\operatorname{Var}_{\mathcal G^2}/q^{d+1}$ 在 $d=6$ 时 $\asymp q^2$（$q=8,16,32$ 时 $/q^2$ 为 $1.24,0.88,1.08$），$/q^4\to0$。提示中的 $q^{2(m-1)}=q^4$ **被否定**。比率 $\operatorname{Var}_{\mathcal G^2}/\operatorname{Var}_{\mathcal G}$ 在 $d=6$ 为 $q^2/3.84,\ q^2/4.94,\ q^2/3.91$（$q=8,16,32$），与 $d=4$ 的精确 $q^2/3$ 同阶。$\operatorname{Var}_{\mathcal G}/q^7$ 在 $q=16,32$ 为 $4.35,4.21$，趋向 KR 值 4（与 $d=4$ 的 $3\ne2$ 不同：$d=6$ 有 4 个零点，单值群看来是满的）。

**机制（`hayes_d6_profile.txt`，$q=2,\dots,32$ 精确）**：$d=6$ 时 $\mathcal G^2$（$q^2$ 个类 $(0,c_2,0,c_4,0)$）上 $\varphi=\Psi-q^7$ **只取三个值**：

| $q$ | $c_2=0$ 的 $q$ 个类 | 份额 | 其余 $(q^2-q)/2$ 个类 | 其余 $(q^2-q)/2$ 个类 |
|---|---|---|---|---|
| 2 | $+32$ | — | $+36$ | $-68$ |
| 4 | $-768$ | — | $-16$ | $+16$ |
| 8 | $+35840$ | 96.8% | $-3264$ | $+1216$ |
| 16 | $-983040$ | 99.78% | $-12032$ | $+12032$ |
| 32 | $+34537472$ | 98.55% | $-816128$ | $+685056$ |

$c_2=0$ 类（即 $f=t^6+f_5t^5+\cdots$ 中 $f_5=0$ 的 Legendre 区间）的偏差对 $k=1..5$ 恰为
$$\varphi=\begin{cases}+q^5+q^4-2q^3&k\text{ 奇}\\ -q^5+q^4&k\text{ 偶}\end{cases}$$
（$32=32+16-16$，$-768=-1024+256$，$35840=32768+4096-1024$，$-983040=-16^5+16^4$，$34537472=32^5+32^4-2\cdot32^3$），相对偏差 $\pm q^{-2}$，符号随 $k$ 的奇偶交替——与 $d=4$ 的 (T3) 同型（那里 $\pm q^3$ 的符号也是 $(-1)^k$）。于是
$$\frac{\operatorname{Var}_{\mathcal G^2}}{q^{7}}=\frac{q\cdot(q^5\pm\cdots)^2+O(q^2\cdot q^{8})}{q^2\cdot q^7}=q^2\,(1+O(1/q)),$$
$d=6$ 的 $q^2$ 律由此解释：$d=4$ 是**一个**类的相对亏损 $1/q$，$d=6$ 是 **$q$ 个**类的相对偏差 $1/q^2$，两者都给出 $\operatorname{Var}_{\mathcal G^2}/q^{d+1}\asymp q^2$（$q^8/(q\cdot q^5)$ 与 $q\cdot q^{10}/(q^2\cdot q^7)$）。单位类本身在 $d=6$ 只占 $3\%\sim12\%$。

**猜想（基于 $k\le5$ 的精确数据）**：对偶数 $d$，$\lim_{q\to\infty}\operatorname{Var}_{\mathcal G^2}(\Psi)/q^{d+3}=:C_d$ 存在，$C_4=1$（定理 2），$C_6=1$（上表：$1.24,0.88,1.08\to1$ 的修正项为 $\pm2/q+O(q^{-2})$）；反常质量集中在 $f_{d-1}=0$（$c_2=0$）的 $q^{d/2-1}$ 个 Legendre 类上。奇数 $d$ 时比率有界（$q=8$：$d=5$ 为 4.3，$d=7$ 为 2.0）。证明 $d=6$ 的闭式需要把 §2 的塔方法推广到 $\mathbb F_{q^{12}}\supset\mathbb F_{q^6}$（6 个共轭的二次式，条件涉及 $e_1,\dots,e_5$）；由于 $12=4\cdot3$，第二层需处理 $\mathbb F_{q^6}\supset\mathbb F_{q^3}$ 的三次结构，原理相同但工作量更大。

### 6.1 $d=6$ 的 Legendre 剖面与 $\rho$ 分布（$q=16$）
$\rho(\varepsilon)$ 按 $(p_1,p_3,p_5)$-支撑分组（$q=16$，共 4095 个非平凡二次特征，$1+\sum\rho=51.82$）：全部 $\rho>0$，单个最大 $0.081$（15 个仅依赖 $p_1$ 的特征），其余 $0.005\sim0.025$，支撑为 $\{1,3,5\}$ 的 3375 个特征贡献 $41.7$。与 $d=4$ 的"$q-1$ 个 $\rho=1$ + 常数"结构完全不同：$d=6$ 的相关性确实弥散，但方向一致（全为正），这正是它仍能叠加出 $q^2$ 量级比率的原因。剖面明细见 `hayes_d6_profile.txt`（$q=8,16,32$ 的全部偏差值及其在 $\sum\varphi^2$ 中的份额）。
