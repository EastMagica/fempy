# 有限元方法

## 从差分方法到变分形式

一般的偏微分方程
$$
\begin{cases}
-u''=\pi\sin(\pi x),\quad x\in(0, 1) \\
u(0)=u(1)=0
\end{cases}
$$

### 差分方法

差分方法是通过**离散**的方法逼近，计算有限个点处的数值

### 变分方法

与差分方法的思想不同，通过**高阶函数**逼近原函数，**维尔斯特拉斯定理**（连续函数必定存在高阶函数逼近）。

平方误差逼近，使得积分误差（$\int_a^b (u-u_h)^2\mathrm{d}x$）最小

定义 n-1 阶逼近函数 $\tilde{u}=\sum^{n-1}_{i=0}a_i x^i$

带入方程中，得到 $-\tilde{u}''=\pi\sin(\pi x)$，而为求解 $\{a_i\}^n_{i=0}$，则需构造 n 个方程，则有如下两种方法：

1. 求导格式

$$
\begin{cases}
-\tilde{u}''=f \\
-\tilde{u}^{(3)}=f' \\
\cdots \\
-\tilde{u}^{(n+1)}=f^{(n-1)}
\end{cases}
$$

2. 积分格式

$$
\begin{cases}
-\tilde{u}''=f \\
-\int\tilde{u}''\mathrm{d}x=\int f\mathrm{d}x \\
\cdots \\
-\underbrace{\idotsint}_{\text{n-1}}\tilde{u}''\mathrm{d}x^n=\underbrace{\idotsint}_{\text{n-1}} f\mathrm{d}x
\end{cases}
$$

3. 谱方法

   定义 n-1 个光滑函数 $\{{\varphi_i(x)}\}^n_{i=1}$，则可定义

$$
\begin{cases}
\cdots \\
-\int_a^b\tilde{u}''(x)\varphi_i(x)\mathrm{d}x=\int_a^b f(x)\varphi_i(x)\mathrm{d}x \\
\end{cases}
$$

4. 谱配置法

   选择 n 个点，拉格朗日方法，构造$\tilde{u}_n=\sum uL(x)$

### 有限元方法

由于高阶函数逼近会出现**龙格现象**，故而引导出了**有限元方法**，构造分段线性函数逼近
$$
V(x)=\{v\in C(0,1)|\ V|_{I_i}=P_1(x),i=1.\cdots,n\} \\
\tilde{u}(x)=\sum u_i \phi_i(x), \quad \phi_i(x)\in V(x)
$$










