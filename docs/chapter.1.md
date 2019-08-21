# 有限元方法计算流程（一维情形）

考虑一维空间中的函数
$$
\begin{cases}
-u''(x)=\pi\sin(\pi x),\quad x\in(0,1)\\
u(0)=u(1)=0
\end{cases}
$$

## 1. 构造变分形式

### 弱形式

$$
-\int u''(x)v(x)\mathrm{d}x=\int f(x)v(x)\mathrm{d}x
$$

分部积分方法，可化为
$$
\int u'(x)v'(x)\mathrm{d}x=\int f(x)v(x)\mathrm{d}x
$$


### 离散形式


$$
\int u'(x)v_h'(x)\mathrm{d}x=\int f(x)v_h(x)\mathrm{d}x,\quad v_h\in H_0^1(0,1)
$$

## 2. 构造函数空间

$V_h \subset H_0^1,\quad \dim(V_h)=N$

$u_h \in V_h$ 

$(u_h', v_h')=(f,v_h'), \quad v_h\in V_h$

$V_h=\mathrm{span}\{\phi_i\}_{i=1}^N$

$u_h=\sum_{i=1}^N u_i \varphi_i$

$\sum^N_{i=1}u_i(\phi_i', \phi_j')=(f, \phi_j)$

构造基函数
$$
\begin{cases}
\varphi_L=\frac{x-x_R}{x_L-x_R} \\
\varphi_R=\frac{x-x_L}{x_R-x_L}
\end{cases}
$$

## 3. 划分网格（Split Mesh）

分割方法略，一般使用一致网格

## 4. 组装矩阵（Assembly Matrix）

$$
AU=F \\
\left[
\begin{matrix}
a_{1,1} & \cdots & a_{1,N} \\
\vdots & \ddots & \vdots \\
a_{N,1} & \cdots & a_{N,N}
\end{matrix}
\right]_{N\times N}
\left[
\begin{matrix}
u_1 \\
\vdots \\
u_N
\end{matrix}
\right]_{N\times 1}
=\left[
\begin{matrix}
F_1 \\
\vdots \\
F_N
\end{matrix}
\right]_{N\times 1}
$$



总刚度矩阵
$$
\left[
\begin{matrix}
({\phi^i_L}', {\phi^i_L}') & ({\phi^i_L}', {\phi^i_R}') \\
({\phi^i_R}', {\phi^i_L}') & ({\phi^i_R}', {\phi^i_R}')
\end{matrix}
\right]
$$


单元刚度矩阵

