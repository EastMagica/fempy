
## basic.py

### function

affine_tri(space_s, space_e, point_s)

$$
\left[
\begin{matrix}
\lambda^0_1 & \lambda^1_1 & \cdots & \lambda^n_1 \\
\lambda^0_2 & \lambda^1_2 & \cdots & \lambda^n_2 \\
\lambda^0_3 & \lambda^1_3 & \cdots & \lambda^n_3 \\
\end{matrix}
\right]
\cdot
\left[
\begin{matrix}
x_0 & y_0 \\
x_1 & y_1 \\
x_2 & y_2
\end{matrix}
\right]
=
\left[
\begin{matrix}
\lambda^0_1 x_0 + \lambda^0_2 x_1 + \lambda^0_3 x_2 & \lambda^0_1 y_0 + \lambda^0_2 y_1 + \lambda^0_3 y_2 \\
\vdots & \vdots \\
\lambda^n_1 x_0 + \lambda^n_2 x_1 + \lambda^n_3 x_2 & \lambda^n_1 y_0 + \lambda^n_2 y_1 + \lambda^n_3 y_2 \\
\end{matrix}
\right]
$$