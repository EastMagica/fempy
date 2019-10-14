#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 18:03
# @file    : 02.2d_laplace.py
# @project : fempy
# software : PyCharm

import numpy as np
import matplotlib.pyplot as plt
from fempy.mesh import Mesh2D
from fempy.fem2d import FEM2D
from fempy.femplot.plot2d import tri_mesh, tri_tripcolor, tri_contour, tri_contourf

r"""
\nabla^2 u(x, y) + \nabla u(x, y) + u(x, y) = f(x, y)


\nabla^2 u
----------

(\nabla^2 u, \phi_j)
= -(\nabla^2 u, \nabla \phi_j)
= -(\sum^n_{i=1} u_i \nabla \phi_i, \nabla \phi_j)
= -\sum^n_{i=1} \left(u_i (\nabla \phi_i, \nabla \phi_j) \right)


\partial u
----------

(\nabla u, \phi_j)
= (\sum^n_{i=1} u_i \partial \phi_i, \phi_j)
= \sum^n_{i=1} \left(u_i (\partial \phi_i, \phi_j) \right)


u
----------

(u, \phi_j)
= (\sum^n_{i=1} u_i \phi_i, \phi_j)
= \sum^n_{i=1} \left(u_i (\phi_i, \phi_j) \right)


matrix compute
--------------

2 dim, triangular mesh, linear basis function

.. Tips: basis grad is constant in linear basis function.

gauss_n = int: N
basis_v = ndarray, (N, 3)
-------------------------
[[basis_0(p[0]), basis_1(p[0]), basis_2(p[0])],
 .............., ............., ..............,
 [basis_0(p[n]), basis_1(p[n]), basis_2(p[n])],]

basis_g = ndarray, (2, 3)
-------------------------
[[basis_0_x(p[0]), basis_1_x(p[0]), basis_2_x(p[0])],
 [basis_0_y(p[0]), basis_1_y(p[0]), basis_2_y(p[0])]]

basis_g = ndarray, (2, N, 3)
----------------------------
[[[basis_0_x(p[0]), basis_1_x(p[0]), basis_2_x(p[0])],
  ................, ..............., ................,
  [basis_0_x(p[n]), basis_1_x(p[n]), basis_2_x(p[n])],],
 [[basis_0_y(p[0]), basis_1_y(p[0]), basis_2_y(p[0])],
  ................, ..............., ................,
  [basis_0_y(p[n]), basis_1_y(p[n]), basis_2_y(p[n])]]]


\nabla^2 u
----------
-\sum^n_{i=1} \left(u_i (\nabla \phi_i, \nabla \phi_j) \right)

if grad is constant
^^^^^^^^^^^^^^^^^^^

> np.dot(basis_g.T, basis_g) * np.sum(gauss_w)

[[u_0_x, u_0_y],
 [u_1_x, u_1_y],
 [u_2_x, u_2_y],]

[[u_0_x, u_1_x, u_2_x, ],
 [u_0_y, u_1_y, u_2_y, ]]


if grad is not constant
^^^^^^^^^^^^^^^^^^^^^^^

> basis_g_x, basis_g_y = basis_g
> np.dot(basis_g_x.T * gauss_w, basis_g_x) + np.dot(basis_g_y.T * gauss_w, basis_g_y)

.. (3, N), (3, N) = (2, 3, N)
.. (3, N) .* (N, ) * (N, 3) + (3, N) .* (N, ) * (N, 3) -> (3, 3)


\partial u
--------
\sum^n_{i=1} \left(u_i (\partial \phi_i, \phi_j) \right)

if grad is constant
^^^^^^^^^^^^^^^^^^^

> basis_g_x, basis_g_y = basis_g.transpose()
> np.dot(basis_v.T, gauss_w) * basis_g_x

.. (N, 3).T * (N,  ) -> (3, 1) * (1, 3) -> (3, 3)

if grad is not constant
^^^^^^^^^^^^^^^^^^^^^^^

> basis_g_x, basis_g_y = basis_g.transpose((0, 1, 2))
> np.dot(basis_v.T * gauss_w, basis_g_x)

u
---

\sum^n_{i=1} \left(u_i (\phi_i, \phi_j) \right)

> np.dot(basis_v.T, basis_v) * gauss_w

.. (3, N) * (N, 3) -> (3, 3)


"""


r"""
A simple example of **2D Laplace equation**. 

.. math::
    \begin{cases}
    \nabla^2 u = 2\pi^2 \sin(\pi x) \sin(\pi y) \\
    u(x, 0) = u(0, y) = u(x, 1) = u(1, y) = 0
    \end{cases}
    
and the variational formulation is

.. math::
    \begin{align}
    (-\nabla^2 u, v) &= (f, v) \\
    (\nabla u, \nabla v) &= (f, v) \\
    (\sum^n_{i=0} u_i \nabla \phi_i, \nabla \phi_j) &= (f, \phi_j),\quad j=1,\cdots,n \\
    \sum^n_{i=0} u_i (\nabla \phi_i, \nabla phi_j) &= (f, \phi_j),\quad j=1,\cdots,n \\
    \end{align}
    
and the *exact solution* is :math:`u = \sin(\pi x) \sin(\pi y)`.

The default is the **linear basis functions**.
"""

__test_name__ = '02.2d_laplace'


# Init condition
# --------------

def f(x, y):
    return 2 * np.pi ** 2 * np.sin(np.pi * x) * np.sin(np.pi * y)


def u_true(x, y):
    return np.sin(np.pi * x) * np.sin(np.pi * y)


def variation(basis_v, basis_g, gauss_p, gauss_w):
    f_v = f(*gauss_p.T)
    f_elem = np.dot(f_v * gauss_w, basis_v)
    a_elem = np.dot(basis_g.T, basis_g) * np.sum(gauss_w)
    return a_elem, f_elem


# compute numerical solution
# --------------------------

mesh = Mesh2D('rect_tri', [(0, 1, 16)])

fem = FEM2D(variation, mesh)
fem.run()

# print and save solution
# -----------------------

# fem.save_data(file_name='02.2d_laplace\data', form='npy')

# print('> Matrix A=\n', fem.a_mat)
# print('> Matrix F=\n', fem.f_lst)
print('> Solution U=\n', fem.u_lst)
print('> L2Error =', fem.error_l2(u_true))

# process data
# ------------

z = fem.u_lst
mtri = fem.mesh.mtri
xp, yp = fem.mesh.points.T

# interpolation
r_mesh, r_u = fem.mesh.refine(z, subdiv=3)

# visualization
# -------------

cm = plt.cm.get_cmap('RdBu_r')

fig, axs = plt.subplots(figsize=(12, 10), nrows=2, ncols=2)
axs = axs.flatten()

axs[0].set_aspect('equal')
tri_mesh(axs[0], mtri, mfc='tab:green')

axs[1].set_aspect('equal')
tcf = tri_contourf(axs[1], mtri, z, levels=10, cmap=cm)
fig.colorbar(tcf, ax=axs[1])

axs[2].set_aspect('equal')
tri_tripcolor(axs[2], mtri, z, edgecolor='k', cmap=cm)

axs[3].set_aspect('equal')
tcf = tri_contour(axs[3], mtri, z, levels=10, cmap=cm)
fig.colorbar(tcf, ax=axs[3])

# plt.savefig(__test_name__ + 'png', dpi=200)
plt.show()
