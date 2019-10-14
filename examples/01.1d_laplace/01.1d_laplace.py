#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/12 17:45
# @file    : 01.1d_laplace.py
# @project : fempy
# software : PyCharm

import numpy as np
import matplotlib.pyplot as plt
from fempy.mesh import Mesh1D
from fempy.fem1d import FEM1D

r"""
A simple example of 1D Laplace equation. 

.. math::
    \begin{cases}
    \nabla^2 u = \pi^2 \sin(\pi x) \\
    u(0) = u(1) = 0
    \end{cases}

and the exact solution is :math:`u = \sin(\pi x)`.
"""

__test_name__ = '01.1d_laplace'


def f(x):
    return np.pi ** 2 * np.sin(np.pi * x)


def u_true(x):
    return np.sin(np.pi * x)


def variation(basis_v, basis_g, gauss_p, gauss_w):
    f_v = f(*gauss_p.T)
    f_elem = np.dot(f_v * gauss_w, basis_v)
    a_elem = np.dot(basis_g.T, basis_g) * np.sum(gauss_w)
    return a_elem, f_elem


mesh = Mesh1D('linspace', (0, 1, 17))

fem = FEM1D(variation, mesh)
fem.run()

u_error_l2 = fem.error_l2(u_true)

# print('> Matrix A=\n', fem.a_mat)
# print('> Matrix F=\n', fem.f_lst)
print('> Solution U=\n', fem.u_lst)
print('> L2Error =', fem.error_l2(u_true))

fem.save_data(file_name='01.1d_laplace', form='npy')

# plt.style.use('grayscale')
fig = plt.figure(figsize=(12, 5))
ax0 = plt.subplot(1, 2, 1)
plt.title('u(x)')
plt.xlim(-0.05, 1.05)
plt.ylim(-0.05, 1.35)
plt.plot(np.linspace(0, 1, 100), np.sin(np.pi * np.linspace(0, 1, 100)), color='b', linestyle='-')
plt.plot(fem.mesh.points, fem.u_lst, color='r', linestyle='--')
ax1 = plt.subplot(1, 2, 2)
plt.title('error')
plt.plot(fem.mesh.points, np.sin(np.pi * fem.mesh.points) - fem.u_lst)
plt.tight_layout()

plt.savefig(__test_name__ + '.png', dpi=200)
