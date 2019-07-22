#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/12 17:45
# @file    : main.py
# @project : fempy
# software : PyCharm

import numpy as np
import matplotlib.pyplot as plt
from fempy.fem1d import FEM1D, Mesh1D, Adaptive1D


eq = {'bnd': np.array([0, 0]),
      'f': lambda x: np.pi**2 * np.sin(np.pi * x)}

mesh = Mesh1D({'domain': np.array([0, 1]), 'opt': ['step', 8]})

fem = FEM1D(eq, mesh)
fem.run()

u_lst = fem.u_lst
u_true = np.sin(np.pi * fem.p_mat)
u_error = u_lst - u_true
u_error_l2 = fem.error_l2(lambda x: np.sin(np.pi * x))


print('u_lst:\n', u_lst)
print('u_true:\n', u_true)
print('u_error:\n', u_error)
print('u_error_l2=', u_error_l2)


# plt.style.use('grayscale')
fig = plt.figure(figsize=(12, 5))
ax0 = plt.subplot(1, 2, 1)
plt.title('u(x)')
plt.plot(np.linspace(0, 1, 100), np.sin(np.pi * np.linspace(0, 1, 100)), color='b', linestyle='-')
plt.plot(fem.p_mat, fem.u_lst, color='r', linestyle='--')
ax1 = plt.subplot(1, 2, 2)
plt.title('error')
plt.plot(fem.p_mat, np.sin(np.pi * fem.p_mat) - fem.u_lst)
plt.tight_layout()
plt.show()
