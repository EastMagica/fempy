#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 18:03
# @file    : test_2d_01.py
# @project : fem
# software : PyCharm

import numpy as np
import matplotlib.pyplot as plt
from fem.fem2d import FEM2D, Mesh2D, plot_tri


eq = {'f': lambda x, y: 2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y),
      'u': lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)}
mesh = Mesh2D(0, 1, 0, 1, 16, 16)

fem = FEM2D(eq, mesh)
fem.run()

print('> Matrix A=\n', fem.a_mat)
print('> Matrix F=\n', fem.f_lst)
print('> Solution U=\n', fem.u_lst)
print('> Exact Solution U=\n', [eq['u'](i[0], i[1]) for i in fem.mesh.p_mat])
print('> L2Error =', fem.error_l2())


fig1 = plt.figure()
ax1 = fig1.gca(projection='3d')
ax1.set_title('Numerical Solution')
plot_tri(ax1, fem.u_lst, mesh)

fig2 = plt.figure()
ax2 = fig2.gca(projection='3d')
ax2.set_title('Exact Solution')
plot_tri(ax2, [eq['u'](i[0], i[1]) for i in fem.mesh.p_mat], mesh)

fig3 = plt.figure()
ax3 = fig3.gca(projection='3d')
ax3.set_title('Error')
plot_tri(ax3, np.array([eq['u'](i[0], i[1]) for i in fem.mesh.p_mat]) - fem.u_lst, mesh)
plt.show()
