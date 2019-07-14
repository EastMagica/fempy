#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 19:52
# @file    : fem_2d_visual.py
# @project : fem
# software : PyCharm

from pyecharts.charts import Bar

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection='3d')

x = [0, 1, 2, 0, 1, 2, 0, 1, 2]
y = [0, 0, 0, 1, 1, 1, 2, 2, 2]
z = [0, 1, 2, 1, 3, 1, 2, 1, 0]

e = [[3, 0, 4], [1, 4, 0], [4, 1, 5], [2, 5, 1],
     [5, 8, 4], [7, 4, 8], [4, 7, 3], [6, 3, 7]]

tri = mtri.Triangulation(x, y, triangles=e)

ax.plot_trisurf(x, y, z, triangles=tri.triangles, cmap=plt.cm.Spectral)

plt.show()
