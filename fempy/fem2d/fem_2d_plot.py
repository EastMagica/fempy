#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 19:52
# @file    : fem_2d_visual.py
# @project : fempy
# software : PyCharm

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D


def plot_tri(ax, u_lst, mesh, color=plt.cm.YlGnBu_r):
    tri = mtri.Triangulation(mesh.p_mat[:, 0], mesh.p_mat[:, 1], triangles=mesh.u_mat)
    ax.plot_trisurf(mesh.p_mat[:, 0], mesh.p_mat[:, 1], u_lst, triangles=tri.triangles, cmap=color)
