#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 19:52
# @file    : plot2d.py
# @project : fempy
# software : PyCharm

"""
转换 ``scipy.spatial.Delaunay`` 类型
为 ``matplotlib.tri.Triangulation``.

plt.contour  : 绘制等高线
plt.clabel   : 标注等高线数据
plt.contourf : 填充色彩
plt.xticks(()); plt.xticks(()) : 去除坐标轴

Notes
-----
等高线图中 ``colors`` 和 ``linewidths`` 代表等高线,
且等高线从 ``level[0]`` 开始, 依次交替,到 ``level[-1]``.

.. code::

    cm = plt.cm.get_cmap('viridis')
    fig, axs = plt.subplots(figsize=(8, 4.5), nrows=2, ncols=2)
    axs = axs.flatten()

"""

from scipy.spatial import Delaunay
from matplotlib.tri import Triangulation


def get_stri(stri):
    x, y = stri.points.T
    triangles = stri.simplices
    return x, y, triangles


def set_mtri(x, y, stri):
    mtri = Triangulation(x, y, triangles=stri)
    return mtri


def is_mtri(tri):
    if isinstance(tri, Delaunay):
        tri = set_mtri(*get_stri(tri))
    return tri


def cube_refine(mtriang, z):
    """
    cube = mtri.CubicTriInterpolator(triang, z)
    cube = mtri.LinearTriInterpolator(triang, z, kind={'min_E', 'geom', 'user'})

    x_mesh, y_mesh = np.meshgrid(np.linspace(0, 1, 16), np.linspace(0, 1, 16))

    z_ref = cube(x_mesh, y_mesh)
    dz_ref = cube.gradient(x_mesh, y_mesh)

    # norm=mpl.color.LogNorm(vmin=z.min(), vmax=z.max())
    # vmin=z.min(), vmax=z.max(), linewidths=4
    c = ax.pcolor(x_mesh, y_mesh, z_ref, edgecolor='k', cmap='RdBu')
    fig.colorbar(c, ax=ax, shrink=1.)

    # automatic layout subplots, to show colorbar
    fig, ax = plt.subplots(constrained_layout=True)
    cs1 = ax.contourf(x_mesh, y_mesh, z_ref, levels={int, array}, cmap=plt.cm.bone)
    cs2 = ax.contour(cs1, levels=cs1.levels[::2], color='red', origin=origin)
    fig.colorbar(c, ax=ax)

    Parameters
    ----------
    mtriang
    z

    Returns
    -------

    """


def tri_mesh(ax, triang, **kwargs):
    mtriang = is_mtri(triang)
    kwargs.setdefault('marker', 'o')
    kwargs['color'] = kwargs.pop('c', 'tab:blue')
    kwargs['linestyle'] = kwargs.pop('ls', '-')
    kwargs['markerfacecolor'] = kwargs.pop('mfc', kwargs['color'])
    ax.triplot(mtriang, **kwargs)


def tri_contourf(ax, triang, z, **kwargs):
    mtriang = is_mtri(triang)
    tcf = ax.tricontourf(mtriang, z, **kwargs)
    return tcf


def tri_contour(ax, triang, z, **kwargs):
    levels = kwargs.pop('levels', 5)
    contour = kwargs.pop('contour', dict())
    contour['colors'] = kwargs.get('colors', ['0.25', '0.25'])
    contour['linewidths'] = kwargs.get('linewidths', [0.5, 0.5])
    mtriang = is_mtri(triang)
    countour = ax.tricontour(mtriang, z, levels=levels, **contour)
    return countour


def tri_tripcolor(ax, triang, z, **kwargs):
    """
    shading='flat', edgecolors='k', cmap=cm

    Parameters
    ----------
    ax
    triang
    z
    shading : optional
    kwargs

    Returns
    -------

    Notes
    -----
    shading in ['flat', 'gouraud']
    edgecolor
    """
    mtriang = is_mtri(triang)
    mtriang = is_mtri(triang)
    tpc = ax.tripcolor(mtriang, z, **kwargs)
    return tpc
