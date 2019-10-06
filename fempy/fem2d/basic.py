#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/18 10:57
# @file    : fempy/fem2d/basic.py
# @project : fempy
# software : PyCharm

import os
import numpy as np
from fempy.basic import _is_ndarray

pnt = np.array([[0, 0], [1, 0], [0, 1]])


def extend_points(points, opts='global', axis=1):
    r"""扩充矩阵

    扩充全局(global)或自然(natural)坐标矩阵.

    .. math:
        \left[\begin{matrix}
            x_1, x_2, \cdots, x_n \\
            y_1, y_2, \cdots, y_n \\
        \end{matrix}\right] \rightarrow
        \left[\begin{matrix}
            x_1, x_2, \cdots, x_n \\
            y_1, y_2, \cdots, y_n \\
            1  , 1  , \cdots, 1   \\
        \end{matrix}\right] \\
        \left[\begin{matrix}
            \lambda_{11}, \cdots, \lambda_{1n} \\
            \lambda_{21}, \cdots, \lambda_{2n} \\
        \end{matrix}\right] \rightarrow
        \left[\begin{matrix}
            \lambda_{11}, \cdots, \lambda_{n1} \\
            \lambda_{12}, \cdots, \lambda_{n2} \\
            \lambda_{13}, \cdots, \lambda_{n3} \\
        \end{matrix}\right]

    Parameters
    ----------
    points : array_like, (N, 2)
    opts   : {'global', 'natural'}, optional
        扩充矩阵类型(全局或自然坐标), 默认全局坐标.
    axis   : {0, 1}, optional
        坐标矩阵方向, 默认为列坐标(axis=1).

    Returns
    -------

    """
    points = _is_ndarray(points)
    # 坐标矩阵转置为列向量
    if axis == 0:
        points = points.transpose()
    elif axis != 1:
        raise ValueError('Axis values can only be 0 or 1.')
    # 构造不同的扩展分量
    if opts.lower() == 'golbal':
        extend_vector = np.zeros(points.shape[1])
    elif opts.lower() == 'natural':
        extend_vector = 1 - np.sum(points, axis=0)
    else:
        raise ValueError(opts + ' coordinates are not supported.')
    # 合并分量
    point_extend = np.vstack((points, extend_vector))
    return point_extend


def area(*points):
    r"""计算三角形面积.

    Parameters
    ----------
    points : array_like (N, 3, 2)
        三角单元顶点坐标.

    Returns
    -------
    s : array_like, (N, )
        三角单元面积

    Note
    ----
    若``points``为(3,2)数组, 则可使用``area(*points)``
    直接计算, 无需为对应参数而拆开数组.

    """
    [p0, p1, p2] = _is_ndarray(points)
    if (p0.shape[1] != 2) or (p1.shape[1] != 2) or (p2.shape[1] != 2):
        raise ValueError('Input a needs to be a (N, 2) Matrix.')
    elif p0.shape[0] != p1.shape[0] != p2.shape[0]:
        raise ValueError('Input matrices need to have same raw.')
    s = ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0])) / 2
    return s


def local_to_global(local_p, global_v):
    r"""
    局部(Local)坐标变换到全局(Global)坐标.

    .. math::
        \vec{x} = X \vec{\lambda} \\
        \left[\begin{matrix}
            x \\ y \\ 1 \\
        \end{matrix}\right] =
        \left[\begin{matrix}
            x_1 & x_2 & x_3 \\
            y_1 & y_2 & y_3 \\
            1   & 1   & 1   \\
        \end{matrix}\right]
        \left[\begin{matrix}
            \lambda_1 \\ \lambda_2 \\ \lambda_3 \\
        \end{matrix}\right]

    Parameters
    ----------
    local_p  : array_like, (N, 2)
        局部坐标中的点坐标.
    global_v : array_like, (3, 2)
        全局坐标中的三角形顶点坐标.

    Returns
    -------
    point_e
        Local Point in Global Area
    """
    local_pex = extend_points(local_p, 'natural', axis=0)
    global_p = np.dot(global_v.transpose(), local_pex)
    return global_p


def global_to_local(global_p, global_v):
    r"""局部(Local)坐标变换到全局(Global)坐标.

    .. math::
        \vec{\lambda} = X^{-1} \vec{x} \\
        \left[\begin{matrix}
            \lambda_1 \\ \lambda_2 \\ \lambda_3 \\
        \end{matrix}\right] =
        \left[\begin{matrix}
            x_1 & x_2 & x_3 \\
            y_1 & y_2 & y_3 \\
            1   & 1   & 1   \\
        \end{matrix}\right]^{-1}
        \left[\begin{matrix}
            x \\ y \\ 1 \\
        \end{matrix}\right]

    Parameters
    ----------
    global_p : array_like, (N, 2)
        全局坐标中的点坐标.
    global_v : array_like, (3, 2)
        全局坐标中的三角形顶点坐标.

    Returns
    -------

    Global Point in Local Area
    """
    global_vex = extend_points(global_v, 'global', axis=0)
    global_pex = extend_points(global_p, 'global', axis=0)
    global_inv = np.linalg.inv(global_vex)
    local_p = np.dot(global_inv, global_pex)
    return local_p


def basis_value(p, v):
    r"""单元基函数.

    由单元基函数的定义, 可以得到单元基函数在点:math:`p`处取值.

    .. math::
        \varphi_1 = \frac{\Delta_{p23}}{\Delta} \\
        \varphi_2 = \frac{\Delta_{p31}}{\Delta} \\
        \varphi_3 = \frac{\Delta_{p12}}{\Delta} \\

    Parameters
    ----------
    p
    v

    Returns
    -------

    """
    area_v = area(v[0], v[1], v[2])
    val = np.array([area(p, v[1], v[2]),
                    area(p, v[2], v[0]),
                    area(p, v[0], v[1])]) / area_v
    return val


def basis_grid(p, v):
    r"""单元基函数的导数.

    由单元基函数的定义, 可以得到单元基函数的导数,
    且三角单元基函数的导数为常数.

    .. math::
        \nabla \varphi_1 =
            \left[\begin{matrix}
            (y_2 - y_3) / \Delta \\ (x_3 - x_2) / \Delta
            \end{matrix}\right],\quad
        \nabla \varphi_1 =
            \left[\begin{matrix}
            (y_3 - y_1) / \Delta \\ (x_1 - x_3) / \Delta
            \end{matrix}\right],\quad
        \nabla \varphi_1 =
            \left[\begin{matrix}
            (y_1 - y_2) / \Delta \\ (x_2 - x_1) / \Delta
            \end{matrix}\right]

    Parameters
    ----------
    p
    v

    Returns
    -------

    """
    area_v = area(v[0], v[1], v[2])
    val = np.array([[v[1, 1] - v[2, 1], v[2, 0] - v[1, 0]],
                    [v[2, 1] - v[0, 1], v[0, 0] - v[2, 0]],
                    [v[0, 1] - v[1, 1], v[1, 0] - v[0, 0]]]) / area_v
    return val


class GaussQuad(object):
    # gauss_temp = TempElement(3)

    @staticmethod
    def gauss(n):
        r"""更改静态变量Gauss模板"""
        GaussQuad.gauss_temp = TempElement(n)

    @staticmethod
    def triquad(tri, *f):
        r"""三角单元上的Gauss求积.

        .. math::
            \iint_{\Delta}f(x,y)dS=\sum_{i=1}^nf(x_i)w_i

        Parameters
        ----------
        f   : function
            被积函数(integrand).
        tri : array_like, (3, 2)
            三角单元顶点坐标.

        Returns
        -------
        gauss_quad : float
            Gauss求积结果.
        """
        gauss = GaussQuad.gauss_temp
        gauss_local_p = gauss.gauss_point()
        gauss_local_w = gauss.gauss_weight()
        gauss_global_p = local_to_global(gauss_local_p, tri)
        gauss_global_w = area(*tri) / gauss.get_area() * gauss_local_w
        gauss_quad = np.sum(f(gauss_global_p) * gauss_global_w)
        return gauss_quad
