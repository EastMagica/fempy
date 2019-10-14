#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 18:00
# @file    : fempy/fem2d/fem2d.py
# @project : fempy
# software : PyCharm

import types
import numpy as np
from fempy.basic import FEM
from .basic import Gaussian, area


class FEM2D(FEM):
    def __init__(self, variation, mesh):
        super().__init__(variation, mesh)
        self.gaussian = Gaussian(6, self.mesh.ndim)
        print('> points:', self.mesh.npoints, ', units:', self.mesh.nsimplices)

    def boundary(self):
        r"""边界条件

        Returns
        -------

        Note
        ----
        此处边界条件仅为第一类(Dirichlet)边界条件.
        """
        bnd = self.mesh.boundary
        self.f_lst[bnd] = 0
        self.a_mat[:, bnd] = 0
        self.a_mat[bnd, :] = 0
        self.a_mat[bnd, bnd] = 1

    @staticmethod
    def basis_value(p, v):
        r"""单元基函数.

        由单元基函数的定义, 可以得到单元基函数在点:math:`p`处取值.

        .. math::
            \varphi_1 = \frac{\Delta_{p23}}{\Delta} \\
            \varphi_2 = \frac{\Delta_{p31}}{\Delta} \\
            \varphi_3 = \frac{\Delta_{p12}}{\Delta} \\

        Parameters
        ----------
        p : array_like, (N, 2)
        v : array_like, (3, 2)

        Returns
        -------
        basis_v : ndarray, (N, 3)
        """
        area_v = area(v[0], v[1], v[2])
        basis_v = np.array([area(p, v[1], v[2]), area(p, v[2], v[0]), area(p, v[0], v[1])]) / area_v
        return basis_v.T

    @staticmethod
    def basis_grid(p, v):
        r"""单元基函数的导数.

        由单元基函数的定义, 可以得到单元基函数的导数,
        且三角单元基函数的导数为常数.

        .. math::
            \nabla \varphi_1 =
                \left[\begin{matrix}
                (y_2 - y_3) / \Delta \\ (x_3 - x_2) / 2\Delta
                \end{matrix}\right],\quad
            \nabla \varphi_1 =
                \left[\begin{matrix}
                (y_3 - y_1) / \Delta \\ (x_1 - x_3) / 2\Delta
                \end{matrix}\right],\quad
            \nabla \varphi_1 =
                \left[\begin{matrix}
                (y_1 - y_2) / \Delta \\ (x_2 - x_1) / 2\Delta
                \end{matrix}\right]

        Parameters
        ----------
        p : array_like, (3, 2)
        v : array_like, (N, 2)

        Returns
        -------
        basis_g : ndarray, (2, 3)

        """
        area_v = area(v[0], v[1], v[2])
        # TODO: Fix Bug, commmit to BLog
        basis_g = np.array([v[[1, 2, 0], 1] - v[[2, 0, 1], 1],
                            v[[2, 0, 1], 0] - v[[1, 2, 0], 0]])
        basis_g = basis_g / 2 / area_v
        return basis_g

