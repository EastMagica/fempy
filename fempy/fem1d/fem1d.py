#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/12 17:04
# @file    : fempy/fem1d/fem1d.py
# @project : fempy
# software : PyCharm

import numpy as np
from scipy.integrate import quad
from fempy.basic import FEM, _is_ndarray
from .basic import Gaussian, area


class FEM1D(FEM):
    def __init__(self, variation, mesh):
        super().__init__(variation, mesh)
        self.gaussian = Gaussian(3, '1')
        print('> points:', self.mesh.npoints, ', units:', self.mesh.nsimplices)

    def boundary(self):
        r"""边界条件

        Returns
        -------

        Note
        ----
        此处边界条件仅为第一类(Dirichlet)边界条件.
        """
        p_bnd = self.mesh.boundary
        self.f_lst[p_bnd] = 0
        self.a_mat[p_bnd, :] = 0
        self.a_mat[:, p_bnd] = 0
        self.a_mat[p_bnd, p_bnd] = 1

    @staticmethod
    def basis_value(p, v):
        """

        Parameters
        ----------
        p
        v

        Returns
        -------
        basis_v : ndarray, (N, 2)
        """
        p = _is_ndarray(p).squeeze()
        value = np.array([v[1] - p, p - v[0]]) / area(*v)
        return value.T

    @staticmethod
    def basis_grid(p, v):
        """

        Parameters
        ----------
        p
        v

        Returns
        -------
        basis_g : ndarray, (1, 2)

        """
        grid = np.array([[1], [-1]]) / area(*v)
        return grid.T


# TODO: refactor adaptive method.
class Adaptive1D(object):
    def __init__(self, variation, mesh, niter):
        self.niter = niter
        self.mesh = mesh
        self.variation = variation
        self.fem = FEM1D(self.variation, self.mesh)
        self.fem.run()

    def estimate_error(self):
        r"""
        误差估计.

        Returns
        -------

        """
        eta_h = np.zeros(self.mesh.nsimplices)
        for k, v in enumerate(self.mesh.simplices):
            l_cor, r_cor = self.mesh.points[v]
            eta_h[k] = (r_cor - l_cor) * np.sqrt(quad(lambda x: self.f(x)**2, l_cor, r_cor)[0])
        return eta_h

    def fem_iter(self, n=1):
        for i in range(self.niter):
            eta_h = self.estimate_error()
            self.adaptive_mesh(eta_h)
            self.fem = FEM1D(self.variation, self.mesh)
            self.fem.run()

    def adaptive_mesh(self, eta_h):
        inx = 0
        eta_max = np.max(eta_h)
        eta_num = np.where(eta_h > eta_max * 0.75)[0]
        lst = np.zeros((self.mesh.e_n + len(eta_num), 2), dtype=np.int)
        print('> eta_num:', eta_num)
        self.mesh.points = np.append(self.mesh.points, np.sum(self.mesh.points[self.mesh.simplices[eta_num]], axis=1) / 2)
        for k in range(self.mesh.e_n):
            lst[k+inx, 0] = self.mesh.simplices[k, 0]
            if k in eta_num:
                lst[k+inx, 1] = self.mesh.npoints + inx
                lst[k+inx+1, 0] = self.mesh.npoints + inx
                inx += 1
            lst[k+inx, 1] = self.mesh.simplices[k, 1]
        self.mesh.simplices = np.array(lst)
        self.mesh.npoints = len(self.mesh.points)
        self.mesh.e_n = len(self.mesh.simplices)
        self.mesh.p_end = np.append(self.mesh.p_end, [False]*len(eta_num))
