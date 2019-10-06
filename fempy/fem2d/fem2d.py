#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 18:00
# @file    : fempy/fem2d/fem2d.py
# @project : fempy
# software : PyCharm

from itertools import product
from fempy.basic import fslove, Gaussian
from .basic import *


class FEM2D(object):
    def __init__(self, eq, mesh):
        self.eq = eq
        self.f = eq['f']
        print("> Constructing triangle mesh...")
        self.mesh = mesh
        print('> point:', self.mesh.p_n, ', edge:', self.mesh.e_n, ', tri:', self.mesh.u_n)
        print("> Load Gauss point...")
        self.temp = Gaussian('2', 3)
        self.g_p = self.temp.gauss_p
        self.g_w = self.temp.gauss_w
        self.g_a = self.temp.gauss_a
        print("> Init Matrix A, F, U...")
        self.u_lst = np.zeros(self.mesh.p_n)
        self.f_lst = np.zeros(self.mesh.p_n)
        self.a_mat = np.zeros((self.mesh.p_n, self.mesh.p_n))

    def assembly_af(self):
        r"""组装总矩阵

        Returns
        -------

        """
        for k, v in enumerate(self.mesh.u_mat):
            # constructing temporary matrix b and f
            a_ind = np.array(list(product(v, repeat=2))).T
            b_mat, f_mat = self.construct_af(self.mesh.p_mat[v])
            self.f_lst[v] += f_mat
            self.a_mat[a_ind[0], a_ind[1]] += b_mat.flatten()

    def construct_af(self, unit_v):
        r"""构造单元矩阵

        .. math::
            \iint_K f\cdot v\mathrm{d}x\mathrm{d}y =
            \sum_{i=1}^n f(p_i) \cdot \phi(p_i) w_i =
            (f, \phi_i)

         .. math::
            \iint_K \nabla u\cdot\nabla v\mathrm{d}x\mathrm{d}y =
            \sum_{i=1}^n (\nabla u(p_i), \nabla u(p_i)) u_j

        Parameters
        ----------
        unit_v : array_like, (3, 2)
            三角单元顶点坐标.

        Returns
        -------
        a_elem : ndarray, (3, 3)
            单元刚度矩阵.
        f_elem : ndarray, (3, )
            单元载荷矩阵.

        Note
        ----
        此处:math:`A_i`的计算方法只针对
        右端项为:math:`\Delta u`的情形.
        """
        f_elem = np.zeros(3)
        a_elem = np.zeros((3, 3))
        gauss_global_p = local_to_global(self.g_p, unit_v)
        gauss_global_w = area(*unit_v) / self.g_a * self.g_w
        basis_v = basis_value(gauss_global_p, unit_v)
        basis_g = basis_grid(gauss_global_p, unit_v)
        f_elem += np.dot(self.f(*gauss_global_p.T), basis_v.T) * gauss_global_w
        a_elem += np.dot(basis_g, basis_g.T) * np.sum(gauss_global_w)
        return a_elem, f_elem

    def singular_condition(self):
        r"""边界条件

        Returns
        -------

        Note
        ----
        此处边界条件仅为第一类(Dirichlet)边界条件.
        """
        for k, v in enumerate(self.mesh.p_bnd):
            if v:
                self.f_lst[k] = 0
                self.a_mat[k] = 0
                self.a_mat[k, k] = 1

    def error_l2(self, u_true):
        r"""计算L2误差

        Returns
        -------

        """
        error2 = 0
        for k, v in enumerate(self.mesh.u_mat):
            u_p = self.mesh.p_mat[v]
            gauss_global_p = local_to_global(self.g_p, u_p)
            gauss_global_w = area(*u_p) / self.g_a * self.g_w
            basis_v = basis_value(gauss_global_p, u_p)
            value_true = u_true(*gauss_global_p.T)
            value_calc = np.dot(self.u_lst[v], basis_v)
            error2 += np.sum((value_true - value_calc)**2 * gauss_global_w)
        return np.sqrt(error2)

    def run(self):
        """
        2 Dimension Finite Element Method.
        """
        print("> Assembly Matrix F...")
        # self.assembly_f()
        print("> Assembly Matrix A...")
        self.assembly_af()
        print("> Apply Boundary Conditions...")
        self.singular_condition()
        print("> Solve Matrix U...")
        self.u_lst = fslove(self.a_mat, self.f_lst)


