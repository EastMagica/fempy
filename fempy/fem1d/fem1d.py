#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/12 17:04
# @file    : fempy/fem1d/fem1d.py
# @project : fempy
# software : PyCharm

from itertools import product
from scipy.integrate import quad
from fempy.basic import fslove, Gaussian
from .basic import *


class FEM1D(object):
    def __init__(self, eq, mesh):
        self.f = eq['f']
        self.bnd = eq['bnd']
        self.mesh = mesh
        self.temp = Gaussian('1', 3)
        self.g_p = self.temp.gauss_p
        self.g_w = self.temp.gauss_w
        self.g_a = self.temp.gauss_a
        self.a_mat = np.zeros((self.mesh.p_n, self.mesh.p_n))
        self.f_lst = np.zeros(self.mesh.p_n)
        self.u_lst = np.zeros(self.mesh.p_n)

    def assembly_af(self):
        for k, v in enumerate(self.mesh.e_mat):
            a_ind = np.array(list(product(v, repeat=2))).T
            a_elem, f_elem = self.construct_af(self.mesh.p_mat[v])
            self.f_lst[v] += f_elem
            self.a_mat[a_ind] += a_elem

    def construct_af(self, unit_v):
        f_elem = np.zeros(2)
        a_elem = np.zeros((2, 2))
        gauss_global_p = local_to_global(self.g_p, unit_v)
        gauss_global_w = (unit_v[1] - unit_v[0]) / self.g_a * self.g_w
        basis_v = basis_value(gauss_global_p, unit_v)
        basis_g = basis_grid(gauss_global_p, unit_v)
        f_elem += np.dot(self.f(gauss_global_p), basis_v.T) * gauss_global_w
        a_elem += np.dot(basis_g, basis_g.T) * np.sum(gauss_global_w)
        return a_elem, f_elem

    def singular_condition(self):
        """
        添加边界条件，修改A、F矩阵
        TODO: 添加第二、第三边界条件的适应方法
        """
        p_end = self.mesh.p_end
        print(f'> p_end: {p_end}')
        self.f_lst[p_end] = self.bnd
        self.a_mat[p_end, :] = 0
        self.a_mat[p_end, p_end] = 1

    def error_l2(self, u_true):
        """
        计算 $L^2$ 误差
        :param u_true:
        :return:
        """
        error2 = 0
        for k, v in enumerate(self.mesh.u_mat):
            unit_v = self.mesh.p_mat[v]
            gauss_global_p = local_to_global(self.g_p, unit_v)
            gauss_global_w = (unit_v[1] - unit_v[0]) / self.g_a * self.g_w
            basis_v = basis_value(gauss_global_p, unit_v)
            value_true = u_true(*gauss_global_p.T)
            value_calc = np.dot(self.u_lst[v], basis_v)
            error2 += np.sum((value_true - value_calc) ** 2 * gauss_global_w)
        return np.sqrt(error2)

    def run(self):
        """
        2 Dimension Finite Element Method.
        """
        print("> Assembly Matrix A and F...")
        self.assembly_af()
        print(f'> a_mat0:\n{self.a_mat}\n> f_lst0:\n{self.f_lst}')
        print("> Apply Boundary Conditions...")
        self.singular_condition()
        print(f'> a_mat:\n{self.a_mat}\n> f_lst:\n{self.f_lst}')
        print("> Solve Matrix U...")
        self.u_lst = fslove(self.a_mat, self.f_lst)


class Adaptive1D(object):
    def __init__(self, eq, mesh):
        self.eq = eq
        self.f = eq['f']
        self.mesh = mesh
        print('> points:', self.mesh.p_mat)
        print('> edges:', self.mesh.e_mat)
        self.fem = FEM1D(eq, mesh)
        self.fem.run()
        print('> u_h:', self.fem.u_lst)
        # print('> el2:', self.fem.error_l2(self.eq['u_true']))
        # plot_1d(self.fem.u_lst, self.mesh, self.eq['u_true'], 'start plot')

    def fem_iter(self, n=1):
        for i in range(n):
            print(f"\n\n{'-'*32}\n> {i+1}'st Split\n{'-'*32}")
            eta_h = self.estimate_error()
            self.adaptive_mesh(eta_h)
            print('> points:', self.mesh.p_mat)
            print('> edges:', self.mesh.e_mat)
            self.fem = FEM1D(self.eq, self.mesh)
            self.fem.run()
            print('> u_h:', self.fem.u_lst)
            # print('> el2:', self.fem.error_l2(self.eq['u_true']))
            # plot_1d(self.fem.u_lst, self.mesh, self.eq['u_true'], str(i+1)+"'st Split")

    def estimate_error(self):
        eta_h = np.zeros(self.mesh.e_n)
        for k, v in enumerate(self.mesh.e_mat):
            l_cor, r_cor = self.mesh.p_mat[v]
            eta_h[k] = (r_cor - l_cor) * np.sqrt(quad(lambda x: self.f(x)**2, l_cor, r_cor)[0])
        return eta_h

    def adaptive_mesh(self, eta_h):
        inx = 0
        eta_max = np.max(eta_h)
        eta_num = np.where(eta_h > eta_max * 0.75)[0]
        lst = np.zeros((self.mesh.e_n + len(eta_num), 2), dtype=np.int)
        print('> eta_num:', eta_num)
        self.mesh.p_mat = np.append(self.mesh.p_mat, np.sum(self.mesh.p_mat[self.mesh.e_mat[eta_num]], axis=1)/2)
        for k in range(self.mesh.e_n):
            lst[k+inx, 0] = self.mesh.e_mat[k, 0]
            if k in eta_num:
                lst[k+inx, 1] = self.mesh.p_n + inx
                lst[k+inx+1, 0] = self.mesh.p_n + inx
                inx += 1
            lst[k+inx, 1] = self.mesh.e_mat[k, 1]
        self.mesh.e_mat = np.array(lst)
        self.mesh.p_n = len(self.mesh.p_mat)
        self.mesh.e_n = len(self.mesh.e_mat)
        self.mesh.p_end = np.append(self.mesh.p_end, [False]*len(eta_num))
