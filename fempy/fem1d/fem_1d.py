#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/12 17:04
# @file    : fempy/fem1d/fem_1d.py
# @project : fempy
# software : PyCharm

import numpy as np
from scipy.integrate import quad
from fempy.basic import fslove
from fempy.fem1d.fem_1d_plot import plot_1d


def phi_l(l, r):
    """
    Left $\phi$ function
    :param l:
    :param r:
    :return:
    """
    return lambda x: (x - r) / (l - r)


def phi_r(l, r):
    """
    Right $\phi$ function
    :param l:
    :param r:
    :return:
    """
    return lambda x: (x - l) / (r - l)


def phi_l_d1(l, r):
    """
    The first erivation of left $\phi$
    :param l:
    :param r:
    :return:
    """
    return lambda x: 1 / (l - r)


def phi_r_d1(l, r):
    """
    The first erivation of right $\phi$
    :param l:
    :param r:
    :return:
    """
    return lambda x: 1 / (r - l)


def inner_product_1d(f0, f1, a, b):
    """
    Inner product of 2 function
    $$\int^b_a f_1(x)f_2(x) \mathrm{d}x$$
    :param f0:
    :param f1:
    :param a:
    :param b:
    :return:
    """
    inner_solve = quad(lambda x: f0(x) * f1(x), a, b)
    return inner_solve[0]


class FEM1D(object):
    def __init__(self, eq, mesh):
        self.f = eq['f']
        self.bnd = eq['bnd']
        self.mesh = mesh
        self.a_mat = np.zeros((self.mesh.p_n, self.mesh.p_n))
        self.f_lst = np.zeros(self.mesh.p_n)
        self.u_lst = np.zeros(self.mesh.p_n)

    def assembly_f(self):
        """
        组装矩阵 $F$
        """
        for k, v in enumerate(self.mesh.e_mat):
            l_inx, r_inx = v
            l_cor, r_cor = self.mesh.p_mat[[l_inx, r_inx]]
            self.f_lst[l_inx] += inner_product_1d(self.f, phi_l(l_cor, r_cor), l_cor, r_cor)
            self.f_lst[r_inx] += inner_product_1d(self.f, phi_r(l_cor, r_cor), l_cor, r_cor)

    def assembly_a(self):
        """
        组装矩阵 $A$
        """
        for k, v in enumerate(self.mesh.e_mat):
            l_inx, r_inx = v
            l_cor, r_cor = self.mesh.p_mat[[l_inx, r_inx]]
            # TODO: package other forms function
            self.a_mat[l_inx][l_inx] += inner_product_1d(phi_l_d1(l_cor, r_cor),
                                                         phi_l_d1(l_cor, r_cor), l_cor, r_cor)
            self.a_mat[l_inx][r_inx] += inner_product_1d(phi_r_d1(l_cor, r_cor),
                                                         phi_l_d1(l_cor, r_cor), l_cor, r_cor)
            self.a_mat[r_inx][l_inx] += inner_product_1d(phi_l_d1(l_cor, r_cor),
                                                         phi_r_d1(l_cor, r_cor), l_cor, r_cor)
            self.a_mat[r_inx][r_inx] += inner_product_1d(phi_r_d1(l_cor, r_cor),
                                                         phi_r_d1(l_cor, r_cor), l_cor, r_cor)

    def singular_condition(self):
        """
        添加边界条件，修改A、F矩阵
        TODO: 添加第二、第三边界条件的适应方法
        """
        p_end = self.mesh.p_end
        print(f'> p_end: {p_end}')
        self.f_lst[p_end] = self.bnd
        self.a_mat[p_end, :] = 0
        # self.a_mat[:, p_end] = 0
        self.a_mat[p_end, p_end] = 1
        # self.f_lst[[1, -2]] -= self.a_mat[[1, -2], [0, -1]] * self.bnd

    def error_l2(self, u_true):
        """
        计算 $L^2$ 误差
        :param u_true:
        :return:
        """
        e2 = 0
        for k, v in enumerate(self.mesh.e_mat):
            l_inx, r_inx = v
            l_cor, r_cor = self.mesh.p_mat[v]
            e2 += quad(lambda x: (u_true(x) - self.u_lst[l_inx] * phi_l(l_cor, r_cor)(x)
                                  - self.u_lst[r_inx] * phi_r(l_cor, r_cor)(x)) ** 2, l_cor, r_cor)[0]
        return np.sqrt(e2)

    def run(self):
        """
        2 Dimension Finite Element Method.
        """
        print("> Assembly Matrix F...")
        self.assembly_f()
        print("> Assembly Matrix A...")
        self.assembly_a()
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
        print('> u_h:', self.fem.u_lst)
        self.fem = FEM1D(eq, mesh)
        self.fem.run()
        print('> u_h:', self.fem.u_lst)
        print('> el2:', self.fem.error_l2(self.eq['u_true']))
        plot_1d(self.fem.u_lst, self.mesh, self.eq['u_true'], 'start plot')

    def fem_iter(self, n=1):
        for i in range(n):
            eta_h = self.estimate_error()
            self.adaptive_mesh(eta_h)
            print('> points:', self.mesh.p_mat)
            print('> edges:', self.mesh.e_mat)
            self.fem = FEM1D(self.eq, self.mesh)
            self.fem.run()
            print('> u_h:', self.fem.u_lst)
            print('> el2:', self.fem.error_l2(self.eq['u_true']))
            plot_1d(self.fem.u_lst, self.mesh, self.eq['u_true'], str(i+1)+"'st Split")

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
        # for k, v in enumerate(eta_num):
        #     lst.append([k + self.mesh.p_n, self.mesh.e_mat[v, 1]])
        #     self.mesh.e_mat[v, 1] = k + self.mesh.p_n
        # self.mesh.e_mat = np.vstack((self.mesh.e_mat, lst))
        for k in range(self.mesh.e_n):
            print(k, inx)
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
