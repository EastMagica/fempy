#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/12 17:04
# @file    : fem_1d.py
# @project : fem
# software : PyCharm

import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve


def phi_l(l, r):
    """
    左 $\phi$ 函数生成

    :param l:
    :param r:
    :return:
    """
    return lambda x: (x - r) / (l - r)


def phi_r(l, r):
    """
    右 $\phi$ 函数生成

    :param l:
    :param r:
    :return:
    """
    return lambda x: (x - l) / (r - l)


def phi_l_d1(l, r):
    """
    左 $\phi$ 函数一阶导数生成
    :param l:
    :param r:
    :return:
    """
    return lambda x: 1 / (l - r)


def phi_r_d1(l, r):
    """
    右 $\phi$ 函数一阶导数生成
    :param l:
    :param r:
    :return:
    """
    return lambda x: 1 / (r - l)


def _interval_split(domain, opt):
    """
    剖分计算区间，提供4种不同的划分方式：
        1. 'step': 分隔区间数量（均匀剖分）
        2. 'h'   : 分隔区间步长（均匀剖分）
        3. 'f'   : 函数划分区间（非均匀剖分）
        4. 'lst' : 给定分隔区间（非均匀剖分）

    :param domain:
    :param opt:
    :return:
    """
    if opt[0] == 'step':
        x_lst = np.linspace(domain[0], domain[1], opt[1] + 1)
    elif opt[0] == 'h':
        x_lst = np.arange(domain[0], domain[1], opt[1])
    elif opt[0] == 'f':
        x_lst = opt[1](domain)
    elif opt[0] == 'lst':
        x_lst = np.array(opt[1])
    else:
        raise ValueError
    return x_lst


def inner_product(f0, f1, a, b):
    """
    计算两函数的内积
    $$\int^b_a f_1(x)f_2(x) \mathrm{d}x$$

    :param f0:
    :param f1:
    :param a:
    :param b:
    :return:
    """
    inner_solve = quad(lambda x: f0(x) * f1(x), a, b)
    return inner_solve[0]


def fslove(a_mat, f_lst):
    """
    求解线性方程组 $AU=F$

    :param a_mat:
    :param f_lst:
    :return:
    """
    return solve(a_mat, f_lst)


def construct_cor(domain, opt):
    """
    构造坐标矩阵、位置矩阵

    :return x_lst:
    :return e_mat:
    """
    x_lst = _interval_split(domain, opt)
    n = len(x_lst)
    e_mat = np.transpose(np.vstack((np.arange(0, n - 1), np.arange(1, n))))
    return x_lst, e_mat, n


def error_l2(u_lst, u_true, x_lst, e_mat):
    e2 = 0
    for k, v in enumerate(e_mat):
        l_inx, r_inx = v
        l_cor, r_cor = x_lst[[l_inx, r_inx]]
        e2 += quad(lambda x: (u_true(x) - u_lst[l_inx] * phi_l(l_cor, r_cor)(x)
                              - u_lst[r_inx] * phi_r(l_cor, r_cor)(x)) ** 2, l_cor, r_cor)[0]
    return np.sqrt(e2)


class FEM1D(object):
    def __init__(self, eq):
        self.f = eq['f']
        self.bnd = eq['bnd']
        self.x_lst, self.e_mat, self.n = construct_cor(eq['domain'], eq['split'])
        self.a_mat = np.zeros((self.n, self.n))
        self.f_lst = np.zeros(self.n)
        self.u_lst = np.zeros(self.n)
        print(self.x_lst, '\n', self.e_mat)

    def assembly_f(self):
        """
        组装矩阵 $F$
        """
        for k, v in enumerate(self.e_mat):
            l_inx, r_inx = v
            l_cor, r_cor = self.x_lst[[l_inx, r_inx]]
            self.f_lst[l_inx] += inner_product(self.f, phi_l(l_cor, r_cor), l_cor, r_cor)
            self.f_lst[r_inx] += inner_product(self.f, phi_r(l_cor, r_cor), l_cor, r_cor)

    def assembly_a(self):
        """
        组装矩阵 $A$
        """
        for k, v in enumerate(self.e_mat):
            l_inx, r_inx = v
            l_cor, r_cor = self.x_lst[[l_inx, r_inx]]
            # TODO: package other forms function
            self.a_mat[l_inx][l_inx] += inner_product(phi_l_d1(l_cor, r_cor),
                                                      phi_l_d1(l_cor, r_cor), l_cor, r_cor)
            self.a_mat[l_inx][r_inx] += inner_product(phi_r_d1(l_cor, r_cor),
                                                      phi_l_d1(l_cor, r_cor), l_cor, r_cor)
            self.a_mat[r_inx][l_inx] += inner_product(phi_l_d1(l_cor, r_cor),
                                                      phi_r_d1(l_cor, r_cor), l_cor, r_cor)
            self.a_mat[r_inx][r_inx] += inner_product(phi_r_d1(l_cor, r_cor),
                                                      phi_r_d1(l_cor, r_cor), l_cor, r_cor)

    def singular_condition(self):
        """
        添加边界条件，修改A、F矩阵
        """
        self.a_mat[[0, -1], [0, -1]] = 1
        self.a_mat[[0, -1, 1, -2], [1, -2, 0, -1]] = 0
        self.f_lst[[0, -1]] = self.bnd
        self.f_lst[[1, -2]] -= self.a_mat[[1, -2], [0, -1]] * self.bnd

    def run(self):
        """
        一维有限元方法
        """
        print("> Assembly Matrix F...")
        self.assembly_f()
        print("> Assembly Matrix A...")
        self.assembly_a()
        print("> Apply Boundary Conditions...")
        self.singular_condition()
        print("> Solve Matrix U...")
        self.u_lst = fslove(self.a_mat, self.f_lst)
