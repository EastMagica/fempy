#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/12 17:04
# @file    : fem_1d.py
# @project : fem
# software : PyCharm

import numpy as np
import scipy as sp
from scipy.integrate import quad


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


def interval_split():
    return True


def inner_product(f0, f1, a, b):
    """
    计算 $L2$ 空间上的内积
    $$\int^1_0 f_1(x)f_2(x) \mathrm{d}x$$
    :param f0:
    :param f1:
    :param a:
    :param b:
    :return:
    """
    inner_solve = quad(lambda x: f0(x) * f1(x), a, b)
    return inner_solve[0]


def assembly_f(x_lst, e_mat, f, n):
    """
    组装矩阵 $F$
    :param x_lst:
    :param e_mat:
    :param f:
    :param n:
    :return:
    """
    f_lst = np.zeros(n)
    for i in range(n):
        l_inx, r_inx = e_mat[i]
        l_cor = x_lst[l_inx]
        r_cor = x_lst[r_inx]
        f_lst[l_inx] += inner_product(f, phi_l(l_cor, r_cor), l_cor, r_cor)
        f_lst[r_inx] += inner_product(f, phi_r(l_cor, r_cor), l_cor, r_cor)
    return f_lst


def assembly_a(x_lst, e_mat, n):
    """
    组装矩阵 $A$
    :param x_lst:
    :param e_mat:
    :param n:
    :return:
    """
    a_mat = np.zeros((n, n))
    for i in range(n):
        l_inx, r_inx = e_mat[i]
        l_cor = x_lst[l_inx]
        r_cor = x_lst[r_inx]
        a_mat[l_inx][l_inx] += inner_product(phi_l(l_cor, r_cor),
                                             phi_l(l_cor, r_cor), l_cor, r_cor)
        a_mat[l_inx][r_inx] += inner_product(phi_r(l_cor, r_cor),
                                             phi_l(l_cor, r_cor), l_cor, r_cor)
        a_mat[r_inx][l_inx] += inner_product(phi_l(l_cor, r_cor),
                                             phi_r(l_cor, r_cor), l_cor, r_cor)
        a_mat[r_inx][r_inx] += inner_product(phi_r(l_cor, r_cor),
                                             phi_r(l_cor, r_cor), l_cor, r_cor)
    return a_mat


def fslove(a_mat, f_lst):
    """
    求解线性方程组 $AU=F$
    :param a_mat:
    :param f_lst:
    :return:
    """
    pass
