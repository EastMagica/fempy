#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 18:01
# @file    : fempy/basic.py
# @project : fempy
# software : PyCharm

from scipy.linalg import solve


def fslove(a_mat, f_lst):
    """
    求解线性方程组 $AU=F$

    :param a_mat:
    :param f_lst:
    :return:
    """
    return solve(a_mat, f_lst)