#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 18:01
# @file    : fempy/basic.py
# @project : fempy
# software : PyCharm

import os
import numpy as np
from scipy.linalg import solve


def _is_ndarray(*args):
    """类型转换: ndarray类型"""
    ndarray_list = [np.asarray(array) for array in args]
    if len(ndarray_list) == 1:
        return ndarray_list[0]
    return ndarray_list


def fslove(a_mat, f_lst):
    """线性方程组求解器

    Parameters
    ----------
    a_mat
    f_lst

    Returns
    -------

    Notes
    -----
    针对不同矩阵选择最优解法.

    """
    return solve(a_mat, f_lst)


class Gaussian(object):
    r"""高斯积分点

    获取n次高斯积分点坐标及权重.

    Attributes
    ----------
    gauss : array_like
        高斯积分点数据.

    Note
    ----
    数据来源.

    """

    def __init__(self, dim, n):
        gauss = self.get_gauss(dim, n)
        self.gauss_a = 2 if dim == 1 else 2
        self.gauss_p = np.transpose(gauss[:-1])
        self.gauss_w = gauss[-1]

    @staticmethod
    def get_gauss(dim, n):
        dim = str(dim)
        if dim in ['1', '2']:
            file_name = 'res/gauss/gauss'+dim+'d/gauss' + str(n) + '.npy'
        else:
            raise ValueError(dim+'-D Gaussian is not supported.')
        print('> load', file_name)
        if os.path.exists(file_name):
            gauss = np.load(file_name)
            return gauss
        else:
            raise ValueError(str(n) + ' Gauss Point does not exists!')

    @staticmethod
    def get_area(dim):
        if str(dim) == '1':
            return 2
        elif str(dim) == '2':
            return 1

    @staticmethod
    def get_gauss_list():
        return [1, 3, 4, 6, 7, 12, 13, 25, 33]