#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : east
# @time    : 2019/10/10 14:02
# @file    : gaussian.py
# @project : fem
# software : PyCharm

import os
import abc
from numpy import load

__project_name__ = "fempy"

_dim_list = ['1', '2']
_gaussian_list = {'1': ['1', '2', '3', '4', '5'],
                  '2': ['1', '3', '4', '6', '7', '12', '13', '25', '33']}
_gaussian_area = {'1': 2, '2': 1 / 2}


def find_gaussian(n, dim):
    r"""
    解析高斯积分参数文件路径.

    Parameters
    ----------
    n
    dim

    Returns
    -------

    """
    project_path = os.path.abspath(__file__).split(__project_name__)[0]
    relative_path = 'res/gauss/gauss' + dim + 'd/gauss' + str(n) + '.npy'
    absolute_path = os.path.join(project_path, relative_path)
    print('> load gaussian: ', absolute_path)
    return absolute_path


def load_gaussian(n, dim):
    r"""
    加载高斯积分参数.

    Parameters
    ----------
    n : int, str
    dim : int, str

    Returns
    -------
    points : ndarray, (N, )
    weight : ndarray, (N, )
    area : float

    """
    n = str(n) if not isinstance(n, str) else n
    dim = str(dim) if not isinstance(dim, str) else dim
    if dim not in _dim_list:
        raise ValueError(dim + '-D Gaussian is not supported.')
    if n not in _gaussian_list[dim]:
        raise ValueError(str(n) + ' Gauss Point does not exists!')
    file_path = find_gaussian(n, dim)
    if not os.path.exists(file_path):
        raise ValueError('gauss' + str(n) + '.npy is not exist.')
    gaussian = load(file_path)
    points = gaussian[:-1]
    weight = gaussian[-1]
    area = _gaussian_area[dim]
    points = points.T if dim == '2' else points.reshape(-1, 1)
    return points, weight, area


class GaussianTemp(abc.ABC):
    """

    Methods
    ----------
    points : ndarray, (N, dim)
    weight : ndarray, (N, )
    area : float
    """
    def __init__(self, n, dim):
        self.__dim = dim
        self.__gauss_n = n
        gaussian = load_gaussian(n, dim)
        self.__points, self.__weight, self.__area = gaussian

    @property
    def gauss_n(self):
        return self.__gauss_n

    @property
    def dim(self):
        return self.__dim

    @property
    def points(self):
        return self.__points

    @property
    def weight(self):
        return self.__weight

    @property
    def area(self):
        return self.__area

    @abc.abstractmethod
    def local_to_global(self, global_v):
        raise NotImplementedError
