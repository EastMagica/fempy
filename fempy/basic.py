#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 18:01
# @file    : fempy/basic.py
# @project : fempy
# software : PyCharm

import os
import abc
import numpy as np
from scipy.linalg import solve
from .gaussian import GaussianTemp

__project_name__ = "fempy"


def _is_ndarray(*args):
    """类型转换: ndarray类型"""
    ndarray_list = [np.array(arr) for arr in args]
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


class FEM(abc.ABC):
    @abc.abstractmethod
    def __init__(self, variation, mesh):
        self.mesh = mesh
        self.variation = variation
        self.gaussian = None
        self.u_lst = np.zeros(self.mesh.npoints)
        self.f_lst = np.zeros(self.mesh.npoints)
        self.a_mat = np.zeros((self.mesh.npoints, self.mesh.npoints))

    @staticmethod
    @abc.abstractmethod
    def basis_value(p, v):
        raise NotImplementedError

    @staticmethod
    @abc.abstractmethod
    def basis_grid(p, v):
        raise NotImplementedError

    @abc.abstractmethod
    def boundary(self):
        raise NotImplementedError

    def assembly_af(self):
        r"""组装总矩阵"""
        for k, v in enumerate(self.mesh.simplices):
            xm, ym = np.meshgrid(v, v)
            a_elem, f_elem = self.construct_af(self.mesh.points[v])
            self.f_lst[v] += f_elem
            self.a_mat[xm, ym] += a_elem

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
        gauss_p, gauss_w = self.gaussian.local_to_global(unit_v)
        basis_v = self.basis_value(gauss_p, unit_v)
        basis_g = self.basis_grid(gauss_p, unit_v)
        a_elem, f_elem = self.variation(basis_v, basis_g, gauss_p, gauss_w)
        return a_elem, f_elem

    def simplices_error(self, u_true):
        error_lst = np.zeros(self.mesh.nsimplices)
        for k, v in enumerate(self.mesh.simplices):
            unit_v = self.mesh.points[v]
            gauss_p, gauss_w = self.gaussian.local_to_global(unit_v)
            basis_v = self.basis_value(gauss_p, unit_v)
            value_true = u_true(*gauss_p.T)
            value_calc = np.dot(self.u_lst[v], basis_v.T)
            error_lst[k] = np.sum((value_true - value_calc) ** 2 * gauss_w)
        return error_lst

    def error_l2(self, u_true):
        r"""计算L2误差

        Parameters
        ----------
        u_true

        Returns
        -------

        """
        error2 = np.sum(self.simplices_error(u_true))
        return np.sqrt(error2)

    def save_data(self, file_name='data', form='npy', path=''):
        r"""
        存储计算数据.

        Parameters
        ----------
        file_name
        form
        path

        Returns
        -------

        """
        path = path if path else os.getcwd()
        if form.lower() == 'npy':
            np.save(file_name + '_a_mat', self.a_mat)
            np.save(file_name + '_f_lst', self.f_lst)
            np.save(file_name + '_u_lst', self.u_lst)
        elif form.lower() == 'txt':
            np.savetxt(file_name + '_a_mat.txt', self.a_mat, delimiter=',', newline=',\n')
            np.savetxt(file_name + '_f_lst.txt', self.f_lst, delimiter=',', newline=',\n')
            np.savetxt(file_name + '_u_lst.txt', self.u_lst, delimiter=',', newline=',\n')
        else:
            raise ValueError('Format ' + form + ' is not supported.')

    def run(self):
        r"""
        计算.

        Returns
        -------

        """
        print("> Assembly Matrix A and F...")
        self.assembly_af()
        # print(f'> a_mat0:\n{self.a_mat}\n> f_lst0:\n{self.f_lst}')
        print("> Apply Boundary Conditions...")
        self.boundary()
        # print(f'> a_mat:\n{self.a_mat}\n> f_lst:\n{self.f_lst}')
        print("> Solve Matrix U...")
        self.u_lst = fslove(self.a_mat, self.f_lst)
