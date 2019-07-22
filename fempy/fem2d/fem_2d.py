#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 18:00
# @file    : fempy/fem2d/fem_2d.py
# @project : fempy
# software : PyCharm

import numpy as np
from scipy.integrate import quad, dblquad
from itertools import product as prod
from fempy.basic import fslove
from fempy.fem2d.mesh_2d import Mesh2D
from fempy.fem2d.basic import *


def inner_product_2d(f, a0, b0, a1, b1):
    return dblquad(f, a0, b0, a1, b1)


def construct_b(p):
    # TODO: Complete function
    b_mat = np.zeros(9)
    sa = (p[0, 0] - p[2, 0]) * (p[1, 1] - p[2, 1]) - (p[1, 0] - p[2, 0]) * (p[0, 1] - p[2, 1]) / 2
    b_mat[0] = ((p[1, 1] - p[2, 1]) ** 2 + (p[2, 0] - p[1, 0]) ** 2) / (2 * sa)
    b_mat[1] = ((p[2, 1] - p[0, 1]) * (p[1, 1] - p[2, 1]) +
                (p[0, 0] - p[2, 0]) * (p[2, 0] - p[1, 0])) / (2 * sa)
    b_mat[2] = ((p[0, 1] - p[1, 1]) * (p[1, 1] - p[2, 1]) +
                (p[1, 0] - p[0, 0]) * (p[2, 0] - p[1, 0])) / (2 * sa)
    b_mat[3] = b_mat[1]
    b_mat[4] = ((p[2, 1] - p[0, 1]) ** 2 + (p[0, 0] - p[2, 0]) ** 2) / (2 * sa)
    b_mat[5] = ((p[0, 1] - p[1, 1]) * (p[2, 1] - p[0, 1]) +
                (p[1, 0] - p[0, 0]) * (p[0, 0] - p[2, 0])) / (2 * sa)
    b_mat[6] = b_mat[2]
    b_mat[7] = b_mat[5]
    b_mat[8] = ((p[0, 1] - p[1, 1]) ** 2 + (p[1, 0] - p[0, 0]) ** 2) / (2 * sa)
    return b_mat


def construct_f(f, p):
    # TODO: Complete function
    f_lst = np.zeros(3)
    det_j = (p[0, 0] - p[2, 0]) * (p[1, 1] - p[2, 1]) - (p[1, 0] - p[2, 0]) * (p[0, 1] - p[2, 1])
    fun = lambda x, y: f(p[2, 0] + x * (p[0, 0] - p[2, 0]) + y * (p[1, 0] - p[2, 0]),
                         p[2, 1] + x * (p[0, 1] - p[2, 1]) + y * (p[1, 1] - p[2, 1]))
    f_lst[0] = dblquad(lambda x, y: fun(x, y) * x * det_j, 0, 1, lambda x: 0, lambda x: 1 - x)[0]
    f_lst[1] = dblquad(lambda x, y: fun(x, y) * y * det_j, 0, 1, lambda x: 0, lambda x: 1 - x)[0]
    f_lst[2] = dblquad(lambda x, y: fun(x, y) * (1 - x - y) * det_j, 0, 1, lambda x: 0, lambda x: 1 - x)[0]
    return f_lst


class FEM2D(object):
    def __init__(self, eq, mesh):
        self.eq = eq
        self.f = eq['f']
        # generate triangular mesh
        print("> Constructing triangle mesh...")
        self.mesh = mesh
        print('> point:', self.mesh.p_n, ', edge:', self.mesh.e_n, ', tri:', self.mesh.u_n)
        # select of Gauss points
        # TODO: Jacobian < 0 by Gass Point Problem
        print("> Load Gauss point...")
        self.temp = TempElement(3)
        self.g_s = self.temp.get_area()
        self.g_p = self.temp.get_gauss_point()
        self.g_w = self.temp.get_gauss_weight()
        # init matrix a, u, f
        print("> Init Matrix A, F, U...")
        self.u_lst = np.zeros(self.mesh.p_n)
        self.f_lst = np.zeros(self.mesh.p_n)
        self.a_mat = np.zeros((self.mesh.p_n, self.mesh.p_n))

    def assembly_af_0(self):
        """
        Assembly matrix A and F by general method.
        :return:
        """
        for k, v in enumerate(self.mesh.u_mat):
            # assembly matrix F
            self.f_lst[v] += construct_f(self.f, self.mesh.p_mat[v])
            # assembly matrix A
            b_mat = construct_b(self.mesh.p_mat[v])
            a_ind = np.transpose(np.array(list(prod(v, repeat=2))))
            self.a_mat[a_ind[0], a_ind[1]] += b_mat.flatten()

    def assembly_af_1(self):
        """
        Assembly matrix A and F by Gauss Points.
        :return:
        """
        for k, v in enumerate(self.mesh.triangles):
            # constructing temporary matrix b and f
            b_mat, f_mat = self.construct_af(self.mesh.points[v])
            # assembly matrix F
            self.f_lst[v] += f_mat
            # assembly matrix A
            a_ind = np.transpose(np.array(list(prod(v, repeat=2))))
            self.a_mat[a_ind[0], a_ind[1]] += b_mat.flatten()

    def construct_af(self, p):
        """
        Constructing temporary matrix b and f
        :param p:
        :return:
        """
        a_mat = np.zeros((3, 3))
        f_lst = np.zeros(3)
        # affine gauss point into xoy space
        q_point = local_to_global(self.g_p, p)
        basis_v = basis_value(q_point, p)
        basis_g = basis_grid(q_point, p)
        jacobian = local_to_global_jacobian(p)
        # f mat element (1x3)
        f_lst += np.dot(self.f(q_point[:, 0], q_point[:, 1]), basis_v.T) * self.g_w * jacobian * self.g_s
        # a mat element (3x3)
        for k in range(3):
            for j in range(3):
                a_mat[j, k] += basis_g[j, 0] * basis_g[k, 0] + basis_g[j, 1] * basis_g[k, 1]
        a_mat *= np.sum(self.g_w * jacobian * self.g_s)
        return a_mat, f_lst

    def singular_condition(self):
        """
        Applying boundary conditions
        :return:
        """
        for k, v in enumerate(self.mesh.p_bnd):
            if v:
                self.f_lst[k] = 0
                self.a_mat[k] = 0
                self.a_mat[k, k] = 1

    def error_l2(self):
        error2 = 0
        for k, v in enumerate(self.mesh.u_mat):
            p = self.mesh.p_mat[v]
            q_point = local_to_global(self.g_p, p)
            basis_v = basis_value(q_point, p)
            jacobian = local_to_global_jacobian(p)
            error2 += np.sum((self.eq['u'](q_point[:, 0], q_point[:, 1]) - np.dot(self.u_lst[v], basis_v))**2 * self.g_w * jacobian * self.g_s)
        return np.sqrt(error2)

    def run(self):
        """
        2 Dimension Finite Element Method.
        """
        print("> Assembly Matrix F...")
        # self.assembly_f()
        print("> Assembly Matrix A...")
        self.assembly_af_1()
        print("> Apply Boundary Conditions...")
        self.singular_condition()
        print("> Solve Matrix U...")
        self.u_lst = fslove(self.a_mat, self.f_lst)


