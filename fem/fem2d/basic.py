#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/18 10:57
# @file    : fem/fem2d/basic.py
# @project : fem
# software : PyCharm

import os
import numpy as np

pnt = np.array([[0, 0], [1, 0], [0, 1]])


def area(p0, p1, p2):
    """

    Tips: if p are array(n, 2), p need use p.T
    :param p0: array_like
    :param p1: array_like
    :param p2:
    :return:
    """
    p0 = np.array(p0).T
    return (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0])


def __jacobian(p, q):
    """

    :param p:
    :param q:
    :return:
    """
    return area(p[0], p[1], p[2]) / area(q[0], q[1], q[2])


def affine_tri(space_s, space_e, point_s):
    """

    :param space_s: Start Affine Space
    :param space_e: Ended Affine Space
    :param point_s: Point in Start Affine Space
    :return point_e: Point in Ended Affine Space
    """
    area_s = area(space_s[0], space_s[1], space_s[2])
    lambda_se = np.array([area(point_s, space_s[1], space_s[2]),
                          area(point_s, space_s[2], space_s[0]),
                          area(point_s, space_s[0], space_s[1])]) / area_s
    point_e = np.dot(lambda_se.T, space_e)
    return point_e


def local_to_global(local_p, global_v):
    """

    :param local_p: Local Point
    :param global_v: Triangle point
    :return : Local Point in Global Area
    """
    return affine_tri(pnt, global_v, local_p)


def global_to_local(global_p, global_v):
    """

    :param global_p: Global Point
    :param global_v: Triangle point
    :return : Global Point in Local Area
    """
    return affine_tri(global_v, pnt, global_p)


def local_to_global_jacobian(global_v):
    """
    :param global_v:
    :return:
    """
    return __jacobian(global_v, pnt)


def global_to_local_jacobian(global_v):
    """
    :param global_v: Triangle point
    :return g_p:
    """
    return __jacobian(pnt, global_v)


def basis_value(p, v):
    """

    :param p:
    :param v:
    :return:
    """
    area_v = area(v[0], v[1], v[2])
    val = np.array([area(p, v[1], v[2]), area(p, v[2], v[0]), area(p, v[0], v[1])]) / area_v
    return val


def basis_grid(p, v):
    """

    :param p:
    :param v:
    :return:
    """
    area_v = area(v[0], v[1], v[2])
    val = np.array([[(v[1, 1] - v[2, 1]), (v[2, 0] - v[1, 0])],
                    [(v[2, 1] - v[0, 1]), (v[0, 0] - v[2, 0])],
                    [(v[0, 1] - v[1, 1]), (v[1, 0] - v[0, 0])]]) / area_v
    return val


class TempElement(object):
    def __init__(self, n):
        file_name = 'res/gauss/gauss' + str(n) + '.npy'
        print('> load', file_name)
        if os.path.exists(file_name):
            self.gauss = np.load(file_name)
        else:
            raise ValueError(str(n) + ' Gauss Point does not exists!')

    def get_area(self):
        return area(pnt[0], pnt[1], pnt[2]) / 2

    def get_gauss_point(self):
        return self.gauss[:2].T

    def get_gauss_weight(self):
        return self.gauss[2]



