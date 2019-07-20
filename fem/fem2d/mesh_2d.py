#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/16 21:04
# @file    : mesh_2d.py
# @project : fem
# software : PyCharm

import numpy as np
from itertools import product


class Mesh2D(object):
    def __init__(self, xs, xe, ys, ye, nx, ny):
        # TODO: Complete function "split triangle mesh by some point"
        self.p_n = (nx + 1) * (ny + 1)
        self.e_n = 3 * nx * ny + nx * ny
        self.u_n = nx * ny * 2
        # 点、边、元
        self.p_mat = np.zeros((self.p_n, 2))
        self.e_mat = np.zeros((self.e_n, 2), dtype=np.int)
        self.u_mat = np.zeros((self.u_n, 3), dtype=np.int)
        self.u_edg = np.zeros((self.u_n, 3), dtype=np.int)
        # 节点坐标
        self.p_mat = np.array(list(product(np.linspace(xs, xe, nx + 1), np.linspace(ys, ye, ny + 1))))
        self.p_bnd = (self.p_mat[:, 0] == xs) | (self.p_mat[:, 0] == xe) | (self.p_mat[:, 1] == ys) | (self.p_mat[:, 1] == ye)
        # 边坐标
        self.e_bnd = np.zeros(self.e_n, dtype=np.bool)
        for i in range(nx):
            # (ny + 1) * nx
            __ei = [(ny + 1) * i, (ny + 1) * (i + 1)]
            self.e_mat[__ei[0]:__ei[1], :] = np.transpose([np.arange(ny + 1) + __ei[0], np.arange(ny + 1) + __ei[1]])
            self.e_bnd[[0 + (ny + 1) * i, ny + (ny + 1) * i]] = True
        for i in range(ny):
            # (nx + 1) * ny
            __ei = [(nx + 1) * i + (ny + 1) * nx, (nx + 1) * (i + 1) + (ny + 1) * nx]
            self.e_mat[__ei[0]:__ei[1], :] = np.transpose([np.arange(nx + 1) * (ny + 1) + i,
                                                      np.arange(nx + 1) * (ny + 1) + (i + 1)])
            self.e_bnd[[(ny + 1) * nx + (nx + 1) * i, (ny + 1) * nx + nx + (nx + 1) * i]] = True
        for i in range(nx):
            # nx * ny
            __ei = [ny * i + 2 * ny * nx + nx + ny, ny * (i + 1) + 2 * ny * nx + nx + ny]
            self.e_mat[__ei[0]:__ei[1], :] = np.transpose([np.arange(ny) + (ny + 1) * i,
                                                      np.arange(ny) + (ny + 1) * (i + 1) + 1])
        # 单元坐标
        for i in range(1, nx + 1):
            self.u_mat[ny * (i - 1):ny * i, :] = np.transpose([np.arange(ny) + i * (ny + 1),
                                                          np.arange(ny) + 1 + i * (ny + 1),
                                                          np.arange(ny) + (i - 1) * (ny + 1)])
            self.u_edg[ny * (i - 1):ny * i, :] = np.transpose([np.arange(ny) + (i - 1) * ny + (nx + 1) * ny + nx * (ny + 1),
                                                          np.arange(ny) + (i - 1) * (ny + 1),
                                                          np.arange(ny) * (nx + 1) + i + (ny + 1) * nx])
        for i in range(1, nx + 1):
            self.u_mat[nx * ny + ny * (i - 1):nx * ny + ny * i, :] = np.transpose([np.arange(ny) + 1 + (i - 1) * (ny + 1),
                                                                              np.arange(ny) + (i - 1) * (ny + 1),
                                                                              np.arange(ny) + 1 + i * (ny + 1)])
            self.u_edg[nx * ny + ny * (i - 1):nx * ny + ny * i, :] = np.transpose(
                [np.arange(ny) + (i - 1) * ny + (nx + 1) * ny + nx * (ny + 1),
                 np.arange(ny) + 1 + (i - 1) * (ny + 1),
                 np.arange(ny) * (nx + 1) + i - 1 + (ny + 1) * nx])