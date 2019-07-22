#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/21 16:38
# @file    : mesh_1d.py
# @project : fempy
# software : PyCharm

import numpy as np


def interval_split(domain, opt):
    """
    Partition calculation interval by 4 ways:
        1. 'step': Number of subdivisions
        2. 'h'   : 分隔区间步长（均匀剖分）
        3. 'f'   : 函数划分区间（非均匀剖分）
        4. 'lst' : 给定分隔区间（非均匀剖分）

    :param domain:
    :param opt:
    :return:
    """
    if opt[0] == 'step':
        p_mat = np.linspace(domain[0], domain[1], opt[1] + 1)
    elif opt[0] == 'h':
        p_mat = np.arange(domain[0], domain[1], opt[1])
    elif opt[0] == 'f':
        p_mat = opt[1](domain)
    elif opt[0] == 'lst':
        p_mat = np.array(opt[1])
    else:
        raise ValueError
    return p_mat


class Mesh1D(object):
    def __init__(self, mesh_opt):
        """
        generate 1-D mesh
        :return p_mat: coordinate point matrix
        :return e_mat:
        """
        self.p_mat = interval_split(mesh_opt['domain'], mesh_opt['opt'])
        self.p_n = len(self.p_mat)
        self.e_mat = np.transpose(np.vstack((np.arange(0, self.p_n - 1), np.arange(1, self.p_n))))
        self.e_n = len(self.e_mat)
        self.p_end = np.zeros(self.p_n, dtype=np.bool)
        self.p_end[[0, -1]] = True
