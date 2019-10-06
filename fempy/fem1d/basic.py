#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/9/29 12:23
# @file    : basic.py
# @project : fem
# software : PyCharm

import numpy as np


def basis_value(p, v):
    """

    Parameters
    ----------
    p
    v

    Returns
    -------

    """
    return np.array([p - v[1], v[0] - p]) / (v[0] - v[1])


def basis_grid(p, v):
    """

    Parameters
    ----------
    p
    v

    Returns
    -------

    """
    return np.tile([[1], [-1]], (1, p.shape[0])) / (v[0] - v[1])


def local_to_global(local_p, global_v):
    """

    Parameters
    ----------
    local_p : array_like, (N, )
    global_v : array_like, (2, )

    Returns
    -------

    """
    return local_p * (global_v[1] - global_v[0]) + global_v[0]
