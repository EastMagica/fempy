#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/21 16:38
# @file    : mesh1d.py
# @project : fempy
# software : PyCharm

import numpy as np


class Mesh1D(object):
    def __init__(self, method='linspace', opt=None):
        if method == 'linspace':
            self.__points = np.linspace(*opt)
            self.__p_bnd = (self.__points == 0) | (self.__points == 1)
            self.__simplices = np.array([[i, i+1] for i in range(len(self.__points)-1)])
        else:
            raise ValueError('Other methods are not supported.')

    @property
    def points(self):
        return self.__points

    @property
    def simplices(self):
        return self.__simplices

    @property
    def npoints(self):
        return len(self.__points)

    @property
    def nsimplices(self):
        return len(self.__simplices)

    @property
    def boundary(self):
        return self.__p_bnd
