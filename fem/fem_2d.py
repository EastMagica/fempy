#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/14 18:00
# @file    : fem_2d.py
# @project : fem
# software : PyCharm

import numpy as np
from scipy.integrate import quad
from fem.basic import fslove


def phi(p0, p1, p2):
    delta = p0[0] * (p1[1] - p2[1]) + \
            p1[0] * (p2[1] - p0[1]) + \
            p2[0] * (p0[1] - p1[1])
    return delta


def phi_2():
    pass


def phi_3():
    pass


def phi_1_d1():
    pass


def phi_2_d1():
    pass


def phi_3_d1():
    pass


def _interval_split():
    pass


def inner_product_2d(f0, f1, a, b, c):
    pass


def construct_cor():
    pass


def error_l2():
    pass


class FEM1D(object):
    def __init__(self):
        pass

    def assembly_f(self):
        pass

    def assembly_a(self):
        pass

    def singular_condition(self):
        pass

    def run(self):
        """
        二维有限元方法
        """
        print("> Assembly Matrix F...")
        self.assembly_f()
        print("> Assembly Matrix A...")
        self.assembly_a()
        print("> Apply Boundary Conditions...")
        self.singular_condition()
        print("> Solve Matrix U...")



