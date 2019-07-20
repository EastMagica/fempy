#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/12 17:49
# @file    : test_1d.py
# @project : fem
# software : PyCharm

import unittest

import numpy as np
import scipy as sp

from scipy.integrate import quad

import fem.fem_1d as fem1d


class Test(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_inner_product(self):
        result0 = fem1d.inner_product_1d(lambda x: x, lambda x: sp.sin(x), 0, 1)
        result1 = quad(lambda x: x * sp.sin(x), 0, 1)
        self.assertEqual(result0, result1, 'Test inner_product fail')
