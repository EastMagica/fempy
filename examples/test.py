#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/20 16:32
# @file    : test.py
# @project : fempy
# software : PyCharm

import os
from fempy.fem1d import *
from fempy.fem2d import *

print(os.getcwd())

# Test 1: class TempElement
# print('>>> Test 1: class TempElement')
# temp = TempElement(3)
# print('get_area:', temp.get_area())
# print('get_gauss_point:', temp.get_gauss_point())
# print('get_gauss_weight', temp.get_gauss_weight())


# Test 2: def area
# print('>>> Test 2: def area')
# print(area([0, 0], [1, 0], [0, 1.0]))
# print(area([0, 0], [1, 0], [0.5, 1]))
# print(area(np.array([[0, 1], [0, 2], [0.5, 1.5]]).T, [0, 0], [1, 0]))
# print(area(np.array([0, 1]).T, [0, 0], [1, 0]))
# print(area(np.array([0, 2]).T, [0, 0], [1, 0]))
# print(area(np.array([0.5, 1.5]).T, [0, 0], [1, 0]))


# Test 3: class Mesh2D
print('>>> Test 3: class Mesh2D')
mesh2d = Mesh2D(0, 1, 0, 1, 4, 4)
print('pn, en, un:', mesh2d.p_n, mesh2d.e_n, mesh2d.u_n)
print('p_mat:', mesh2d.points)
print('p_bnd:', mesh2d.p_bnd)
print('e_mat:', mesh2d.edges)
print('e_bnd:', mesh2d.e_bnd)
print('u_mat:', mesh2d.triangles)
print('u_edg:', mesh2d.u_edg)
