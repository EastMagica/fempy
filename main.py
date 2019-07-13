#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/12 17:45
# @file    : main.py
# @project : fem
# software : PyCharm

from fem.fem_1d import *


print(inner_product(lambda x: x, lambda x: x, 0, 1))

print(quad(lambda x: x**2, 0, 1))
