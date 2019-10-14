#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : east
# @time    : 2019/10/8 23:41
# @file    : test_basic.py
# @project : fem
# software : PyCharm

from fempy.basic import GaussianTemp

gauss1d = GaussianTemp('1', 3)
print(f'gauss1d.points:\n{gauss1d.points}')
print(f'gauss1d.points:\n{gauss1d.weight}')

print('----------------------------------')

gauss2d = GaussianTemp('2', 3)
print(f'gauss2d.points:\n{gauss2d.points}')
print(f'gauss2d.points:\n{gauss1d.weight}')
