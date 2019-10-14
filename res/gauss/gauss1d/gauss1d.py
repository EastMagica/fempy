#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : east
# @time    : 2019/10/7 0:23
# @file    : gauss1d.py
# @project : fem
# software : PyCharm

import numpy as np

"""
Data from Wikipedia [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature).
"""

n = [1, 2, 3, 4, 5]
name = ['gauss1', 'gauss2', 'gauss3', 'gauss4', 'gauss5']

gauss1 = np.array([[0], [2]])

gauss2 = np.array([[-1/np.sqrt(3), 1/np.sqrt(3)], [1, 1]], dtype=np.float64)

gauss3 = np.array([[-np.sqrt(3/5), 0, np.sqrt(3/5)],
                   [5 / 9, 8 / 9, 5 / 9]], dtype=np.float64)

gauss4 = np.array([[-np.sqrt(3/7+2/7*np.sqrt(6/5)), -np.sqrt(3/7-2/7*np.sqrt(6/5)),
                    np.sqrt(3/7-2/7*np.sqrt(6/5)), np.sqrt(3/7-2/7*np.sqrt(6/5))],
                   [(18-np.sqrt(30))/36, (18+np.sqrt(30))/36,
                    (18+np.sqrt(30))/36, (18-np.sqrt(30))/36]], dtype=np.float64)

gauss5 = np.array([[-np.sqrt(5+2*np.sqrt(10/7))/3, -np.sqrt(5-2*np.sqrt(10/7))/3,
                    0, np.sqrt(5-2*np.sqrt(10/7))/3, np.sqrt(5+2*np.sqrt(10/7))/3],
                   [(322-13*np.sqrt(70))/900, (322+13*np.sqrt(70))/900, 128/225,
                    (322+13*np.sqrt(70))/900, (322-13*np.sqrt(70))/900]], dtype=np.float64)

for k, v in enumerate([gauss1, gauss2, gauss3, gauss4, gauss5]):
    np.save(name[k], v)
