#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/21 16:24
# @file    : test_1d_02.py
# @project : fempy
# software : PyCharm

import numpy as np
import matplotlib.pyplot as plt
from fempy.fem1d import FEM1D, Mesh1D, Adaptive1D

mesh = Mesh1D({'domain': np.array([-5, 5]), 'opt': ['step', 4]})
eq = {'bnd': np.array([np.e ** (-5**2), np.e ** (-5**2)]),
      'f': lambda x: (2 - 4 * x**2) * np.e ** (-x**2),
      'u_true': lambda x: np.e ** (-x**2)}

# mesh = Mesh1D({'domain': np.array([-1.5, 1.5]), 'opt': ['step', 4]})
# eq = {'bnd': np.array([np.e ** (-16 * 1.5**2), np.e ** (-16 * 1.5**2)]),
#       'f': lambda x: (2 * 16 - 4 * 16**2 * x**2) * np.e ** (-16 * x**2),
#       'u_true': lambda x: np.e ** (-16 * x**2)}

# mesh = Mesh1D({'domain': np.array([-np.sqrt(2), np.sqrt(2)]), 'opt': ['step', 4]})
# eq = {'bnd': np.array([0, 0]),
#       'f': lambda x: 12 * x**2 - 4,
#       'u_true': lambda x: 2 * x**2 - x**4}

ada = Adaptive1D(eq, mesh)
ada.fem_iter(5)
