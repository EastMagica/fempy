#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/21 22:02
# @file    : plot1d.py
# @project : fem
# software : PyCharm

import numpy as np
import matplotlib.pyplot as plt


def plot_1d(u_lst, mesh, u_true, title):
    # plt.style.use('grayscale')
    fig = plt.figure(figsize=(12, 5))
    ax0 = plt.subplot(1, 2, 1)
    ax0.set_title('u(x)')
    # ax0.plot(np.linspace(mesh.p_mat[mesh.p_end][0], mesh.p_mat[mesh.p_end][1], 500),
    #          u_true(np.linspace(mesh.p_mat[mesh.p_end][0], mesh.p_mat[mesh.p_end][1], 500)), color='b', linestyle='-')
    ax0.plot(np.sort(mesh.p_mat), u_lst[np.argsort(mesh.p_mat)], color='r', linestyle='--')
    ax0.scatter(mesh.p_mat, u_lst, color='black')
    # ax0.set_xticks(np.sort(mesh.p_mat))
    ax0.set_xlim(mesh.p_mat[mesh.p_end][0], mesh.p_mat[mesh.p_end][1])
    # ax0.set_ylim(-0.05, 1.2)
    ax1 = plt.subplot(1, 2, 2)
    ax1.set_title('split points')
    # ax1.plot(np.sort(mesh.p_mat), (u_true(mesh.p_mat) - u_lst)[np.argsort(mesh.p_mat)])
    # ax1.scatter(mesh.p_mat, u_true(mesh.p_mat) - u_lst)
    ax1.scatter(mesh.p_mat, [0]*mesh.p_n)
    ax1.set_xlim(mesh.p_mat[mesh.p_end][0], mesh.p_mat[mesh.p_end][1])
    fig.suptitle(title)
    plt.show()
