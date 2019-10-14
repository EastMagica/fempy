#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/16 21:04
# @file    : mesh2d.py
# @project : fempy
# software : PyCharm

# imports
import numpy as np
from scipy.spatial import Delaunay
from matplotlib.tri import Triangulation, UniformTriRefiner


# functions

def rect_tri(x_opt, y_opt=None):
    """

    Parameters
    ----------
    x_opt :
    y_opt :

    Returns
    -------

    """
    y_opt = x_opt if y_opt is None else y_opt
    print(x_opt, y_opt)
    x = x_opt if isinstance(x_opt, np.ndarray) else np.linspace(*x_opt)
    y = y_opt if isinstance(y_opt, np.ndarray) else np.linspace(*y_opt)
    print(x, y)
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError("x and y must be equal-length 1-D arrays")
    x, y = np.meshgrid(x, y)
    points = np.transpose([x.flatten(), y.flatten()])
    striang = Delaunay(points)
    # c = np.unique(self.stri.convex_hull)
    x0 = points[:, 0] == 0
    xn = points[:, 0] == 1
    y0 = points[:, 1] == 0
    yn = points[:, 1] == 1
    cond = x0 | xn | y0 | yn
    return striang, cond


# classes

class Mesh2D(object):
    """
    Construct Mesh
    """

    def __init__(self, method='rect_tri', opt=None):
        """

        Parameters
        ----------
        method
        opt
        """
        if method == 'rect_tri':
            self.__stri, self.__p_bnd = rect_tri(*opt)
        else:
            raise ValueError('Other methods are not supported.')

    @property
    def simplices(self):
        """

        Returns
        -------
        simplices : np.ndarray, (Ntri, 3)
        """
        return self.__stri.simplices

    @property
    def points(self):
        """

        Returns
        -------
        points : np.ndarray, (N, 2)
        """
        return self.__stri.points

    @property
    def npoints(self):
        return self.__stri.npoints

    @property
    def nsimplices(self):
        return self.__stri.nsimplex

    @property
    def ndim(self):
        return self.__stri.ndim

    @property
    def boundary(self):
        return self.__p_bnd

    @property
    def mtri(self):
        triangles = self.simplices
        x, y = np.transpose(self.points)
        return Triangulation(x, y, triangles)

    def refine(self, z, subdiv=3, triinterpolator=None):
        r"""
        对三角形网格重新采样.

        Parameters
        ----------
        z : array_like, (N, )
        subdiv : int, optional
            采样深度, 将三角形单元分为 4**subdiv.
        triinterpolator : TriInterpolator, optional
            插值方法, 默认使用 CubicTriInterpolator.

        Returns
        -------

        """
        refiner = UniformTriRefiner(self.mtri)
        mtri_ref, z_ref = refiner.refine_field(z, subdiv=subdiv, triinterpolator=triinterpolator)
        return mtri_ref, z_ref
