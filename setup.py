#!/usr/bin/python
# -*- coding:utf-8 -*-
# @author  : East
# @time    : 2019/7/12 16:54
# @file    : setup.py
# @project : fem
# software : PyCharm

from setuptools import setup, find_packages


setup(
    name='fem',
    version='0.0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['numpy', 'scipy', 'matplotlib']
)
