# !/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   This file is part of the tspex package, available at:
#   https://github.com/apcamargo/tspex
#
#   Tspex is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   Contact: antoniop.camargo@gmail.com

"""The setup script."""

from setuptools import setup, find_packages

setup(
    name='tspex',
    version='0.1.0',
    packages=find_packages(include=['tspex']),
    license='GNU General Public License v3.0',
    description='A Python package for calculating tissue-specificity metrics for gene expression.',
    long_description=open('README.md').read(),
    install_requires=['numpy', 'pandas >= 0.23'],
    python_requires= '>=3',
    url='https://github.com/apcamargo/tspex',
    keywords=['gene expression', 'tissue-specificity', 'transcriptomics'],
    author='Antonio Pedro Camargo',
    author_email='antoniop.camargo@gmail.com',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Programming Language :: Python :: 3',
    ],
)
