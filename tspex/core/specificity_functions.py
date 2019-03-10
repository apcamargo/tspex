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

"""

"""

import numpy as np

from tspex.core.auxiliary_functions import dpm, entropy, js_distance, roku


def counts(vector, **kwargs):
    threshold = kwargs.pop('threshold', 0)
    n = len(vector)
    if n <= 1:
        return 0.0
    else:
        cts = np.sum(vector >= threshold)
        if cts == 0:
            return 0.0
        else:
            cts_transformed = (1-(cts/n)) * (n/(n-1))
            return cts_transformed

def tsi(vector, **kwargs):
    if not np.any(vector):
        return 0.0
    else:
        tissue_specificity_index = max(vector) / np.sum(vector)
        return tissue_specificity_index

def tau(vector, **kwargs):
    if not np.any(vector):
        return 0.0
    else:
        n = len(vector)
        vector_r = vector / np.max(vector)
        tau_index = np.sum(1-vector_r) / (n-1)
        return tau_index

def gini(vector, **kwargs):
    transform = kwargs.pop('transform', True)
    if not np.any(vector):
        return 0.0
    else:
        vector = np.sort(vector)
        n = len(vector)
        index = np.arange(1, n+1)
        gini_coefficient = (np.sum((2*index-n-1)*vector)) / (n*np.sum(vector))
        if transform:
            transformed_gini_coefficient = gini_coefficient * (n/(n-1))
            return transformed_gini_coefficient
        else:
            return gini_coefficient

def simpson(vector, **kwargs):
    transform = kwargs.pop('transform', True)
    if not np.any(vector):
        return 0.0
    else:
        p = vector / np.sum(vector)
        simpson_index = np.sum(p**2)
        if transform:
            min_simpson = 1 / len(vector)
            transformed_simpson_index = (simpson_index-min_simpson) / (1-min_simpson)
            return transformed_simpson_index
        else:
            return simpson_index

def shannon_specificity(vector, **kwargs):
    transform = kwargs.pop('transform', True)
    if not np.any(vector):
        return 0.0
    else:
        n = len(vector)
        ss = np.log2(n) - entropy(vector)
        if transform:
            ss_transformed = ss / np.log2(n)
            return ss_transformed
        else:
            return ss

def roku_specificity(vector, **kwargs):
    transform = kwargs.pop('transform', True)
    if not np.any(vector):
        return 0.0
    else:
        n = len(vector)
        rs = np.log2(n) - roku(vector)
        if transform:
            rs_transformed = rs / np.log2(n)
            return rs_transformed
        else:
            return rs

def zscore(vector, **kwargs):
    transform = kwargs.pop('transform', True)
    n = len(vector)
    std = np.std(vector, ddof=1)
    if std == 0:
        return np.zeros(n)
    else:
        zs = (vector-np.mean(vector)) / std
        if transform:
            max_zs = (n-1) / np.sqrt(n)
            zs_transformed = (zs+max_zs) / (2*max_zs)
            return zs_transformed
        else:
            return zs

def spm(vector, **kwargs):
    n = len(vector)
    spm_vector = []
    for i in range(n):
        if vector[i] == 0:
            spm_vector.append(0)
        else:
            vector_i = np.zeros(n)
            vector_i[i] = vector[i]
            spm_i = (vector[i]**2) / (np.linalg.norm(vector)*vector[i])
            spm_vector.append(spm_i)
    return np.array(spm_vector)

def spm_dpm(vector, **kwargs):
    spm_vector = spm(vector)
    spm_dispersion = dpm(spm_vector)
    return spm_dispersion

def js_specificity(vector, **kwargs):
    n = len(vector)
    js_vector = []
    for i in range(n):
        if vector[i] == 0:
            js_vector.append(0)
        else:
            vector_i = np.zeros(n)
            vector_i[i] = vector[i]
            js_i = 1 - js_distance(vector, vector_i)
            js_vector.append(js_i)
    return np.array(js_vector)

def js_specificity_dpm(vector, **kwargs):
    js_vector = js_specificity(vector)
    js_specificity_dispersion = dpm(js_vector)
    return js_specificity_dispersion
