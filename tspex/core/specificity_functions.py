import numpy as np

from tspex.core.auxiliary_functions import dpm, entropy, js_distance, roku


def counts(vector, **kwargs):
    threshold = kwargs.pop('threshold', 0)
    n = len(vector)
    if n <= 1:
        return 0.0
    else:
        counts = np.sum(vector>=threshold)
        if counts == 0:
            return 0.0
        else:
            counts_transformed = (1-(counts/n)) * (n/(n-1))
            return counts_transformed

def tsi(vector, **kwargs):
    if not np.any(vector):
        return 0.0
    else:
        tsi = max(vector) / np.sum(vector)
        return tsi

def tau(vector, **kwargs):
    if not np.any(vector):
        return 0.0
    else:
        n = len(vector)
        vector_r = vector / np.max(vector)
        tau = np.sum(1-vector_r) / (n-1)
        return tau

def gini(vector, **kwargs):
    transform = kwargs.pop('transform', True)
    if not np.any(vector):
        return 0.0
    else:
        vector = np.sort(vector)
        n = len(vector)
        index = np.arange(1, n+1)
        gini = (np.sum((2*index-n-1)*vector)) / (n*np.sum(vector))
        if transform:
            transformed_gini = gini * (n/(n-1))
            return transformed_gini
        else:
            return gini

def simpson(vector, **kwargs):
    transform = kwargs.pop('transform', True)
    if not np.any(vector):
        return 0.0
    else:
        p = vector / np.sum(vector)
        simpson = np.sum(p**2)
        if transform:
            min_simpson = 1 / len(vector)
            transformed_simpson = (simpson-min_simpson) / (1-min_simpson)
            return transformed_simpson
        else:
            return simpson

def shannon_specificity(vector, **kwargs):
    transform = kwargs.pop('transform', True)
    if not np.any(vector):
        return 0.0
    else:
        n = len(vector)
        shannon_specificity = np.log2(n) - entropy(vector)
        if transform:
            shannon_specificity_transformed = shannon_specificity / np.log2(n)
            return shannon_specificity_transformed
        else:
            return shannon_specificity

def roku_specificity(vector, **kwargs):
    transform = kwargs.pop('transform', True)
    if not np.any(vector):
        return 0.0
    else:
        n = len(vector)
        roku_specificity = np.log2(n) - roku(vector)
        if transform:
            roku_specificity_transformed = roku_specificity / np.log2(n)
            return roku_specificity_transformed
        else:
            return roku_specificity

def zscore(vector, **kwargs):
    transform = kwargs.pop('transform', True)
    n = len(vector)
    std = np.std(vector, ddof=1)
    if std == 0:
        return np.zeros(n)
    else:
        zscore = (vector-np.mean(vector)) / std
        if transform:
            max_zscore = (n-1) / np.sqrt(n)
            zscore_transformed = (zscore+max_zscore) / (2*max_zscore)
            return zscore_transformed
        else:
            return zscore

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
    spm_dpm = dpm(spm_vector)
    return spm_dpm

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
    js_specificity_dpm = dpm(js_vector)
    return js_specificity_dpm
