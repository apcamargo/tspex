# -*- coding: utf-8 -*-
#
#   This file is part of the tspex package, available at:
#   https://github.com/apcamargo/tspex
#
#   Tspex is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
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
Functions to compute tissue-specificity metrics from expression arrays.
"""

import numpy as np

from tspex.core.auxiliary_functions import dpm, entropy, js_distance, roku


def counts(vector, **kwargs):
    """
    Quantify tissue-specificity as the proportion of tissues above an
    expression threshold.

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.
    threshold : int or float, default 0
        Value above which the gene is considered to be expressed. By default,
        any positive expression value is considered.

    Returns
    -------
    float
        Single summary of the tissue-specificity. Ranges from 0 (ubiquitous
        expression) to 1 (specific expression).

    References
    ----------
    .. [1] Duret, Laurent, and Dominique Mouchiroud. "Determinants of
           substitution rates in mammalian genes: expression pattern affects
           selection intensity but not mutation rate." Molecular biology and
           evolution 17.1 (2000)
    .. [2] Ponger, Loïc, Laurent Duret, and Dominique Mouchiroud. "Determinants
           of CpG islands: expression in early embryo and isochore structure."
           Genome research 11.11 (2001)
    .. [3] Lercher, Martin J., Araxi O. Urrutia, and Laurence D. Hurst.
           "Clustering of housekeeping genes provides a unified model of gene
           order in the human genome." Nature genetics 31.2 (2002)
    .. [4] Vinogradov, Alexander E. "Isochores and tissue‐specificity." Nucleic
           acids research 31.17 (2003)
    .. [5] Subramanian, Sankar, and Sudhir Kumar. "Gene expression intensity
           shapes evolutionary rates of the proteins encoded by the vertebrate
           genome." Genetics 168.1 (2004)
    .. [6] Park, Seung Gu, and Sun Shim Choi. "Expression breadth and expression
           abundance behave differently in correlations with evolutionary
           rates." BMC evolutionary biology 10.1 (2010)
    """

    threshold = kwargs.pop('threshold', 0)
    n = len(vector)
    if n <= 1:
        return 0.0
    else:
        cts = np.sum(vector > threshold)
        if cts == 0:
            return 0.0
        else:
            cts_transformed = (1 - (cts / n)) * (n / (n - 1))
            return cts_transformed


def tau(vector, **kwargs):
    """
    Quantify tissue-specificity as the Tau index [1].

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.

    Returns
    -------
    float
        Single summary of the tissue-specificity. Ranges from 0 (ubiquitous
        expression) to 1 (specific expression).

    References
    ----------
    .. [1] Yanai, Itai, et al. "Genome-wide midrange transcription profiles
           reveal expression level relationships in human tissue specification."
           Bioinformatics 21.5 (2004)
    """

    if not np.any(vector):
        return 0.0
    else:
        n = len(vector)
        vector_r = vector / np.max(vector)
        tau_index = np.sum(1 - vector_r) / (n - 1)
        return tau_index


def gini(vector, **kwargs):
    """
    Quantify tissue-specificity as the Gini coefficient [1].

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.
    transform : bool, default True
        Transform the computed value so it lies in the [0,1] range. By
        default, the value is transformed.

    Returns
    -------
    float
        Single summary of the tissue-specificity. If transformed, it ranges
        from 0 (ubiquitous expression) to 1 (specific expression).

    References
    ----------
    .. [1] Ceriani, Lidia, and Paolo Verme. "The origins of the Gini index:
           extracts from Variabilità e Mutabilità (1912) by Corrado Gini." The
           Journal of Economic Inequality 10.3 (2012)
    """

    transform = kwargs.pop('transform', True)
    if not np.any(vector):
        return 0.0
    else:
        vector = np.sort(vector)
        n = len(vector)
        index = np.arange(1, n + 1)
        gini_coefficient = (np.sum((2 * index - n - 1) * vector)) / (n * np.sum(vector))
        if transform:
            transformed_gini_coefficient = gini_coefficient * (n / (n - 1))
            return transformed_gini_coefficient
        else:
            return gini_coefficient


def simpson(vector, **kwargs):
    """
    Quantify tissue-specificity as the Simpson index [1].

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.
    transform : bool, default True
        Transform the computed value so it lies in the [0,1] range. By
        default, the value is transformed.

    Returns
    -------
    float
        Single summary of the tissue-specificity. If transformed, it ranges
        from 0 (ubiquitous expression) to 1 (specific expression).

    References
    ----------
    .. [1] Simpson, Edward H. "Measurement of diversity." Nature 163.4148 (1949)
    """

    transform = kwargs.pop('transform', True)
    if not np.any(vector):
        return 0.0
    else:
        p = vector / np.sum(vector)
        simpson_index = np.sum(p ** 2)
        if transform:
            min_simpson = 1 / len(vector)
            transformed_simpson_index = (simpson_index - min_simpson) / (1 - min_simpson)
            return transformed_simpson_index
        else:
            return simpson_index


def shannon_specificity(vector, **kwargs):
    """
    Quantify tissue-specificity using Shannon entropy [1]. Specifically, this
    metric corresponds to the difference between the maximum and the observed
    vector entropy.

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.
    transform : bool, default True
        Transform the computed value so it lies in the [0,1] range. By
        default, the value is transformed.

    Returns
    -------
    float
        Single summary of the tissue-specificity. If transformed, it ranges
        from 0 (ubiquitous expression) to 1 (specific expression).

    References
    ----------
    .. [1] Shannon, Claude Elwood. "A mathematical theory of communication."
           Bell system technical journal 27.3 (1948)
    """

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
    """
    Quantify tissue-specificity using the ROKU method [1], in which the
    expression vector is processed and used to compute the Shannon entropy [2].
    Data processing is done by subtracting the one-step Tukey's biweight and
    taking the absolute value. This metric corresponds to the difference
    between the maximum and the observed processed vector entropy.

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.
    transform : bool, default True
        Transform the computed value so it lies in the [0,1] range. By
        default, the value is transformed.

    Returns
    -------
    float
        Single summary of the tissue-specificity. If transformed, it ranges
        from 0 (ubiquitous expression) to 1 (specific expression).

    References
    ----------
    .. [1] Kadota, Koji, et al. "ROKU: a novel method for identification of
           tissue-specific genes." BMC bioinformatics 7.1 (2006)
    .. [2] Shannon, Claude Elwood. "A mathematical theory of communication."
           Bell system technical journal 27.3 (1948)
    """

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


def tsi(vector, **kwargs):
    """
    Quantify tissue-specificity as the ratio between the expression vector
    and the sum of the expression values in all tissues.

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.

    Returns
    -------
    numpy.array
        Summary of the tissue-specificity level in each tissue. It ranges from
        0 (ubiquitous expression) to 1 (specific expression).

    References
    ----------
    .. [1] Julien, Philippe, et al. "Mechanisms and evolutionary patterns of
           mammalian and avian dosage compensation." PLoS biology 10.5 (2012)
    """

    if not np.any(vector):
        return 0.0
    else:
        tissue_specificity_index = vector / np.sum(vector)
        return tissue_specificity_index


def zscore(vector, **kwargs):
    """
    Quantify tissue-specificity as z-scores.

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.
    transform : bool, default True
        Transform the computed value so it lies in the [0,1] range. By
        default, the value is transformed.

    Returns
    -------
    numpy.array Summary of the tissue-specificity level in each tissue. If
    transformed, it ranges from 0 (ubiquitous expression) to 1 (specific
    expression).

    References
    ----------
    .. [1] Kryuchkova-Mostacci, Nadezda, and Marc Robinson-Rechavi. "A benchmark
           of gene expression tissue-specificity metrics." Briefings in
           bioinformatics 18.2 (2017)
    """

    transform = kwargs.pop('transform', True)
    n = len(vector)
    std = np.std(vector, ddof=1)
    if std == 0:
        return np.zeros(n)
    else:
        zs = (vector - np.mean(vector)) / std
        if transform:
            max_zs = (n - 1) / np.sqrt(n)
            zs_transformed = (zs + max_zs) / (2 * max_zs)
            return zs_transformed
        else:
            return zs


def spm(vector, **kwargs):
    """
    Quantify tissue-specificity as the Specificity Measure (SPM) [1,2].

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.

    Returns
    -------
    numpy.array
        Summary of the tissue-specificity level in each tissue. It ranges from
        0 (ubiquitous expression) to 1 (specific expression).

    References
    ----------
    .. [1] Xiao, Sheng-Jian, et al. "TiSGeD: a database for tissue-specific
           genes." Bioinformatics 26.9 (2010)
    .. [2] Pan, Jian-Bo, et al. "PaGeFinder: quantitative identification of
           spatiotemporal pattern genes." Bioinformatics 28.11 (2012)
    """

    n = len(vector)
    spm_vector = []
    for i in range(n):
        if vector[i] == 0:
            spm_vector.append(0)
        else:
            spm_i = (vector[i] ** 2) / (np.linalg.norm(vector) * vector[i])
            spm_vector.append(spm_i)
    return np.array(spm_vector)


def spm_dpm(vector, **kwargs):
    """
    Quantify tissue-specificity as the Dispersion Measure (DPM) [1] computed
    with SPM values.

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.

    Returns
    -------
    float
        Single summary of the tissue-specificity. It ranges from 0 (ubiquitous
        expression) to 1 (specific expression).

    References
    ----------
    .. [1] Pan, Jian-Bo, et al. "PaGeFinder: quantitative identification of
           spatiotemporal pattern genes." Bioinformatics 28.11 (2012)
    """

    spm_vector = spm(vector)
    spm_dispersion = dpm(spm_vector)
    return spm_dispersion


def js_specificity(vector, **kwargs):
    """
    Quantify tissue-specificity using the Jensen-Shannon distance, as
    described in [1].

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.

    Returns
    -------
    numpy.array
        Summary of the tissue-specificity level in each tissue. It ranges from
        0 (ubiquitous expression) to 1 (specific expression).

    References
    ----------
    .. [1] Cabili, Moran N., et al. "Integrative annotation of human large
           intergenic noncoding RNAs reveals global properties and specific
           subclasses." Genes & development 25.18 (2011)
    """

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
    """
    Quantify tissue-specificity as the Dispersion Measure (DPM) [1] computed
    with Jensen-Shannon distance-based specificity [2] values.

    Parameters
    ----------
    vector : numpy.array
        Gene expression vector. Each value corresponds to the gene expression
        in a given tissue.

    Returns
    -------
    float
        Single summary of the tissue-specificity. It ranges from 0 (ubiquitous
        expression) to 1 (specific expression).

    References
    ----------
    .. [1] Pan, Jian-Bo, et al. "PaGeFinder: quantitative identification of
           spatiotemporal pattern genes." Bioinformatics 28.11 (2012)
    .. [2] Cabili, Moran N., et al. "Integrative annotation of human large
           intergenic noncoding RNAs reveals global properties and specific
           subclasses." Genes & development 25.18 (2011)
    """

    js_vector = js_specificity(vector)
    js_specificity_dispersion = dpm(js_vector)
    return js_specificity_dispersion
