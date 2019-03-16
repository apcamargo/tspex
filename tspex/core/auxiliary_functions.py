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
Auxiliary functions.
"""

import numpy as np


def dpm(vector):
    """
    Compute the Dispersion Measure (DPM) [1] of a vector.

    Parameters
    ----------
    vector : numpy.array
        Input vector.

    Returns
    -------
    float
        Dispersion Measure (DPM) of the vector.

    References
    ----------
    .. [1] Pan, Jian-Bo, et al. "PaGeFinder: quantitative identification of
           spatiotemporal pattern genes." Bioinformatics 28.11 (2012)
    """

    n = len(vector)
    dispersion_measure = np.std(vector, ddof=1) * np.sqrt(n)
    return dispersion_measure


def tukey_biweight(vector, c=5, epsilon=1e-4):
    """
    Compute the one-step Tukey's biweight of a vector. Taken from the affy R
    package.

    Parameters
    ----------
    vector : numpy.array
        Input vector.
    c : int or float, default 5
        Tuning constant.
    epsilon : int or float, default 1e-4
        Fuzzy value to avoid division by zero.

    Returns
    -------
    float
        One-step Tukey's biweight of the vector.
    """

    m = np.median(vector)
    s = np.median(np.abs(vector-m))
    u = (vector-m) / ((c*s)+epsilon)
    i = (np.abs(u) > 1)
    w = (1-u**2) ** 2
    w[i] = 0
    tbi = np.sum(w*vector) / np.sum(w)
    return tbi


def entropy(vector):
    """
    Compute the Shannon entropy [1] of a vector.

    Parameters
    ----------
    vector : numpy.array
        Input vector.

    Returns
    -------
    float
        Shannon entropy of the vector.

    References
    ----------
    .. [1] Shannon, Claude Elwood. "A mathematical theory of communication."
           Bell system technical journal 27.3 (1948)
    """

    n = len(vector)
    if not np.any(vector):
        return np.log2(n)
    else:
        p = vector / np.sum(vector)
        p = p[np.nonzero(p)[0]]
        h = -1 * np.dot(p, np.log2(p))
        return h


def roku(vector):
    """
    Compute the Shannon entropy of a vector processed using the ROKU method
    [2]. Data processing is done by subtracting the one-step Tukey's biweight
    and taking the absolute value.

    Parameters
    ----------
    vector : numpy.array
        Input vector.

    Returns
    -------
    float
        Shannon entropy of the ROKU-processed vector.

    References
    ----------
    .. [1] Shannon, Claude Elwood. "A mathematical theory of communication."
           Bell system technical journal 27.3 (1948)
    .. [2] Kadota, Koji, et al. "ROKU: a novel method for identification of
           tissue-specific genes." BMC bioinformatics 7.1 (2006)
    """

    tbi = tukey_biweight(vector)
    vector_p = np.abs(vector-tbi)
    h = entropy(vector_p)
    return h


def js_distance(p, q):
    """
    Compute the Jensen-Shannon distance [1] between two vectors.

    Parameters
    ----------
    p : numpy.array
        First vector.
    q : numpy.array
        Second vector.

    Returns
    -------
    float
        Jensen-Shannon distance between the two vectors.

    References
    ----------
    .. [1] Endres, Dominik Maria, and Johannes E. Schindelin. "A new metric for
           probability distributions." IEEE Transactions on Information theory
           (2003).
    """

    p = p / np.sum(p)
    q = q / np.sum(q)
    left = entropy((p+q)/2)
    right = (entropy(p)+entropy(q)) / 2
    js = left - right
    jsd = np.sqrt(js)
    return jsd
