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
TissueSpecificity class of the tspex library.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tspex.core.specificity_functions import (counts, gini, js_specificity,
                                              js_specificity_dpm,
                                              roku_specificity,
                                              shannon_specificity, simpson,
                                              spm, spm_dpm, tau, tsi, zscore)


class TissueSpecificity:
    """
    Create an object of the TissueSpecificity class.

    Parameters
    ----------
    expression_data : pandas.core.frame.DataFrame
        Pandas DataFrame containing the expression matrix, with rows
        corresponding to genes and columns to tissues/conditions.
    method : str
        A string representing which tissue-expression metric should be
        calculated. One of: 'counts', 'tsi', 'tau', 'gini', 'simpson',
        'shannon_specificity', 'roku_specificity', 'zscore', 'spm', 'spm_dpm',
        'js_specificity', 'js_specificity_dpm'.
    log : bool, default False
        Log-transform the expression matrix before computing tissue-expression
        by taking the base-2 logarithm of one plus the expression values. By
        default, no transformation is performed.
    transform : bool, default True
        Transform the computed value so it lies in the [0,1] range. By
        default, the value is transformed. The following metrics are affected
        by changes this parameter: 'gini', 'simpson', 'shannon_specificity',
        'roku_specificity', 'zscore'.
    threshold : int or float, default 0
        Value above which the gene is considered to be expressed. By default,
        any positive expression value is considered. Only the 'counts' metric
        is affected by changes in this parameter.

    Attributes
    ----------
    expression_data : pandas.DataFrame
        Expression matrix used to compute the tissue-specificity values. If the
        log parameter was set to True, the values will be log-transformed.
    tissue_specificity : pandas.Series or pandas.DataFrame
        Tissue-specificity values computed from the input expression matrix.
    """

    def __init__(self, expression_data, method, log=False, **kwargs):
        self._function_dictionary = {
            'counts': counts,
            'tsi': tsi,
            'tau': tau,
            'gini': gini,
            'simpson': simpson,
            'shannon_specificity': shannon_specificity,
            'roku_specificity': roku_specificity,
            'zscore': zscore,
            'spm': spm,
            'spm_dpm': spm_dpm,
            'js_specificity': js_specificity,
            'js_specificity_dpm': js_specificity_dpm
        }
        if log:
            self.expression_data = expression_data.astype('float')
            self.expression_data = self.expression_data.apply(lambda x: np.log(x+1))
        else:
            self.expression_data = expression_data.astype('float')
        self._method = str(method)
        self._transform = kwargs.pop('transform', True)
        self._threshold = kwargs.pop('threshold', 0)
        self.tissue_specificity = self._compute_tissue_specificity()

    def _compute_tissue_specificity(self):
        func = self._function_dictionary[self._method]
        if self._method in ['zscore', 'spm', 'js_specificity']:
            tissue_specificity = self.expression_data.apply(func, axis=1, result_type='broadcast',
                                                            transform=self._transform)
        else:
            tissue_specificity = self.expression_data.apply(func, axis=1, result_type='reduce',
                                                            transform=self._transform,
                                                            threshold=self._threshold)
        tissue_specificity = tissue_specificity.round(4)
        return tissue_specificity


    def histogram(self, bins=50, size=(7,4), dpi=100):
        """
        Plot a histogram of the tissue-specificity values. If the chosen metric
        is one of 'zscore', 'spm' or 'js_specificity', the maximum row value is used
        as a representative of the gene tissue-specificity.

        Parameters
        ----------
        bins : int, default 50
            Number of bins in the histogram.
        size : tuple, default (7,4)
            Size of the figure.
        dpi : int, default 100
            The resolution in dots per inch.
        """

        with plt.style.context('seaborn-whitegrid'):
            if self._method in ['zscore', 'spm', 'js_specificity']:
                data = self.tissue_specificity.max(axis=1).values
            else:
                data = self.tissue_specificity.values
            fig, ax = plt.subplots(figsize=size, dpi=dpi)
            ax.hist(data, bins=bins, alpha=0.85, color='#262626')
            ax.set_ylabel('Number of genes')
            ax.set_xlabel(self._method)
            ax.set_title('Histogram of {} values'.format(self._method), loc='left')


    def to_file(self, filename):
        """
        Write the tissue-specificity values into a tab-separated values (tsv)
        file.

        Parameters
        ----------
        filename : str
            A string containing a path to a filename.
        """

        self.tissue_specificity.to_csv(filename, sep='\t')

    def _repr_html_(self):
        if isinstance(self.tissue_specificity, pd.core.frame.DataFrame):
            return self.tissue_specificity._repr_html_()

    def __repr__(self):
        return self.tissue_specificity.__repr__()
