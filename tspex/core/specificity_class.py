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

"""

import numpy as np
import pandas as pd

from tspex.core.specificity_functions import (counts, gini, js_specificity,
                                              js_specificity_dpm,
                                              roku_specificity,
                                              shannon_specificity, simpson,
                                              spm, spm_dpm, tau, tsi, zscore)

func_dictionary = {
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

class TissueSpecificity:
    def __init__(self, expression_data, method, log=False, **kwargs):
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
        func = func_dictionary[self._method]
        if self._method in ['zscore', 'spm', 'js_specificity']:
            tissue_specificity = self.expression_data.apply(func, axis=1, result_type='broadcast',
                                                            transform=self._transform)
        else:
            tissue_specificity = self.expression_data.apply(func, axis=1, result_type='reduce',
                                                            transform=self._transform,
                                                            threshold=self._threshold)
        tissue_specificity = tissue_specificity.round(4)
        return tissue_specificity

    def to_file(self, filename):
        self.tissue_specificity.to_csv(filename, sep='\t')

    def _repr_html_(self):
        if isinstance(self.tissue_specificity, pd.core.frame.DataFrame):
            return self.tissue_specificity._repr_html_()

    def __repr__(self):
        return self.tissue_specificity.__repr__()
