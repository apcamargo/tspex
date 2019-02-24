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
        if log == True:
            self.expression_data = expression_data.astype('float').apply(lambda x: np.log(x+1))
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
        return tissue_specificity

    def to_file(self, filename):
        self.tissue_specificity.to_csv(filename, sep='\t')

    def _repr_html_(self):
        if isinstance(self.tissue_specificity, pd.core.frame.DataFrame):
            return self.tissue_specificity._repr_html_()

    def __repr__(self):
        return self.tissue_specificity.__repr__()
