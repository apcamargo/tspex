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

import numpy as np
import pandas as pd
import pytest

from tspex import TissueSpecificity

test_data = pd.read_csv(
    'tests/test_data.tsv', index_col=0, header=0, sep=None, thousands=',', engine='python'
)
duplicated_test_data = pd.concat([test_data, test_data])
negative_test_data = -1 * test_data
non_numerical_test_data = test_data
non_numerical_test_data.loc[:, 'non_numerical'] = [
    'a',
    'b',
    'c',
    'd',
    'e',
    'f',
    'g',
    'h',
    'i',
    'j',
]


def test_specificity_class_negative_values():
    pytest.raises(ValueError, TissueSpecificity, negative_test_data, method='gini')

def test_specificity_class_duplicated_genes():
    pytest.raises(ValueError, TissueSpecificity, duplicated_test_data, method='gini')

def test_specificity_class_plot_histogram():
    # General scoring metric
    assert (
        TissueSpecificity(non_numerical_test_data, method='gini').plot_histogram() is None
    )
    # Individualized scoring metric
    assert (
        TissueSpecificity(non_numerical_test_data, method='spm').plot_histogram() is None
    )

def test_specificity_class_plot_heatmap():
    # General scoring metric
    assert (
        TissueSpecificity(non_numerical_test_data, method='gini').plot_heatmap(
            threshold=0.9
        )
        is None
    )
    # Individualized scoring metric
    assert (
        TissueSpecificity(non_numerical_test_data, method='spm').plot_heatmap(
            threshold=0.9
        )
        is None
    )
    # sort_genes=True
    assert (
        TissueSpecificity(non_numerical_test_data, method='gini').plot_heatmap(
            threshold=0.9, sort_genes=True
        )
        is None
    )
    # use_zscore=True
    assert (
        TissueSpecificity(non_numerical_test_data, method='gini').plot_heatmap(
            threshold=0.9, use_zscore=True
        )
        is None
    )
    # gene_names=False
    assert (
        TissueSpecificity(non_numerical_test_data, method='gini').plot_heatmap(
            threshold=0.9, gene_names=False
        )
        is None
    )
    # tissue_names=False
    assert (
        TissueSpecificity(non_numerical_test_data, method='gini').plot_heatmap(
            threshold=0.9, tissue_names=False
        )
        is None
    )
    # All genes below threshold
    assert (
        TissueSpecificity(non_numerical_test_data, method='gini').plot_heatmap(
            threshold=1
        )
        is None
    )

def test_specificity_class_non_numerical_column():
    assert np.all(
        TissueSpecificity(
            non_numerical_test_data, method='gini', log=False, transform=False
        ).tissue_specificity.values
        == np.array(
            [
                0.5299,
                0.1861,
                0.0844,
                0.3474,
                0.2539,
                0.2295,
                0.2259,
                0.8214,
                0.2342,
                0.2582,
            ]
        )
    )

def test_specificity_class_general_scoring():
    assert np.all(
        TissueSpecificity(
            test_data, method='gini', log=False, transform=False
        ).tissue_specificity.values
        == np.array(
            [
                0.5299,
                0.1861,
                0.0844,
                0.3474,
                0.2539,
                0.2295,
                0.2259,
                0.8214,
                0.2342,
                0.2582,
            ]
        )
    )
    assert np.all(
        TissueSpecificity(
            test_data, method='gini', log=True, transform=True
        ).tissue_specificity.values
        == np.array(
            [
                0.2575,
                0.1067,
                0.0444,
                0.2259,
                0.185,
                0.1458,
                0.1306,
                0.7753,
                0.1307,
                0.1573,
            ]
        )
    )

def test_specificity_class_individualized_scoring():
    # transform=False and log=False
    assert np.all(
        TissueSpecificity(
            test_data, method='zscore', log=False, transform=False
        ).tissue_specificity.values
        == np.array(
            [
                [-0.5635, -0.2404, -0.5129, -0.3526, 2.0271, -0.3577],
                [-1.5805, 0.0456, 0.8062, 0.2104, 1.1733, -0.655],
                [0.6634, 0.246, 1.1194, -0.009, -1.7866, -0.2331],
                [-0.5055, 1.3064, 1.1814, -0.2528, -1.0678, -0.6617],
                [0.1812, -0.8203, 0.3897, -1.548, 1.1542, 0.6432],
                [-1.3225, 0.1824, 0.7338, -1.0406, 0.199, 1.2478],
                [-0.8264, -1.2235, 0.018, 1.2556, -0.2982, 1.0744],
                [-0.4111, -0.4086, -0.4122, -0.4049, 2.0412, -0.4045],
                [1.0544, -1.0087, 1.3507, -0.3583, -0.9545, -0.0837],
                [-1.1506, 1.2017, -1.0321, -0.3351, 0.3514, 0.9648],
            ]
        )
    )
    # transform=True and log=True
    assert np.all(
        TissueSpecificity(
            test_data, method='zscore', log=True, transform=True
        ).tissue_specificity.values
        == np.array(
            [
                [0.2792, 0.5161, 0.3322, 0.4542, 0.9673, 0.4511],
                [0.0708, 0.5439, 0.6849, 0.5771, 0.7434, 0.3799],
                [0.6584, 0.5705, 0.7488, 0.514, 0.0457, 0.4626],
                [0.4182, 0.7904, 0.7733, 0.4957, 0.1608, 0.3616],
                [0.5794, 0.3545, 0.615, 0.0706, 0.7257, 0.6548],
                [0.1439, 0.5787, 0.68, 0.255, 0.582, 0.7604],
                [0.3072, 0.1549, 0.5432, 0.7795, 0.465, 0.7501],
                [0.3876, 0.3995, 0.3818, 0.4151, 0.9992, 0.4169],
                [0.7555, 0.226, 0.8044, 0.4464, 0.2477, 0.52],
                [0.1893, 0.7626, 0.2393, 0.4639, 0.6186, 0.7263],
            ]
        )
    )