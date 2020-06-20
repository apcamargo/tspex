import numpy as np
import pytest
from tspex.core.specificity_functions import (
    counts,
    gini,
    js_specificity,
    js_specificity_dpm,
    roku_specificity,
    shannon_specificity,
    simpson,
    spm,
    spm_dpm,
    tau,
    tsi,
    zscore,
)

standard_vector = np.array([0, 1.2, 0, 0.6, 7, 0.8])
all_zero_vector = np.array([0, 0, 0, 0, 0])
single_element_vector = np.array([1])

def test_counts():
    # No threshold
    assert round(counts(standard_vector), 4) == 0.4
    assert round(counts(all_zero_vector), 4) == 0.0
    assert round(counts(single_element_vector), 4) == 0.0
    # Threshold above some of the values
    assert round(counts(standard_vector, threshold=1), 4) == 0.8
    assert round(counts(all_zero_vector, threshold=1), 4) == 0.0
    assert round(counts(single_element_vector, threshold=1), 4) == 0.0

def test_gini():
    # With transform=False
    assert round(gini(standard_vector, transform=False), 4) == 0.6736
    assert round(gini(all_zero_vector, transform=False), 4) == 0.0
    assert round(gini(single_element_vector, transform=False), 4) == 0.0
    # With transform=True
    assert round(gini(standard_vector, transform=True), 4) == 0.8083
    assert round(gini(all_zero_vector, transform=True), 4) == 0.0
    assert round(gini(single_element_vector, transform=True), 4) == 0.0

def test_js_specificity():
    # With transform=False
    assert all(js_specificity(standard_vector, transform=False).round(4) == np.array([0, 0.1533, 0, 0.0898, 0.6117, 0.1123]))
    assert all(js_specificity(all_zero_vector, transform=False).round(4) == np.array([0, 0, 0, 0, 0]))
    assert js_specificity(single_element_vector, transform=False).round(4) == np.array([0])
    # With transform=True
    assert all(js_specificity(standard_vector, transform=True).round(4) == np.array([0, 0.1533, 0, 0.0898, 0.6117, 0.1123]))
    assert all(js_specificity(all_zero_vector, transform=True).round(4) == np.array([0, 0, 0, 0, 0]))
    assert js_specificity(single_element_vector, transform=True).round(4) == np.array([0])

def test_js_specificity_dpm():
    assert round(js_specificity_dpm(standard_vector), 4) == 0.5612
    assert round(js_specificity_dpm(all_zero_vector), 4) == 0.0
    assert round(js_specificity_dpm(single_element_vector), 4) == 0.0

def test_roku_specificity():
    # With transform=False
    assert round(roku_specificity(standard_vector, transform=False), 4) == 1.2847
    assert round(roku_specificity(all_zero_vector, transform=False), 4) == 0.0
    assert round(roku_specificity(single_element_vector, transform=False), 4) == 0.0
    # With transform=True
    assert round(roku_specificity(standard_vector, transform=True), 4) == 0.497
    assert round(roku_specificity(all_zero_vector, transform=True), 4) == 0.0
    assert round(roku_specificity(single_element_vector, transform=True), 4) == 0.0

def test_shannon_specificity():
    # With transform=False
    assert round(shannon_specificity(standard_vector, transform=False), 4) == 1.3289
    assert round(shannon_specificity(all_zero_vector, transform=False), 4) == 0.0
    assert round(shannon_specificity(single_element_vector, transform=False), 4) == 0.0
    # With transform=True
    assert round(shannon_specificity(standard_vector, transform=True), 4) == 0.5141
    assert round(shannon_specificity(all_zero_vector, transform=True), 4) == 0.0
    assert round(shannon_specificity(single_element_vector, transform=True), 4) == 0.0

def test_simpson():
    # With transform=False
    assert round(simpson(standard_vector, transform=False), 4) == 0.5582
    assert round(simpson(all_zero_vector, transform=False), 4) == 0.0
    assert round(simpson(single_element_vector, transform=False), 4) == 0.0
    # With transform=True
    assert round(simpson(standard_vector, transform=True), 4) == 0.4698
    assert round(simpson(all_zero_vector, transform=True), 4) == 0.0
    assert round(simpson(single_element_vector, transform=True), 4) == 0.0

def test_spm():
    # With transform=False
    assert all(spm(standard_vector, transform=False).round(4) == np.array([0, 0.1673, 0, 0.0837, 0.976 , 0.1115]))
    assert all(spm(all_zero_vector, transform=False).round(4) == np.array([0, 0, 0, 0, 0]))
    assert spm(single_element_vector, transform=False).round(4) == np.array([0])
    # With transform=True
    assert all(spm(standard_vector, transform=True).round(4) == np.array([0, 0.1673, 0, 0.0837, 0.976 , 0.1115]))
    assert all(spm(all_zero_vector, transform=True).round(4) == np.array([0, 0, 0, 0, 0]))
    assert spm(single_element_vector, transform=True).round(4) == np.array([0])

def test_spm_dpm():
    # With transform=False
    assert round(spm_dpm(standard_vector, transform=False), 4) == 0.9174
    assert round(spm_dpm(all_zero_vector, transform=False), 4) == 0.0
    assert round(spm_dpm(single_element_vector, transform=False), 4) == 0.0
    # With transform=True
    assert round(spm_dpm(standard_vector, transform=True), 4) == 0.9174
    assert round(spm_dpm(all_zero_vector, transform=True), 4) == 0.0
    assert round(spm_dpm(single_element_vector, transform=True), 4) == 0.0

def test_tau():
    assert round(tau(standard_vector), 4) == 0.9257
    assert round(tau(all_zero_vector), 4) == 0.0
    assert round(tau(single_element_vector), 4) == 0.0

def test_tsi():
    assert all(tsi(standard_vector).round(4) == np.array([0, 0.125 , 0, 0.0625, 0.7292, 0.0833]))
    assert all(tsi(all_zero_vector).round(4) == np.array([0, 0, 0, 0, 0]))
    assert tsi(single_element_vector).round(4) == np.array([0])

def test_zscore():
    # With transform=False
    assert all(zscore(standard_vector, transform=False).round(4) == np.array([-0.5956, -0.1489, -0.5956, -0.3723, 2.0102, -0.2978]))
    assert all(zscore(all_zero_vector, transform=False).round(4) == np.array([0, 0, 0, 0, 0]))
    assert zscore(single_element_vector, transform=False).round(4) == np.array([0])
    # With transform=True
    assert all(zscore(standard_vector, transform=True).round(4) == np.array([0.3541, 0.4635, 0.3541, 0.4088, 0.9924, 0.4271]))
    assert all(zscore(all_zero_vector, transform=True).round(4) == np.array([0, 0, 0, 0, 0]))
    assert zscore(single_element_vector, transform=True).round(4) == np.array([0])