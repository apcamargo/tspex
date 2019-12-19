# tspex


- [Overview](#overview)
- [Documentation](#documentation)
- [Installation](#installation)
- [Python API tutorial](#python-api-tutorial)
- [Command-line interface](#command-line-interface)
- [Examples](#examples)

## Overview

tspex is a tissue-specificity calculator tool. It provides both an easy-to-use object-oriented Python API and a command-line interface (CLI) for calculating a variety of tissue-specificity metrics from gene expression data.

tspex features include:
  - Twelve different tissue-specificity metrics.
  - Integration with popular data analysis libraries, such as NumPy, SciPy, and pandas.
  - Visualization functions.
  - Support for Jupyter notebooks.

## Documentation

A complete documentation for tspex can be found at [https://apcamargo.github.io/tspex/](https://apcamargo.github.io/tspex/).

## Installation

There are two ways to install tspex:

- Using pip:

```
pip install tspex
```

- Using conda:

```
conda install -c bioconda tspex
```


## Python API tutorial

For a detailed guide on how to use the Python API, please check the [Jupyter notebook tutorial](https://github.com/apcamargo/tspex/blob/master/docs/python_api.ipynb).


## Command-line interface

tspex can be executed from the command line using the `tspex` command. It takes an expression matrix file as input and outputs the computed tissue-specificity values.


```
usage: tspex [-h] [-l] [-d] [-t THRESHOLD] input_file output_file method

Compute gene tissue-specificity from an expression matrix and save the output.

positional arguments:
  input_file            Expression matrix file in the TSV, CSV or Excel
                        formats.
  output_file           Output TSV file containing tissue-specificity values.
  method                Tissue-specificity metric. Allowed values are:
                        "counts", "tau", "gini", "simpson",
                        "shannon_specificity", "roku_specificity", "tsi",
                        "zscore", "spm", "spm_dpm", "js_specificity",
                        "js_specificity_dpm".

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -l, --log             Log-transform expression values. (default: False)
  -d, --disable_transformation
                        By default, tissue-specificity values are transformed
                        so that they range from 0 (perfectly ubiquitous) to 1
                        (perfectly tissue-specific). If this parameter is
                        used, transformation will be disabled and each metric
                        will have have a diferent range of possible values.
                        (default: False)
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold to be used with the "counts" metric. If
                        another method is chosen, this parameter will be
                        ignored. (default: 0)
```

### Examples

- Using the `spm` metric to compute tissue-specificity values from a log-transformed expression matrix:

```
tspex --log gene_expression.tsv tspex_spm.tsv spm
```

- Using the `counts` method to compute tissue-specificity by counting the number of tissues in which the gene expression is greater than 10:

```
tspex --threshold 10 gene_expression.tsv tspex_counts.tsv counts
```

- Using the `zscore` without transformation to quantify tissue-specificity as the number of standard deviations away from the mean gene expression:

```
tspex --disable_transformation gene_expression.tsv tspex_zscore.tsv zscore
```


## Cite tspex

If you use tspex in your research, it would be appreciated if you could cite it.

[![DOI](https://zenodo.org/badge/172384291.svg)](https://zenodo.org/badge/latestdoi/172384291)

> Camargo, Antonio P., "*tspex: tissue-specificity calculator*" (2019). doi:10.5281/zenodo.3558708
