# tspex

`tspex` is a Python package for calculating tissue-specificity metrics from gene expression data.

`tspex` features include:
  - An easy-to-use object-oriented interface.
  - A command-line interface.
  - Twelve different tissue-specificity metrics.
  - Integration with pandas.
  - Graphing functions.
  - Support for Jupyter notebooks.

## Installation

```
# Using pip
$ pip install tspex

# using conda
$ conda install -c bioconda tspex
```

## Command-line interface

`tspex` can be run from the command line to rapidly process a expression matrix file and generate an
output containing the computed tissue-specificity values. Here is the usage of this interface:

```
usage: tspex [-h] [-l] [-d] [-t THRESHOLD] input_file output_file method

Compute gene tissue-specificity from an expression matrix and save the output.

positional arguments:
  input_file                           Expression matrix file in the TSV or CSV formats.
  output_file                          Output TSV file containing tissue-specificity
                                       values.
  method                               Tissue-specificity metric. Allowed values are:
                                       "counts", "tsi", "tau", "gini", "simpson",
                                       "shannon_specificity", "roku_specificity",
                                       "zscore", "spm", "spm_dpm", "js_specificity",
                                       "js_specificity_dpm".

optional arguments:
  -h, --help                           show this help message and exit
  -l, --log                            Log-transform expression values.
  -d, --disable_transformation         By default, tissue-specificity values are
                                       transformed so that they range from 0 (perfectly
                                       ubiquitous) to 1 (perfectly tissue-specific). If
                                       this parameter is used, transformation will be
                                       disabled and each metric will have have a diferent
                                       range of possible values.
  -t THRESHOLD, --threshold THRESHOLD  Threshold to be used with the "counts" metric. If
                                       another method is chosen, this parameter will be
                                       ignored.
```

To compute the SPM values of a log-transformed expression matrix through the command-line interface:

```
tspex --log gene_expression.tsv gene_spm.tsv spm
```
