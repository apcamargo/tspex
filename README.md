# tspex

tspex is a Python package for calculating tissue-specificity metrics from gene expression data. It provides both an easy-to-use object-oriented interface and a command-line interface.

tspex features include:
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

tspex can be executed from the command line to process an expression matrix and output the computed tissue-specificity values into a file. Here is the usage of this interface:

```
usage: tspex [-h] [-l] [-d] [-t THRESHOLD] input_file output_file method

Compute gene tissue-specificity from an expression matrix and save the output.

positional arguments:
  input_file            Expression matrix file in the TSV or CSV formats.
  output_file           Output TSV file containing tissue-specificity values.
  method                Tissue-specificity metric. Allowed values are:
                        "counts", "tsi", "tau", "gini", "simpson",
                        "shannon_specificity", "roku_specificity", "zscore",
                        "spm", "spm_dpm", "js_specificity",
                        "js_specificity_dpm".

optional arguments:
  -h, --help            show this help message and exit
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