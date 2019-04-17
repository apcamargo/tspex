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
Command-line interface for tspex.
"""

import argparse
import sys

import pandas as pd
import tspex


def execute_tspex(input_file, output_file, method, log, disable_transformation, threshold):
    """Compute gene tissue-specificity from a expression matrix file."""
    transform = not disable_transformation
    expression_matrix = pd.read_csv(input_file, index_col=0, header=0, sep=None, engine='python')
    tissue_specificity = tspex.TissueSpecificity(
        expression_matrix, method, log, transform=transform, threshold=threshold)
    tissue_specificity.tissue_specificity.to_csv(output_file, sep='\t')


def main():
    method_choices = ['counts', 'tsi', 'tau', 'gini', 'simpson', 'shannon_specificity',
                      'roku_specificity', 'zscore', 'spm', 'spm_dpm', 'js_specificity',
                      'js_specificity_dpm']
    parser = argparse.ArgumentParser(
        description='Compute gene tissue-specificity from a expression matrix file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_file',
                        help='Expression matrix file in the TSV or CSV format.')
    parser.add_argument('output_file',
                        help='Output TSV file containing tissue-specificity values.')
    parser.add_argument(
        'method', choices=method_choices, metavar='method',
        help='Tissue-specificity metric. Allowed values are: "' + '", "'.join(method_choices) + '".')
    parser.add_argument('-l', '--log',
                        action='store_true',
                        help='Log-transform expression values.')
    parser.add_argument(
        '-d', '--disable_transformation', action='store_true',
        help=('By default, tissue-specificity values are transformed so that they range from 0 '
              '(perfectly ubiquitous) to 1 (perfectly tissue-specific). If this parameter is used, '
              'transformation will be disabled and each metric will have have a diferent range of '
              'values.'))
    parser.add_argument(
        '-t', '--threshold', default=0, type=int,
        help=('Threshold to use with the "counts" metric. If another method if chosen, this '
              'parameter will be ignored.'))
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    execute_tspex(**vars(args))


if __name__ == '__main__':
    main()
