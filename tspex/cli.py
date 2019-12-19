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


def tspex_cli(input_file, output_file, method, log, disable_transformation, threshold):
    """Compute gene tissue-specificity from a expression matrix file and save an output file."""
    transform = not disable_transformation
    if input_file.rsplit('.', 1)[1].lower() in ['xls', 'xlsx']:
        expression_matrix = pd.read_excel(
            input_file, index_col=0, header=0, thousands=','
        )
    else:
        expression_matrix = pd.read_csv(
            input_file, index_col=0, header=0, sep=None, thousands=',', engine='python'
        )
    tissue_specificity = tspex.TissueSpecificity(
        expression_matrix, method, log, transform=transform, threshold=threshold
    )
    tissue_specificity.tissue_specificity.to_csv(output_file, sep='\t')


def main():
    method_choices = [
        'counts',
        'tau',
        'gini',
        'simpson',
        'shannon_specificity',
        'roku_specificity',
        'tsi',
        'zscore',
        'spm',
        'spm_dpm',
        'js_specificity',
        'js_specificity_dpm',
    ]
    parser = argparse.ArgumentParser(
        description='Compute gene tissue-specificity from an expression matrix and save the output.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--version', action='version', version='%(prog)s 0.6.1')
    parser.add_argument(
        'input_file', help='Expression matrix file in the TSV, CSV or Excel formats.'
    )
    parser.add_argument(
        'output_file', help='Output TSV file containing tissue-specificity values.'
    )
    parser.add_argument(
        'method',
        choices=method_choices,
        metavar='method',
        help='Tissue-specificity metric. Allowed values are: "'
        + '", "'.join(method_choices)
        + '".',
    )
    parser.add_argument(
        '-l', '--log', action='store_true', help='Log-transform expression values.'
    )
    parser.add_argument(
        '-d',
        '--disable_transformation',
        action='store_true',
        help=(
            'By default, tissue-specificity values are transformed so that they range from 0 '
            '(perfectly ubiquitous) to 1 (perfectly tissue-specific). If this parameter is used, '
            'transformation will be disabled and each metric will have have a diferent range of '
            'possible values.'
        ),
    )
    parser.add_argument(
        '-t',
        '--threshold',
        default=0,
        type=int,
        help=(
            'Threshold to be used with the "counts" metric. If another method is chosen, this '
            'parameter will be ignored.'
        ),
    )
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    tspex_cli(**vars(args))
