import argparse
import sys

import pandas as pd
import tspex


def main(input_file, output_file, method, log, transform, threshold):
    """Compute gene tissue-specificity from a expression matrix file."""
    expression_matrix = pd.read_csv(input_file, index_col=0, header=0, sep=None, engine='python')
    tissue_specificity = tspex.TissueSpecificity(expression_matrix, method, log, method, threshold)
    tissue_specificity.tissue_specificity.to_csv(output_file, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Classify sequences from a input FASTA file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_file',
                        help='expression matrix file in the TSV or CSV format.')
    parser.add_argument('output_file',
                        help='output TSV file containing the tissue-specificity values.')
    parser.add_argument('method',
                        choices=['counts','tsi', 'tau', 'gini', 'simpson', 'shannon_specificity',
                                 'roku_specificity', 'zscore', 'spm', 'spm_dpm', 'js_specificity',
                                 'js_specificity_dpm'],
                        help='tissue-specificity metric.')
    parser.add_argument('--log',
                        action='store_true',
                        help='foo.')
    parser.add_argument('--disable_transformation',
                        action='store_false', dest='transform',
                        help='foo.')
    parser.add_argument('--threshold',
                        default=0, type=int,
                        help='foo.')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    main(**vars(args))
