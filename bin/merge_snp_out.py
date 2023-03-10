#!/usr/bin/env python3

import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infiles', nargs='+', required=True)
    parser.add_argument('--amr_csv', required=True)
    parser.add_argument('--outfile', required=False)
    return parser.parse_args()

def update_sample_column(f, amr):
    split_string = '_AMR_analytic_matrix_with_SNP_confirmation.csv'
    sample = f.split(split_string)[0].split('/')[-1]
    temp = pd.read_csv(f)
    amr[sample] = temp[sample]

def main():
    args = parse_args()
    amr = pd.read_csv(args.amr_csv)
    [update_sample_column(f, amr) for f in args.infiles]
    amr.to_csv(args.outfile, index=False)


if __name__ == '__main__':
    main()