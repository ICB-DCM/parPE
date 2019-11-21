#!/usr/bin/env python3
"""Read final parameter vectors from parPE (multi-start) optimization
result file and write to TSV file."""

import h5py
import argparse
import pandas as pd
import numpy as np


def final_parameters_to_df(hdf5_result_file: str,
                           df: pd.DataFrame = None) -> pd.DataFrame:
    """Read final parameter vectors from parPE (multi-start) optimization
    result file

    Row names are parameter Ids. Parameter vectors as rows.

    Arguments:
        hdf5_result_file: Filename to read from
        df: Dataframe to append to
    """
    rows = []
    indices = []

    with h5py.File(hdf5_result_file, 'r') as f:
        parameter_names = f['/inputData/parameters/parameterNames'][:]
        for mspath in f['/multistarts/']:
            # for potentially multiple cost function evaluations
            parameters_path = '/multistarts/%s/iterCostFunParameters' % mspath
            if parameters_path not in f:
                print('Skipping due to missing parameter trajectory.')
                return df

            parameters = f[parameters_path][:]

            for par_vec_idx in range(parameters.shape[1] - 1, 0, -1):
                if not np.all(np.isnan(parameters[:, par_vec_idx])):
                    d = dict(zip(parameter_names, parameters[:, par_vec_idx]))
                    indices.append(f'{hdf5_result_file}-start-{mspath}')
                    rows.append(d)
                    break

    if not rows:
        return df

    if df is None:
        df = pd.DataFrame(rows, index=indices)
    else:
        df = df.append(pd.DataFrame(rows, index=indices), sort=True)

    return df


def parse_cli_args():
    """Parse command line arguments

    Returns:
        parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Read final parameter vectors from parPE (multi-start) '
                    'optimization result file and write to tsv.')

    parser.add_argument('infiles', nargs='*',
                        help='Result files to read from')

    parser.add_argument('outfile', help='Output file')

    args = parser.parse_args()

    if not args.infiles:
        raise ValueError('No input files specified.')

    if not args.outfile:
        raise ValueError('No output file specified.')

    return args


def main():
    args = parse_cli_args()

    df = None
    for infile in args.infiles:
        print('Processing', infile)
        df = final_parameters_to_df(hdf5_result_file=infile, df=df)

    print('Writing to', args.outfile)
    df.to_csv(args.outfile, sep='\t')


if __name__ == '__main__':
    main()
