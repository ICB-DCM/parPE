#!/usr/bin/env python3
"""
Merge parPE result files.

2017 Daniel Weindl <daniel.weindl@helmholtz-muenchen.de>

"""

import h5py
import glob
import argparse


def merge_files(in_files, out_file, save_iteration, same_input):
    """
    Read final multistart optimization results from HDF5 in_files and write to HDF5 out_file.
    """
    if isinstance(in_files, list):
        in_files_list = in_files
    else:
        in_files_list = glob.glob(in_files)

    f_out = h5py.File(out_file, "w")

    start_idx = 0
    if same_input:
        f_out.create_group('multistarts')
        for in_file in in_files_list:
            print(in_file)
            with h5py.File(in_file, "r") as fIn:
                if '/inputData' not in f_out:
                    fIn.copy('/inputData', f_out)
                copy_multistarts(fIn, f_out, start_idx, save_iteration)
                f_out['multistarts/' + str(start_idx)].attrs['source_file'] = in_file.split(sep='/')[-1]
                start_idx += 1
    else:
        for in_file in in_files_list:
            print(in_file)
            with h5py.File(in_file, "r") as fIn:
                fIn.copy('/', f_out, in_file.split(sep='/')[-1])


def copy_multistarts(f_in, f_out, start_idx, save_iteration):
    """
    Copy all starts in fIn to fOut.
    """
    for start in f_in['/multistarts']:
        if save_iteration:
            f_in.copy('/multistarts/' + start, f_out, 'multistarts/' + str(start_idx))
        else:
            for key in f_in['/multistarts/' + start].keys():
                if not (key == 'iteration'):
                    f_in.copy('/multistarts/' + start + '/' + key, f_out, 'multistarts/' + str(start_idx) + '/' + key)


def parse_cli_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        description='Merge h5 resultfiles.')

    parser.add_argument('-f', '--in-files', dest='in_files',
                        required=True,
                        help='Path for files to be merged. E.g. results/*-rank00000.h5. Or list with file names.')

    parser.add_argument('-o', '--out-file', dest='out_file',
                        required=True,
                        help='Name of HDF5 file to be generated')

    parser.add_argument('-i', '--save-iteration', dest='save_iteration', default=False,
                        help='Save optimizer information for all iterations')

    parser.add_argument('-s', '--same-input', dest='same_input', default=True,
                        help='Merges multistarts with the same inputData.')
    args = parser.parse_args()

    return args


def main():

    args = parse_cli_args()
    
    merge_files(args.in_files, args.out_file, args.save_iteration, args.same_input)


if __name__ == "__main__":
    main()
