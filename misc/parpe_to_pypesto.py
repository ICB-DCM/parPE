#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Union

import h5py
from pypesto.result.optimize import OptimizerResult as PypestoOptimizerResult
from pypesto.result.result import Result as PypestoResult
from pypesto.store import write_result as write_pypesto_result
from pypesto.C import HISTORY, TRACE, N_ITERATIONS


def write_pypesto_history(
        fn_parpe: Union[str, Path],
        fn_pypesto: Union[str, Path],
        i_ms: Union[int, str]
):
    """
    Write the History from parPE file directly to pypesto file.
    """
    with h5py.File(fn_parpe, 'r') as f:
        with h5py.File(fn_pypesto, 'a') as g:
            trace_grp = g.require_group(f'/{HISTORY}/{i_ms}/{TRACE}')

            iterations = len(f[f'multistarts/{i_ms}/iteration'])

            for i_iter in range(iterations):
                start_group_in = f[f'multistarts/{i_ms}']
                iteration_grp = trace_grp.require_group(str(i_iter))

                iteration_grp['fval'] = \
                    start_group_in['iterCostFunCost'][:, i_iter]
                iteration_grp['grad'] = \
                    start_group_in['iterCostFunGradient'][:, i_iter]
                iteration_grp['x'] = \
                    start_group_in['iterCostFunParameters'][:, i_iter]
                iteration_grp['time'] = \
                    start_group_in['iterCostFunWallSec'][:, i_iter].squeeze()

            trace_grp.attrs[N_ITERATIONS] = iterations


def pypesto_optimizer_result_from_parpe(
        fn_parpe: Union[Path, str],
        i_ms: Union[int, str]
) -> PypestoOptimizerResult:
    """
    Fill in a `pypesto.OptimizerResult` object from parpe history data.
    """
    result = PypestoOptimizerResult()

    with h5py.File(fn_parpe, 'r') as f:
        result.id = str(i_ms)
        ms_grp = f[f'multistarts/{i_ms}/']
        result.x = ms_grp['finalParameters'][()]
        result.grad = ms_grp['iterCostFunGradient'][:, -1]
        result.fval = ms_grp['finalCost'][()]
        result.x0 = ms_grp['initialParameters'][()]
        result.fval0 = ms_grp['iterCostFunParameters'][0, :]
        result.exitflag = ms_grp['exitStatus'][()]
        result.n_fval = 0
        result.n_grad = 0
        for i_iter in ms_grp['iteration']:
            iter_grp = ms_grp[f'iteration/{i_iter}']
            result.n_grad += iter_grp['costFunGradient'].shape[1]
            result.n_fval += iter_grp['costFunCost'].shape[1]

    return result


def parpe_to_pypesto_history(
        fn_parpe: Union[str, Path],
        fn_pypesto: Union[str, Path],
):
    """
    Convert a parPE history file to a pypesto History file.
    """
    result = PypestoResult()
    with h5py.File(fn_parpe, 'r') as f:
        for i_ms in f['multistarts']:
            write_pypesto_history(
                fn_parpe=fn_parpe,
                fn_pypesto=fn_pypesto,
                i_ms=i_ms
            )
            result.optimize_result.append(
                pypesto_optimizer_result_from_parpe(fn_parpe, i_ms)
            )
    write_pypesto_result(
        result, fn_pypesto,
        optimize=True, problem=False, sample=False, profile=False
    )


def parse_args():
    """
    Parse command line arguments

    :return:
        Parsed CLI arguments from :mod:`argparse`.
    """
    parser = argparse.ArgumentParser(
        description='Convert parPE result file to pypesto HDF5 result file')

    parser.add_argument(
        '-i', dest='parpe_file',
        help='parPE result file to convert',
        required=True
    )
    parser.add_argument(
        '-o', dest='pypesto_file',
        help='pypesto file to write',
        required=True
    )

    return parser.parse_args()


def main():
    args = parse_args()
    parpe_to_pypesto_history(fn_pypesto=args.pypesto_file,
                             fn_parpe=args.parpe_file)


if __name__ == '__main__':
    main()
