import h5py
import pypesto.optimize
import pypesto.objective
from pypesto.store.save_to_hdf5 import get_or_create_group
from pypesto.result import Result
from pypesto.store import write_result


def write_pypesto_history(fn_parpe,
                          fn_pypesto,
                          i_ms):
    """
    Write the History from parPE file directly to pypesto file.
    """
    f = h5py.File(fn_parpe, 'r')
    g = h5py.File(fn_pypesto, 'w')

    history_grp = get_or_create_group(g, 'history')
    ms_grp = get_or_create_group(history_grp, str(i_ms))
    trace_grp = get_or_create_group(ms_grp, 'trace')
    iterations = len(f[f'multistarts/{i_ms}/iteration'])

    for i_iter in range(iterations):
        iteration = get_or_create_group(trace_grp, str(i_iter))
        iteration['fval'] = f[f'multistarts/{i_ms}/iterCostFunCost'][:, i_iter]
        iteration['grad'] = f[f'multistarts/{i_ms}/iterCostFunGradient'][:, i_iter]
        iteration['x'] = f[f'multistarts/{i_ms}/iterCostFunParameters'][:, i_iter]
        iteration['time'] = f[f'multistarts/{i_ms}/iterCostFunWallSec'][:, i_iter]

    f.close()
    g.close()


def get_optimizer_result_from_parPE(fn_parpe,
                                    i_ms):
    """
    Fill in a `OptimizerResult` object from parpe history data.
    """

    result = pypesto.optimize.OptimizerResult()

    with h5py.File(fn_parpe, 'r') as f:
        result.id = str(i_ms)
        result.x = f[f'multistarts/{i_ms}/finalParameters'][()]
        result.grad = f[f'multistarts/{i_ms}/iterCostFunGradient'][:, -1]
        result.fval = f[f'multistarts/{i_ms}/finalCost'][()]
        result.x0 = f[f'multistarts/{i_ms}/initialParameters'][()]
        result.fval0 = f[f'multistarts/{i_ms}/iterCostFunParameters'][0, :]
        result.n_fval = 0
        result.n_grad = 0
        for i_iter in f[f'multistarts/{i_ms}/iteration']:
            result.n_grad += f[f'multistarts/{i_ms}/iteration/{i_iter}/costFunGradient'].shape[1]
            result.n_fval += f[f'multistarts/{i_ms}/iteration/{i_iter}/costFunCost'].shape[1]
        result.exitflag = f[f'multistarts/{i_ms}/exitStatus'][()]

    return result


def parpe_to_pypesto_history(fn_parpe,
                             fn_pypesto):
    """
    Convert a parPE history file to a pypesto History file.
    """

    with h5py.File(fn_parpe, 'r') as f:
        result = Result()
        for i_ms in f['multistarts']:
            write_pypesto_history(fn_parpe=fn_parpe,
                                  fn_pypesto=fn_pypesto,
                                  i_ms=i_ms)
            result.optimize_result.append(get_optimizer_result_from_parPE(fn_parpe, i_ms))
        write_result(result, fn_pypesto, problem=False, sample=False, profile=False)
