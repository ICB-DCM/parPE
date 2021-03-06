#!/usr/bin/env python3
"""Run steadystate example for testing

Expects to be run from within the steadystate example build directory
"""

import contextlib
import os
import shutil
import subprocess

import h5py
import pytest

# General setup
script_path = os.path.dirname(os.path.abspath(__name__))
# test executables are expected in working directory
#  (i.e. working dir is build/examples/parpeamici/steadystate)
cwd = os.getcwd()
HDF5_FILE = os.path.join(
    cwd, 'steadystate_scaled-prefix/src/steadystate_scaled/',
    'example_data.h5')
HDF5_FILE_TEST = os.path.join(
    cwd, 'steadystate_scaled-prefix/src/steadystate_scaled/',
    'example_data-testset.h5')
# build directory is expected to be a direct subdirectory of the parPE root dir
parpe_root = os.path.join(script_path, "../../../..")
OPTIM_OPT = os.path.join(parpe_root, "misc/optimizationOptions.py")
MPIEXEC = os.environ.get('PARPE_TESTS_MPIEXEC',
                         "mpiexec -n 5 --oversubscribe").split(" ")
optim_exe = './example_steadystate_multi'
sim_exe = './example_steadystate_multi_simulator'
print('Files:', HDF5_FILE, HDF5_FILE_TEST, MPIEXEC)

result_filename = '_rank00000.h5'


@pytest.mark.parametrize("hierarchical", [False, True])
def test_nompi_gradient_check(hierarchical):
    """Test gradient check without MPI"""

    outdir = f'example_steadystate_multi-test-gradient-{int(hierarchical)}'
    shutil.rmtree(outdir, ignore_errors=True)

    subprocess.run([OPTIM_OPT, HDF5_FILE, '-s', 'hierarchicalOptimization',
                    str(int(hierarchical))],
                   capture_output=True, check=True, encoding="utf-8")
    # Gradient check may fail in certain parameter regimes. Therefore, try
    # multiple times.
    max_tries = 4
    for _ in range(max_tries):
        ret = subprocess.run(f'{optim_exe} -t gradient_check'
                             f' -o {outdir}/ {HDF5_FILE}'.split(' '),
                             capture_output=True,
                             check=True, encoding="utf-8")
        if '[ERR]' not in ret.stdout:
            # one good pass suffices
            return

    pytest.fail(f"Gradient-check did not pass once in {max_tries} tries.")


def test_nompi_optimization():
    """Run optimization and simulation without MPI with default settings"""

    outdir = 'example_steadystate_multi-test-optimize'
    shutil.rmtree(outdir, ignore_errors=True)

    ret = subprocess.run([optim_exe, '-o', outdir + '/', HDF5_FILE],
                         capture_output=True,
                         check=True, encoding="utf-8")
    assert '[ERR]' not in ret.stdout
    check_optimization_results(
        filename=os.path.join(outdir, result_filename))

    # Simulate at optimum
    sim_file = os.path.abspath('simulate1.h5')
    with contextlib.suppress(FileNotFoundError):
        os.remove(sim_file)

    cmd = [sim_exe,
           os.path.join(outdir, result_filename), '/inputData',
           os.path.join(outdir, result_filename), '/',
           sim_file, '/',
           '--at-optimum', '--nompi', '--compute-inner']
    ret = subprocess.run(cmd, capture_output=True,
                         check=True, encoding="utf-8")
    assert os.path.isfile(sim_file)
    assert '[ERR]' not in ret.stdout
    assert '[WRN]' not in ret.stdout
    assert 'exception' not in ret.stdout

    check_simulation_results(filename=sim_file, sim_type="at-optimum")


def test_mpi_optimization():
    """Test optimization and simulation with MPI"""

    # Run optimization with default settings
    outdir = 'example_steadystate_multi-test-optimize'
    shutil.rmtree(outdir, ignore_errors=True)
    cmd = [*MPIEXEC, optim_exe, '--mpi', '-o', outdir + '/', HDF5_FILE]
    ret = subprocess.run(cmd, capture_output=True,
                         check=True, encoding="utf-8")
    assert 'Maximum Number of Iterations Exceeded' in ret.stdout \
           or 'Solved To Acceptable Level' in ret.stdout

    check_optimization_results(
        filename=os.path.join(outdir, result_filename))

    # Simulate along trajectory
    sim_file = os.path.abspath('simulate2.h5')
    with contextlib.suppress(FileNotFoundError):
        os.remove(sim_file)

    cmd = [*MPIEXEC, sim_exe,
           os.path.join(outdir, result_filename), '/inputData',
           os.path.join(outdir, result_filename), '/',
           sim_file, '/',
           '--along-trajectory', '--mpi', '--nocompute-inner']
    ret = subprocess.run(cmd, capture_output=True,
                         check=True, encoding="utf-8")

    assert os.path.isfile(sim_file), ret
    assert '[ERR]' not in ret.stdout
    assert '[WRN]' not in ret.stdout
    assert 'exception' not in ret.stdout

    check_simulation_results(filename=sim_file, sim_type="along-trajectory")

    # Simulate on test set
    sim_file = os.path.abspath('simulate3.h5')
    with contextlib.suppress(FileNotFoundError):
        os.remove(sim_file)

    cmd = [*MPIEXEC, sim_exe,
           HDF5_FILE_TEST, '/',
           os.path.join(outdir, result_filename), '/',
           sim_file, '/',
           '--at-optimum', '--mpi', '--compute-inner']
    ret = subprocess.run(cmd, capture_output=True,
                         check=True, encoding="utf-8")
    assert os.path.isfile(sim_file), ret
    assert '[ERR]' not in ret.stdout
    assert '[WRN]' not in ret.stdout
    assert 'exception' not in ret.stdout

    check_simulation_results(filename=sim_file, sim_type="at-optimum")


def check_optimization_results(filename: str) -> None:
    """Check for presence of optimization results"""

    try:
        with h5py.File(filename, "r") as f:
            # TODO: extend checks
            assert f["/multistarts/0/iterCostFunParameters"].size > 0
    except Exception as e:
        import sys
        raise type(e)(str(e) + f' occurred in {filename}') \
            .with_traceback(sys.exc_info()[2])


def check_simulation_results(filename: str, sim_type: str) -> None:
    """Check for presence of simulation results"""

    try:
        with h5py.File(filename, "r") as f:
            # TODO: extend checks
            if sim_type == "at-optimum":
                assert len(f["/multistarts/0/yMes"])\
                       == len(f["/multistarts/0/ySim"])
            elif sim_type == "along-trajectory":
                assert len(f["/multistarts/0/iter/0/yMes"])\
                       == len(f["/multistarts/0/iter/0/ySim"])
            else:
                raise ValueError(f"Unknown sim_type: {sim_type}")

    except Exception as e:
        import sys
        raise type(e)(str(e) + f' occurred in {filename}') \
            .with_traceback(sys.exc_info()[2])
