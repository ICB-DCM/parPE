#!/usr/bin/env python3
"""Run steadystate example for testing

Expects to be run from within the steadystate example build directory
"""

import contextlib
import os
import shutil
import subprocess


# General setup
script_path = os.path.dirname(os.path.abspath(__name__))
# test executables are expected in working directory
cwd = os.getcwd()
HDF5_FILE = os.path.join(
    cwd, 'steadystate_scaled-prefix/src/steadystate_scaled/',
    'example_data.h5')
HDF5_FILE_TEST = os.path.join(
    cwd, 'steadystate_scaled-prefix/src/steadystate_scaled/',
    'example_data-testset.h5')

MPIEXEC = ['mpiexec', '--oversubscribe', '-n', '5']
optim_exe = './example_steadystate_multi'
sim_exe = './example_steadystate_multi_simulator'
print('Files:', HDF5_FILE, HDF5_FILE_TEST, MPIEXEC)

result_filename = '_rank00000.h5'

with contextlib.suppress(subprocess.CalledProcessError):
    # Allow running as root in docker
    subprocess.run('grep docker /proc/1/cgroup -qa', shell=True, check=True)
    MPIEXEC.append("--allow-run-as-root")

    # If we are running in docker, we generally don't have SYS_PTRACE
    #  permissions  and thus, cannot use vader. Also disable Infiniband.
    subprocess.run('mpiexec --version | grep open-mpi', shell=True, check=True)
    MPIEXEC.extend(["--oversubscribe",
                    "--mca", "btl_vader_single_copy_mechanism", "none",
                    "--mca", "btl", "^openib",
                    "--mca", "oob_tcp_if_include", "lo",
                    "--mca", "btl_tcp_if_include", "lo",
                    "--mca", "orte_base_help_aggregate", "0"])


def test_nompi_gradient_check():
    """Test gradient check without MPI"""

    outdir = 'example_steadystate_multi-test-gradient'
    shutil.rmtree(outdir, ignore_errors=True)
    ret = subprocess.run(f'{optim_exe} -t gradient_check'
                         f' -o {outdir}/ {HDF5_FILE}'.split(' '),
                         capture_output=True,
                         check=True, encoding="utf-8")
    assert '[ERR]' not in ret.stdout


def test_nompi_optimization():
    """Run optimization and simulation without MPI with default settings"""

    outdir = 'example_steadystate_multi-test-optimize'
    shutil.rmtree(outdir, ignore_errors=True)

    ret = subprocess.run([optim_exe, '-o', outdir + '/', HDF5_FILE],
                         capture_output=True,
                         check=True, encoding="utf-8")
    assert '[ERR]' not in ret.stdout

    # Simulate at optimum
    sim_file = 'simulate1.h5'
    with contextlib.suppress(FileNotFoundError):
        os.remove(sim_file)

    cmd = [sim_exe, os.path.join(outdir, result_filename), '/',
           sim_file, '/',
           '--at-optimum', '--nompi']
    ret = subprocess.run(cmd, capture_output=True,
                         check=True, encoding="utf-8")
    assert os.path.isfile(sim_file)
    assert '[ERR]' not in ret.stdout
    assert '[WRN]' not in ret.stdout
    assert 'exception' not in ret.stdout
    # test dataset exists
    subprocess.run(['h5dump', '-d', '/multistarts/0/yMes/3', sim_file],
                   check=True)


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

    # Simulate along trajectory
    sim_file = 'simulate2.h5'
    with contextlib.suppress(FileNotFoundError):
        os.remove(sim_file)

    cmd = [*MPIEXEC, sim_exe, os.path.join(outdir, result_filename), '/',
           sim_file, '/',
           '--along-trajectory', '--mpi']
    ret = subprocess.run(cmd, capture_output=True,
                         check=True, encoding="utf-8")

    assert os.path.isfile(sim_file)
    assert '[ERR]' not in ret.stdout
    assert '[WRN]' not in ret.stdout
    assert 'exception' not in ret.stdout
    # test dataset exists
    subprocess.run(['h5dump', '-d', '/multistarts/0/iter/1/yMes/3', sim_file],
                   check=True)

    # Simulate on test set
    sim_file = 'simulate3.h5'
    with contextlib.suppress(FileNotFoundError):
        os.remove(sim_file)

    cmd = [*MPIEXEC, sim_exe, HDF5_FILE_TEST, '/',
           os.path.join(outdir, result_filename), '/',
           sim_file, '/',
           '--at-optimum', '--mpi']
    ret = subprocess.run(cmd, capture_output=True,
                         check=True, encoding="utf-8")
    assert os.path.isfile(sim_file)
    assert '[ERR]' not in ret.stdout
    assert '[WRN]' not in ret.stdout
    assert 'exception' not in ret.stdout
    # test dataset exists
    subprocess.run(['h5dump', '-d', '/multistarts/0/yMes/3', sim_file],
                   check=True)
