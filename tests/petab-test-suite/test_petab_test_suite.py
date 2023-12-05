#!/usr/bin/env python3

"""Run PEtab test suite (https://github.com/PEtab-dev/petab_test_suite)"""
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Union

import petabtests
import pytest
from _pytest.outcomes import Skipped

from amici.logging import get_logger, set_log_level

logger = get_logger(__name__, logging.DEBUG)
set_log_level(get_logger("amici.petab_import"), logging.DEBUG)
stream_handler = logging.StreamHandler()
logger.addHandler(stream_handler)


def test_case(case: Union[int, str]) -> None:
    """Wrapper for _test_case for handling test outcomes"""
    try:
        _test_case(case)
    except Exception as e:
        if isinstance(e, NotImplementedError) \
                or "models with timepoint specific mappings" in str(e):
            logger.info(f"Case {case} expectedly failed. "
                        "Required functionality is not yet "
                        f"implemented: {e}")
            pytest.skip(str(e))
        else:
            raise e


def check_run(cmd: List[str]) -> subprocess.CompletedProcess:
    """Run a given external command, verify 0-return code, print output on
    failure and raise or return CompletedProcess"""

    ret = subprocess.run(cmd, check=False, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, encoding='utf-8')
    if ret.returncode != 0:
        raise AssertionError(f"{' '.join(map(str, cmd))} exited with status "
                             f"{ret.returncode}:\n"
                             f"{ret.stdout}")

    return ret


def _test_case(case: Union[int, str]) -> None:
    """Run a single PEtab test suite case"""
    model_type = 'sbml'
    version = "v1.0.0"

    case = petabtests.test_id_str(case)
    logger.debug(f"Case {case}")

    # Ignore some that are not supported by AMICI
    #  0012 compartments
    #  0009 newton failure preeq
    #  0018 requires fixing import code for RateRules
    if case in ['0012', '0018']:
        raise NotImplementedError(case)

    case_dir = petabtests.get_case_dir(case, model_type, version)
    yaml_file = case_dir / petabtests.problem_yaml_name(case)
    yaml_file = case_dir / petabtests.problem_yaml_name(case)
    test_dir = Path(__file__).parent.absolute()
    parpe_dir = test_dir.parents[1]
    model_name = f'model_{case}'
    amici_model_dir = Path('petab_test_models') / f'model_{case}'
    parpe_model_dir = Path('petab_test_models') / f'model_{case}_parpe'
    hdf5_input = parpe_model_dir / 'input.h5'
    hdf5_output = parpe_model_dir / 'out.h5'

    # create amici model from PEtab
    cmd = ['amici_import_petab', '-y', yaml_file, '-o', amici_model_dir,
           '-n', model_name]
    if case == '0006':
        cmd.append('--flatten')
    print(" ".join(map(str, cmd)))
    check_run(cmd)

    # must not exist when calling setup_amici_model.sh
    if os.path.isdir(parpe_model_dir):
        shutil.rmtree(parpe_model_dir)

    # set up for parPE
    cmd = [parpe_dir / 'misc' / 'setup_amici_model.sh',
           amici_model_dir, parpe_model_dir]
    print(" ".join(map(str, cmd)))
    check_run(cmd)

    # create input hdf5 file
    cmd = ['parpe_petab_to_hdf5',
           '-y', yaml_file, '-d', amici_model_dir,
           '-n', model_name, '-o', hdf5_input]
    if case == '0006':
        cmd.append('--flatten')
    print(" ".join(map(str, cmd)))
    check_run(cmd)

    # simulate model using nominal parameters
    cmd = [parpe_model_dir / 'build' / f'simulateNominal_{model_name}',
           hdf5_input, hdf5_output]
    print(" ".join(map(str, cmd)))
    ret = check_run(cmd)

    # check output
    g = re.search(r'Likelihood: (\d+\.\d+)', ret.stdout)[0]
    llh_actual = - float(g.split(' ')[1])
    print("Actual llh:", llh_actual)
    solution = petabtests.load_solution(case, model_type, version=version)

    gt_llh = solution[petabtests.LLH]
    assert llh_actual == pytest.approx(gt_llh,
                                       rel=solution[petabtests.TOL_LLH])

    # FIXME
    #  0011 init conc condition table
    #  0013 parametric init conc condition table


def run() -> None:
    """Run the full PEtab test suite"""

    n_success = 0
    n_skipped = 0
    all_cases = list(petabtests.get_cases('sbml', version="v1.0.0"))
    for case in all_cases:
        try:
            test_case(case)
            n_success += 1
        except Skipped:
            n_skipped += 1
        except Exception as e:
            # run all despite failures
            logger.error(f"Case {case} failed.")
            logger.error(e)

    logger.info(f"{n_success} / {len(all_cases)} successful, "
                f"{n_skipped} skipped")
    if n_success != len(all_cases):
        sys.exit(1)


if __name__ == '__main__':
    run()
