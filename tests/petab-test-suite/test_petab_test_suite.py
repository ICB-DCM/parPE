#!/usr/bin/env python3

"""Run PEtab test suite (https://github.com/PEtab-dev/petab_test_suite)"""
import logging
import os
import re
import shutil
import subprocess
import sys
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
                or "Timepoint-specific parameter overrides" in str(e):
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
        raise AssertionError(f"{' '.join(cmd)} exited with status "
                             f"{ret.returncode}:\n"
                             f"{ret.stdout}")

    return ret


def _test_case(case: Union[int, str]) -> None:
    """Run a single PEtab test suite case"""

    case = petabtests.test_id_str(case)
    logger.debug(f"Case {case}")

    # Ignore some
    if case in ['0012']:
        raise NotImplementedError(case)

    case_dir = os.path.join(petabtests.CASES_DIR, case)
    yaml_file = os.path.join(case_dir, petabtests.problem_yaml_name(case))
    test_dir = os.path.dirname(os.path.abspath(__file__))
    parpe_dir = os.path.dirname(os.path.dirname(test_dir))
    model_name = f'model_{case}'
    amici_model_dir = f'petab_test_models/model_{case}'
    parpe_model_dir = f'petab_test_models/model_{case}_parpe'
    hdf5_input = os.path.join(parpe_model_dir, 'input.h5')
    hdf5_output = os.path.join(parpe_model_dir, 'out.h5')

    # create amici model from PEtab
    cmd = ['amici_import_petab', '-y', yaml_file, '-o', amici_model_dir,
           '-n', model_name]
    check_run(cmd)

    # must not exist when calling setup_amici_model.sh
    if os.path.isdir(parpe_model_dir):
        shutil.rmtree(parpe_model_dir)

    # set up for parPE
    cmd = [os.path.join(parpe_dir, 'misc', 'setup_amici_model.sh'),
           amici_model_dir, parpe_model_dir]
    check_run(cmd)

    # create input hdf5 file
    cmd = ['parpe_petab_to_hdf5',
           '-y', yaml_file, '-d', amici_model_dir,
           '-n', model_name, '-o', hdf5_input]
    check_run(cmd)

    # simulate model using nominal parameters
    cmd = [os.path.join(parpe_model_dir, 'build',
                        f'simulateNominal_{model_name}'),
           hdf5_input, hdf5_output]
    ret = check_run(cmd)

    # check output
    print(ret.stdout)
    g = re.search(r'Likelihood: (\d+\.\d+)', ret.stdout).group(0)
    llh_actual = - float(g.split(' ')[1])
    print(llh_actual)
    solution = petabtests.load_solution(case)

    gt_llh = solution[petabtests.LLH]
    assert llh_actual == pytest.approx(gt_llh)

    # FIXME 0007, 0009, 0011, 0010, 0013


def run() -> None:
    """Run the full PEtab test suite"""

    n_success = 0
    n_skipped = 0
    for case in petabtests.CASES_LIST:
        try:
            test_case(case)
            n_success += 1
        except Skipped:
            n_skipped += 1
        except Exception as e:
            # run all despite failures
            logger.error(f"Case {case} failed.")
            logger.error(e)

    logger.info(f"{n_success} / {len(petabtests.CASES_LIST)} successful, "
                f"{n_skipped} skipped")
    if n_success != len(petabtests.CASES_LIST):
        sys.exit(1)


if __name__ == '__main__':
    run()
