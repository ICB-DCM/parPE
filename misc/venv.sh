#!/bin/bash
#
# Setup virtual environment for building/testing parPE
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
PARPE_ROOT=$(cd ${SCRIPT_PATH}/.. && pwd)
BUILD_DIR=$1

# set default build dir if not provided
if [[ -z "${BUILD_DIR}" ]]; then BUILD_DIR="${PARPE_ROOT}/build"; fi

# Save the time for installing AMICI if the venv already exists
# NOTE: Must remove folder if AMICI is updated
if [[ ! -d ${BUILD_DIR}/venv ]]; then
    # create venv
    python3 -m venv ${BUILD_DIR}/venv
    source ${BUILD_DIR}/venv/bin/activate
    pip3 install wheel

    # install AMICI
    cd ${PARPE_ROOT}/deps/AMICI/python/sdist
    pip3 install -e .

    # install parPE
    cd ${PARPE_ROOT}/python
    pip3 install -e .
    pip3 install -U git+https://github.com/ICB-DCM/PEtab.git@develop
fi
