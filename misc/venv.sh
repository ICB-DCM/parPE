#!/bin/bash
#
# Setup virtual environment for building/testing parPE
#

SCRIPT_PATH=$(dirname $BASH_SOURCE)
PARPE_ROOT=$(cd ${SCRIPT_PATH}/.. && pwd)
BUILD_DIR=$1

# set default build dir if not provided
if [[ -z "${BUILD_DIR}" ]]; then BUILD_DIR="${PARPE_ROOT}/build"; fi

# Save the time for installing AMICI if the venv already exists
# NOTE: Must remove folder if AMICI is updated
if [[ ! -d ${BUILD_DIR}/venv ]]; then
    # create venv
    python3 -m venv ${BUILD_DIR}/venv 2>/dev/null
    # in case this fails (usually due to missing ensurepip), try getting pip
    # manually
    if [[ $? ]]; then
        set -e
        python3 -m venv ${BUILD_DIR}/venv --clear --without-pip
        source ${BUILD_DIR}/venv/bin/activate
        curl https://bootstrap.pypa.io/get-pip.py -o ${BUILD_DIR}/get-pip.py
        python3 ${BUILD_DIR}/get-pip.py
    else
        set -e
        source ${BUILD_DIR}/venv/bin/activate
    fi

    pip3 install wheel pytest

    # install AMICI
    cd ${PARPE_ROOT}/deps/AMICI/python/sdist
    pip3 install -e .

    # install parPE
    cd ${PARPE_ROOT}/python
    pip3 install -e .
    #pip3 install https://github.com/ICB-DCM/PEtab/archive/develop.zip
fi
