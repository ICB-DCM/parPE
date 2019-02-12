#!/bin/bash
#
# Setup virtual environment
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
PARPE_ROOT=$(cd $SCRIPT_PATH/.. && pwd)
BUILD_DIR=$1

if [ ! -d ${BUILD_DIR}/venv ]; then
# create venv
python3 -m venv ${BUILD_DIR}/venv
source ${BUILD_DIR}/venv/bin/activate

# install amici
cd ${PARPE_ROOT}/deps/AMICI/python/sdist
pip3 install -e .
pip3 install termcolor

deactivate
fi

source ${BUILD_DIR}/venv/bin/activate

# install parpe
cd ${PARPE_ROOT}/python
pip3 install -e .
