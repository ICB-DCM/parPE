#!/bin/bash
#
# Setup virtual environment
#
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
PARPE_ROOT=$(cd $SCRIPT_PATH/.. && pwd)


python3 -m venv ${PARPE_ROOT}/build/venv
source ${PARPE_ROOT}/build/venv/bin/activate
cd ${PARPE_ROOT}/deps/AMICI/python/sdist
pip3 install -e .
