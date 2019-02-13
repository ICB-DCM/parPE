#!/usr/bin/env bash
# For an AMICI-imported model, set up parPE build

set -e

# invalid options:
if [[ $# -lt 2 ]] || [[ $# -gt 2 ]]; then
    echo "For an AMICI-imported model, set up parPE build"
    echo "USAGE: $(basename "$0") MODEL_DIR OUTPUT_DIR"
    exit 1;
fi

SCRIPT_PATH=$(dirname $BASH_SOURCE)
MODEL_DIR=$1
OUTPUT_DIR=$2
TEMPLATE_DIR=${SCRIPT_PATH}/../templates
MODEL_NAME=$(basename ${MODEL_DIR})

if [[ ! -d ${MODEL_DIR} ]]; then
    echo "ERROR: Model directory ${MODEL_DIR} does not exist."
    exit 1
fi

if [[ -e ${OUTPUT_DIR} ]]; then
    echo "ERROR: Output directory ${OUTPUT_DIR} exists. Change or delete."
    exit 1
fi

echo "Copying model to output directory $OUTPUT_DIR ..."
mkdir -p ${OUTPUT_DIR}
cp -R ${MODEL_DIR} ${OUTPUT_DIR}/model

echo "Adding CMake and source files ..."
cp ${TEMPLATE_DIR}/CMakeLists.template.txt ${OUTPUT_DIR}/CMakeLists.txt
sed -ri "s/mymodel/${MODEL_NAME}/" ${OUTPUT_DIR}/CMakeLists.txt
# TODO change in AMICI to avoid naming conflicts
sed -ri "s/simulate_/amici_/" ${OUTPUT_DIR}/model/CMakeLists.txt
cp ${TEMPLATE_DIR}/main*.cpp ${OUTPUT_DIR}

echo "Setting up build ..."
mkdir ${OUTPUT_DIR}/build
cd ${OUTPUT_DIR}/build
CC=mpicc CXX=mpiCC cmake ..

echo "Building ..."
make ${MAKE_OPTS}
