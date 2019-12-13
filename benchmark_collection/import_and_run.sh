#!/usr/bin/env bash
# Import, build, and run model from the benchmark collection

set -e

get_abs_filename() {
  # $1 : relative filename
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

if [[ $# -lt 1 ]] || [[ $# -gt 1 ]]; then
    echo "Import, build, and run model from the benchmark collection"
    echo "USAGE: $(basename "$0") /path/to/model/dir"
    exit 1;
fi

PETAB_MODEL_DIR=$1
SCRIPT_PATH=$(get_abs_filename $(dirname $BASH_SOURCE))
PARPE_DIR=${SCRIPT_PATH}/..
MODEL_NAME=$(basename ${PETAB_MODEL_DIR})
AMICI_MODEL_DIR=${SCRIPT_PATH}/${MODEL_NAME}

cd ${PETAB_MODEL_DIR}

echo "Running petablint..."
petablint -v -n ${MODEL_NAME}

# import AMICI model
if [[ ! -d ${AMICI_MODEL_DIR} ]]; then
    echo "Importing model..."
    CMD="amici_import_petab.py --verbose -n ${MODEL_NAME} -o ${AMICI_MODEL_DIR}"
    echo "${CMD}"
    ${CMD}
fi

cd ${SCRIPT_PATH}

# parPE build
if [[ ! -d parpe_${MODEL_NAME} ]]; then
    echo "Setting up parPE..."
    ${PARPE_DIR}/misc/setup_amici_model.sh ${AMICI_MODEL_DIR} parpe_${MODEL_NAME}
else
    if [[ -z "${AMICI_ROOT}" ]]
        then AMICI_ROOT=${SCRIPT_PATH}/../deps/AMICI/
    fi
    (cd parpe_${MODEL_NAME}/build && cmake -DAmici_DIR=${AMICI_ROOT}/build .. && make)
fi

# generate data file
echo "Importing data..."
CMD="${PARPE_DIR}/misc/generateHDF5DataFileFromText.py \
    -o parpe_${MODEL_NAME}/${MODEL_NAME}.h5 \
    -s ${PETAB_MODEL_DIR}/model_${MODEL_NAME}.xml \
    -d ${AMICI_MODEL_DIR} \
    -m ${PETAB_MODEL_DIR}/measurementData_${MODEL_NAME}.tsv \
    -c ${PETAB_MODEL_DIR}/experimentalCondition_${MODEL_NAME}.tsv \
    -n ${MODEL_NAME} \
    -p ${PETAB_MODEL_DIR}/parameters_${MODEL_NAME}.tsv"
echo $CMD
$CMD

echo ""
echo "Start parameter estimation by:"
echo "parpe_${MODEL_NAME}/build/estimate_${MODEL_NAME} -o parpe_${MODEL_NAME}_results/ parpe_${MODEL_NAME}/${MODEL_NAME}.h5"

echo ""
echo "Running simulation with nominal parameters..."
parpe_${MODEL_NAME}/build/simulateNominal_${MODEL_NAME} parpe_${MODEL_NAME}/${MODEL_NAME}.h5
