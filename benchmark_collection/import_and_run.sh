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
SCRIPT_PATH=$(get_abs_filename "$(dirname "$BASH_SOURCE")")
PARPE_DIR=${SCRIPT_PATH}/..
MODEL_NAME=$(basename "${PETAB_MODEL_DIR}")
AMICI_MODEL_DIR=${SCRIPT_PATH}/${MODEL_NAME}
petab_yaml=${MODEL_NAME}.yaml
if [[ -z "${AMICI_ROOT}" ]]
    then AMICI_ROOT=${PARPE_DIR}/deps/AMICI/
fi


cd "${PETAB_MODEL_DIR}"

echo "Running petablint on ${petab_yaml}..."
petablint -v -y "${MODEL_NAME}".yaml

# import AMICI model
if [[ ! -d ${AMICI_MODEL_DIR} ]]; then
    echo "Importing model..."
    CMD="amici_import_petab.py --verbose -y ${MODEL_NAME}.yaml -o ${AMICI_MODEL_DIR}"
    echo "${CMD}"
    ${CMD}
fi

cd "${SCRIPT_PATH}"

# parPE build
if [[ ! -d parpe_${MODEL_NAME} ]]; then
    echo "Setting up parPE..."
    "${PARPE_DIR}"/misc/setup_amici_model.sh "${AMICI_MODEL_DIR}" parpe_"${MODEL_NAME}"
else
    (cd parpe_"${MODEL_NAME}"/build && cmake -DAmici_DIR="${AMICI_ROOT}"/build .. && make)
fi

# generate data file
echo "Importing data..."
CMD="parpe_petab_to_hdf5 \
    -o parpe_${MODEL_NAME}/${MODEL_NAME}.h5 \
    -d ${AMICI_MODEL_DIR} \
    -y ${PETAB_MODEL_DIR}/${MODEL_NAME}.yaml \
    -n ${MODEL_NAME}"
echo "$CMD"
$CMD

echo ""
echo "Start parameter estimation by:"
echo "parpe_${MODEL_NAME}/build/estimate_${MODEL_NAME} -o parpe_${MODEL_NAME}_results/ parpe_${MODEL_NAME}/${MODEL_NAME}.h5"

echo ""
echo "Running simulation with nominal parameters..."
cmd="parpe_${MODEL_NAME}/build/simulateNominal_${MODEL_NAME} parpe_${MODEL_NAME}/${MODEL_NAME}.h5"
echo "$cmd"
$cmd  |& tee tmp.out

# Check output
NLLH=$(grep Likelihood tmp.out | tr -cd '[:print:]' | sed -r 's/.*Likelihood: (.*)\[.*/\1/')
#rm tmp.out
REFERENCE_FILE="${AMICI_ROOT}/tests/benchmark-models/benchmark_models.yaml"
REF=$(shyaml get-value "${MODEL_NAME}".llh < "$REFERENCE_FILE")
ABS="define abs(i) {\\nif (i < 0) return (-i) \nreturn (i)\n}\n"

# Do we match within tolerance?
# LLH VS NLLH!
if (( $(echo -e "${ABS}\nabs($NLLH + $REF) < 0.001 * abs($REF)" | bc) )); then
  echo "OKAY: Expected llh $REF, got nllh $NLLH"
else
  echo "FAILED: Expected llh $REF, got nllh $NLLH"
  exit 1
fi
