#!/usr/bin/env bash
# Import, build, and run model from the benchmark collection

set -eou pipefail

get_abs_filename() {
  # $1 : relative filename
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

check_output() {
  # Check output
  local nllh=$(grep Likelihood "${output_file}" | tr -cd '[:print:]' | sed -r 's/.*Likelihood: (.*)\[.*/\1/')
  rm "${output_file}"
  local reference_file="${amici_root}/tests/benchmark-models/benchmark_models.yaml"
  local ref=$(shyaml get-value "${model_name}".llh <"$reference_file")
  local abs="define abs(i) {\\nif (i < 0) return (-i) \nreturn (i)\n}\n"

  # Do we match within tolerance?
  # LLH VS nllh!
  if (($(echo -e "${abs}\nabs($nllh + $ref) < 0.001 * abs($ref)" | bc))); then
    echo "OKAY: Expected llh $ref, got nllh $nllh"
  else
    echo "FAILED: Expected llh $ref, got nllh $nllh"
    exit 1
  fi
}

if [[ $# -lt 1 ]] || [[ $# -gt 1 ]]; then
  echo "Import, build, and run model from the benchmark collection"
  echo "USAGE: $(basename "$0") /path/to/model/dir"
  exit 1
fi

petab_model_dir=$1
script_path=$(get_abs_filename "$(dirname "$BASH_SOURCE")")
parpe_dir=${script_path}/..
model_name=$(basename "${petab_model_dir}")
amici_model_dir=${script_path}/${model_name}
petab_yaml=$(python -c "import benchmark_models_petab; print(benchmark_models_petab.get_problem_yaml_path('${model_name}'))")
output_file=tmp.out
hdf5_infile="parpe_${model_name}/${model_name}.h5"
estimate_exe="parpe_${model_name}/build/estimate_${model_name}"

amici_root="${AMICI_ROOT:-${parpe_dir}/deps/AMICI/}"

cd "${petab_model_dir}"

echo "Running petablint on ${petab_yaml}..."
petablint -v ${petab_yaml}

# problems we need to flatten
to_flatten=(
  "Bruno_JExpBot2016" "Chen_MSB2009" "Crauste_CellSystems2017"
  "Fiedler_BMCSystBiol2016" "Fujita_SciSignal2010" "SalazarCavazos_MBoC2020"
)
flatten=""
for item in "${to_flatten[@]}"; do
  if [[ "$item" == "$model_name" ]]; then
    flatten="--flatten"
    break
  fi
done

# import AMICI model
if [[ ! -d ${amici_model_dir} ]]; then
  echo "Importing model..."
  cmd="amici_import_petab --verbose ${petab_yaml} -n ${model_name} -o ${amici_model_dir} ${flatten}"
  echo "${cmd}"
  ${cmd}
fi

cd "${script_path}"

# parPE build
if [[ ! -d parpe_${model_name} ]]; then
  echo "Setting up parPE..."
  "${parpe_dir}"/misc/setup_amici_model.sh "${amici_model_dir}" parpe_"${model_name}"
else
  (cd parpe_"${model_name}"/build && cmake -DAmici_DIR="${amici_root}"/build .. && make)
fi

# generate data file
echo "Importing data..."
cmd="parpe_petab_to_hdf5 \
    -o ${hdf5_infile} \
    -d ${amici_model_dir} \
    -y ${petab_yaml} \
    -n ${model_name} \
    --ignore-initialization-priors \
    ${flatten}"
echo "$cmd"
$cmd

echo ""
echo "Start parameter estimation by:"
echo "${estimate_exe} -o parpe_${model_name}_results/ ${hdf5_infile}"

echo ""
echo "Running simulation with nominal parameters..."
cmd="parpe_${model_name}/build/simulateNominal_${model_name} ${hdf5_infile}"
echo "$cmd"
$cmd |& tee "${output_file}"

check_output

echo "Checking gradient..."
echo PARPE_NO_DEBUG=1 "${estimate_exe}" -t gradient_check -o "parpe_${model_name}_results/gradient_check" "${hdf5_infile}"
PARPE_NO_DEBUG=1 "${estimate_exe}" -t gradient_check -o "parpe_${model_name}_results/gradient_check" "${hdf5_infile}"
