#!/usr/bin/env bash
# For an AMICI-imported model, set up parPE build

set -euo pipefail

# invalid options:
if [[ $# -lt 2 ]] || [[ $# -gt 2 ]]; then
    echo "For an AMICI-imported model, set up parPE build"
    echo "USAGE: $(basename "$0") model_dir output_dir"
    exit 1;
fi

script_path=$(realpath $(dirname "$BASH_SOURCE"))
model_dir="$1"
output_dir="$2"
template_dir="${script_path}/../templates"
model_name=$(basename "${model_dir}")
make_opts="${MAKE_OPTS:-}"

if [[ ! -d "${model_dir}" ]]; then
    echo "ERROR: Model directory ${model_dir} does not exist."
    exit 1
fi

if [[ -e "${output_dir}" ]]; then
    echo "ERROR: Output directory ${output_dir} exists. Change or delete."
    exit 1
fi

echo "Copying model to output directory ${output_dir} ..."
mkdir -p "${output_dir}"
cp -R "${model_dir}" "${output_dir}/model"

echo "Adding CMake and source files ..."
cp "${template_dir}/CMakeLists.template.txt" "${output_dir}/CMakeLists.txt"
sed -ri "s/mymodel/${model_name}/" "${output_dir}/CMakeLists.txt"
# TODO change in AMICI to avoid naming conflicts
sed -ri "s/simulate_/amici_/" "${output_dir}/model/CMakeLists.txt"
cp "${template_dir}"/main*.cpp "${output_dir}"

echo "Setting up build ..."

amici_root=${AMICI_ROOT:-${script_path}/../deps/AMICI/}

cmake -S "${output_dir}" \
      -B "${output_dir}/build" \
      -DAmici_DIR="${amici_root}/build" \
      -DParPE_DIR="${script_path}/../build"

echo "Building ..."
cmake --build "${output_dir}/build" ${make_opts} -- VERBOSE=1
