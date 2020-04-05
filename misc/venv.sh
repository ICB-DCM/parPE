#!/usr/bin/env bash
#
# Setup virtual environment for building/testing parPE
#

script_path=$(dirname "$BASH_SOURCE")
parpe_root=$(cd "${script_path}/.." && pwd)

build_dir=${1:-"${parpe_root}/build"}
venv_dir="${build_dir}/venv"


# Save the time for installing AMICI if the venv already exists
# NOTE: Must remove folder if AMICI is updated
if [[ ! -d "${venv_dir}" ]]; then
    # create venv
    python3 -m venv "${venv_dir}" 2>/dev/null
    # in case this fails (usually due to missing ensurepip), try getting pip
    # manually
    if [[ $? ]]; then
        set -e
        python3 -m venv "${venv_dir}" --clear --without-pip
        source "${venv_dir}/bin/activate"
        curl https://bootstrap.pypa.io/get-pip.py -o "${build_dir}/get-pip.py"
        python3 "${build_dir}/get-pip.py"
    else
        set -e
        source "${venv_dir}/bin/activate"
    fi

    pip3 install wheel pytest

    # install AMICI
    cd "${parpe_root}/deps/AMICI/python/sdist"
    pip3 install -e .

    # install parPE
    cd "${parpe_root}/python"
    pip3 install -e .
    #pip3 install https://github.com/ICB-DCM/PEtab/archive/develop.zip
fi
