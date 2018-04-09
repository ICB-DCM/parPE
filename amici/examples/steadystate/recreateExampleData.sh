#!/bin/bash
SCRIPT_DIR="`dirname \"$0\"`"
SCRIPT_DIR="`( cd \"$SCRIPT_DIR\" && pwd )`"

# create SBML
SBML_MODEL_NAME="model_steadystate_scaled.sbml"
./createSteadystateExampleSBML.py > $SBML_MODEL_NAME

# create Amici model
AMICI_PATH="${SCRIPT_DIR}/../../../deps/AMICI/"
COMPILED_MODEL_NAME="model_steadystate_scaled"
SYMS_FUNCNAME="model_steadystate_scaled_syms"
matlab -nojvm -nosplash -nodesktop -r "addpath('$AMICI_PATH'); addpath('$SCRIPT_PATH'); installAMICI; SBML2AMICI('$SBML_MODEL_NAME', '$COMPILED_MODEL_NAME'); amiwrap('$COMPILED_MODEL_NAME','$SYMS_FUNCNAME'); quit;"
rm "${COMPILED_MODEL_NAME}"*nom.mat

# copy c++ model
MODEL_SRC="${AMICI_PATH}/models/${COMPILED_MODEL_NAME}"
MODEL_DEST="${SCRIPT_DIR}/${COMPILED_MODEL_NAME}"
mkdir -p "${MODEL_DEST}"
cp ${MODEL_SRC}/*.cpp ${MODEL_DEST} 
cp ${MODEL_SRC}/*.h ${MODEL_DEST}