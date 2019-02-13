#!/bin/bash
SCRIPT_DIR="`dirname \"$0\"`"
SCRIPT_DIR="`( cd \"$SCRIPT_DIR\" && pwd )`"

set -e

BUILD_DIR=$1
VENV=${BUILD_DIR}/venv

source $VENV/bin/activate
cd $SCRIPT_DIR

cd ${BUILD_DIR}/examples/parpeamici/lucarelli/lucarelli_12-prefix/src/lucarelli_12
$SCRIPT_DIR/build_lucarelli_model.py
