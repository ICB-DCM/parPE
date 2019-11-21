#!/usr/bin/env bash
# Run a given command in the given Python virtual environment

set -e

# invalid options:
if [[ $# -lt 2 ]]; then
    echo "Run a given command in the given Python virtual environment"
    echo "USAGE: $(basename "$0") VENV_DIR COMMAND [ARG ...]"
    exit 1;
fi

source $1/bin/activate
"${@:2}"
