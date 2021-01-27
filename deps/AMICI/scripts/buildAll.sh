#!/bin/bash
#
# Build AMICI along with dependencies and test suite
#
set -e

script_path=$(dirname "$BASH_SOURCE")
script_path=$(cd "$script_path" && pwd)

"${script_path}/buildDependencies.sh"
"${script_path}/buildAmici.sh"
