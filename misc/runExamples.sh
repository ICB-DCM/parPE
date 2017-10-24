#!/bin/bash

set -e
set -x

mpiexec -n 5 loadbalancer/examples/loadbalancer/example_loadbalancer

amici/examples/steadystate/example_steadystate

amici/examples/steadystate/example_steadystate_parallel

mpiexec -n 5 amici/examples/steadystate/example_steadystate_parallel

mpiexec -n 5 amici/examples/steadystate/example_steadystate_multi -o deleteme ../amici/examples/steadystate/data.h5

