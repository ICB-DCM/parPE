#!/bin/bash

set -e
set -x

mpirun -np 5 loadbalancer/examples/loadbalancer/example_loadbalancer

mpirun -np 5 amici/examples/steadystate/example_steadystate

mpirun -np 5 amici/examples/steadystate/example_steadystate_parallel

mpirun -np 5 amici/examples/steadystate/example_steadystate_multi
