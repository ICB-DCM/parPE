# *parPE* steadystate example

This directory contains file to generate a simple parameter estimation problem
and some Jupyter notebooks to demonstrate basic parPE usage.

* `parpeExampleSteadystateBasic.ipynb` describes the parameter estimation
  problem and demonstrates basic parPE usage

* `parpeExampleSteadystateHierarchical.ipynb` demonstrates the hierarchical
  optimization approach implemented in parPE

* `parpeExampleSteadystateMinibatch.ipynb` demonstrates the use of mini-batch
  optimization implemented in parPE

The example problem is generated by

* `createSteadystateExampleSBML.py`: Create SBML model from scratch
* `create_steadystate_amici_model.py`: Create example artificial measurements
  and other PEtab files
