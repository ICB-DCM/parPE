# Description of input format for misc/generateHDF5DataFileFromText.py

This document describes the input data for 
`misc/generateHDF5DataFileFromText.py`, which generates HDF5 files for 
parameter estimation. 

For usage of the conversion scripts, run `misc/generateHDF5DataFileFromText.py`
without any arguments.

Required inputs are:
- the SBML model
- the folder of the corresponding python-AMICI generated model
- a text file with condition specific parameters
- a text file the measurements to be used for model training
- a text file describing optimization parameters

The SBML model and text files are expected to adhere to the 
[PEtab format](https://github.com/ICB-DCM/PEtab/blob/master/doc/documentation_data_format.md)
