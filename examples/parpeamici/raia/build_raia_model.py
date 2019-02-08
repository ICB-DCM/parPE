#!/usr/bin/env python3

import subprocess
import os
import sys
import amici
import numpy as np
import pandas as pd
import h5py


def createConditionDataframe(indices, conditions, parameters):
    """Create dataframe with fixed-parameters for each condition to simulate"""

    conditionDf = pd.DataFrame(index=indices)
    conditionDf['ID'] = conditionDf.index
    for icondition, condition in enumerate(conditions):
        conditionDf['condition-%d' % icondition] = condition

    return conditionDf


def getReturnDataForCondition(model, solver, condition, simulationParameters,
                              sigmay):
    model.setParameters(amici.DoubleVector(simulationParameters))

    # simulate without measurements
    edata = amici.ExpData(model.get())
    edata.fixedParameters = amici.DoubleVector(condition)
    edata.my = amici.DoubleVector(
        np.full(shape=model.nt() * model.nytrue, fill_value=np.nan))
    rdata = amici.runAmiciSimulation(model, solver, edata)
    # fixedParametersPreequilibration =

    # confirm gradient is 0 for real measurements and save expected llh
    measurement = rdata['y']
    measurement = np.random.normal(loc=rdata['y'], scale=sigmay)
    measurement = np.random.normal(
        loc=rdata['y'],
        scale=sigmay[0][0])
    # print(measurement)

    edata.my = amici.DoubleVector(measurement.flatten())
    edata.sigmay = amici.DoubleVector(sigmay.flatten())
    model.requireSensitivitiesForAllParameters()
    rdata = amici.runAmiciSimulation(model, solver, edata)
    # return generated noisy measurents
    rdata['y'] = measurement
    return rdata


def main():
    script_dir = os.path.splitext(os.path.abspath(__file__))[0]
    script_dir = os.path.split(script_dir)[0]
    sbml_file = os.path.join(script_dir, 'raia.xml')
    model_name = 'raia'
    model_output_dir = os.path.join(os.getcwd())

    print("Importing model from", sbml_file)
    print("Generating files in", model_output_dir)

    # set observables and constants
    observables_list = ['observable_RecSurf', 'observable_IL13_cell',
                        'observable_pIL4Ra', 'observable_pJAK2',
						'observable_SOCS3mRNA', 'observable_CD274mRNA',
                        'observable_SOCS3', 'observable_pSTAT5']
    fixed_parameters = ['il13_level']

    # wrap the model
    sbmlImporter = amici.SbmlImporter(sbml_file)

    observables = amici.assignmentRules2observables( \
        sbmlImporter.sbml,  # the libsbml model object
        filter_function=lambda variable: variable.getId() in observables_list)

    sbmlImporter.sbml2amici(model_name,
                            output_dir=model_output_dir,
                            observables=observables,
                            constantParameters=fixed_parameters)

    # load model
    sys.path.insert(0, model_output_dir)
    import raia as modelModule

    model = modelModule.getModel()
    solver = model.getSolver()
    rdata = amici.runAmiciSimulation(model, solver)

    expectedLlh = 0.0
    sigma_default = 0.1  # parameters are lin
    timepoints = np.linspace(0, 100, 5)

    model = modelModule.getModel()
    true_parameters = np.array(model.getParameters())

    # setup model & solver
    model = modelModule.getModel()
    model.setTimepoints(amici.DoubleVector(timepoints))
    model.setParameters(amici.DoubleVector(true_parameters))

    solver = model.getSolver()
    solver.setMaxSteps(10000)

    # generate condition-vectors
    nConditions = 10
    init_conditions = [np.array([0.])]
    for i_cond in range(1, nConditions):
    	if (i_cond < 5):
    		init_conditions.append(2. * np.array([float(i_cond)]))
    	else:
    		init_conditions.append(15. * np.array([float(i_cond - 4)]))

    indices = fixed_parameters
    conditionDf = createConditionDataframe(indices, init_conditions,
                                           true_parameters)

    df = pd.DataFrame(data={
        'observable': [],
        'condition': [],
        'conditionRef': [],
        'time': [],
        'measurement': [],
    })

    for icondition, condition in enumerate(init_conditions):

        simulationParameters = true_parameters
        sigmay = np.ones(shape=(model.nt(), model.nytrue)) * sigma_default

        # simulate condition
        rdata = getReturnDataForCondition(model, solver, condition,
                                          simulationParameters, sigmay)

        expectedLlh += rdata['llh']

        conditionName = 'condition-%d' % icondition

        # Append data
        for iy, observableName in enumerate(observables.keys()):
            scalingParameter = ['']
            sigma = sigmay[:, iy]

            df = df.append(pd.DataFrame(
                {'observable': [observableName] * model.nt(),
                 'condition': [conditionName] * model.nt(),
                 'conditionRef': [''] * model.nt(),
                 'scalingParameter': scalingParameter * model.nt(),
                 'time': np.array(model.getTimepoints()),
                 'measurement': rdata['y'][:, iy],
                 'sigma': sigma
                 }), ignore_index=True)

    # write data frames to file
    measurement_file = 'example_data.tsv'
    fixed_parameter_file = 'example_data_fixed.tsv'
    hdf5File = 'example_data.h5'

    df.to_csv(measurement_file, sep='\t', index=False)
    conditionDf.to_csv(fixed_parameter_file, sep='\t', index=False)

    # convert to HDF5
    subprocess.call(["/bin/bash", "-c",
                     "if [[ -f example_data.h5 ]]; then cp example_data.h5 example_data.h5.bak; fi"])

    out = subprocess.run(['%s/generateHDF5DataFileFromText.py' % os.path.join(script_dir, '..', '..', '..', 'misc'),
                          hdf5File,
                          sbml_file,
                          model_output_dir,
                          measurement_file,
                          fixed_parameter_file], stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    print(out.stdout.decode("utf-8"))

    # changes some solver options in the hdf5 file
    f_hdf5 =  h5py.File(hdf5File, 'r+')
    amici_options = f_hdf5['amiciOptions']
    amici_options.attrs['atol'] = 1.0e-8
    amici_options.attrs['rtol'] = 1.0e-6
    amici_options.attrs['sensi_meth'] = 1

if __name__ == '__main__':
    main()
