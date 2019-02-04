#!/usr/bin/env python3

"""Create AMICI model and data for steadystate example"""

import amici
import os
import subprocess
import numpy as np
import pandas as pd
import sys
import argparse
import h5py


def parse_cli_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        description='Create AMICI model and data for steadystate example.')

    parser.add_argument('-s', '--sbml', dest='sbml_file_name',
                        default='model_steadystate_scaled.sbml',
                        help='SBML model filename')
    parser.add_argument('-o', '--model-output-dir', dest='model_output_dir',
                        default='model_steadystate_scaled',
                        help='SBML model filename')
    '''
    parser.add_argument('-m', '--measurements', dest='measurement_file_name',
                        help='Measurement table', required=True)
    parser.add_argument('-c', '--conditions', dest='condition_file_name',
                        help='Conditions table', required=True)
    parser.add_argument('-p', dest='parameter_file_name',
                        help='Parameter table', required=True)
    '''
    args = parser.parse_args()

    return args


def create_module(sbml_model_file, model_name, model_output_dir):
    """Create Python module from SBML model"""

    sbml_importer = amici.SbmlImporter(sbml_model_file)
    sbml = sbml_importer.sbml
    observables = amici.assignmentRules2observables(
        sbml,
        filter_function=lambda variableId:
        variableId.getId().startswith('observable_') and not variableId.getId().endswith('_sigma'))

    fixed_parameters = ['k0']

    print('Observables:', observables)
    print('Fixed parameters', fixed_parameters)

    sbml_importer.sbml2amici(model_name,
                            output_dir=model_output_dir,
                            observables=observables,
                            constantParameters=fixed_parameters,
                            sigmas={'observable_x1withsigma': 'observable_x1withsigma_sigma'})

    return fixed_parameters, observables


def print_model_info(sbml_file):
    """Show model info"""

    import libsbml
    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file)
    sbml_model = sbml_doc.getModel()

    print('Model', sbml_file, f'(pwd: {os.getcwd()})')
    print('Species: ', [s.getId() for s in sbml_model.getListOfSpecies()])

    print('Parameters: ',
          [p.getId() for p in sbml_model.getListOfParameters()])

    print('\nReactions:')
    for reaction in sbml_model.getListOfReactions():
        reactants = ' + '.join(['%s %s' % (
        int(r.getStoichiometry()) if r.getStoichiometry() > 1 else '',
        r.getSpecies()) for r in reaction.getListOfReactants()])
        products = ' + '.join(['%s %s' % (
        int(r.getStoichiometry()) if r.getStoichiometry() > 1 else '',
        r.getSpecies()) for r in reaction.getListOfProducts()])
        reversible = '<' if reaction.getReversible() else ''
        print('%3s: %10s %1s->%10s\t\t[%s]' % (reaction.getId(),
                                               reactants,
                                               reversible,
                                               products,
                                               libsbml.formulaToL3String(
                                                   reaction.getKineticLaw().getMath())))


def create_data(modelModule, fixed_parameters, observables,
                sbml_file,
                model_output_dir
                ):
    """Create in silico experimental data for parameter estimation

    - Simulate time-course for four different conditions
    - Add gaussian noise according to selected sigma parameter
    - Mimic 2 experimental batches: odd-numbered condition indices and even-numbered conditions have different offset parameter
    """

    sigma_default = 0.1  # parameters are lin
    sigma_parameter = 0.2
    offset_batch_1 = 3.0
    offset_batch_2 = 4.0
    offsetted_observable_idx = 4
    sigma_parameter_observable_idx = 5
    model_offset_parameter_idx = 6
    sigma_parameter_idx = 7
    timepoints = np.logspace(-5, 1, 20)

    model = modelModule.getModel()
    default_parameters = np.array(model.getParameters())
    default_parameters[sigma_parameter_idx] = sigma_parameter
    true_parameters = default_parameters.copy()
    true_parameters = np.append(true_parameters,
                                offset_batch_2)  # add second offset parameter
    print('true_parameters:\t%s' % true_parameters)

    expectedLlh = 0.0

    # setup model & solver
    model = modelModule.getModel()
    model.setTimepoints(amici.DoubleVector(timepoints))
    model.setParameters(amici.DoubleVector(default_parameters))
    print('Default parameters:\t%s' % default_parameters)

    solver = model.getSolver()
    solver.setMaxSteps(10000)

    # generate conditon-vectors
    conditions = [np.array(model.getFixedParameters())]
    conditions.append(conditions[0] * 1.1)
    conditions.append(conditions[0] * 1.2)
    conditions.append(conditions[0] * 1.3)

    conditionDf = createConditionDataframe(fixed_parameters, conditions)

    df = pd.DataFrame(data={
        'observable': [],
        'condition': [],
        'conditionRef': [],
        'scalingParameter': [],
        'time': [],
        'measurement': [],
        'sigma': []
    })

    print()

    for icondition, condition in enumerate(conditions):
        print('Condition %d: %s' % (icondition, condition))

        # different offset for two "batches"
        batch_id = icondition % 2
        if batch_id == 0:
            simulationParameters = default_parameters
            simulationParameters[
                model_offset_parameter_idx] = offset_batch_1
        else:
            simulationParameters = default_parameters
            simulationParameters[
                model_offset_parameter_idx] = offset_batch_2

        sigmay = np.ones(shape=(model.nt(), model.nytrue)) * sigma_default
        sigmay[:,
        sigma_parameter_observable_idx] = np.nan  # observable with sigma parameter

        # simulate condition
        rdata = getReturnDataForCondition(model, solver, condition,
                                          simulationParameters, sigmay,
                                          sigma_parameter_observable_idx,
                                          sigma_parameter_idx
                                          )

        print('\tllh: ', rdata['llh'])
        print('\tsllh', rdata['sllh'])

        expectedLlh += rdata['llh']

        conditionName = 'condition-%d' % icondition

        # Append data
        for iy, observableName in enumerate(observables.keys()):
            scalingParameter = ['']
            sigma = sigmay[:, iy]

            if observableName == 'observable_x1_scaled':
                # scalingParameter = ['scaling_x1_%s' % conditionName]
                scalingParameter = ['scaling_x1_common']
            elif observableName == 'observable_x2_offsetted':
                # scalingParameter = ['offset_x2_%s' % conditionName]
                # scalingParameter = ['offset_x2_common']
                scalingParameter = ['offset_x2_batch-%d' % batch_id]
            elif observableName == 'observable_x1withsigma':
                # scalingParameter = ['observable_x1withsigma_sigma_%s' % conditionName]
                scalingParameter = ['observable_x1withsigma_sigma_common']

            df = df.append(pd.DataFrame(
                {'observable': [observableName] * model.nt(),
                 'condition': [conditionName] * model.nt(),
                 'conditionRef': [''] * model.nt(),
                 'scalingParameter': scalingParameter * model.nt(),
                 'time': np.array(model.getTimepoints()),
                 'measurement': rdata['y'][:, iy],
                 'sigma': sigma
                 }), ignore_index=True)
        print()

    print('Expected llh: ', expectedLlh)


    # write data frames to file
    measurement_file = 'example_data.tsv'
    fixed_parameter_file = 'example_data_fixed.tsv'
    hdf5File = 'example_data.h5'

    df.to_csv(measurement_file, sep='\t', index=False)
    conditionDf.to_csv(fixed_parameter_file, sep='\t', index=False)

    cmd = "bash -c 'if [[ -f example_data.h5 ]]; then cp example_data.h5 example_data.h5.bak; fi'"
    out = subprocess.check_output(cmd, shell=True)
    print(out.decode('utf-8'))

    # convert to HDF5
    cmd = ['%s/generateHDF5DataFileFromText.py' % os.path.join(
        os.path.split(os.path.abspath(__file__))[0], '..', '..', '..', 'misc'),
                          hdf5File,
                          sbml_file,
                          model_output_dir,
                          measurement_file,
                          fixed_parameter_file]
    print("Running: ", ' '.join(cmd))
    out = subprocess.run(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    print(out.stdout.decode("utf-8"))
    assert out.returncode == 0

    with h5py.File(hdf5File, 'r+') as f:
        f.require_dataset(
            '/parameters/true_parameters',
            shape=(len(true_parameters),), dtype="f8", data=true_parameters)
        f.require_dataset(
            '/parameters/true_llh',
            shape=(1,), dtype="f8", data=expectedLlh)


    hdf5FileMinibatch = 'example_data_minibatch.h5'
    from shutil import copyfile
    copyfile(hdf5File, hdf5FileMinibatch)




    # write true parameters as first starting point, an perturbed additional points
    # two times the same point to check for reproducibility
    with h5py.File(hdf5File, 'r+') as f:
        pscale = f['/parameters/pscale'][:]
        true_parameters_scaled = true_parameters.copy()
        for i, p in enumerate(pscale):
            if p == 2:
                true_parameters_scaled[i] = np.log10(true_parameters[i])

        for i in range(10):
            parameters = true_parameters_scaled
            parameters = parameters + np.random.normal(0.0, 0.2 + i * 0.1,
                                                       true_parameters.shape)
            # parameters = np.random.uniform(-3, 5, true_parameters.shape)

            # print(parameters)
            f['/optimizationOptions/randomStarts'][:, 2 * i] = parameters
            f['/optimizationOptions/randomStarts'][:, 2 * i + 1] = parameters

    copyfile(hdf5File, hdf5FileMinibatch)



def createConditionDataframe(fixed_parameters, conditions):
    """Create dataframe with fixed-parameters for each condition to simulate"""
    conditionDf = pd.DataFrame(index=fixed_parameters)
    conditionDf['ID'] = conditionDf.index
    for icondition, condition in enumerate(conditions):
        conditionDf['condition-%d' % icondition] = condition

    return conditionDf


def getReturnDataForCondition(model, solver, condition,
                              simulationParameters, sigmay,
                              sigma_parameter_observable_idx,
                              sigma_parameter_idx):
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
    print("\tSigma mean per observable:", sigmay.mean(axis=0))
    measurement[:, sigma_parameter_observable_idx] = np.random.normal(
        loc=rdata['y'][:, sigma_parameter_observable_idx],
        scale=simulationParameters[sigma_parameter_idx])
    print("\tMean abs. relative measurement error per observable:")
    print("\t",
          np.mean(np.abs((measurement - rdata['y']) / rdata['y']), axis=0))

    edata.my = amici.DoubleVector(measurement.flatten())
    edata.sigmay = amici.DoubleVector(sigmay.flatten())
    solver.setSensitivityMethod(amici.SensitivityMethod_forward)
    solver.setSensitivityOrder(amici.SensitivityOrder_first)
    model.requireSensitivitiesForAllParameters()
    rdata = amici.runAmiciSimulation(model, solver, edata)
    # return generated noisy measurents
    rdata['y'] = measurement
    return rdata


def main():
    args = parse_cli_args()

    sbml_file = args.sbml_file_name
    model_output_dir = args.model_output_dir
    script_path = os.path.split(os.path.abspath(__file__))[0]
    model_name = 'model_steadystate_scaled'

    print(f'{__file__ } running in {os.getcwd()}')
    print(f'Creating model {sbml_file}')

    # Create sbml model from scratch
    cmd = f'bash -c "{script_path}/createSteadystateExampleSBML.py > {sbml_file}"'
    print(cmd)
    out = subprocess.check_output(cmd, shell=True)
    print(out.decode('utf-8'))
    print()

    print_model_info(sbml_file)
    print()

    observables = []

    fixed_parameters, observables = create_module(sbml_file, model_name, model_output_dir)
    print()

    # load model
    sys.path.insert(0, model_output_dir)
    import model_steadystate_scaled as modelModule

    print()
    print("--- Creating data ---")
    create_data(modelModule, fixed_parameters, observables, sbml_file,
                          model_output_dir)


if __name__ == '__main__':
    main()
