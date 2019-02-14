#!/usr/bin/env python3

"""Create AMICI model and synthetic data for steadystate example"""

import amici
import os
import subprocess
import numpy as np
import pandas as pd
import sys
import argparse
import h5py
import petab
import libsbml
import importlib


def parse_cli_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        description='Create AMICI model and data for steadystate example.')

    parser.add_argument('-s', '--sbml', dest='sbml_file_name',
                        default='model_steadystate_scaled.sbml',
                        help='SBML model filename')

    parser.add_argument('-o', '--model-output-dir', dest='model_output_dir',
                        default='model_steadystate_scaled',
                        help='Name of the model directory to be created')

    parser.add_argument('-m', '--measurements', dest='measurement_file_name',
                        default='example_data.tsv',
                        help='Name of measurement table to be generated')

    parser.add_argument('-c', '--conditions', dest='condition_file_name',
                        default='example_data_fixed.tsv',
                        help='Name of conditions table to generate')

    parser.add_argument('-f', dest='hdf5_file_name',
                        default='example_data.h5',
                        help='Name of HDF5 file to generate')

    parser.add_argument('-p', dest='parameter_file_name',
                        default='example_data_parameter.tsv',
                        help='Name of parameter table to be created')

    args = parser.parse_args()

    return args


def create_module(sbml_model_file, model_name, model_output_dir):
    """Create Python module from SBML model

    Arguments:
        sbml_model_file: SBML file
        model_name: Name of the model
        model_output_dir: Directory for model code

    Returns:
        list of parameters,
        dictionary with observables (id: formula)
    """

    sbml_importer = amici.SbmlImporter(sbml_model_file)
    sbml_model = sbml_importer.sbml

    observables = petab.get_observables(sbml_model, remove=True)
    sigmas = petab.get_sigmas(sbml_model, remove=True)
    fixed_parameters = ['k0']

    sbml_importer.sbml2amici(
        modelName=model_name,
        output_dir=model_output_dir,
        observables=observables,
        constantParameters=fixed_parameters,
        sigmas=sigmas
    )

    return fixed_parameters, observables


def print_model_info(sbml_file):
    """Show model info"""

    sbml_reader = libsbml.SBMLReader()
    sbml_doc = sbml_reader.readSBML(sbml_file)
    sbml_model = sbml_doc.getModel()

    print('Model', sbml_file, f'(pwd: {os.getcwd()})')
    print('Species: ', [s.getId() for s in sbml_model.getListOfSpecies()])

    print('Parameters: ',
          [p.getId() for p in sbml_model.getListOfParameters()])

    print('\nReactions:')
    for reaction in sbml_model.getListOfReactions():
        reactants = ' + '.join(
            ['%s %s' % (int(r.getStoichiometry())
                        if r.getStoichiometry() > 1
                        else '',
                        r.getSpecies()) for r in reaction.getListOfReactants()]
        )
        products = ' + '.join(
            ['%s %s' % (int(r.getStoichiometry())
                        if r.getStoichiometry() > 1
                        else '',
                        r.getSpecies()) for r in reaction.getListOfProducts()]
        )
        reversible = '<' if reaction.getReversible() else ''
        print('%3s: %10s %1s->%10s\t\t[%s]' %
              (reaction.getId(),
               reactants,
               reversible,
               products,
               libsbml.formulaToL3String(reaction.getKineticLaw().getMath())))


def create_data_tables(model, fixed_parameters,
        measurement_file, fixed_parameter_file):
    """Create synthetic data for parameter estimation

    - Simulate time-course for four different conditions
    - Add gaussian noise according to selected sigma parameter
    - Mimic 2 experimental batches: odd-numbered condition indices
    and even-numbered conditions have different offset parameter

    Arguments:

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

    # set true parameters
    default_parameters = np.array(model.getParameters())
    default_parameters[sigma_parameter_idx] = sigma_parameter
    print('Default model parameters:')
    for p, val in zip(model.getParameterIds(), model.getParameters()):
        print(f'\t{p}: {val}')
    print()

    true_parameters = {pid: val for pid, val in zip(model.getParameterIds(), default_parameters)}
    true_parameters['scaling_x1_common'] = \
        true_parameters['observableParameter1_x1_scaled']
    # extend to optimization parameter vector: add second offset parameter
    true_parameters['offset_x2_batch-0'] = offset_batch_1
    true_parameters['offset_x2_batch-1'] = offset_batch_2
    true_parameters['x1withsigma_sigma'] = sigma_parameter

    print('True parameters:\t%s' % true_parameters)

    # setup model
    model.setTimepoints(timepoints)
    model.setParameters(default_parameters)

    # setup solver
    solver = model.getSolver()
    solver.setMaxSteps(10000)

    # generate four condition-vectors
    condition_df = petab.create_condition_df(fixed_parameters)
    condition_df.loc['condition-0'] = model.getFixedParameters()
    condition_df.loc['condition-1'] = np.array(model.getFixedParameters()) * 1.1
    condition_df.loc['condition-2'] = np.array(model.getFixedParameters()) * 1.2
    condition_df.loc['condition-3'] = np.array(model.getFixedParameters()) * 1.3

    print(condition_df)
    measurement_df = petab.create_measurement_df()

    print()

    # set sigmas
    sigmay = np.ones(shape=(model.nt(), model.nytrue)) * sigma_default
    # observable with sigma parameter
    sigmay[:, sigma_parameter_observable_idx] = np.nan

    # llh for noisy simulated data with true parameters
    expected_llh = 0.0

    for condition_idx, condition_id in enumerate(condition_df.index.values):
        condition_parameters = condition_df.loc[condition_id, :]
        print(f'Condition {condition_id}:  {condition_parameters}')

        # different offset for two "batches"
        batch_id = condition_idx % 2
        model_parameters = default_parameters
        if batch_id == 0:
            model_parameters[
                model_offset_parameter_idx] = offset_batch_1
        else:
            model_parameters[
                model_offset_parameter_idx] = offset_batch_2

        # simulate condition
        rdata = getReturnDataForCondition(
            model, solver, condition_parameters,
            model_parameters, sigmay,
            sigma_parameter_observable_idx,
            sigma_parameter_idx
        )

        print('\tllh: ', rdata['llh'])
        print('\tsllh', rdata['sllh'])

        expected_llh += rdata['llh']

        measurement_df = append_measurements_for_condition(
            model, measurement_df, sigmay, condition_id, batch_id, rdata)
        print()

    print('Expected llh: ', expected_llh)

    # write data frames to file
    measurement_df.to_csv(measurement_file, sep='\t', index=False)
    condition_df.to_csv(fixed_parameter_file, sep='\t', index=True)

    return true_parameters, expected_llh


def getReturnDataForCondition(model, solver, fixed_parameters,
                              dynamic_parameters, sigmay,
                              sigma_parameter_observable_idx,
                              sigma_parameter_idx):

    model.setParameters(amici.DoubleVector(dynamic_parameters))

    # Simulate without measurements for noise-free trajectories
    edata = amici.ExpData(model.get())
    edata.fixedParameters = fixed_parameters
    edata.my = np.full(shape=model.nt() * model.nytrue, fill_value=np.nan)
    rdata = amici.runAmiciSimulation(model, solver, edata)
    # fixedParametersPreequilibration =

    # synthetic_data = rdata['y'] # noise-free

    # Add noise to simulation
    synthetic_data = np.random.normal(loc=rdata['y'], scale=sigmay)
    print("\tSigma mean per observable:", sigmay.mean(axis=0))
    # Apply correct sigma parameter
    synthetic_data[:, sigma_parameter_observable_idx] = \
        np.random.normal(loc=rdata['y'][:, sigma_parameter_observable_idx],
                         scale=dynamic_parameters[sigma_parameter_idx])

    print("\tMean abs. relative measurement error per observable:")
    print("\t",
          np.mean(np.abs((synthetic_data - rdata['y']) / rdata['y']), axis=0))

    # Use synthetic data to get expected llh
    edata.my = synthetic_data.flatten()
    edata.sigmay = sigmay.flatten()
    solver.setSensitivityMethod(amici.SensitivityMethod_forward)
    solver.setSensitivityOrder(amici.SensitivityOrder_first)
    model.requireSensitivitiesForAllParameters()
    rdata = amici.runAmiciSimulation(model, solver, edata)

    # TODO: confirm gradient is 0 for real measurements and save expected llh

    # return generated noisy measurements
    rdata['y'] = synthetic_data

    return rdata


def append_measurements_for_condition(
        model, measurement_df, sigmay, condition_id, batch_id, rdata):
    # Append data to measurement table
    for iy, observable_id in enumerate(model.getObservableIds()):
        scaling_parameter = ['']
        noise_parameter = sigmay[:, iy]

        if observable_id == 'observable_x1_scaled':
            scaling_parameter = ['scaling_x1_common']
        elif observable_id == 'observable_x2_offsetted':
            scaling_parameter = ['offset_x2_batch-%d' % batch_id]
        elif observable_id == 'observable_x1withsigma':
            noise_parameter = ['x1withsigma_sigma'] * model.nt()

        measurement_df = measurement_df.append(pd.DataFrame(
            {'observableId': [observable_id[len('observable_'):]] * model.nt(),
             'simulationConditionId': [condition_id] * model.nt(),
             'time': np.array(model.getTimepoints()),
             'measurement': rdata['y'][:, iy],
             'observableParameters': scaling_parameter * model.nt(),
             'noiseParameters': noise_parameter
             }), ignore_index=True, sort=False)
    return measurement_df


def generate_hdf5_file(sbml_file_name, model_output_dir, measurement_file_name,
                       condition_file_name, hdf5_file_name,
                       parameter_file_name, model_name):

    cmd = f"bash -c 'if [[ -f {hdf5_file_name} ]]; then "\
        f"cp {hdf5_file_name} {hdf5_file_name}.bak; fi'"
    out = subprocess.check_output(cmd, shell=True)
    print(out.decode('utf-8'))

    # convert to HDF5
    script_file = os.path.join(os.path.split(os.path.abspath(__file__))[0],
                               '..', '..', '..', 'misc',
                               'generateHDF5DataFileFromText.py')
    cmd = [script_file,
           '-o', hdf5_file_name,
           '-s', sbml_file_name,
           '-d', model_output_dir,
           '-m', measurement_file_name,
           '-c', condition_file_name,
           '-n', model_name,
           '-p', parameter_file_name]
    print("Running: ", ' '.join(cmd))
    out = subprocess.run(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    print(out.stdout.decode("utf-8"))
    assert out.returncode == 0


def save_expected_results(hdf5_file_name, true_parameters_dict, expected_llh):
    # write true parameters
    with h5py.File(hdf5_file_name, 'r+') as f:
        true_parameters = [true_parameters_dict[p] for p in
                           f['/parameters/parameterNames']]
        f.require_dataset(
            '/parameters/true_parameters',
            shape=(len(true_parameters),), dtype="f8", data=true_parameters)
        f.require_dataset(
            '/parameters/true_llh',
            shape=(1,), dtype="f8", data=expected_llh)


def write_starting_points(hdf5_file_name, true_parameters):
    """Write true parameters as first starting point, an perturbed additional
    points, two times the same point to check for reproducibility"""
    with h5py.File(hdf5_file_name, 'r+') as f:
        pscale = f['/parameters/pscale'][:]
        true_parameters_scaled = np.array([true_parameters[p] for p in
                           f['/parameters/parameterNames']])
        for i, p in enumerate(pscale):
            if p == amici.ParameterScaling_log10:
                true_parameters_scaled[i] = np.log10(true_parameters_scaled[i])

        for i in range(10):
            parameters = true_parameters_scaled
            parameters += np.random.normal(0.0, 0.2 + i * 0.1,
                                           true_parameters_scaled.shape)
            # parameters = np.random.uniform(-3, 5, true_parameters.shape)

            # print(parameters)
            f['/optimizationOptions/randomStarts'][:, 2 * i] = parameters
            f['/optimizationOptions/randomStarts'][:, 2 * i + 1] = parameters


def create_parameter_table(sbml_file, condition_file, measurement_file,
                           parameter_file, nominal_parameters):
    """Create PEtab parameter table"""

    problem = petab.Problem(sbml_file, condition_file, measurement_file)
    df = problem.create_parameter_df(lower_bound=-3,
                                     upper_bound=5)
    # TODO: move to peTAB
    df['hierarchicalOptimization'] = 0
    df.loc['scaling_x1_common', 'hierarchicalOptimization'] = 1
    df.loc['offset_x2_batch-0', 'hierarchicalOptimization'] = 1
    df.loc['offset_x2_batch-1', 'hierarchicalOptimization'] = 1
    df.loc['x1withsigma_sigma', 'hierarchicalOptimization'] = 1
    #df.parameterScale = 'lin'
    #df.estimate = 0

    print(nominal_parameters)

    for pid, val in nominal_parameters.items():
        if pid in df.index:
            df.loc[pid, 'nominalValue'] = val
            df.loc[pid, 'parameterScale'] = 'log10'
            df.loc[pid, 'estimate'] = 1
        elif pid.startswith('noiseParameter') \
            or pid.startswith('observableParameter'):
            continue
        else:
            print("extra parameter", pid, val)
    df.to_csv(parameter_file, sep="\t", index=True)


def main():
    args = parse_cli_args()

    script_path = os.path.split(os.path.abspath(__file__))[0]
    model_name = 'model_steadystate_scaled'

    print(f'{__file__} running in {os.getcwd()}')
    print(f'Processing model {args.sbml_file_name}')

    # Create sbml model from scratch
    cmd = f'bash -c "{script_path}/createSteadystateExampleSBML.py > {args.sbml_file_name}"'
    print(cmd)
    out = subprocess.check_output(cmd, shell=True)
    print(out.decode('utf-8'))
    print()

    print_model_info(args.sbml_file_name)
    print()

    fixed_parameters, observables = create_module(
        args.sbml_file_name, model_name, args.model_output_dir)

    print('Observables:', observables)
    print('Fixed parameters', fixed_parameters)
    print()

    # load model
    sys.path.insert(0, args.model_output_dir)
    model_module = importlib.import_module(model_name)

    print()
    print("--- Creating data ---")
    true_parameters, expected_llh = create_data_tables(
        model=model_module.getModel(),
        measurement_file=args.measurement_file_name,
        fixed_parameter_file=args.condition_file_name,
        fixed_parameters=fixed_parameters
    )

    create_parameter_table(args.sbml_file_name,
                           args.condition_file_name,
                           args.measurement_file_name,
                           args.parameter_file_name,
                           true_parameters
                           )

    # check for valid PEtab
    pp = petab.Problem(args.sbml_file_name,
                       args.condition_file_name,
                       args.measurement_file_name,
                       args.parameter_file_name)
    petab.lint_problem(pp)

    generate_hdf5_file(
        sbml_file_name=args.sbml_file_name,
        model_output_dir=args.model_output_dir,
        measurement_file_name=args.measurement_file_name,
        condition_file_name=args.condition_file_name,
        hdf5_file_name=args.hdf5_file_name,
        parameter_file_name=args.parameter_file_name,
        model_name=model_name
    )

    save_expected_results(args.hdf5_file_name, true_parameters, expected_llh)

    write_starting_points(args.hdf5_file_name, true_parameters)

    # TODO
    #hdf5FileMinibatch = 'example_data_minibatch.h5'
    #from shutil import copyfile
    #copyfile(hdf5_file_name, hdf5FileMinibatch)


if __name__ == '__main__':
    main()
