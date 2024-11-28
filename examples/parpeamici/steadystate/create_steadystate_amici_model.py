#!/usr/bin/env python3

"""Create AMICI model and synthetic data for steadystate example"""

import argparse
import importlib
import os
import subprocess
import sys

import amici
import h5py
import libsbml
import numpy as np
import pandas as pd
import petab.v1 as petab
import petab.C as ptc
import yaml
from petab.v1.models.sbml_model import SbmlModel


def parse_cli_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        description='Create AMICI model PEtab files for steadystate example.')

    parser.add_argument('-o', '--model-output-dir', dest='model_output_dir',
                        default='model_steadystate_scaled',
                        help='Name of the AMICI model directory to be created')

    parser.add_argument('-f', dest='hdf5_file_name',
                        default='example_data.h5',
                        help='Name of HDF5 file to generate')

    parser.add_argument('-p', dest='petab_dir',
                        default='steadystate_petab',
                        help='Directory to write PEtab files to')

    args = parser.parse_args()

    return args


def create_module(petab_problem, model_name, model_output_dir):
    """Create AMICI Python module from SBML model

    Arguments:
        sbml_model_file: SBML file
        model_name: Name of the model
        model_output_dir: Directory for model code
    """

    from amici.petab.petab_import import import_model
    import_model(petab_problem=petab_problem,
                 model_name=model_name, model_output_dir=model_output_dir,
                 verbose=True, )


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


def create_data_tables(model, condition_df):
    """Create synthetic data for parameter estimation

    - Simulate time-course for four different conditions
    - Add gaussian noise according to selected sigma parameter
    - Mimic 2 experimental batches: odd-numbered condition indices
    and even-numbered conditions have different offset parameter

    Arguments:

    """
    timepoints = np.logspace(-5, 1, 20)
    sigma_default = 0.1  # parameters are lin
    sigma_parameter = 0.2
    offset_batch_1 = 3.0
    offset_batch_2 = 4.0

    parameter_ids = list(model.getParameterIds())
    observable_ids = list(model.getObservableIds())

    sigma_parameter_observable_idx = \
        observable_ids.index('obs_x1withsigma')
    model_offset_parameter_idx = \
        parameter_ids.index('observableParameter1_obs_x2_offsetted')
    sigma_parameter_idx = \
        parameter_ids.index('noiseParameter1_obs_x1withsigma')

    # set true parameters
    default_parameters = np.array(model.getParameters())
    default_parameters[sigma_parameter_idx] = sigma_parameter
    print('Default model parameters:')
    for p, val in zip(model.getParameterIds(), model.getParameters()):
        print(f'\t{p}: {val}')
    print()

    true_parameters = {pid: val for pid, val in zip(model.getParameterIds(),
                                                    default_parameters)}
    # output parameters don't have default values from SBML mdoel
    true_parameters['observableParameter1_obs_x1_scaled'] = 2.0
    true_parameters['noiseParameter1_obs_x1withsigma'] = 0.2
    true_parameters['noiseParameter1_obs_x1_scaled'] = 0.2
    true_parameters['noiseParameter1_obs_x2'] = 0.2
    true_parameters['noiseParameter1_obs_x1'] = 0.2
    true_parameters['noiseParameter1_obs_x2_offsetted'] = 0.2
    true_parameters['noiseParameter1_obs_x3'] = 0.2
    true_parameters['observableParameter1_obs_x2_offsetted'] = 3.0

    true_parameters['scaling_x1_common'] = \
        true_parameters['observableParameter1_obs_x1_scaled']
    # extend to optimization parameter vector: add second offset parameter
    true_parameters['offset_x2_batch_0'] = offset_batch_1
    true_parameters['offset_x2_batch_1'] = offset_batch_2
    true_parameters['x1withsigma_sigma'] = sigma_parameter

    print('True parameters:\t%s' % true_parameters)

    # setup model
    model.setTimepoints(timepoints)
    model.setParameters(default_parameters)

    # setup solver
    solver = model.getSolver()
    solver.setMaxSteps(10000)

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
        print(f'Condition {condition_idx} "{condition_id}":  {condition_parameters}')

        # different offset for two "batches"
        batch_id = condition_idx % 2
        model_parameters = default_parameters[:]
        if batch_id == 0:
            model_parameters[
                model_offset_parameter_idx] = offset_batch_1
        else:
            model_parameters[
                model_offset_parameter_idx] = offset_batch_2

        print('Model parameters:', model_parameters)

        # simulate condition
        rdata = get_return_data_for_condition(
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

    return measurement_df, true_parameters, expected_llh


def get_return_data_for_condition(model, solver, fixed_parameters,
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
    # due to noise, there may be negative measurements. we don't want them.
    synthetic_data = np.abs(synthetic_data)
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
    rdata._swigptr.y = amici.DoubleVector(synthetic_data.flatten())

    return rdata


def append_measurements_for_condition(
        model, measurement_df, sigmay, condition_id, batch_id, rdata):
    # Append data to measurement table
    for iy, observable_id in enumerate(model.getObservableIds()):
        scaling_parameter = ['']
        noise_parameter = sigmay[:, iy]

        if observable_id == 'obs_x1_scaled':
            scaling_parameter = ['scaling_x1_common']
        elif observable_id == 'obs_x2_offsetted':
            scaling_parameter = ['offset_x2_batch_%d' % batch_id]
        elif observable_id == 'obs_x1withsigma':
            noise_parameter = ['x1withsigma_sigma'] * model.nt()

        measurement_df = pd.concat([measurement_df, pd.DataFrame(
            {ptc.OBSERVABLE_ID: [observable_id] * model.nt(),
             ptc.SIMULATION_CONDITION_ID: [condition_id] * model.nt(),
             ptc.TIME: np.array(model.getTimepoints()),
             ptc.MEASUREMENT: rdata['y'][:, iy],
             ptc.OBSERVABLE_PARAMETERS: scaling_parameter * model.nt(),
             ptc.NOISE_PARAMETERS: noise_parameter
             })], ignore_index=True, sort=False)
    return measurement_df


def generate_hdf5_file(yaml_file, model_output_dir,
                       hdf5_file_name, model_name):

    cmd = f"bash -c 'if [[ -f {hdf5_file_name} ]]; then "\
        f"cp {hdf5_file_name} {hdf5_file_name}.bak; fi'"
    out = subprocess.check_output(cmd, shell=True)
    print(out.decode('utf-8'))

    # convert to HDF5
    script_file = 'parpe_petab_to_hdf5'
    cmd = [script_file,
           '-o', hdf5_file_name,
           '-d', model_output_dir,
           '-n', model_name,
           '-y', yaml_file]
    print("Running: ", ' '.join(cmd))
    out = subprocess.run(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    print(out.stdout.decode("utf-8"))
    assert out.returncode == 0


def save_expected_results(hdf5_file_name, true_parameters_dict, expected_llh):
    # write true parameters
    with h5py.File(hdf5_file_name, 'r+') as f:
        true_parameters = [true_parameters_dict[p] for p in
                           f['/parameters/parameterNames'].asstr()]
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
        pscale = f['/parameters/pscaleOptimization'][:]
        true_parameters_scaled = np.array([true_parameters[p] for p in
                           f['/parameters/parameterNames'].asstr()])
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


def create_parameter_table(problem: petab.Problem,
                           nominal_parameters):
    """Create PEtab parameter table"""

    df = petab.create_parameter_df(
        condition_df=problem.condition_df,
        observable_df=problem.observable_df,
        measurement_df=problem.measurement_df,
        model=SbmlModel(problem.sbml_model),
        include_optional=True, lower_bound=1e-3, upper_bound=1e5)

    df['hierarchicalOptimization'] = 0
    df.loc['scaling_x1_common', 'hierarchicalOptimization'] = 1
    df.loc['offset_x2_batch_0', 'hierarchicalOptimization'] = 1
    df.loc['offset_x2_batch_1', 'hierarchicalOptimization'] = 1
    df.loc['x1withsigma_sigma', 'hierarchicalOptimization'] = 1

    for pid, val in nominal_parameters.items():
        if pid in df.index:
            df.loc[pid, ptc.NOMINAL_VALUE] = val
            df.loc[pid, ptc.PARAMETER_SCALE] = ptc.LOG10
            df.loc[pid, ptc.ESTIMATE] = 1
        elif pid.startswith('noiseParameter') \
                or pid.startswith('observableParameter'):
            continue
        else:
            print("extra parameter", pid, val)

    # offsets can be negative: adapt scaling and bounds:
    offsets = df.index.str.startswith('offset_')
    df.loc[offsets, ptc.PARAMETER_SCALE] = ptc.LIN

    problem.parameter_df = df


def create_test_data(measurement_file_name, parameter_file_name, yaml_config,
                     yaml_file_name_test, model_output_dir, model_name,
                     hdf5_file_name):
    """Create some synthetic data to emulate a test set"""

    test_measurement_file_name = \
        "-testset".join(os.path.splitext(measurement_file_name))
    test_parameter_file_name = \
        "-testset".join(os.path.splitext(parameter_file_name))

    # measurements
    df = petab.get_measurement_df(measurement_file_name)
    df.loc[df.observableParameters == 'scaling_x1_common', 'measurement'] = \
        df.loc[df.observableParameters == 'scaling_x1_common', 'measurement'] \
        * 2.0
    df.loc[~df.observableParameters.isnull(), 'observableParameters'] = \
        df.loc[~df.observableParameters.isnull(), 'observableParameters'] \
        + "_test"

    petab.write_measurement_df(df, test_measurement_file_name)

    # parameters
    df = petab.get_parameter_df(parameter_file_name)
    df.rename(index={'scaling_x1_common' : 'scaling_x1_common_test',
                     'offset_x2_batch_0': 'offset_x2_batch_0_test',
                     'offset_x2_batch_1': 'offset_x2_batch_1_test'},
              inplace=True)
    petab.write_parameter_df(df, test_parameter_file_name)

    # yaml
    yaml_config[ptc.PARAMETER_FILE] = test_parameter_file_name
    yaml_config[ptc.PROBLEMS][0][ptc.MEASUREMENT_FILES][0] = \
        test_measurement_file_name
    with open(yaml_file_name_test, 'w') as outfile:
        yaml.dump(yaml_config, outfile, default_flow_style=False)

    generate_hdf5_file(
        yaml_file=yaml_file_name_test,
        model_output_dir=model_output_dir,
        hdf5_file_name="-testset".join(os.path.splitext(hdf5_file_name)),
        model_name=model_name
    )


def main():
    args = parse_cli_args()

    script_path = os.path.split(os.path.abspath(__file__))[0]
    model_name = 'model_steadystate_scaled'
    sbml_file_name = "model_steadystate_scaled.sbml"
    measurement_file_name = 'example_data.tsv'
    condition_file_name = 'example_data_fixed.tsv'
    parameter_file_name = 'example_data_parameter.tsv'
    observable_file_name = 'model_steadystate_observables.tsv'
    yaml_file_name = 'model_steadystate.yaml'
    yaml_file_name_test = 'model_steadystate_test.yaml'

    print(f'{__file__} running in {os.getcwd()}')
    print(f'Processing model {sbml_file_name}')

    # Create sbml model from scratch
    cmd = f'bash -c "{script_path}/createSteadystateExampleSBML.py '\
          f'> {sbml_file_name}"'
    print(cmd)
    out = subprocess.check_output(cmd, shell=True)
    print(out.decode('utf-8'))
    print()

    print_model_info(sbml_file_name)
    print()

    # create condition table
    condition_df = pd.DataFrame(data={
        ptc.CONDITION_ID: ["condition_0", "condition_1",
                           "condition_2", "condition_3"],
        "k0": [1, 1.1, 1.2, 1.3]
    })
    condition_df.set_index([ptc.CONDITION_ID], inplace=True)

    # create observables
    observable_df = pd.DataFrame(data={
        ptc.OBSERVABLE_ID:
            ["obs_x1", "obs_x2", "obs_x3",
             "obs_x1_scaled", "obs_x2_offsetted", "obs_x1withsigma"],
        ptc.OBSERVABLE_FORMULA:
            ["x1", "x2", "x3", "observableParameter1_obs_x1_scaled * x1",
             "observableParameter1_obs_x2_offsetted + x2", "x1"],
        ptc.NOISE_FORMULA:
            ["noiseParameter1_obs_x1", "noiseParameter1_obs_x2",
             "noiseParameter1_obs_x3",
             "noiseParameter1_obs_x1_scaled",
             "noiseParameter1_obs_x2_offsetted",
             "noiseParameter1_obs_x1withsigma"],
    })
    observable_df.set_index([ptc.OBSERVABLE_ID], inplace=True)


    # assemble PEtab problem
    pp = petab.Problem(model=SbmlModel.from_file(sbml_file_name))
    pp.observable_df = observable_df
    pp.condition_df = condition_df
    # dummy parameter and measurement df needed for amici-import
    pp.measurement_df = petab.create_measurement_df()
    pp.parameter_df = pp.create_parameter_df(include_optional=True, lower_bound=1e-3, upper_bound=1e5)
    pp.parameter_df[ptc.ESTIMATE] = 1


    create_module(petab_problem=pp, model_name=model_name,
                  model_output_dir=args.model_output_dir)

    # load model
    sys.path.insert(0, args.model_output_dir)
    model_module = importlib.import_module(model_name)

    print()
    print("--- Creating data ---")

    measurement_df, true_parameters, expected_llh = create_data_tables(
        model=model_module.getModel(),
        condition_df=condition_df)

    pp.measurement_df = measurement_df
    create_parameter_table(problem=pp, nominal_parameters=true_parameters)

    # check for valid PEtab
    petab.lint_problem(pp)

    # Save remaining tables
    pp.to_files(measurement_file=measurement_file_name,
                condition_file=condition_file_name,
                observable_file=observable_file_name,
                parameter_file=parameter_file_name)

    # Create PEtab yaml file
    config = {
        'format_version': petab.__format_version__,
        'parameter_file': parameter_file_name,
        'problems': [
            {
                ptc.SBML_FILES: [sbml_file_name],
                ptc.CONDITION_FILES: [condition_file_name],
                ptc.MEASUREMENT_FILES: [measurement_file_name],
                ptc.OBSERVABLE_FILES: [observable_file_name],
            },
        ]
    }
    petab.validate(config)#, path_prefix=model_dir)
    with open(yaml_file_name, 'w') as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

    # create training data
    generate_hdf5_file(yaml_file=yaml_file_name,
                       model_output_dir=args.model_output_dir,
                       hdf5_file_name=args.hdf5_file_name,
                       model_name=model_name)

    create_test_data(measurement_file_name, parameter_file_name, config,
                     yaml_file_name_test, args.model_output_dir, model_name,
                     args.hdf5_file_name)

    save_expected_results(args.hdf5_file_name, true_parameters, expected_llh)

    write_starting_points(args.hdf5_file_name, true_parameters)

    hdf5_file_minibatch = os.path.join(os.path.dirname(args.hdf5_file_name),
                                       'example_data_minibatch.h5')
    from shutil import copyfile
    copyfile(args.hdf5_file_name, hdf5_file_minibatch)


if __name__ == '__main__':
    main()
