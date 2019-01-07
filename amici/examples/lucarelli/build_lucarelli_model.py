# import all necessary modules
import amici
import os
import sys
import numpy as np
import pandas as pd
import libsbml

# set the model specific variables and names
sbml_file = 'lucarelli_12.xml'
model_name = 'lucarelli_12'
model_output_dir = 'lucarelli_12'
SBMLreader = libsbml.SBMLReader()

# import the model from sbml
sbml_doc = SBMLreader.readSBML(sbml_file)
sbml_model = sbml_doc.getModel()
for reaction in sbml_model.getListOfReactions():
    reactants = ' + '.join(['%s %s' % (
        int(r.getStoichiometry()) if r.getStoichiometry() > 1 else '',
        r.getSpecies()) for r in reaction.getListOfReactants()])
    products = ' + '.join(['%s %s' % (
        int(r.getStoichiometry()) if r.getStoichiometry() > 1 else '',
        r.getSpecies()) for r in reaction.getListOfProducts()])
    reversible = '<' if reaction.getReversible() else ''

# set observables and constants
observables_list = ['observable_Ski', 'observable_Skil', 'observable_Dnmt3a',
                    'observable_Sox4', 'observable_Jun', 'observable_Smad7',
                    'observable_Klf10', 'observable_Bmp4', 'observable_Cxcl15',
                    'observable_Dusp5', 'observable_Tgfa',
                    'observable_Pdk4']  #
fixed_parameters = ['init_TGFb', 'init_S2', 'init_S3', 'init_S4']

# wrap the model
sbmlImporter = amici.SbmlImporter(sbml_file)
sbml = sbmlImporter.sbml

observables = amici.assignmentRules2observables( \
    sbmlImporter.sbml,  # the libsbml model object
    filter_function=lambda variable: variable.getId() in observables_list)

sbmlImporter.sbml2amici(model_name,
output_dir = model_output_dir,
observables = observables,
constantParameters = fixed_parameters)

# load model
sys.path.insert(0, model_output_dir)
import lucarelli_12 as modelModule

model = modelModule.getModel()
solver = model.getSolver()
rdata = amici.runAmiciSimulation(model, solver)
default_parameters = np.array(model.getParameters())

# create data as a data frame
from pandas import DataFrame


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

# generate conditon-vectors
sigma = 0.5
mu = 1.5
nConditions = 10
init_TGFb = np.linspace(0, 1, 11)
init_conditions = [10.0 ** (mu + sigma * np.random.randn(4))]
init_conditions[0][0] = init_TGFb[np.random.randint(0, 11)]
for i in range(1, nConditions):
    init_conditions.append(np.array(10.0 ** (mu + sigma * np.random.randn(4))))
    init_conditions[i][0] = init_TGFb[np.random.randint(0, 11)]

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
import subprocess

subprocess.call(["/bin/bash", "-c", "if [[ -f example_data.h5 ]]; then cp example_data.h5 example_data.h5.bak; fi"])

out = subprocess.run(['%s/generateHDF5DataFileFromText.py' % os.path.join(
    os.getcwd(), '..', '..', '..', 'misc'),
                      hdf5File,
                      sbml_file,
                      model_output_dir,
                      measurement_file,
                      fixed_parameter_file], stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT)
print(out.stdout.decode("utf-8"))
assert 0 == out.returncode, "Returned " + str(out.returncode)
