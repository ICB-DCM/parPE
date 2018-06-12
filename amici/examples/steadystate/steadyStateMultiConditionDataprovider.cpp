#include "steadyStateMultiConditionDataprovider.h"
#include <optimizationOptions.h>
#include <cstdio>
#include <misc.h>
#include <multiConditionProblemResultWriter.h>

SteadyStateMultiConditionDataProvider::SteadyStateMultiConditionDataProvider(std::unique_ptr<amici::Model> model, std::string hdf5Filename, std::string rootPath)
    : MultiConditionDataProviderHDF5(std::move(model), hdf5Filename, rootPath),
      solver(this->model->getSolver())
{
    setupModelAndSolver(*this->model, *this->solver);

}

std::unique_ptr<amici::Model> SteadyStateMultiConditionDataProvider::getModel() const
{
    return std::unique_ptr<amici::Model>(model->clone());
}

std::unique_ptr<amici::Solver> SteadyStateMultiConditionDataProvider::getSolver() const
{
    return std::unique_ptr<amici::Solver>(solver->clone());
}


void SteadyStateMultiConditionDataProvider::setupModelAndSolver(amici::Model &model, amici::Solver &solver) const {
    //hsize_t m = 0, n = 0;
    //model.setTimepoints(amici::hdf5::getDoubleDataset2D(fileId, rootPath + "/parameters/t", m, n));
    // set model constants

    // calculate sensitivities for all parameters
    model.requireSensitivitiesForAllParameters();
    //model.setParameterScale(amici::AMICI_SCALING_LOG10);

    solver.setSensitivityOrder(amici::AMICI_SENSI_ORDER_FIRST);
    //solver.setSensitivityMethod(amici::AMICI_SENSI_FSA);
    solver.setSensitivityMethod(amici::AMICI_SENSI_ASA);
    solver.setMaxSteps(1e4);
    solver.setNewtonMaxLinearSteps(100);
    solver.setNewtonMaxSteps(40);
}

