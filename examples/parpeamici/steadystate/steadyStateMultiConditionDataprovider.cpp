#include "steadyStateMultiConditionDataprovider.h"

#include <parpecommon/misc.h>
#include <parpeoptimization/optimizationOptions.h>

#include <cstdio>

SteadyStateMultiConditionDataProvider::SteadyStateMultiConditionDataProvider(
        std::unique_ptr<amici::Model> model,
        std::string const& hdf5Filename,
        std::string const& rootPath)
    : MultiConditionDataProviderHDF5(std::move(model), hdf5Filename, rootPath),
      solver_(this->model_->getSolver())
{
    setupModelAndSolver(*this->model_, *this->solver_);

}

std::unique_ptr<amici::Model> SteadyStateMultiConditionDataProvider::getModel() const
{
    return std::unique_ptr<amici::Model>(model_->clone());
}

std::unique_ptr<amici::Solver> SteadyStateMultiConditionDataProvider::getSolver() const
{
    return std::unique_ptr<amici::Solver>(solver_->clone());
}


void SteadyStateMultiConditionDataProvider::setupModelAndSolver(amici::Model &model, amici::Solver &solver) const {
    //hsize_t m = 0, n = 0;
    //model.setTimepoints(amici::hdf5::getDoubleDataset2D(fileId, rootPath + "/parameters/t", m, n));
    // set model constants

    // calculate sensitivities for all parameters
    model.requireSensitivitiesForAllParameters();
    //model.setParameterScale(amici::AMICI_SCALING_LOG10);

    solver.setSensitivityOrder(amici::SensitivityOrder::first);
    //solver.setSensitivityMethod(amici::AMICI_SENSI_FSA);
    solver.setSensitivityMethod(amici::SensitivityMethod::adjoint);
    solver.setMaxSteps(10000);
    solver.setNewtonMaxLinearSteps(100);
    solver.setNewtonMaxSteps(40);
}

