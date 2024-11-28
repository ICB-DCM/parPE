#include "steadyStateMultiConditionDataprovider.h"

#include <parpecommon/misc.h>
#include <parpeoptimization/optimizationOptions.h>

#include <cstdio>

SteadyStateMultiConditionDataProvider::SteadyStateMultiConditionDataProvider(
        std::unique_ptr<amici::Model> model,
        std::string const& hdf5Filename,
        std::string const& rootPath)
    : MultiConditionDataProviderHDF5(std::move(model), hdf5Filename, rootPath)
{
    solver_ = MultiConditionDataProviderHDF5::getModel()->getSolver();
    setupModelAndSolver();
}


std::unique_ptr<amici::Solver> SteadyStateMultiConditionDataProvider::getSolver() const
{
    return std::unique_ptr<amici::Solver>(solver_->clone());
}


void SteadyStateMultiConditionDataProvider::setupModelAndSolver() const {
    // calculate sensitivities for all parameters
    //model.requireSensitivitiesForAllParameters();

    solver_->setSensitivityOrder(amici::SensitivityOrder::first);
    // solver_->setSensitivityMethod(amici::SensitivityMethod::adjoint);
    // solver_->setMaxSteps(10000);
    // solver_->setNewtonMaxSteps(40);
}

