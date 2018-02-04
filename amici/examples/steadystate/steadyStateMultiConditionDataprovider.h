#ifndef STEADYSTATEMULTICONDITIONPROBLEM_H
#define STEADYSTATEMULTICONDITIONPROBLEM_H

#include "MultiConditionDataProvider.h"
#include "multiConditionProblem.h"
#include "steadystateProblem.h"

#include <amici_hdf5.h>
#include <amici.h>

#include <memory>
#include <assert.h>

#include "wrapfunctions.h"

/**
 * @brief The SteadyStateMultiConditionDataProvider class provides the interface
 * to a HDF5 data file. Some non-default paths within the hdf5 file are set here.
 */
class SteadyStateMultiConditionDataProvider
    : public parpe::MultiConditionDataProvider {

  public:
    SteadyStateMultiConditionDataProvider(std::unique_ptr<amici::Model> model,
                                          std::string hdf5Filename, std::string rootPath = "");

    int getNumConditionSpecificParametersPerSimulation() const override;

    std::unique_ptr<amici::Model> getModel() const override;
    std::unique_ptr<amici::Solver> getSolver() const override;

    ~SteadyStateMultiConditionDataProvider() = default;

private:
    void setupModelAndSolver(amici::Model& model, amici::Solver& solver) const;

    std::unique_ptr<amici::Solver> solver;

};

#endif // STEADYSTATEMULTICONDITIONPROBLEM_H
