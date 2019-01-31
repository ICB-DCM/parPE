#ifndef STEADYSTATEMULTICONDITIONPROBLEM_H
#define STEADYSTATEMULTICONDITIONPROBLEM_H

#include <parpeamici/multiConditionDataProvider.h>
#include <parpeamici/multiConditionProblem.h>
#include "steadystateProblem.h"

#include <amici/hdf5.h>
#include <amici/amici.h>

#include <memory>
#include <assert.h>

#include "wrapfunctions.h"

/**
 * @brief The SteadyStateMultiConditionDataProvider class provides the interface
 * to a HDF5 data file. Some non-default paths within the hdf5 file are set here.
 */
class SteadyStateMultiConditionDataProvider
    : public parpe::MultiConditionDataProviderHDF5 {

  public:
    SteadyStateMultiConditionDataProvider(std::unique_ptr<amici::Model> model,
                                          const std::string &hdf5Filename, const std::string &rootPath = "");

    std::unique_ptr<amici::Model> getModel() const override;
    std::unique_ptr<amici::Solver> getSolver() const override;

    ~SteadyStateMultiConditionDataProvider() override = default;

private:
    void setupModelAndSolver(amici::Model& model, amici::Solver& solver) const;

    std::unique_ptr<amici::Solver> solver;

};

#endif // STEADYSTATEMULTICONDITIONPROBLEM_H
