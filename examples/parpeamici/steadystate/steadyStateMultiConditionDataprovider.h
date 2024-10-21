#ifndef STEADYSTATEMULTICONDITIONPROBLEM_H
#define STEADYSTATEMULTICONDITIONPROBLEM_H

#include "steadystateProblem.h"
#include <parpeamici/multiConditionDataProvider.h>
#include <parpeamici/multiConditionProblem.h>

#include <amici/amici.h>
#include <amici/hdf5.h>

#include <assert.h>
#include <memory>

/**
 * @brief The SteadyStateMultiConditionDataProvider class provides the interface
 * to a HDF5 data file. Some non-default paths within the hdf5 file are set
 * here.
 */
class SteadyStateMultiConditionDataProvider
    : public parpe::MultiConditionDataProviderHDF5 {

  public:
    SteadyStateMultiConditionDataProvider(
        std::unique_ptr<amici::Model> model,
        std::string const& hdf5Filename,
        std::string const& rootPath = "");

    std::unique_ptr<amici::Solver> getSolver() const override;

    ~SteadyStateMultiConditionDataProvider() override = default;

  private:
    void setupModelAndSolver() const;

    std::unique_ptr<amici::Solver> solver_;
};

#endif // STEADYSTATEMULTICONDITIONPROBLEM_H
