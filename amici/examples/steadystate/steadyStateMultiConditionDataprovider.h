#ifndef STEADYSTATEMULTICONDITIONPROBLEM_H
#define STEADYSTATEMULTICONDITIONPROBLEM_H

#include "steadystateProblem.h"

#include <MultiConditionDataProvider.h>
#include <multiConditionProblem.h>

#include <amici_hdf5.h>

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
    SteadyStateMultiConditionDataProvider(Model *model,
                                          std::string hdf5Filename, std::string rootPath = "");

    int getNumConditionSpecificParametersPerSimulation() const override;

    std::unique_ptr<amici::UserData> getUserData() const override;

    ~SteadyStateMultiConditionDataProvider() = default;

private:
    void setupUserData(UserData *udata) const;

    std::unique_ptr<UserData> udata;
};

#endif // STEADYSTATEMULTICONDITIONPROBLEM_H
