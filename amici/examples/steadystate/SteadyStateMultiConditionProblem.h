#ifndef STEADYSTATEMULTICONDITIONPROBLEM_H
#define STEADYSTATEMULTICONDITIONPROBLEM_H

#include "MultiConditionDataProvider.h"
#include "assert.h"
#include "multiConditionProblem.h"
#include "steadystateProblem.h"
#include "wrapfunctions.h"
#include <amici_hdf5.h>
#include <memory>

/**
 * @brief The SteadyStateMultiConditionDataProvider class provides the interface
 * to a HDF5 data file
 */
class SteadyStateMultiConditionDataProvider
    : public parpe::MultiConditionDataProvider {

  public:
    SteadyStateMultiConditionDataProvider(Model *model,
                                          std::string hdf5Filename);

    int getNumberOfConditions() const override { return 12; }

    int getNumConditionSpecificParametersPerSimulation() const override {
        return 0;
    }

    int updateFixedSimulationParameters(int conditionIdx,
                                        UserData *udata) const override;

    UserData *getUserData() const override;

    ~SteadyStateMultiConditionDataProvider() = default;

private:
    void setupUserData(UserData *udata) const;

    std::unique_ptr<UserData> udata;
};

class SteadyStateMultiConditionProblem : public parpe::MultiConditionProblem {
  public:
    SteadyStateMultiConditionProblem(SteadyStateMultiConditionDataProvider *dp,
                                     parpe::LoadBalancerMaster *loadBalancer);

    ~SteadyStateMultiConditionProblem() = default;
};

#endif // STEADYSTATEMULTICONDITIONPROBLEM_H
