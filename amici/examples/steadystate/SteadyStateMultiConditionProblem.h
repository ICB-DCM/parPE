#ifndef STEADYSTATEMULTICONDITIONPROBLEM_H
#define STEADYSTATEMULTICONDITIONPROBLEM_H

#include "MultiConditionDataProvider.h"
#include "assert.h"
#include "multiConditionProblem.h"
#include "steadystateProblem.h"
#include "wrapfunctions.h"
#include <amici_hdf5.h>

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


    void setupUserData(UserData *udata) const;

    UserData *getUserData() const override;

    ~SteadyStateMultiConditionDataProvider();

    UserData *udata;
};

class SteadyStateMultiConditionProblem : public parpe::MultiConditionProblem {
  public:
    SteadyStateMultiConditionProblem(SteadyStateMultiConditionDataProvider *dp,
                                     parpe::LoadBalancerMaster *loadBalancer);

    void setSensitivityOptions(bool sensiRequired) override;

    ~SteadyStateMultiConditionProblem();
};

#endif // STEADYSTATEMULTICONDITIONPROBLEM_H
