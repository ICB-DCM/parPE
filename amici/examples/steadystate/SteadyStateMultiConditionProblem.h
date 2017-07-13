#ifndef STEADYSTATEMULTICONDITIONPROBLEM_H
#define STEADYSTATEMULTICONDITIONPROBLEM_H

#include "multiConditionProblem.h"
#include "MultiConditionDataProvider.h"
#include "steadystateProblem.h"
#include "ami_hdf5.h"
#include "wrapfunctions.h"
#include "assert.h"

class SteadyStateMultiConditionDataProvider : public MultiConditionDataProvider
{
public:
    SteadyStateMultiConditionDataProvider(const char *hdf5Filename);

    int getNumberOfConditions() const { return 12; }

    int getNumConditionSpecificParametersPerSimulation() const { return 0; }

    int updateFixedSimulationParameters(int conditionIdx, UserData *udata) const;

    ExpData *getExperimentalDataForCondition(int conditionIdx) const;

    void setupUserData(UserData *udata) const;

    UserData *getUserData() const;

    UserData *udata;
};


class SteadyStateMultiConditionProblem : public MultiConditionProblem
{
public:
    SteadyStateMultiConditionProblem(SteadyStateMultiConditionDataProvider *dp);

    void setSensitivityOptions(bool sensiRequired);

};

#endif // STEADYSTATEMULTICONDITIONPROBLEM_H
