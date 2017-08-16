#ifndef MULTICONDITIONSTEADYSTATEPROBLEM_H
#define MULTICONDITIONSTEADYSTATEPROBLEM_H

#include "multiConditionProblem.h"

class MultiConditionSteadyStateProblem : public MultiConditionProblem {
  public:
    static ReturnData *
    runAndLogSimulation(UserData *udata,
                        MultiConditionDataProvider *dataProvider,
                        JobIdentifier path, int jobId, int *status);
};

#endif // MULTICONDITIONSTEADYSTATEPROBLEM_H
