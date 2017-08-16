#ifndef STEADYSTATEPROBLEM_PARALLEL_H
#define STEADYSTATEPROBLEM_PARALLEL_H

#include "steadystateProblem.h"

class SteadystateProblemParallel : public SteadystateProblem {
  public:
    SteadystateProblemParallel();

    int evaluateObjectiveFunction(const double *parameters, double *objFunVal,
                                  double *objFunGrad);

    int evaluateParallel(const double *parameters, double *objFunVal,
                         double *objFunGrad);

    int evaluateSerial(const double *parameters, double *objFunVal,
                       double *objFunGrad);

    ~SteadystateProblemParallel();

  protected:
    int commSize;

    int numConditions;
};

#endif // STEADYSTATEPROBLEM_PARALLEL_H
