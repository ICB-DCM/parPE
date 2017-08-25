#ifndef STEADYSTATEPROBLEM_PARALLEL_H
#define STEADYSTATEPROBLEM_PARALLEL_H

#include "steadystateProblem.h"
#include <LoadBalancerMaster.h>

class SteadystateProblemParallel : public ExampleSteadystateProblem {
  public:
    SteadystateProblemParallel(LoadBalancerMaster *loadBalancer);

    int evaluateObjectiveFunction(const double *parameters, double *objFunVal,
                                  double *objFunGrad);

    int evaluateParallel(const double *parameters, double *objFunVal,
                         double *objFunGrad);

    int evaluateSerial(const double *parameters, double *objFunVal,
                       double *objFunGrad);

    ~SteadystateProblemParallel();

    LoadBalancerMaster *loadBalancer;

  protected:
    int commSize;

    int numConditions;

};

#endif // STEADYSTATEPROBLEM_PARALLEL_H
