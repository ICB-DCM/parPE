#ifndef STEADYSTATEPROBLEM_PARALLEL_H
#define STEADYSTATEPROBLEM_PARALLEL_H

#include "steadystateProblem.h"
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>

class SteadystateProblemParallel : public ExampleSteadystateProblem,
                                   public LoadBalancerWorker {
  public:
    SteadystateProblemParallel(LoadBalancerMaster *loadBalancer);

    int evaluateObjectiveFunction(const double *parameters, double *objFunVal,
                                  double *objFunGrad);

    int evaluateParallel(const double *parameters, double *objFunVal,
                         double *objFunGrad);

    int evaluateSerial(const double *parameters, double *objFunVal,
                       double *objFunGrad);

    void messageHandler(char **buffer, int *size, int jobId);

    ~SteadystateProblemParallel();

    LoadBalancerMaster *loadBalancer;

  protected:
    int commSize;

    int numConditions;
};

#endif // STEADYSTATEPROBLEM_PARALLEL_H
