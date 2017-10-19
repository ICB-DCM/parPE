#ifndef STEADYSTATEPROBLEM_PARALLEL_H
#define STEADYSTATEPROBLEM_PARALLEL_H

#include "steadystateProblem.h"
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>

class SteadystateProblemParallel : public ExampleSteadystateProblem,
                                   public parpe::LoadBalancerWorker {
  public:
    SteadystateProblemParallel(parpe::LoadBalancerMaster *loadBalancer);

    int evaluateObjectiveFunction(const double *parameters, double *objFunVal,
                                  double *objFunGrad) override;

    int evaluateParallel(const double *parameters, double *objFunVal,
                         double *objFunGrad);

    int evaluateSerial(const double *parameters, double *objFunVal,
                       double *objFunGrad);

    void messageHandler(std::vector<char> &buffer, int jobId) override;

    ~SteadystateProblemParallel();

    parpe::LoadBalancerMaster *loadBalancer;

  protected:
    int commSize;

    int numConditions;
};

#endif // STEADYSTATEPROBLEM_PARALLEL_H
