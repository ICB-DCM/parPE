#ifndef STEADYSTATEPROBLEM_PARALLEL_H
#define STEADYSTATEPROBLEM_PARALLEL_H

#include "steadystateProblem.h"
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>
#include <memory>
#include "SteadyStateMultiConditionProblem.h"

class SteadystateProblemParallel : public parpe::OptimizationProblem, public parpe::LoadBalancerWorker {
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

    std::unique_ptr<SteadyStateMultiConditionDataProvider> dataProvider;
    std::unique_ptr<Model> model;
    std::unique_ptr<amici::UserData> udata;
};

#endif // STEADYSTATEPROBLEM_PARALLEL_H
