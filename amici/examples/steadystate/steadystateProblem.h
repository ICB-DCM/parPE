#ifndef STEADYSTATEPROBLEM_H
#define STEADYSTATEPROBLEM_H

#include "include/amici_interface_cpp.h"
#include "optimizationProblem.h"
#include <hdf5.h>

class Model;

class ExampleSteadystateProblem : public parPE::OptimizationProblem {
  public:
    ExampleSteadystateProblem();

    int evaluateObjectiveFunction(const double *parameters, double *objFunVal,
                                  double *objFunGrad) override;

    int intermediateFunction(int alg_mod, int iter_count, double obj_value,
                             double inf_pr, double inf_du, double mu,
                             double d_norm, double regularization_size,
                             double alpha_du, double alpha_pr,
                             int ls_trials) override;

    void logObjectiveFunctionEvaluation(const double *parameters,
                                        double objectiveFunctionValue,
                                        const double *objectiveFunctionGradient,
                                        int numFunctionCalls,
                                        double timeElapsed) override;

    void logOptimizerFinished(double optimalCost,
                              const double *optimalParameters,
                              double masterTime, int exitStatus) override;

    ~ExampleSteadystateProblem();

    void requireSensitivities(bool sensitivitiesRequired);
    void readFixedParameters(int conditionIdx);
    void readMeasurement(int conditionIdx);

    UserData *udata = nullptr;
    ExpData *edata = nullptr;
    Model *model = nullptr;

  protected:
    void setupUserData(int conditionIdx);
    void setupExpData(int conditionIdx);

    hid_t fileId = -1;
};

#endif // STEADYSTATEPROBLEM_H
