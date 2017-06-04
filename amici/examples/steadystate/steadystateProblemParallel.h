#ifndef STEADYSTATEPROBLEM_H
#define STEADYSTATEPROBLEM_H

#include "optimizationProblem.h"
#include "include/amici_interface_cpp.h"

class SteadystateProblemParallel : public OptimizationProblem
{
public:
    SteadystateProblemParallel(int numConditions);

    int evaluateObjectiveFunction(const double *parameters, double *objFunVal, double *objFunGrad);

    int evaluateParallel(const double *parameters, double *objFunVal, double *objFunGrad);

    int evaluateSerial(const double *parameters, double *objFunVal, double *objFunGrad);

    int intermediateFunction(int alg_mod,
                             int iter_count,
                             double obj_value,
                             double inf_pr, double inf_du,
                             double mu,
                             double d_norm,
                             double regularization_size,
                             double alpha_du, double alpha_pr,
                             int ls_trials);


    void logObjectiveFunctionEvaluation(const double *parameters,
                                                double objectiveFunctionValue,
                                                const double *objectiveFunctionGradient,
                                                int numFunctionCalls,
                                                double timeElapsed);

    void logOptimizerFinished(double optimalCost,
                                const double *optimalParameters,
                                double masterTime,
                                int exitStatus);

    ~SteadystateProblemParallel();


    void setupUserData();
    void setupExpData();

    UserData *udata;
    ExpData *edata;
    int commSize;

    int numConditions;

    double *fixedParameters;
};

#endif // STEADYSTATEPROBLEM_H
