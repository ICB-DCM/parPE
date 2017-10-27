#include <localOptimizationToms611.h>
#include "logging.h"
#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include <toms611.h>

namespace parpe {

void calcf(integer &n, doublereal *x, integer &nf, doublereal &f,
           OptimizationProblem *problem, doublereal *urparm, void *ufparm) {
    static_assert(sizeof(doublereal) == sizeof(double), "Float size mismatch");

    problem->evaluateObjectiveFunction(x, &f, nullptr);
}

void calcg(integer &n, doublereal *x, integer &nf, doublereal *g,
           OptimizationProblem *problem, doublereal *urparm, void *ufparm) {

    static_assert(sizeof(doublereal) == sizeof(double), "Float size mismatch");
    double unusedFVal;
    problem->evaluateObjectiveFunction(x, &unusedFVal, g);

}


int OptimizerToms611TrustRegionSumsl::optimize(OptimizationProblem *problem)
{
    integer numOptimizationVariables = problem->getNumOptimizationParameters();

    integer liv = toms611_sumsl_iv_min_length;
    integer iv[liv];
    iv[0] = 0; // fresh start, make sumsl_ call deflt_

    integer lv = toms611_sumsl_v_min_length(numOptimizationVariables);
    doublereal v[lv];

    doublereal scaling[numOptimizationVariables];
    std::fill(scaling, scaling + numOptimizationVariables, 1.0);

    doublereal parameters[numOptimizationVariables];
    double const* startingPoint = problem->getInitialParameters();
    std::copy(startingPoint, startingPoint + numOptimizationVariables, parameters);

    sumsl_(numOptimizationVariables,
           scaling,
           parameters,
           reinterpret_cast<S_fp>(calcf), (S_fp)calcg,
           iv, liv,
           lv, v,
           reinterpret_cast<integer *>(problem), // sumsl_ only lets us pass integer, real or function...
           nullptr, nullptr);

    return iv[0] >= first_error_code;
}


} // namespace parpe
