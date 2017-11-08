#include <localOptimizationToms611.h>
#include "logging.h"
#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include <cmath>
#include <toms611.h>

namespace parpe {

static_assert(sizeof(doublereal) == sizeof(double), "Float size mismatch");

void calcf(integer &n, doublereal *x, integer &nf, doublereal &f,
           OptimizationProblem *problem, doublereal *urparm, void *ufparm) {

    problem->evaluateObjectiveFunction(x, &f, nullptr);
    *urparm = f;
}

void calcg(integer &n, doublereal *x, integer &nf, doublereal *g,
           OptimizationProblem *problem, doublereal *urparm, void *ufparm) {
    static int __thread numFunctionCalls = 0;
    problem->evaluateObjectiveFunction(x, urparm, g);
    problem->intermediateFunction(0, numFunctionCalls, *urparm, 0, 0, 0, 0, 0, 0, 0, 0);
    ++numFunctionCalls;
}


int OptimizerToms611TrustRegionSumsl::optimize(OptimizationProblem *problem)
{
    integer numOptimizationVariables = problem->getNumOptimizationParameters();

    // allocate toms 611 memory and set options
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

    double fval = NAN; // the last computed cost function value; is this necessarily the one for the final parameters?
    clock_t timeBegin = clock();

    sumsl_(numOptimizationVariables,
           scaling,
           parameters,
           reinterpret_cast<S_fp>(calcf), (S_fp)calcg,
           iv, liv,
           lv, v,
           reinterpret_cast<integer *>(problem), // sumsl_ only lets us pass integer, real or function...
           &fval, nullptr);

    clock_t timeEnd = clock();
    double wallTime = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;

    problem->logOptimizerFinished(fval, parameters, wallTime, iv[0]);

    return iv[0] >= first_error_code;
}


} // namespace parpe
