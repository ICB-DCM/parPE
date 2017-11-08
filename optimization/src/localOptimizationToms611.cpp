#include <localOptimizationToms611.h>
#include "logging.h"
#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include <cmath>
#include <toms611.h>
#include <misc.h>

namespace parpe {

static_assert(sizeof(doublereal) == sizeof(double), "Float size mismatch");

void calcf(integer const &n, doublereal const *x, integer &nf, doublereal &f,
           OptimizationProblem *problem, doublereal *urparm, void *ufparm) {

    if(!withinBounds(n, x, problem->getParametersMin(), problem->getParametersMax())) {
        nf = 0; // tells optimizer to choose a shorter step
        return;
    }

    problem->evaluateObjectiveFunction(x, &f, nullptr);
    *urparm = f;

    if(std::isnan(f)) {
        nf = 0; // tells optimizer to choose a shorter step
        return;
    }
}

void calcg(integer const &n, doublereal const *x, integer &nf, doublereal *g,
           OptimizationProblem *problem, doublereal *urparm, void *ufparm) {
    static __thread int numFunctionCalls = 0;

    if(!withinBounds(n, x, problem->getParametersMin(), problem->getParametersMax())) {
        nf = 0; // tells optimizer to choose a shorter step
        return;
    }

    problem->evaluateObjectiveFunction(x, urparm, g);
    problem->intermediateFunction(0, numFunctionCalls, *urparm, 0, 0, 0, 0, 0, 0, 0, 0);
    ++numFunctionCalls;

    if(std::isnan(*urparm)) {
        nf = 0; // tells optimizer to choose a shorter step
        return;
    }
}


int OptimizerToms611TrustRegionSumsl::optimize(OptimizationProblem *problem)
{
    integer numOptimizationVariables = problem->getNumOptimizationParameters();

    // allocate toms 611 memory and set options
    integer liv = toms611_sumsl_iv_min_length;
    integer iv[liv];
    iv[0] = 0; // fresh start, make sumsl_ call deflt_

    integer lv = toms611_sumsl_v_min_length(numOptimizationVariables)*10;
    std::vector<doublereal> v(lv);

    doublereal scaling[numOptimizationVariables];
    std::fill(scaling, scaling + numOptimizationVariables, 1.0);

    // Initialize optimizer options and memory
    integer deflt_algorithm = 2; /* general unconstrained optimization constants */
    deflt_(deflt_algorithm, iv, liv, lv, v.data());

    if(iv[0] != 12) // dflt_ success
        throw std::exception();

    // change default options
    auto o = problem->getOptimizationOptions();
    iv[17] = o.maxOptimizerIterations; // mxiter
    iv[23] = 0; // x0prt: don't print x0 and scaling

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
           lv, v.data(),
           reinterpret_cast<integer *>(problem), // sumsl_ only lets us pass integer, real or function...
           &fval, nullptr);

    clock_t timeEnd = clock();
    double wallTime = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;

    problem->logOptimizerFinished(fval, parameters, wallTime, iv[0]);

    return iv[0] >= first_error_code;
}


} // namespace parpe
