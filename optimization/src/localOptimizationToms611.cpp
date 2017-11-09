#include <localOptimizationToms611.h>
#include "logging.h"
#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include <cmath>
#include <toms611.h>
#include <misc.h>

namespace parpe {

static_assert(sizeof(doublereal) == sizeof(double), "Float size mismatch");

// For C -> 0-based! (iv)
enum toms611intOptionsIndices {
    inith = 24,
    mxfcal = 16,
    mxiter = 17,
    outlev = 18,
    solprt = 21,
    x0prt = 23
};

// For C -> 0-based! (v)
enum toms611realOptionsIndices {
    bias = 42,
    afctol = 30,
    dinit = 37,
    lmax0 = 34,
    lmaxs = 35,
    rfctol = 31,
    sctol = 36,
    tuner1 = 25,
    xctol = 32,
    xftol = 33,
};


void setToms611Option(const std::pair<const std::string, const std::string> &pair, std::pair<integer *, doublereal *> options) {
    // for iterating over OptimizationOptions

    const std::string &key = pair.first;
    const std::string &val = pair.second;

    integer *iv = options.first;
    doublereal *v = options.second;

    if(key == "mxfcal") {
        iv[mxfcal] = std::stoi(val);
    } else if(key == "mxiter") {
        iv[mxiter] = std::stoi(val);
    } else if(key == "outlev") {
        iv[outlev] = std::stoi(val);
    } else if(key == "solprt") {
        iv[solprt] = std::stoi(val);
    } else if(key == "x0prt") {
        iv[x0prt] = std::stoi(val);
    } else if(key == "bias") {
        v[bias] = std::stod(val);
    } else if(key == "afctol") {
        v[afctol] = std::stod(val);
    } else if(key == "dinit") {
        v[dinit] = std::stod(val);
    } else if(key == "lmax0") {
        v[lmax0] = std::stod(val);
    } else if(key == "lmaxs") {
        v[lmaxs] = std::stod(val);
    } else if(key == "rfctol") {
        v[rfctol] = std::stod(val);
    } else if(key == "sctol") {
        v[sctol] = std::stod(val);
    } else if(key == "tuner1") {
        v[tuner1] = std::stod(val);
    } else if(key == "bias") {
        v[bias] = std::stod(val);
    } else if(key == "xftol") {
        v[xftol] = std::stod(val);
    } else {
        logmessage(LOGLVL_WARNING, "Ignoring unknown optimization option %s.", key.c_str());
        return;
    }

    logmessage(LOGLVL_WARNING, "Set optimization option %s to %s.", key.c_str(), val.c_str());
}


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

    integer lv = toms611_sumsl_v_min_length(numOptimizationVariables);
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

    problem->getOptimizationOptions().for_each< std::pair<integer *, doublereal *> >(setToms611Option, std::pair<integer *, doublereal *>(iv, v.data()));


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
