#include "localOptimizationIpoptTNLP.h"
#include "optimizationProblem.h"
#include <cassert>
#include <cstring>
#include <IpIpoptData.hpp>
#include <IpIpoptCalculatedQuantities.hpp>

namespace parpe {


LocalOptimizationIpoptTNLP::LocalOptimizationIpoptTNLP(
    OptimizationProblem *problem, double& finalCost, std::vector<double>& finalParameters)
    : problem(problem), reporter(problem->getReporter()), finalCost(finalCost), finalParameters(finalParameters) {
    cachedGradient.resize(problem->costFun->numParameters());
}


bool LocalOptimizationIpoptTNLP::get_nlp_info(Index &n, Index &m,
                                              Index &nnz_jac_g,
                                              Index &nnz_h_lag,
                                              IndexStyleEnum &index_style) {

    n = problem->costFun->numParameters();
    m = 0;                       // number of constrants
    nnz_jac_g = 0;               // numNonZeroElementsConstraintJacobian
    nnz_h_lag = 0;               // numNonZeroElementsLagrangianHessian
    index_style = TNLP::C_STYLE; // array layout for sparse matrices

    return true;
}

bool LocalOptimizationIpoptTNLP::get_bounds_info(Index n, Number *x_l,
                                                 Number *x_u, Index m,
                                                 Number *g_l, Number *g_u) {
    // parameter bounds
    problem->fillParametersMin(x_l);
    problem->fillParametersMax(x_u);

    // no constraints supported for now -> no constraint bounds

    return true;
}

bool LocalOptimizationIpoptTNLP::get_starting_point(Index n, bool init_x,
                                                    Number *x, bool init_z,
                                                    Number *z_L, Number *z_U,
                                                    Index m, bool init_lambda,
                                                    Number *lambda) {
    /* this function is called twice by IpOpt which is a problem if problem->fillInitialParameters provides random parameters, therefore initial point needs to be stored */

    if (init_x) {
        if(initialParameters.size() == 0) {
            initialParameters.resize(n);
            problem->fillInitialParameters(initialParameters.data());
            if(reporter && reporter->starting(n, initialParameters.data()))
                return false;
        }
        std::copy(initialParameters.begin(), initialParameters.end(), x);
    }

    assert(init_z == false);
    assert(init_lambda == false);

    return true;
}

bool LocalOptimizationIpoptTNLP::eval_f(Index n, const Number *x, bool new_x,
                                        Number &obj_value) {
    auto unlockIpOpt = ipOptReleaseLock();

    // update cached parameters
    finalParameters.resize(n);
    std::copy(x, x + n, finalParameters.begin());

    int errors = 0;
    if(reporter && reporter->beforeCostFunctionCall(n, x) != 0)
        return true;

    if (new_x || !haveCachedCost) {
        // Have to compute anew
        errors = problem->costFun->evaluate(x, obj_value, nullptr);
        cachedErrors = errors;
        cachedCost = obj_value;
        haveCachedCost = true;
        haveCachedGradient = false;
    } else {
        // recycle old result
        errors = cachedErrors;
        obj_value = cachedCost;
    }

    if(reporter && reporter->afterCostFunctionCall(n, x, cachedCost, nullptr))
        return true;

    return errors == 0;
}

bool LocalOptimizationIpoptTNLP::eval_grad_f(Index n, const Number *x,
                                             bool new_x, Number *grad_f) {
    auto unlockIpOpt = ipOptReleaseLock();

    // update cached parameters
    finalParameters.resize(n);
    std::copy(x, x + n, finalParameters.begin());

    if(reporter && reporter->beforeCostFunctionCall(n, x) != 0)
        return true;

    int errors = 0;

    if (new_x || !haveCachedGradient) {
        // Have to compute anew
        errors = problem->costFun->evaluate(x, cachedCost, grad_f);
        std::copy(grad_f, grad_f + n, cachedGradient.begin());
        cachedErrors = errors;
        haveCachedCost = true;
        haveCachedGradient = true;
    } else {
        // recycle old result
        std::copy(cachedGradient.begin(), cachedGradient.end(), grad_f);
        errors = cachedErrors;
    }

    if(reporter)
        return (errors == 0) && !reporter->afterCostFunctionCall(n, x, cachedCost, cachedGradient.data());

    return errors == 0;
}

bool LocalOptimizationIpoptTNLP::eval_g(Index n, const Number *x, bool new_x,
                                        Index m, Number *g) {

    assert(false && "no constraints, should never get here");

    return false;
}

bool LocalOptimizationIpoptTNLP::eval_jac_g(Index n, const Number *x,
                                            bool new_x, Index m, Index nele_jac,
                                            Index *iRow, Index *jCol,
                                            Number *values) {
    // no constraints, nothing to do here, but will be called once

    if (new_x) {
        // because next function will be called with
        // new_x==false, but we didn't prepare anything
        haveCachedCost = false;
        haveCachedGradient = false;
    }
    return true;
}

bool LocalOptimizationIpoptTNLP::intermediate_callback(
    AlgorithmMode mode, Index iter, Number obj_value, Number inf_pr,
    Number inf_du, Number mu, Number d_norm, Number regularization_size,
    Number alpha_du, Number alpha_pr, Index ls_trials, const IpoptData *ip_data,
    IpoptCalculatedQuantities *ip_cq) {

    auto unlockIpOpt = ipOptReleaseLock();

    // better: find x in ip_data->curr()->x();
    int status = reporter->iterationFinished(problem->costFun->numParameters(),
                                             finalParameters.data(),
                                             obj_value,
                                             haveCachedGradient?cachedGradient.data():nullptr);

#ifdef INSTALL_SIGNAL_HANDLER
    if (caughtTerminationSignal) {
        logmessage(LOGLVL_CRITICAL, "CAUGHT SIGTERM... EXITING.");
        return false;
    }
#endif

    return status == 0;
}

void LocalOptimizationIpoptTNLP::finalize_solution(
    SolverReturn status, Index n, const Number *x, const Number *z_L,
    const Number *z_U, Index m, const Number *g, const Number *lambda,
    Number obj_value, const IpoptData *ip_data,
    IpoptCalculatedQuantities *ip_cq) {

    auto unlockIpOpt = ipOptReleaseLock();

    finalCost = obj_value;
    finalParameters.resize(n);
    std::copy(x, x + n, finalParameters.begin());

    reporter->finished(obj_value, x, status);
}

InverseUniqueLock ipOptReleaseLock()
{
    return InverseUniqueLock(&mutexIpOpt);
}

InverseUniqueLock::InverseUniqueLock(mutexIpOptType *mutex)
    : mutex(mutex)
{
    mutex->unlock();
}

InverseUniqueLock::~InverseUniqueLock()
{
    mutex->lock();
}

std::unique_lock<mutexIpOptType> ipOptGetLock()
{
    return std::unique_lock<mutexIpOptType>(mutexIpOpt);
}

} // namespace parpe
