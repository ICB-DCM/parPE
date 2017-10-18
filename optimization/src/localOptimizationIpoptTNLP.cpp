#include "localOptimizationIpoptTNLP.h"
#include "optimizationProblem.h"
#include <cassert>
#include <cstring>

LocalOptimizationIpoptTNLP::LocalOptimizationIpoptTNLP(
    OptimizationProblem *problem, pthread_mutex_t *ipoptMutex)
    : problem(problem), ipoptMutex(ipoptMutex) {
    timeBegin = clock();
}

LocalOptimizationIpoptTNLP::~LocalOptimizationIpoptTNLP() {
    if (lastGradient)
        delete[] lastGradient;
}

bool LocalOptimizationIpoptTNLP::get_nlp_info(Index &n, Index &m,
                                              Index &nnz_jac_g,
                                              Index &nnz_h_lag,
                                              IndexStyleEnum &index_style) {
    pthread_mutex_unlock(ipoptMutex);

    n = problem->getNumOptimizationParameters();
    m = 0;                       // number of constrants
    nnz_jac_g = 0;               // numNonZeroElementsConstraintJacobian
    nnz_h_lag = 0;               // numNonZeroElementsLagrangianHessian
    index_style = TNLP::C_STYLE; // array layout for sparse matrices

    pthread_mutex_lock(ipoptMutex);

    return true;
}

bool LocalOptimizationIpoptTNLP::get_bounds_info(Index n, Number *x_l,
                                                 Number *x_u, Index m,
                                                 Number *g_l, Number *g_u) {
    // parameter bounds
    memcpy(x_l, problem->getParametersMin(), sizeof(Number) * n);
    memcpy(x_u, problem->getParametersMax(), sizeof(Number) * n);

    // no constraints -> no constraint bounds

    return true;
}

bool LocalOptimizationIpoptTNLP::get_starting_point(Index n, bool init_x,
                                                    Number *x, bool init_z,
                                                    Number *z_L, Number *z_U,
                                                    Index m, bool init_lambda,
                                                    Number *lambda) {
    if (init_x) {
        const double *startingPoint = problem->getInitialParameters();
        if (startingPoint) {
            memcpy(x, startingPoint, sizeof(Number) * n);
        } else {
            problem->fillInitialParameters(x);
        }
    }

    assert(init_z == false);
    assert(init_lambda == false);

    return true;
}

bool LocalOptimizationIpoptTNLP::eval_f(Index n, const Number *x, bool new_x,
                                        Number &obj_value) {
    static __thread int numFunctionCalls = 0;
    ++numFunctionCalls;
    // logmessage(LOGLVL_DEBUG, "Eval_F (%d) #%d.", new_x, numFunctionCalls);

    pthread_mutex_unlock(ipoptMutex);

    int errors = 0;

    clock_t timeBegin = clock();

    if (new_x || !lastCostP) {
        errors = problem->evaluateObjectiveFunction(x, &obj_value, NULL);
        if (lastGradient) // invalidate
            delete[] lastGradient;
        lastGradient = NULL;
        lastErrors = errors;
        lastCost = obj_value;
        lastCostP = &lastCost;
    } else {
        errors = lastErrors;
        obj_value = lastCost;
    }

    clock_t timeEnd = clock();
    double wallTime = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;

    problem->logObjectiveFunctionEvaluation(x, obj_value, NULL,
                                            numFunctionCalls, wallTime);

    pthread_mutex_lock(ipoptMutex);

    return errors == 0;
}

bool LocalOptimizationIpoptTNLP::eval_grad_f(Index n, const Number *x,
                                             bool new_x, Number *grad_f) {
    static __thread int numFunctionCalls = 0;
    ++numFunctionCalls;
    // logmessage(LOGLVL_DEBUG, "Eval_Grad_F (%d) #%d", new_x,
    // numFunctionCalls);

    pthread_mutex_unlock(ipoptMutex);

    int errors = 0;

    clock_t timeBegin = clock();

    if (new_x || !lastCostP || !lastGradient) {
        errors = problem->evaluateObjectiveFunction(x, &lastCost, grad_f);

        if (!lastGradient)
            lastGradient = new Number[problem->getNumOptimizationParameters()];
        memcpy(lastGradient, grad_f, sizeof(Number) * n);
        lastCostP = &lastCost;
        lastErrors = errors;
    } else {
        memcpy(grad_f, lastGradient, sizeof(Number) * n);
        errors = lastErrors;
    }

    clock_t timeEnd = clock();
    double timeElapsed = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;

    problem->logObjectiveFunctionEvaluation(x, lastCost, grad_f,
                                            numFunctionCalls, timeElapsed);

    pthread_mutex_lock(ipoptMutex);

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

    if (new_x)
        lastCostP = NULL; // because next function will be called with
                          // new_x==false, but we didn't prepare anything

    return true;
}

bool LocalOptimizationIpoptTNLP::intermediate_callback(
    AlgorithmMode mode, Index iter, Number obj_value, Number inf_pr,
    Number inf_du, Number mu, Number d_norm, Number regularization_size,
    Number alpha_du, Number alpha_pr, Index ls_trials, const IpoptData *ip_data,
    IpoptCalculatedQuantities *ip_cq) {

    pthread_mutex_unlock(ipoptMutex);

    int status = true;
    status = problem->intermediateFunction(
        (int)mode, iter, obj_value, inf_pr, inf_du, mu, d_norm,
        regularization_size, alpha_du, alpha_pr, ls_trials);

#ifdef INSTALL_SIGNAL_HANDLER
    if (caughtTerminationSignal) {
        logmessage(LOGLVL_CRITICAL, "CAUGHT SIGTERM... EXITING.");
        return false;
    }
#endif

    pthread_mutex_lock(ipoptMutex);

    return status == 0;
}

void LocalOptimizationIpoptTNLP::finalize_solution(
    SolverReturn status, Index n, const Number *x, const Number *z_L,
    const Number *z_U, Index m, const Number *g, const Number *lambda,
    Number obj_value, const IpoptData *ip_data,
    IpoptCalculatedQuantities *ip_cq) {
    pthread_mutex_unlock(ipoptMutex);

    clock_t timeEnd = clock();
    double timeElapsed = (double)(timeEnd - timeBegin) / CLOCKS_PER_SEC;

    problem->logOptimizerFinished(obj_value, x, timeElapsed, status);

    pthread_mutex_lock(ipoptMutex);
}
