#include "localOptimizationIpopt.h"
#include "optimizationProblem.h"

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <signal.h>
#include <alloca.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// This should fix `error "don't have header file for stddef"' with some IpOpt versions
#define HAVE_CSTDDEF
#include <IpTNLP.hpp>
#undef HAVE_CSTDDEF

#include <IpIpoptApplication.hpp>

#include <pthread.h>
#include "logging.h"

#ifdef INSTALL_SIGNAL_HANDLER
extern volatile sig_atomic_t caughtTerminationSignal;
#endif

/**
 * @brief ipoptMutex Ipopt seems not to be thread safe. Lock this mutex every time
 * when control is passed to ipopt functions.
 */
static pthread_mutex_t ipoptMutex = PTHREAD_MUTEX_INITIALIZER;

using namespace Ipopt;

class MyNLP : public Ipopt::TNLP
{
public:
    MyNLP(OptimizationProblem *problem) : problem(problem)
    {
        timeBegin = clock();
    }

    virtual ~MyNLP()
    {
        if(lastGradient)
            delete[] lastGradient;
    }

    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, IndexStyleEnum& index_style)
    {
        n = problem->numOptimizationParameters;
        m = 0; // number of constrants
        nnz_jac_g = 0; // numNonZeroElementsConstraintJacobian
        nnz_h_lag = 0; // numNonZeroElementsLagrangianHessian
        index_style = TNLP::C_STYLE; // array layout for sparse matrices

        return true;
    }


    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u)
    {
        // parameter bounds
        memcpy(x_l, problem->parametersMin, sizeof(Number) * n);
        memcpy(x_u, problem->parametersMax, sizeof(Number) * n);

        // no constraints -> no constraint bounds

        return true;
    }



    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda)
    {
        if(init_x) {
            if(problem->initialParameters) {
                memcpy(x, problem->initialParameters, sizeof(Number) * n);
            } else {
                getRandomStartingpoint(problem->parametersMin, problem->parametersMax, problem->numOptimizationParameters, x);
            }
        }

        assert(init_z == false);
        assert(init_lambda == false);

        return true;
    }


    virtual bool eval_f(Index n, const Number* x, bool new_x,
                        Number& obj_value)
    {
        static __thread int numFunctionCalls = 0;
        ++numFunctionCalls;
        // logmessage(LOGLVL_DEBUG, "Eval_F (%d) #%d.", new_x, numFunctionCalls);

        pthread_mutex_unlock(&ipoptMutex);

        int errors = 0;

        clock_t timeBegin = clock();

        if(new_x || !lastCostP) {
            errors = problem->evaluateObjectiveFunction(x, &obj_value, NULL);
            if(lastGradient) // invalidate
                delete[] lastGradient;
            lastGradient = NULL;
            lastErrors = errors;
            lastCost = obj_value;
        } else {
            errors = lastErrors;
            obj_value = lastCost;
        }

        clock_t timeEnd = clock();
        double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;

        problem->logObjectiveFunctionEvaluation(x, obj_value, NULL, numFunctionCalls, timeElapsed);

        pthread_mutex_lock(&ipoptMutex);

        return errors == 0;
    }


    virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
                             Number* grad_f)
    {
        static __thread int numFunctionCalls = 0;
        ++numFunctionCalls;
        // logmessage(LOGLVL_DEBUG, "Eval_Grad_F (%d) #%d", new_x, numFunctionCalls);

        pthread_mutex_unlock(&ipoptMutex);

        int errors = 0;

        clock_t timeBegin = clock();

        if(new_x || !lastCostP || !lastGradient) {
            errors = problem->evaluateObjectiveFunction(x, &lastCost, grad_f);

            if(lastGradient) delete[] lastGradient;
            lastGradient = new Number[problem->numOptimizationParameters];
            memcpy(lastGradient, grad_f, sizeof(Number) * n);
            lastErrors = errors;
        } else {
            memcpy(grad_f, lastGradient, sizeof(Number) * n);
            errors = lastErrors;
        }

        clock_t timeEnd = clock();
        double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;

        problem->logObjectiveFunctionEvaluation(x, lastCost, grad_f, numFunctionCalls, timeElapsed);

        pthread_mutex_lock(&ipoptMutex);

        return errors == 0;
    }


    virtual bool eval_g(Index n, const Number* x, bool new_x,
                        Index m, Number* g)
    {
        // no constraints, should never get here
        assert(false);
        return true;

    }

    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow,
                            Index *jCol, Number* values)
    {
        // no constraints, nothing to do here, but will be called once

        if(new_x)
            lastCostP = NULL; // because next function will be called with new_x==false, but we didn't prepare anything

        return true;
    }


    virtual bool intermediate_callback(AlgorithmMode mode,
                                       Index iter, Number obj_value,
                                       Number inf_pr, Number inf_du,
                                       Number mu, Number d_norm,
                                       Number regularization_size,
                                       Number alpha_du, Number alpha_pr,
                                       Index ls_trials,
                                       const IpoptData* ip_data,
                                       IpoptCalculatedQuantities* ip_cq)
    {
        int status = true;
        status = problem->intermediateFunction((int)mode,
                                                   iter, obj_value,
                                                   inf_pr,  inf_du,
                                                   mu, d_norm,
                                                   regularization_size,
                                                   alpha_du,  alpha_pr,
                                                   ls_trials);

    #ifdef INSTALL_SIGNAL_HANDLER
        if(caughtTerminationSignal) {
            logmessage(LOGLVL_CRITICAL, "CAUGHT SIGTERM... EXITING.");
            return false;
        }
    #endif

        return status == 0;

    }

    virtual void finalize_solution(SolverReturn status,
                                   Index n, const Number* x, const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value,
                                   const IpoptData* ip_data,
                                   IpoptCalculatedQuantities* ip_cq)
    {
        clock_t timeEnd = clock();
        double timeElapsed = (double) (timeEnd - timeBegin) / CLOCKS_PER_SEC;

        problem->logOptimizerFinished(obj_value, x, timeElapsed, status);
    }

    OptimizationProblem *problem;
    clock_t timeBegin;

    // for caching
    Number *lastGradient = NULL;
    Number lastCost = INFINITY;
    Number *lastCostP = &lastCost;
    int lastErrors = 0;
};



int getLocalOptimumIpopt(OptimizationProblem *problem) {
    pthread_mutex_lock(&ipoptMutex);

    assert(sizeof(double) == sizeof(Number));

    SmartPtr<TNLP> mynlp = new MyNLP(problem);
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();


    if(problem->optimizationOptions->printToStdout) {
        app->Options()->SetIntegerValue("print_level", 5);
        app->Options()->SetStringValue("print_user_options", "yes");
    } else {
        app->Options()->SetIntegerValue("print_level", 0);
        app->Options()->SetStringValue("print_user_options", "no");
        app->Options()->SetStringValue("sb", "yes"); // suppress copyright message
    }

    //    AddIpoptStrOption(nlp, "derivative_test", "first-order");
    //    AddIpoptIntOption(nlp, "derivative_test_first_index", 4130);

    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetStringValue("limited_memory_update_type", "bfgs");

    app->Options()->SetIntegerValue("max_iter", problem->optimizationOptions->maxOptimizerIterations);
    app->Options()->SetNumericValue("tol", 1e-9);

    app->Options()->SetIntegerValue("acceptable_iter", 1);
    app->Options()->SetNumericValue("acceptable_tol", 1e20); // set ridiculously high, so only the acceptable_* options below matter
    //AddIpoptNumOption(nlp, "acceptable_constr_viol_tol", 1e20);
    //AddIpoptNumOption(nlp, "acceptable_dual_inf_tol", 1e20);
    //AddIpoptNumOption(nlp, "acceptable_compl_inf_tol", 1e20);
    app->Options()->SetNumericValue("acceptable_obj_change_tol", problem->optimizationOptions->functionTolerance);

    // TODO check further limited memory options http://www.coin-or.org/Ipopt/documentation/node53.html#opt:hessian_approximation

    ApplicationReturnStatus status;
    status = app->Initialize();
    status = app->OptimizeTNLP(mynlp);

    pthread_mutex_unlock(&ipoptMutex);

    return (int)status < Maximum_Iterations_Exceeded;
}


