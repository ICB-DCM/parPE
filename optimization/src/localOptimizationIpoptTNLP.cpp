#include "localOptimizationIpoptTNLP.h"
#include "optimizationProblem.h"
#include <cassert>
#include <cstring>
#include <IpIpoptData.hpp>
#include <IpDenseVector.hpp>
//#include <IpIpoptCalculatedQuantities.hpp>

namespace parpe {


LocalOptimizationIpoptTNLP::LocalOptimizationIpoptTNLP(OptimizationProblem &problem, OptimizationReporter &reporter)
    : problem(problem), reporter(reporter)
{

}


bool LocalOptimizationIpoptTNLP::get_nlp_info(Index &n, Index &m,
                                              Index &nnz_jac_g,
                                              Index &nnz_h_lag,
                                              IndexStyleEnum &index_style) {

    n = reporter.numParameters();
    m = 0;                       // number of constrants
    nnz_jac_g = 0;               // numNonZeroElementsConstraintJacobian
    nnz_h_lag = 0;               // numNonZeroElementsLagrangianHessian
    index_style = TNLP::C_STYLE; // array layout for sparse matrices

    return true;
}

bool LocalOptimizationIpoptTNLP::get_bounds_info(Index n, Number *x_l,
                                                 Number *x_u, Index  /*m*/,
                                                 Number * /*g_l*/, Number * /*g_u*/) {
    // parameter bounds
    problem.fillParametersMin(gsl::make_span<double>(x_l, n));
    problem.fillParametersMax(gsl::make_span<double>(x_u, n));

    // no constraints supported for now -> no constraint bounds

    return true;
}

bool LocalOptimizationIpoptTNLP::get_starting_point(Index n, bool init_x,
                                                    Number *x, bool init_z,
                                                    Number * /*z_L*/, Number * /*z_U*/,
                                                    Index  /*m*/, bool init_lambda,
                                                    Number * /*lambda*/) {
    /* this function is called twice by IpOpt which is a problem if problem->fillInitialParameters provides random parameters, therefore initial point needs to be stored */

    if (init_x) {
        if(initialParameters.empty()) {
            initialParameters.resize(n);
            problem.fillInitialParameters(initialParameters);
            if(reporter.starting(initialParameters))
                return false;
        }
        std::copy(initialParameters.begin(), initialParameters.end(), x);
    }

    assert(init_z == false);
    assert(init_lambda == false);

    return true;
}

bool LocalOptimizationIpoptTNLP::eval_f(Index n, const Number *x, bool  /*new_x*/,
                                        Number &obj_value) {
    auto unlockIpOpt = ipOptReleaseLock();

    return reporter.evaluate(gsl::make_span<double const>(x, n), obj_value, gsl::span<double>()) == functionEvaluationSuccess;
}

bool LocalOptimizationIpoptTNLP::eval_grad_f(Index n, const Number *x,
                                             bool  /*new_x*/, Number *grad_f) {
    auto unlockIpOpt = ipOptReleaseLock();

    double obj_value;
    return reporter.evaluate(gsl::make_span<double const>(x, n), obj_value, gsl::make_span<double>(grad_f, n)) == functionEvaluationSuccess;
}

bool LocalOptimizationIpoptTNLP::eval_g(Index  /*n*/, const Number * /*x*/, bool  /*new_x*/,
                                        Index  /*m*/, Number * /*g*/) {

    assert(false && "no constraints, should never get here");

    return false;
}

bool LocalOptimizationIpoptTNLP::eval_jac_g(Index  /*n*/, const Number * /*x*/,
                                            bool  /*new_x*/, Index  /*m*/, Index  /*nele_jac*/,
                                            Index * /*iRow*/, Index * /*jCol*/,
                                            Number * /*values*/) {
    // no constraints, nothing to do here, but will be called once
    return true;
}

bool LocalOptimizationIpoptTNLP::intermediate_callback(
    AlgorithmMode  /*mode*/, Index  /*iter*/, Number obj_value, Number  /*inf_pr*/,
    Number  /*inf_du*/, Number  /*mu*/, Number  /*d_norm*/, Number  /*regularization_size*/,
    Number  /*alpha_du*/, Number  /*alpha_pr*/, Index  /*ls_trials*/, const IpoptData *ip_data,
    IpoptCalculatedQuantities * /*ip_cq*/) {

    auto unlockIpOpt = ipOptReleaseLock();

    // get current parameters from IpOpt which are not available directly
    gsl::span<double const> parameters;
    auto x = ip_data->curr()->x();
    auto xx = dynamic_cast<const Ipopt::DenseVector*>(Ipopt::GetRawPtr(x));
    if(xx)
        parameters = gsl::span<double const>(xx->Values(), xx->Dim());
    else
        logmessage(LOGLVL_WARNING, "Not Ipopt::DenseVector in LocalOptimizationIpoptTNLP::intermediate_callback");

    // is always the last step accepted?
    int status = reporter.iterationFinished(parameters, obj_value, gsl::span<double>());

#ifdef INSTALL_SIGNAL_HANDLER
    if (caughtTerminationSignal) {
        logmessage(LOGLVL_CRITICAL, "CAUGHT SIGTERM... EXITING.");
        return false;
    }
#endif

    return status == 0;
}

void LocalOptimizationIpoptTNLP::finalize_solution(
    SolverReturn status, Index n, const Number *x, const Number * /*z_L*/,
    const Number * /*z_U*/, Index  /*m*/, const Number * /*g*/, const Number * /*lambda*/,
    Number obj_value, const IpoptData * /*ip_data*/,
    IpoptCalculatedQuantities * /*ip_cq*/) {

    auto unlockIpOpt = ipOptReleaseLock();

    reporter.finished(obj_value, gsl::span<double const>(x, n), status);
}


std::unique_lock<mutexIpOptType> ipOptGetLock()
{
    return std::unique_lock<mutexIpOptType>(mutexIpOpt);
}

InverseUniqueLock<mutexIpOptType> ipOptReleaseLock()
{
    return InverseUniqueLock<mutexIpOptType>(&mutexIpOpt);
}

} // namespace parpe
