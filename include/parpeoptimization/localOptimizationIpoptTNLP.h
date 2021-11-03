#ifndef LOCALOPTIMIZATIONIPOPTTNLP_H
#define LOCALOPTIMIZATIONIPOPTTNLP_H

// This should fix `error "don't have header file for stddef"' with some IpOpt
// versions
#ifndef HAVE_CSTDDEF
#define HAVE_CSTDDEF
#include <IpIpoptApplication.hpp>
#undef HAVE_CSTDDEF
#else
#include <IpIpoptApplication.hpp>
#endif

#include <cmath>
#include <memory>
#include <mutex>

#include <parpecommon/misc.h>

namespace parpe {

/** Mutex for managing access to IpOpt routines which are not thread-safe */
using mutexIpOptType = std::recursive_mutex;

/**
 * @brief ipoptMutex Ipopt seems not to be thread safe. Lock this mutex every
 * time that control is passed to ipopt functions.
 */
static mutexIpOptType mutexIpOpt {};


InverseUniqueLock<mutexIpOptType> ipOptReleaseLock();

std::unique_lock<mutexIpOptType> ipOptGetLock();


class OptimizationProblem;
class OptimizationReporter;

class LocalOptimizationIpoptTNLP : public Ipopt::TNLP {
  public:

    LocalOptimizationIpoptTNLP(OptimizationProblem &problem, OptimizationReporter &reporter);

    ~LocalOptimizationIpoptTNLP() override = default;

    bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m,
                              Ipopt::Index &nnz_jac_g,
                              Ipopt::Index &nnz_h_lag,
                              IndexStyleEnum &index_style) override;

    bool get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l,
                                 Ipopt::Number *x_u, Ipopt::Index m,
                                 Ipopt::Number *g_l, Ipopt::Number *g_u) override;

    bool get_starting_point(Ipopt::Index n, bool init_x,
                                    Ipopt::Number *x,
                                    bool init_z, Ipopt::Number *z_L,
                                    Ipopt::Number *z_U,
                                    Ipopt::Index m, bool init_lambda,
                                    Ipopt::Number *lambda) override;

    bool eval_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
                        Ipopt::Number &obj_value) override;

    /**
     * @brief See Ipopt::TNLP::eval_grad_f.
     *
     * Note: Failure in eval_f (i.e. returning non-finite value or false) will make IpOpt try a new
     * step, unless this happens at the starting point. However, if eval_f succeeds, but eval_grad_f
     * fails, then IpOpt will terminate.
     *
     * @param n
     * @param x
     * @param new_x
     * @param grad_f
     * @return
     */
    bool eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
                             Ipopt::Number *grad_f) override;

    bool eval_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
                        Ipopt::Index m,
                        Ipopt::Number *g) override;

    bool eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
                            Ipopt::Index m,
                            Ipopt::Index nele_jac, Ipopt::Index *iRow,
                            Ipopt::Index *jCol,
                            Ipopt::Number *values) override;

    bool intermediate_callback(
        Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
        Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
        Ipopt::Number d_norm, Ipopt::Number regularization_size,
        Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials,
        const Ipopt::IpoptData *ip_data, Ipopt::IpoptCalculatedQuantities *ip_cq) override;

    void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                                   const Ipopt::Number *x,
                                   const Ipopt::Number *z_L,
                                   const Ipopt::Number *z_U, Ipopt::Index m,
                                   const Ipopt::Number *g,
                                   const Ipopt::Number *lambda,
                                   Ipopt::Number obj_value,
                                   const Ipopt::IpoptData *ip_data,
                                   Ipopt::IpoptCalculatedQuantities *ip_cq) override;
private:
    OptimizationProblem &problem;
    OptimizationReporter &reporter;

    // need to store initial parameters, because IpOpt asks twice
    std::vector<double> initialParameters;

};

void setIpOptOption(const std::pair<const std::string, const std::string> &pair,
                    Ipopt::SmartPtr<Ipopt::OptionsList>* o);

} // namespace parpe

#endif // LOCALOPTIMIZATIONIPOPTTNLP_H
