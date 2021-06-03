#include <parpeoptimization/localOptimizationFides.h>

#include "parpecommon/logging.h"
#include "parpeoptimization/optimizationOptions.h"
#include "parpeoptimization/optimizationProblem.h"

#include <blaze/math/CustomVector.h>
#include <fides/minimize.hpp>

using blaze::columnVector;
using blaze::CustomVector;
using blaze::DynamicVector;
using blaze::unaligned;
using blaze::unpadded;
using UnalignedUnpadded = CustomVector<int, unaligned, unpadded, columnVector>;

namespace gsl {

template<typename Type, bool TF, typename Alloc, typename Tag>
auto
make_span(DynamicVector<Type, TF, Alloc, Tag>& dv)
{
    return gsl::span<Type>(dv.data(), dv.size());
}

template<typename Type, bool TF, typename Alloc, typename Tag>
auto
make_span(DynamicVector<Type, TF, Alloc, Tag> const& dv)
{
    return gsl::span<const Type>(dv.data(), dv.size());
}

} // namespace gsl

namespace parpe {

fides::Options
get_optimization_options(OptimizationOptions const& parpe_options)
{
    fides::Options fides_options;
    fides_options.maxiter = parpe_options.maxOptimizerIterations;

    parpe_options.for_each<int>(
      [&fides_options](const std::pair<const std::string, const std::string>& pair, int) {
          const std::string& key = pair.first;
          const std::string& val = pair.second;
          auto &options = fides_options;
          if (key == "maxiter") {
              options.maxiter = std::stoi(val);
          } else if (key == "maxtime") {
              options.maxtime = std::chrono::seconds(std::stoi(val));
          } else if (key == "fatol") {
              options.fatol = std::stod(val);
          } else if (key == "frtol") {
              options.frtol = std::stod(val);
          } else if (key == "xtol") {
              options.xtol = std::stod(val);
          } else if (key == "gatol") {
              options.gatol = std::stod(val);
          } else if (key == "grtol") {
              options.grtol = std::stod(val);
          } else if (key == "subspace_solver") {
              auto result = std::find_if(
                fides::subspace_dim_to_str.cbegin(),
                fides::subspace_dim_to_str.cend(),
                [key](const auto& kv) { return kv.second == key; });
              if (result != fides::subspace_dim_to_str.cend())
                  options.subspace_solver = result->first;
              else
                  logmessage(LOGLVL_WARNING,
                             "Invalid value %s provided for option "
                             "'subspace_solver'. Ignoring.",
                             val.c_str());
          } else if (key == "stepback_strategy") {
              auto result = std::find_if(
                fides::step_back_strategy_str.cbegin(),
                fides::step_back_strategy_str.cend(),
                [key](const auto& kv) { return kv.second == key; });
              if (result != fides::step_back_strategy_str.cend())
                  options.stepback_strategy = result->first;
              else
                  logmessage(LOGLVL_WARNING,
                             "Invalid value %s provided for option "
                             "'stepback_strategy'. Ignoring.",
                             val.c_str());
          } else if (key == "theta_max") {
              options.theta_max = std::stoi(val);
          } else if (key == "delta_init") {
              options.delta_init = std::stoi(val);
          } else if (key == "mu") {
              options.mu = std::stoi(val);
          } else if (key == "eta") {
              options.eta = std::stoi(val);
          } else if (key == "gamma1") {
              options.gamma1 = std::stoi(val);
          } else if (key == "gamma2") {
              options.gamma2 = std::stoi(val);
          } else if (key == "refine_stepback") {
              options.refine_stepback = std::stoi(val);
          } else {
              logmessage(LOGLVL_WARNING,
                         "Ignoring unknown optimization option %s.",
                         key.c_str());
              return;
          }

          logmessage(LOGLVL_DEBUG,
                     "Set optimization option %s to %s.",
                     key.c_str(),
                     val.c_str());
      },
      0);

    return fides_options;
}

std::tuple<int, double, std::vector<double>>
OptimizerFides::optimize(OptimizationProblem* problem)
{
    auto reporter = problem->getReporter();
    auto numParams =
      static_cast<std::size_t>(problem->cost_fun_->numParameters());

    DynamicVector<double> x0(numParams);
    problem->fillInitialParameters(x0);

    DynamicVector<double> lb(numParams);
    problem->fillParametersMin(lb);

    DynamicVector<double> ub(numParams);
    problem->fillParametersMax(ub);

    WallTimer optimization_timer;

    auto fun = [&problem](DynamicVector<double> x) {
        static __thread int numFunctionCalls = 0;
        ++numFunctionCalls;
        DynamicVector<double> g(x.size(), NAN);
        double fval = NAN;
        problem->cost_fun_->evaluate(
          gsl::make_span(x), fval, gsl::make_span(g), nullptr, nullptr);

        return std::make_tuple(fval, g, blaze::DynamicMatrix<double>());
    };

    auto fides_options =
      get_optimization_options(problem->getOptimizationOptions());
    // TODO to config
    fides::BFGS hessian_approximation;
    fides::Optimizer optimizer(
      fun, lb, ub, fides_options, &hessian_approximation);

    reporter->starting(x0);
    auto [fval, x, grad, hess] = optimizer.minimize(x0);
    reporter->finished(fval, x, static_cast<int>(optimizer.exit_flag_));

    return std::make_tuple(static_cast<int>(optimizer.exit_flag_) <= 0,
                           fval,
                           std::vector<double>(x.begin(), x.end()));
}

} // namespace parpe
