#include <parpeoptimization/localOptimizationFides.h>

#include "parpecommon/logging.h"
#include "parpeoptimization/optimizationOptions.h"
#include "parpeoptimization/optimizationProblem.h"

#include <fides/minimize.hpp>
#include <blaze/math/CustomVector.h>

using blaze::unaligned;
using blaze::unpadded;
using blaze::CustomVector;
using blaze::columnVector;
using blaze::DynamicVector;
using UnalignedUnpadded = CustomVector<int,unaligned,unpadded,columnVector>;

namespace gsl {
template< typename Type, bool TF, typename Alloc, typename Tag >
auto make_span(DynamicVector<Type,TF,Alloc,Tag> & dv) {
    return gsl::span<Type>(dv.data(), dv.size());
}

template< typename Type, bool TF, typename Alloc, typename Tag >
auto make_span(DynamicVector<Type,TF,Alloc,Tag> const& dv) {
    return gsl::span<const Type>(dv.data(), dv.size());
}

}
namespace parpe {



std::tuple<int, double, std::vector<double>> OptimizerFides::optimize(OptimizationProblem *problem)
{
    int numParams = problem->cost_fun_->numParameters();

    DynamicVector<double> x0(static_cast<std::size_t>(numParams));
    problem->fillInitialParameters(gsl::make_span(x0));

    DynamicVector<double> lb(static_cast<std::size_t>(numParams));
    problem->fillParametersMin(gsl::make_span(lb));

    DynamicVector<double> ub(static_cast<std::size_t>(numParams));
    problem->fillParametersMax(gsl::make_span(ub));

    WallTimer optimization_timer;

    auto fun = [&problem](DynamicVector<double> x) {
        static __thread int numFunctionCalls = 0;
        ++numFunctionCalls;
        DynamicVector<double> g(x.size(), NAN);
        double fval = NAN;
        problem->cost_fun_->evaluate(gsl::make_span(x), fval,
                                     gsl::make_span(g), nullptr, nullptr);

        return std::make_tuple(fval, g, blaze::DynamicMatrix<double>());
    };

    auto parpe_options = problem->getOptimizationOptions();
    fides::Options fides_options;
    fides_options.maxiter = parpe_options.maxOptimizerIterations;

    fides::BFGS hessian_approximation;
    fides::Optimizer optimizer(fun, lb, ub, 10, fides_options,
                               &hessian_approximation);
    auto [fval, x, grad, hess] = optimizer.minimize(x0);
    return std::make_tuple(
        static_cast<int>(optimizer.exit_flag_) > 0,
        fval, std::vector<double>(x.begin(), x.end()));
}

} // namespace parpe
