#include <parpecommon/functions.h>

namespace parpe {

FunctionEvaluationStatus GradientFunction::evaluate(
    gsl::span<double const> parameters,
    double& fval,
    gsl::span<double> gradient) const {
    return evaluate(parameters, fval, gradient, nullptr, nullptr);
}

} // namespace parpe
