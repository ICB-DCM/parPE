#ifndef LOCAL_OPTIMIZATION_DLIB_H
#define LOCAL_OPTIMIZATION_DLIB_H

#include <parpeoptimization/optimizationProblem.h>
#include <parpeoptimization/optimizer.h>

namespace parpe {

class OptimizerDlibLineSearch : public Optimizer {
  public:
    OptimizerDlibLineSearch() = default;

    /**
     * @brief Minimize an objective function given as OptimizationProblem using
     * dlib line search algorithm
     *
     * TODO: no options are specifiable for the moment
     *
     * @param problem
     * @return Returns 0 on success.
     */

    std::tuple<int, double, std::vector<double> >
    optimize(OptimizationProblem *problem) override;
};

} // namespace parpe

#endif
