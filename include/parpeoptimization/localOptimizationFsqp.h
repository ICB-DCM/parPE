#ifndef LOCAL_OPTIMIZATION_FSQP_H
#define LOCAL_OPTIMIZATION_FSQP_H

#include <parpeoptimization/optimizationProblem.h>
#include <parpeoptimization/optimizer.h>

#include <vector>

namespace parpe {

/**
 * @brief Interface to the FSQP solver. (Tested with FFSQP Version 3.7b)
 *
 * This solver is not included in the parPE repository. A license must be obtained separately.
 */

class OptimizerFsqp : public Optimizer {
public:
    OptimizerFsqp() = default;

    std::tuple<int, double, std::vector<double> >
    optimize(parpe::OptimizationProblem *problem) override;
};


} // namespace parpe

#endif
