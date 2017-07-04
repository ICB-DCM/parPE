#ifndef OPTIMIZATIONOPTIONS_H
#define OPTIMIZATIONOPTIONS_H

#include <hdf5Misc.h>

typedef enum optimizer_tag {OPTIMIZER_IPOPT, OPTIMIZER_CERES} optimizerEnum;

/** Type to describe an optimization (minimization) problem */

class OptimizationOptions
{
public:
    OptimizationOptions();

    /** Optimizer to use */
    optimizerEnum optimizer;

    /** Optimizer log file */
    char *logFile;

    /** Print progress to stdout */
    bool printToStdout;

    /** Maximum number of optimizer iterations*/
    int maxOptimizerIterations;

    static OptimizationOptions* fromHDF5(const char* fileName);

    static OptimizationOptions* fromHDF5(hid_t fileId);

    /** Number of starts for local optimization (only used for multi-start optimization */
    int numStarts;

    /** Retry optimization in case of infeasibility (only used for multi-start optimization */
    int retryOptimization;

    /** Convergence criterium for relative change in subsequent objective function value change */
    double functionTolerance;
};

#endif // OPTIMIZATIONOPTIONS_H
