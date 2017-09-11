#ifndef OPTIMIZATIONOPTIONS_H
#define OPTIMIZATIONOPTIONS_H

#include <hdf5Misc.h>

#include <string>
class Optimizer;

typedef enum optimizer_tag { OPTIMIZER_IPOPT, OPTIMIZER_CERES } optimizerEnum;

/** Type to describe an optimization (minimization) problem */

class OptimizationOptions {
  public:
    OptimizationOptions() = default;

    /** Optimizer factory method depending on OptimizationOptions::optimizer */
    Optimizer *createOptimizer();

    /** Optimizer to use */
    optimizerEnum optimizer = OPTIMIZER_IPOPT;

    /** Optimizer log file */
    char *logFile = nullptr;

    /** Print progress to stdout */
    bool printToStdout = true;

    /** Maximum number of optimizer iterations*/
    int maxOptimizerIterations = 1000;

    static OptimizationOptions *fromHDF5(const char *fileName);

    static OptimizationOptions *fromHDF5(hid_t fileId);

    static double *getStartingPoint(hid_t fileId, int index);

    /** Number of starts for local optimization (only used for multi-start
     * optimization */
    int numStarts = 1;

    /** Retry optimization in case of infeasibility (only used for multi-start
     * optimization */
    int retryOptimization = false;

    /** Convergence criterion for relative change in subsequent objective
     * function value change */
    double functionTolerance = 1e-18;

    /** see IpOpt */
    int watchdog_shortened_iter_trigger = 10;

    /** see IpOpt */
    bool accept_every_trial_step = false;

    std::string toString();
};

#endif // OPTIMIZATIONOPTIONS_H
