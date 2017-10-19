#ifndef OPTIMIZATIONRESULTWRITER_H
#define OPTIMIZATIONRESULTWRITER_H

#include "optimizationProblem.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <string>

namespace parpe {

class OptimizationProblem;

/**
 * @brief The OptimizationResultWriter class receives results during an
 * optimizer run
 */
class OptimizationResultWriter {
  public:
    /**
     * @brief Default constructor, for testing only
     */
    OptimizationResultWriter();

    /**
     * @brief Write to pre-opened HDF5 file
     * @param problem
     * @param file_id
     */
    OptimizationResultWriter(hid_t file_id);

    /**
     * @brief Open HDF5 file and write there
     * @param problem
     * @param filename Name of the result file
     * @param overwrite Overwrite output file if already exists
     */
    OptimizationResultWriter(const std::string &filename, bool overwrite);

    virtual std::string getOptimizationPath();

    virtual std::string getIterationPath(int iterationIdx);

    virtual void logParPEVersion();

    /**
     * @brief Function to be called after each objective function f(x) or
     * gradient f'(x) evaluation
     * @param parameters Function parameters x
     * @param objectiveFunctionValue  f(x)
     * @param objectiveFunctionGradient f'(x) or NULL
     * @param numFunctionCalls Number of times the objective function has been
     * called (f(x) and f'(x) are counter individually (?))
     * @param timeElapsedInSeconds CPU time for the last objective function
     * evaluation (wall time)
     */
    virtual void logLocalOptimizerObjectiveFunctionEvaluation(
        const double *parameters, int numParameters,
        double objectiveFunctionValue, const double *objectiveFunctionGradient,
        int numFunctionCalls, double timeElapsedInSeconds);

    /**
     * @brief Function to be called after each optimizer iteration. (For
     * parameters, see above or IpOpt intermediate function)
     * @param numIterations
     * @param theta
     * @param objectiveFunctionValue
     * @param gradient
     * @param timeElapsedInSeconds
     * @param alg_mod
     * @param inf_pr
     * @param inf_du
     * @param mu
     * @param d_norm
     * @param regularization_size
     * @param alpha_du
     * @param alpha_pr
     * @param ls_trials
     */
    virtual void logLocalOptimizerIteration(
        int numIterations, double *theta, int numParameters,
        double objectiveFunctionValue, const double *gradient,
        double timeElapsedInSeconds, int alg_mod, double inf_pr, double inf_du,
        double mu, double d_norm, double regularization_size, double alpha_du,
        double alpha_pr, int ls_trials);

    /**
     * @brief Write buffered output to file
     */
    virtual void flushResultWriter();

    /**
     * @brief CPU time for whole application run
     * @param timeInSeconds
     */
    virtual void saveTotalCpuTime(const double timeInSeconds);

    /**
     * @brief Function to be called when local optimization is finished.
     * @param finalNegLogLikelihood Final cost f(x)
     * @param optimalParameters Final parameters x
     * @param masterTime Wall time for this optimization
     * @param exitStatus Exit status (cause of optimizer termination)
     */
    virtual void saveLocalOptimizerResults(double finalNegLogLikelihood,
                                           const double *optimalParameters,
                                           int numParameters, double masterTime,
                                           int exitStatus);

    void setRootPath(std::string const& path);

    virtual ~OptimizationResultWriter();

    hid_t file_id;

  protected:
    int initResultHDFFile(const char *filename, bool overwrite);

    void closeResultHDFFile();

private:
    std::string rootPath = "/";

};

} // namespace parpe

#endif // OPTIMIZATIONRESULTWRITER_H
