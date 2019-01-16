#ifndef OPTIMIZATIONRESULTWRITER_H
#define OPTIMIZATIONRESULTWRITER_H

#include <string>

#include <hdf5.h>
#include <gsl/gsl-lite.hpp>

namespace parpe {

/**
 * @brief The OptimizationResultWriter class receives results during an
 * optimizer run. A new instance is to be created for each run.
 *
 * TODO: change into interface; add OptimizationResultWriterHDF5
 */

class OptimizationResultWriter {
public:
    /**
     * @brief Default constructor, for testing only
     */
    OptimizationResultWriter() = default;

    /**
     * @brief Write to pre-opened HDF5 file (will be re-opened)
     * @param problem
     * @param file_id
     */
    OptimizationResultWriter(hid_t file_id,
                             std::string rootPath);

    /**
     * @brief Open HDF5 file and write there
     * @param problem
     * @param filename Name of the result file
     * @param overwrite Overwrite output file if already exists
     */
    OptimizationResultWriter(const std::string &filename,
                             bool overwrite,
                             std::string rootPath);

    OptimizationResultWriter(OptimizationResultWriter const& other);
    
    virtual ~OptimizationResultWriter();

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
    virtual void logObjectiveFunctionEvaluation(gsl::span<const double> parameters,
                                                double objectiveFunctionValue,
                                                gsl::span<const double> objectiveFunctionGradient,
                                                int numIterations,
                                                int numFunctionCalls,
                                                double timeElapsedInSeconds);

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
    virtual void logOptimizerIteration(int numIterations,
                                       gsl::span<const double> parameters,
                                       double objectiveFunctionValue,
                                       gsl::span<const double> gradient,
                                       double wallSeconds,
                                       double cpuSeconds);
    
    void setLoggingEachIteration(bool logGradient);

    void setLoggingEachFunctionEvaluation(bool logGradient,
                                          bool logParameters);
    
    /*, int alg_mod, double inf_pr, double inf_du,
     double mu, double d_norm, double regularization_size, double alpha_du,
     double alpha_pr, int ls_trials*/

    virtual void starting(gsl::span<const double> initialParameters);

    /**
     * @brief Function to be called when local optimization is finished.
     * @param finalNegLogLikelihood Final cost f(x)
     * @param optimalParameters Final parameters x
     * @param masterTime Wall time for this optimization
     * @param exitStatus Exit status (cause of optimizer termination)
     */
    virtual void saveOptimizerResults(double finalNegLogLikelihood,
                                      gsl::span<const double> optimalParameters,
                                      double wallSec,
                                      double cpuSec,
                                      int exitStatus) const;

    hid_t getFileId() const;

    virtual std::string const& getRootPath() const;
    
    bool logParametersEachFunctionEvaluation = true;
    
    bool logGradientEachFunctionEvaluation = true;

    bool logGradientEachIteration = true;
    
    void setRootPath(std::string const& path);

protected:
    /**
     * @brief Write buffered output to file
     */
    virtual void flushResultWriter() const;

private:
    virtual std::string getIterationPath(int iterationIdx) const;

    hid_t file_id = 0;

    std::string rootPath = "/";

};

} // namespace parpe

#endif // OPTIMIZATIONRESULTWRITER_H
