#ifndef OPTIMIZATIONRESULTWRITER_H
#define OPTIMIZATIONRESULTWRITER_H

#include <string>

#include <H5Cpp.h>
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
     * @brief Write to pre-opened HDF5 file (will be re-opened)
     * @param file
     * @param rootPath
     */
    OptimizationResultWriter(const H5::H5File &file,
                             std::string rootPath);

    /**
     * @brief Open HDF5 file and write there
     * @param filename Name of the result file
     * @param overwrite Overwrite output file if already exists
     * @param rootPath
     */
    OptimizationResultWriter(const std::string &filename,
                             bool overwrite,
                             std::string rootPath);

    OptimizationResultWriter& operator=(const OptimizationResultWriter& other) = delete;

    OptimizationResultWriter(OptimizationResultWriter const& other);

    virtual ~OptimizationResultWriter();

    /**
     * @brief Function to be called after each objective function f(x) or
     * gradient f'(x) evaluation
     * @param parameters Function parameters x
     * @param objectiveFunctionValue  f(x)
     * @param objectiveFunctionGradient f'(x) or NULL
     * @param numFunctionCalls Number of times the objective function has been
     * called (f(x) and f'(x) are counted individually (?))
     * @param timeElapsedInSeconds CPU time for the last objective function
     * evaluation (wall time)
     */
    virtual void logObjectiveFunctionEvaluation(
            gsl::span<const double> parameters,
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

    /**
     * @brief Log optimizer start
     * @param initialParameters
     */
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

    H5::H5File const& getH5File() const;

    virtual std::string const& getRootPath() const;

    bool logParametersEachFunctionEvaluation = true;

    bool logGradientEachFunctionEvaluation = true;

    bool logGradientEachIteration = true;

    /**
     * @brief Set root path in HDF5 file and create the respective group.
     * @param path
     */
    void setRootPath(std::string const& path);

protected:
    /**
     * @brief Write buffered output to file
     */
    virtual void flushResultWriter() const;

private:
    virtual std::string getIterationPath(int iterationIdx) const;

    H5::H5File file = 0;

    /** Root path within HDF5 file */
    std::string rootPath = "/";

};

} // namespace parpe

#endif // OPTIMIZATIONRESULTWRITER_H
