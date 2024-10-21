#ifndef OPTIMIZATIONOPTIONS_H
#define OPTIMIZATIONOPTIONS_H

#include <parpecommon/hdf5Misc.h>

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace parpe {

class Optimizer;

enum class optimizerName {
    OPTIMIZER_IPOPT,
    OPTIMIZER_CERES,
    OPTIMIZER_DLIB,
    OPTIMIZER_TOMS611,
    OPTIMIZER_FSQP,
    OPTIMIZER_FIDES,
    OPTIMIZER_MINIBATCH_1 = 10
};

/** Type to describe an optimization (minimization) problem */

class OptimizationOptions {
  public:
    OptimizationOptions() = default;

    /** Optimizer factory method depending on OptimizationOptions::optimizer */
    std::unique_ptr<Optimizer> createOptimizer() const;

    /** Optimizer to use */
    optimizerName optimizer = optimizerName::OPTIMIZER_IPOPT;

    /** Optimizer log file */
    char* logFile = nullptr;

    /** Print progress to stdout */
    bool printToStdout = true;

    /** Maximum number of optimizer iterations*/
    int maxOptimizerIterations = 100;

    static std::unique_ptr<OptimizationOptions>
    fromHDF5(std::string const& fileName);
    static std::unique_ptr<OptimizationOptions> fromHDF5(
        const H5::H5File& file,
        std::string const& path = "/optimizationOptions");

    static std::vector<double>
    getStartingPoint(const H5::H5File& file, int index);

    /** Number of starts for local optimization (only used for multi-start
     * optimization */
    int numStarts = 1;

    /** Retry optimization in case of infeasibility (only used for multi-start
     * optimization */
    int retryOptimization = false;

    /** use hierarchical optimization if respective configuration is available
     * see hierarchicalOptimization.cpp */
    int hierarchicalOptimization = true;

    int multistartsInParallel = true;

    std::string toString();

    int getIntOption(std::string const& key);
    double getDoubleOption(std::string const& key);
    std::string getStringOption(std::string const& key);

    void setOption(std::string const& key, int value);
    void setOption(std::string const& key, double value);
    void setOption(std::string const& key, std::string value);

    template <typename T>
    void for_each(
        std::function<
            void(std::pair<std::string const, std::string const> const, T)> f,
        T arg) const {
        std::for_each(
            options.cbegin(),
            options.cend(),
            std::bind(f, std::placeholders::_1, arg));
    }

  private:
    std::map<std::string, std::string> options;
};

std::unique_ptr<Optimizer> optimizerFactory(optimizerName optimizer);

/**
 * @brief Print list of supported optimizers
 */
void printAvailableOptimizers(std::string const& prefix = "");
} // namespace parpe

#endif // OPTIMIZATIONOPTIONS_H
