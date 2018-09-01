#ifndef PARPE_OPTIMIZATION_MINIBATCH_OPTIMIZATION_H
#define PARPE_OPTIMIZATION_MINIBATCH_OPTIMIZATION_H

#include <optimizationOptions.h>
#include <optimizationProblem.h>
#include <misc.h>
#include <model.h>

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include <sstream>

namespace parpe {

/**
 * @brief Return status for minibatch optimizer
 */
enum class minibatchExitStatus {
    gradientNormConvergence,
    maxEpochsExceeded,
    invalidNumber
};


/**
 * Problem definition for minibatch optimization
 */
template<typename T>
class MinibatchOptimizationProblem : public OptimizationProblem {
public:
    MinibatchOptimizationProblem() = default;
    MinibatchOptimizationProblem(std::unique_ptr<SummedGradientFunction<T>> costFun,
                                 std::unique_ptr<Logger> logger)
        : OptimizationProblem(costFun, logger)
    {
    }

    MinibatchOptimizationProblem(MinibatchOptimizationProblem const& other) = delete;

    virtual ~MinibatchOptimizationProblem() override = default;

    virtual std::vector<T> getTrainingData() const = 0;

    SummedGradientFunction<T>* getGradientFunction() const {
        auto summedGradientFunction = dynamic_cast<SummedGradientFunction<T>*>(costFun.get());
        RELEASE_ASSERT(summedGradientFunction, "");
        return summedGradientFunction;
    }
};


/**
 * @brief Split a vector into batches of the given size.
 *
 * All but possibly the last one will be of Size batchSize.
 *
 * @param data Data to be split into batches
 * @param batchSize Number of elements in each batch
 * @return Vector batches of data elements
 */
template<typename T>
std::vector<std::vector<T>> getBatches(gsl::span<const T> data, int batchSize) {
    int numBatches = ceil(static_cast<double>(data.size()) / batchSize);

    std::vector<std::vector<T>> batches(numBatches, std::vector<T>());
    for(int i = 0, batchIdx = -1; static_cast<typename decltype (data)::index_type>(i) < data.size(); ++i) {
        if(i % batchSize == 0) {
            ++batchIdx;
            int remaining = data.size() - i;
            batches[batchIdx].resize(std::min(batchSize, remaining));
        }
        batches[batchIdx][i % batchSize] = data[i];
    }

    return batches;
}


/**
 * @brief Get Euclidean (l2) norm of vector.
 * @param v
 * @return the norm
 */
double getVectorNorm(gsl::span<const double> v);

/**
 * @brief Interface for parameter updaters for minibatch optimizers
 */
class ParameterUpdater {
public:
    /**
     * @brief Update parameter vector
     * @param gradient Cost function gradient at parameters
     * @param parameters In: Current parameters, Out: Updated parameters
     */
    virtual void updateParameters(gsl::span<const double> gradient,
                                  gsl::span<double> parameters,
                                  gsl::span<const double> lowerBounds = gsl::span<const double>(),
                                  gsl::span<const double> upperBounds = gsl::span<const double>()) = 0;

    virtual ~ParameterUpdater() = default;
};


class ParameterUpdaterVanilla : public ParameterUpdater {
public:
    ParameterUpdaterVanilla() = default;
    ParameterUpdaterVanilla(double learningRate);

    void updateParameters(gsl::span<const double> gradient, gsl::span<double> parameters,
                          gsl::span<const double> lowerBounds = gsl::span<const double>(),
                          gsl::span<const double> upperBounds = gsl::span<const double>());

    double learningRate = 0.1;
};


class ParameterUpdaterRmsProp : public ParameterUpdater {
public:
    ParameterUpdaterRmsProp() = default;
    ParameterUpdaterRmsProp(double learningRate);

    void updateParameters(gsl::span<const double> gradient, gsl::span<double> parameters,
                          gsl::span<const double> lowerBounds = gsl::span<const double>(),
                          gsl::span<const double> upperBounds = gsl::span<const double>());

    int updates = 0;
    double learningRate = 1.0;
    double decayRate = 0.9;
    double eps = 1e-4;
    std::vector<double> gradientCache;
};


template<typename BATCH_ELEMENT>
class MinibatchOptimizer {
public:
    /**
     * @brief Minimize the given function using mini-batch gradient descent.
     *
     * @param f Function to minize
     * @param data Full data set on which f will be evaluated
     * @param initialParameters Starting point for optimization
     * @param reporter OptimizationReporter instance for tracking progress
     * @param logger Logger instance for status messages
     * @return Tuple (exit code, final cost, final parameters)
     */
    std::tuple<int, double, std::vector<double> > optimize(
            SummedGradientFunction<BATCH_ELEMENT> const& f,
            gsl::span<const BATCH_ELEMENT> data,
            gsl::span<const double> initialParameters,
            gsl::span<const double> lowerParameterBounds,
            gsl::span<const double> upperParameterBounds,
            OptimizationReporter *reporter,
            Logger *logger_)
    {
        RELEASE_ASSERT((unsigned) f.numParameters() == initialParameters.size(), "");
        Logger logger = logger_ ? *logger_ : Logger();

        // We don't change the user inputs but work with copies
        std::vector<double> parameters(initialParameters.begin(), initialParameters.end());
        std::vector<BATCH_ELEMENT> shuffledData(data.begin(), data.end());

        std::vector<double> gradient(parameters.size(), 0.0);
        std::random_device rd;
        std::mt19937 rng(rd());
        minibatchExitStatus status = minibatchExitStatus::gradientNormConvergence;
        double cost = NAN;

        if(reporter) reporter->starting(initialParameters);

        for(int epoch = 0; epoch < maxEpochs; ++epoch) {
            auto epochLogger = logger.getChild(std::string("e") + std::to_string(epoch));

            // Create randomized batches
            std::shuffle(shuffledData.begin(), shuffledData.end(), rng);
            auto batches = getBatches<BATCH_ELEMENT>(shuffledData, batchSize);

            for(int batchIdx = 0; (unsigned) batchIdx < batches.size(); ++batchIdx) {
                auto batchLogger = epochLogger->getChild(std::string("b") + std::to_string(batchIdx));

                if(reporter) {
                    reporter->beforeCostFunctionCall(parameters);
                    reporter->logger->setPrefix(batchLogger->getPrefix());
                }
                double cpuTime = 0.0;
                auto status = f.evaluate(parameters, batches[batchIdx], cost, gradient, batchLogger.get(), &cpuTime);
                if(reporter) {
                    reporter->cpuTimeIterationSec += cpuTime;
                    reporter->cpuTimeTotalSec += cpuTime;
                    reporter->afterCostFunctionCall(parameters, cost, gradient);
                }

                std::stringstream ss;
                ss<<": p: "<<parameters<<" Cost: "<<cost<<" Gradient:"<<gradient
                 << " |g|2: "<<getVectorNorm(gradient)
                 << " Batch: "<<batches[batchIdx]<<std::endl;
                batchLogger->logmessage(LOGLVL_DEBUG, ss.str().c_str());

                if(status == functionEvaluationFailure) {
                    // TODO: do something smarter
                    batchLogger->logmessage(LOGLVL_ERROR, "Minibatch cost function evaluation failed.");
                    reporter->finished(cost, parameters, (int)minibatchExitStatus::invalidNumber);
                    return std::tuple<int, double, std::vector<double> >((int)minibatchExitStatus::invalidNumber, cost, parameters);
                }

                if(reporter) reporter->iterationFinished(parameters, cost, gradient);

                parameterUpdater->updateParameters(gradient, parameters,
                                                   lowerParameterBounds, upperParameterBounds);
            }

            if(getVectorNorm(gradient) <= gradientNormThreshold) {
                epochLogger->logmessage(LOGLVL_INFO, "Convergence: gradientNormThreshold reached.");
                break;
            }
        }

        reporter->finished(cost, parameters, (int)status);
        return std::tuple<int, double, std::vector<double> >((int)status, cost, parameters);
    }

    std::unique_ptr<ParameterUpdater> parameterUpdater = std::make_unique<ParameterUpdaterVanilla>();

    int batchSize = 1;
    int maxEpochs = 10;
    double gradientNormThreshold = 0.0;
};


/**
 * @brief Apply the given key,value option to optimizer
 * @param pair
 * @param optimizer
 */
void setMinibatchOption(const std::pair<const std::string, const std::string> &pair,
                        MinibatchOptimizer<int>* optimizer);


/**
 * @brief Create and setup a minibatch optimizer according to the given options
 * @param options
 * @return
 */
template<typename BATCH_ELEMENT>
std::unique_ptr<MinibatchOptimizer<BATCH_ELEMENT>> getMinibatchOptimizer(OptimizationOptions const& options)
{
    auto optim = std::make_unique<MinibatchOptimizer<BATCH_ELEMENT>>();

    options.for_each<MinibatchOptimizer<BATCH_ELEMENT>*>(setMinibatchOption, optim.get());

    return optim;
}

std::tuple<int, double, std::vector<double> > runMinibatchOptimization(MinibatchOptimizationProblem<int> *problem);

} // namespace parpe

/**
 * @brief Clip values to given element-wise bounds.
 * @param lowerBounds
 * @param upperBounds
 * @param x
 */
template<typename T>
void clipToBounds(gsl::span<const T> lowerBounds, gsl::span<const T> upperBounds, gsl::span<T> x) {
    if(lowerBounds.empty() && upperBounds.empty())
        return;

    RELEASE_ASSERT(lowerBounds.size() == upperBounds.size(), "");
    RELEASE_ASSERT(lowerBounds.size() == x.size(), "");

    for(int i = 0; static_cast<typename gsl::span<const T>::index_type>(i) < x.size(); ++i)
        x[i] = std::min(std::max(lowerBounds[i], x[i]), upperBounds[i]);
}

#endif
