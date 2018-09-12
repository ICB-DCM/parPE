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

enum class learningRateInterp {
    linear,
    inverseLinear,
    logarithmic
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
 * @brief learning rate updaters for minibatch optimizers
 */
class LearningRateUpdater {
public:
    /**
     * @brief Update the learning rate
     * @param startLearningRate Learning rate the the beginning of optimization
	 * @param endLearningRate Learning rate the the end of optimization
     * @param LearningRateadaptionMode Type of interpolation between startLearningRate and endLearningRate 
     */

    LearningRateUpdater(int maxEpochs, learningRateInterp learningRateInterpMode){};
    // LearningRateUpdater(int maxEpochs, learningRateInterp learningRateInterpMode, double startLearningRate, double endLearningRate);

	// Update function, to be called in every epoch or optimization iteration
    void updateLearningRate(int iteration);

	// Update function, to be called if parameter update did not work well
    void reduceLearningRate();

	// Update function, to be called if parameter update worked well
    void increaseLearningRate();

	// Function to retrieve the learning rate
	double getCurrentLearningRate();

	// Function to set the reduction factor directly
	void setReductionFactor(double newReductionFactor);

	// Function to set the new maximum epoch number
	void setMaxEpochs(int newMaxEpochs);	
	
	int maxEpochs = 0;
	double currentLearningRate = 0;
	double reductionFactor = 1;
    double startLearningRate = 0.1;
    double endLearningRate = 0.001;
	learningRateInterp learningRateInterpMode = learningRateInterp::linear;
};



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
    virtual void updateParameters(double learningRate,
                                  gsl::span<const double> gradient,
                                  gsl::span<double> parameters,
                                  gsl::span<const double> lowerBounds = gsl::span<const double>(),
                                  gsl::span<const double> upperBounds = gsl::span<const double>()) = 0;
    
	virtual void undoLastStep() = 0;
	
	virtual void clearCache() = 0;
	
    virtual ~ParameterUpdater() = default;

};


class ParameterUpdaterVanilla : public ParameterUpdater {
public:
    ParameterUpdaterVanilla() = default;
//    ParameterUpdaterVanilla(learningRateAdaption learningRateAdaptionMode, double startLearningRate, double endLearningRate);

    void updateParameters(double learningRate,
                          gsl::span<const double> gradient, 
                          gsl::span<double> parameters,
                          gsl::span<const double> lowerBounds = gsl::span<const double>(),
                          gsl::span<const double> upperBounds = gsl::span<const double>());

    void undoLastStep();
    
	void clearCache();
};


class ParameterUpdaterRmsProp : public ParameterUpdater {
public:
    ParameterUpdaterRmsProp() = default;
//    ParameterUpdaterRmsProp(learningRateAdaption learningRateAdaptionMode, double startLearningRate, double endLearningRate);

    void updateParameters(double learningRate,
                          gsl::span<const double> gradient, 
                          gsl::span<double> parameters,
                          gsl::span<const double> lowerBounds = gsl::span<const double>(),
                          gsl::span<const double> upperBounds = gsl::span<const double>());

    void undoLastStep();
    
	void clearCache();

    int updates = 0;
    double decayRate = 0.9;
    double delta = 1e-7;
    std::vector<double> gradientNormCache;
    std::vector<double> oldGradientNormCache;
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
        std::vector<double> oldParameters(parameters.size(), NAN);
        std::vector<BATCH_ELEMENT> shuffledData(data.begin(), data.end());

        std::vector<double> gradient(parameters.size(), 0.0);
        std::vector<double> oldGradient(parameters.size(), NAN);
        std::random_device rd;
        std::mt19937 rng(rd());
        double cost = NAN;
		int subsequentFails = 0;
		int maxSubsequentFails = 10;
		bool finalFail = false;
		bool coldRestartActive = false;
		learningRateUpdater->setMaxEpochs(maxEpochs);

        if(reporter) reporter->starting(initialParameters);

        for(int epoch = 0; epoch < maxEpochs; ++epoch) {
            auto epochLogger = logger.getChild(std::string("e") + std::to_string(epoch));

            // Create randomized batches
            std::shuffle(shuffledData.begin(), shuffledData.end(), rng);
            auto batches = getBatches<BATCH_ELEMENT>(shuffledData, batchSize);

			// Update learning rate according to epoch
			learningRateUpdater->updateLearningRate(epoch);

            for(int batchIdx = 0; (unsigned) batchIdx < batches.size(); ++batchIdx) {
                auto batchLogger = epochLogger->getChild(std::string("b") + std::to_string(batchIdx));

                auto status = evaluate(f, parameters, batches[batchIdx], cost, gradient, batchLogger.get(), reporter);
                
                // Give some output
                learningRate = learningRateUpdater->getCurrentLearningRate();
                std::stringstream ss;
                ss << ": Cost: " << cost
                 << " |g|2: " << getVectorNorm(gradient)
                 << " Batch: " << batches[batchIdx]
				 << " LearningRate: " << learningRate << std::endl;
                batchLogger->logmessage(LOGLVL_DEBUG, ss.str().c_str());

                if(reporter) reporter->iterationFinished(parameters, cost, gradient);

                if(status == functionEvaluationFailure) {
                	// Check, if the interceptor should be used (should alwayss be the case, except for study purpose...
                	if(interceptor > 0) {
						// Cost function evaluation failed: We need to intercept
						while(status == functionEvaluationFailure) {
							// If the objective function evaluation failed, we want to undo the step
							subsequentFails++;
							parameterUpdater->undoLastStep();
							gradient = oldGradient;
							parameters = oldParameters;
	
							// Check if there are NaNs in the parameter vector now (e.g., fail at first iteration)						
							for(int ip = 0; ip < parameters.size(); ip++) {
								if(std::isnan(parameters[ip])) {
									finalFail = true;
									break;
								}
							}
	
							// If too many fails: cancel optimization 	
							if(subsequentFails >= maxSubsequentFails || finalFail) {
								if(interceptor > 1 && !coldRestartActive) {
									/* Reducing step size did not work. Yet, a small step in descent direction
									 * should actually do the job. So clear all cached gradients and retry with 
									 * a very small step size (e.g. 1e-5)
									 */ 
									subsequentFails = 0;
									parameterUpdater->clearCache();
									learningRateUpdater->setReductionFactor(1e-5);
								} else {
									// Really everything failed, there is no hope for this run any more
									return finish(cost, parameters, minibatchExitStatus::invalidNumber, reporter, batchLogger.get());
								}
							} else {
								// If we did not fail too often, we reduce the step size and try to redo the step
								learningRateUpdater->reduceLearningRate();
							}
							// Do the next step
							learningRate = learningRateUpdater->getCurrentLearningRate();
							parameterUpdater->updateParameters(learningRate, gradient, parameters,
															   lowerParameterBounds, upperParameterBounds);
	
							// Re-evaluate the cost function and hope for the best
							auto status = evaluate(f, parameters, batches[batchIdx], cost, gradient, batchLogger.get(), reporter);
						}
					} else {
	                	// Directly stop optimization
	                	return finish(cost, parameters, minibatchExitStatus::invalidNumber, reporter, batchLogger.get());
					}
                } else {
					// Cost function evaluation was succeful, so we can increse the step size
					// (if was reduced at some earlier point)
					subsequentFails = std::max(subsequentFails - 1, 0);
					learningRateUpdater->increaseLearningRate();

					// Overwrite old parameters and old gradient, since they won't be needed any more
					oldGradient = gradient;
					oldParameters = parameters;
                }

				learningRate = learningRateUpdater->getCurrentLearningRate();
                parameterUpdater->updateParameters(learningRate, gradient, parameters,
                                                   lowerParameterBounds, upperParameterBounds);
            }

            if(getVectorNorm(gradient) <= gradientNormThreshold) {
                // evaluate on full data set
                auto dataSpan = std::vector<BATCH_ELEMENT>(data.cbegin(), data.cend());
                auto status = evaluate(f, parameters, dataSpan, cost, gradient, epochLogger.get(), reporter);
                return finish(cost, parameters, minibatchExitStatus::gradientNormConvergence, reporter, epochLogger.get());
            }
        }

        // evaluate on full data set
        auto dataSpan = std::vector<BATCH_ELEMENT>(data.cbegin(), data.cend());
        auto status = evaluate(f, parameters, dataSpan, cost, gradient, &logger, reporter);

        return finish(cost, parameters, minibatchExitStatus::maxEpochsExceeded, reporter, &logger);
    }

    FunctionEvaluationStatus evaluate(
            SummedGradientFunction<BATCH_ELEMENT> const& f,
            gsl::span<const double> parameters,
            std::vector<BATCH_ELEMENT> datasets,
            double &cost,
            gsl::span<double> gradient,
            Logger *logger,
            OptimizationReporter *reporter
            ) const
    {
        if(reporter) {
            reporter->beforeCostFunctionCall(parameters);
            reporter->logger->setPrefix(logger->getPrefix());
        }

        double cpuTime = 0.0;
        auto status = f.evaluate(parameters, datasets, cost, gradient, logger, &cpuTime);

        if(reporter) {
            reporter->cpuTimeIterationSec += cpuTime;
            reporter->cpuTimeTotalSec += cpuTime;
            reporter->afterCostFunctionCall(parameters, cost, gradient);
        }

        // Normalize to batch size
        double batchSize = datasets.size();
        cost /= batchSize;
        for(auto &g: gradient)
            g /= batchSize;

        return status;
    }

    std::tuple<int, double, std::vector<double> > finish(double cost,
                std::vector<double> const& parameters,
                minibatchExitStatus status,
                OptimizationReporter *reporter,
                Logger *logger)
    {
        if(logger) {
            switch(status) {
            case minibatchExitStatus::invalidNumber:
                logger->logmessage(LOGLVL_ERROR, "Minibatch cost function evaluation failed.");
                break;
            case minibatchExitStatus::gradientNormConvergence:
                logger->logmessage(LOGLVL_INFO, "Convergence: gradientNormThreshold reached.");
                break;
            case minibatchExitStatus::maxEpochsExceeded:
                logger->logmessage(LOGLVL_INFO, "Number of epochs exceeded.");
            }
        }

        if(reporter) reporter->finished(cost, parameters, (int)status);

        return std::tuple<int, double, std::vector<double> >((int)status, cost, parameters);
    }

    std::unique_ptr<ParameterUpdater> parameterUpdater = std::make_unique<ParameterUpdaterVanilla>();

    // Set some default values
    int interceptor = 2;
    int batchSize = 1;
    int maxEpochs = 3;
    double gradientNormThreshold = 0.0;
	double learningRate = 0.001;

    std::unique_ptr<LearningRateUpdater> learningRateUpdater = 
			std::make_unique<LearningRateUpdater>(maxEpochs, learningRateInterp::linear);
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
