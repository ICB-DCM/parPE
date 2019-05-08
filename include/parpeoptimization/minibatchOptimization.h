#ifndef PARPE_OPTIMIZATION_MINIBATCH_OPTIMIZATION_H
#define PARPE_OPTIMIZATION_MINIBATCH_OPTIMIZATION_H

#include <parpeoptimization/optimizationOptions.h>
#include <parpeoptimization/optimizationProblem.h>
#include <parpeoptimization/optimizationResultWriter.h>
#include <parpecommon/misc.h>
#include <parpecommon/model.h>

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <functional>

namespace parpe {

/**
 * @brief Return status for minibatch optimizer
 */
enum class minibatchExitStatus {
    gradientNormConvergence, maxEpochsExceeded, invalidNumber
};

/**
 * @brief Shape of learning rate interpolation
 */
enum class learningRateInterp {
    linear, inverseLinear, logarithmic
};

/**
 * @brief Reaction upon ODE solver crashes
 */
enum class interceptType {
    none, reduceStep, reduceStepAndRestart
};

/**
 * @brief Problem definition for mini-batch optimization.
 *
 * This class provides cost function and training data for a mini-batch optimizer.
 * Data maybe be the actual data or just an index list referencing the data, the
 * cost function will operate on.
 */
template<typename T>
class MinibatchOptimizationProblem: public OptimizationProblem {
public:
    MinibatchOptimizationProblem() = default;

    MinibatchOptimizationProblem(std::unique_ptr<SummedGradientFunction<T>> costFun,
                                 std::unique_ptr<Logger> logger) :
            OptimizationProblem(costFun, logger) {
    }

    MinibatchOptimizationProblem(MinibatchOptimizationProblem const& other) = delete;

    virtual ~MinibatchOptimizationProblem() override = default;

    /** vector of training data */
    virtual std::vector<T> getTrainingData() const = 0;

    /** mini batch cost function */
    SummedGradientFunction<T>* getGradientFunction() const {
        auto summedGradientFunction = dynamic_cast<SummedGradientFunction<T>*>(costFun.get());
        RELEASE_ASSERT(summedGradientFunction, "");
        return summedGradientFunction;
    }
};

/**
 * @brief learning rate updaters for minibatch optimizers
 *
 * The LearningRateUpdater provides the possibility to reduce the learning rate per epoch
 * and makes it possible to adapt the learning rate according to success or failure of
 * the ODE solver.
 */
class LearningRateUpdater {
public:
    /**
     * @brief Update the learning rate
     * @param maxEpochs Maximum number of epochs in optimization
     * @param LearningRateadaptionMode Type of interpolation between startLearningRate and endLearningRate
     */
    LearningRateUpdater(int maxEpochs,
                        learningRateInterp learningRateInterpMode);

    /** Update function, to be called in every epoch or optimization iteration */
    void updateLearningRate(int iteration);

    /** Update function, to be called if parameter update did not work well */
    void reduceLearningRate();

    /** Update function, to be called if parameter update worked well */
    void increaseLearningRate();

    /** Function to retrieve the learning rate */
    double getCurrentLearningRate();

    /** Function to set the reduction factor directly */
    void setReductionFactor(double newReductionFactor);

    /** Function to set the new maximum epoch number */
    void setMaxEpochs(int newMaxEpochs);

    /** Function to set start learning rate */
    void setStartLearningRate(double learningRate);

    /** Function to set end learning rate */
    void setEndLearningRate(double learningRate);

private:
    /** Maximum number of epochs, will be set after creation of problem instance */
    int maxEpochs = 0;

    /** Learning rate, i.e. step size, at the moment of optimization */
    double currentLearningRate = 0.0;

    /** If an optimization step is not succesful, the learning rate, i.e., step size, will be reduced by this factor */
    double reductionFactor = 4.0;

    /** Learning rate, i.e. step size, at the beginning of optimization */
    double startLearningRate = 0.1;

    /** Learning rate, i.e. step size, at the end of optimization */
    double endLearningRate = 0.001;

    /** Mode of interpolation between the beginning and the end of optimization */
    learningRateInterp learningRateInterpMode = learningRateInterp::linear;
};

/**
 * @brief Interface for parameter updaters for minibatch optimizers
 */
class ParameterUpdater {
public:
    /**
     * @brief Update parameter vector
     * @param learningRate Current learning rate, i.e., step-size
     * @param iteration Current iteration, i.e., epoch
     * @param gradient Cost function gradient at parameters
     * @param parameters In: Current parameters, Out: Updated parameters
     */
    virtual void updateParameters(double learningRate,
                                  int iteration,
                                  gsl::span<const double> gradient,
                                  gsl::span<double> parameters,
                                  gsl::span<const double> lowerBounds = gsl::span<const double>(),
                                  gsl::span<const double> upperBounds = gsl::span<const double>()) = 0;

    /** If ODE becomes non-integrable, the last step must be undone using this method */
    virtual void undoLastStep() = 0;

    /** If the ODE is repeatedly non-integrable, a cold restart is performed using this method */
    virtual void clearCache() = 0;

    /** Intialize the parameter updater */
    virtual void initialize(unsigned int numParameters) = 0;

    virtual ~ParameterUpdater() = default;

};

/** 
 * @brief Minibatch optimizer: Vanilla SGD Updater
 * The simplest mini batch algorithm.
 */
class ParameterUpdaterVanilla: public ParameterUpdater {
public:
    ParameterUpdaterVanilla() = default;

    void updateParameters(double learningRate,
                          int iteration,
                          gsl::span<const double> gradient,
                          gsl::span<double> parameters,
                          gsl::span<const double> lowerBounds = gsl::span<const double>(),
                          gsl::span<const double> upperBounds = gsl::span<const double>()) override;

    void undoLastStep() override;

    void clearCache() override;

    void initialize(unsigned int numParameters) override;
};

/**
 * @brief Minibatch optimizer: RMSProp Updater
 * A so-called adaptive mini batching algorithm without momentum
 */
class ParameterUpdaterRmsProp: public ParameterUpdater {
public:
    ParameterUpdaterRmsProp() = default;

    void updateParameters(double learningRate,
                          int iteration,
                          gsl::span<const double> gradient,
                          gsl::span<double> parameters,
                          gsl::span<const double> lowerBounds = gsl::span<const double>(),
                          gsl::span<const double> upperBounds = gsl::span<const double>()) override;

    void undoLastStep() override;

    void clearCache() override;

    void initialize(unsigned int numParameters) override;

private:

    /** Rate for memorizing gradient norms (between 0 and 1, high rates mean long memory) */
    double decayRate = 0.9;

    /** Stabilization factor for gradient normalization (avoid deviding by 0) */
    double delta = 1e-7;

    /** Memorized gradient norms (decaying average) from last steps */
    std::vector<double> gradientNormCache;

    /** Memorized gradient norms (decaying average), one step back (if one step must be undone) */
    std::vector<double> oldGradientNormCache;
};

 /**
 * @brief Minibatch optimizer: Adam Updater
 * A momentum-based and so-called adaptive mini batching algorithm
 */
class ParameterUpdaterAdam: public ParameterUpdater {
public:
    ParameterUpdaterAdam() = default;

    void updateParameters(double learningRate,
                          int iteration,
                          gsl::span<const double> gradient,
                          gsl::span<double> parameters,
                          gsl::span<const double> lowerBounds = gsl::span<const double>(),
                          gsl::span<const double> upperBounds = gsl::span<const double>()) override;

    void undoLastStep() override;
    
    void clearCache() override;

    void initialize(unsigned int numParameters) override;

private:

    /** Rate for memorizing gradients (between 0 and 1, high rates mean long memory) */
    double decayRateGradient = 0.9;

    /** Rate for memorizing gradient norms (between 0 and 1, high rates mean long memory) */
    double decayRateGradientNorm = 0.9;

    /** Stabilization factor for gradient normalization (avoid deviding by 0) */
    double delta = 1e-7;

    /** Memorized gradient norms (decaying average) from last steps */
    std::vector<double> gradientNormCache;

    /** Memorized gradient norms (decaying average), one step back (if one step must be undone) */
    std::vector<double> oldGradientNormCache;

    /** Memorized gradients (decaying average) from last steps */
    std::vector<double> gradientCache;

    /** Memorized gradients (decaying average), one step back (if one step must be undone) */
    std::vector<double> oldGradientCache;
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
std::vector<std::vector<T>> getBatches(gsl::span<const T> data,
                                       int batchSize) {
    int numBatches = ceil(static_cast<double>(data.size()) / batchSize);

    std::vector < std::vector < T >> batches(numBatches, std::vector<T>());
    for (int i = 0, batchIdx = -1; static_cast<typename decltype (data)::index_type>(i) < data.size(); ++i) {
        if (i % batchSize == 0) {
            ++batchIdx;
            int remaining = data.size() - i;
            batches[batchIdx].resize(std::min(batchSize, remaining));
        }
        batches[batchIdx][i % batchSize] = data[i];
    }

    return batches;
}

/**
 * @brief Get scalar product of two vectors.
 * @param v
 * @param w
 * @return the scalar product
 */
double getScalarProduct(gsl::span<const double> v,
                        gsl::span<const double> w);

/**
 * @brief Get Euclidean (l2) norm of vector.
 * @param v
 * @return the norm
 */
double getVectorNorm(gsl::span<const double> v);

/**
 * @brief Get difference of two vectors.
 * @param v
 * @param w
 * @return the norm
 */
std::vector<double> getVectorDifference(gsl::span<const double> v,
                                        gsl::span<const double> w);

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
    std::tuple<int, double, std::vector<double> > optimize(SummedGradientFunction<BATCH_ELEMENT> const& f,
                                                           gsl::span<const BATCH_ELEMENT> data,
                                                           gsl::span<const double> initialParameters,
                                                           gsl::span<const double> lowerParameterBounds,
                                                           gsl::span<const double> upperParameterBounds,
                                                           OptimizationReporter *reporter,
                                                           Logger *logger_) {
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
        int iteration = 0;
        int subsequentFails = 0;
        learningRateUpdater->setMaxEpochs(maxEpochs);
        parameterUpdater->initialize(initialParameters.size());

        if (reporter) {
            reporter->starting(initialParameters);
            reporter->resultWriter->setLoggingEachFunctionEvaluation(false, true);
            reporter->resultWriter->setLoggingEachIteration(false);
        }

        for (int epoch = 0; epoch < maxEpochs; ++epoch) {
            auto epochLogger = logger.getChild(std::string("e") + std::to_string(epoch));

            // Create randomized batches
            std::shuffle(shuffledData.begin(), shuffledData.end(), rng);
            auto batches = getBatches<BATCH_ELEMENT>(shuffledData, batchSize);

            // Update learning rate according to epoch
            learningRateUpdater->updateLearningRate(epoch);

            for (int batchIdx = 0; (unsigned) batchIdx < batches.size(); ++batchIdx) {
                auto batchLogger = epochLogger->getChild(std::string("b") + std::to_string(batchIdx));
                iteration++;

                auto status = evaluate(f, parameters, batches[batchIdx], cost, gradient, batchLogger.get(), reporter);

                // Give some output
                learningRate = learningRateUpdater->getCurrentLearningRate();
                std::stringstream ss;
                ss << ": Cost: " << cost << " |g|2: " << getVectorNorm(gradient) << " Batch: " << batches[batchIdx]
                        << " LearningRate: " << learningRate << std::endl;
                batchLogger->logmessage(LOGLVL_DEBUG, ss.str().c_str());

                if (status == functionEvaluationFailure) {
                    // Check, if the interceptor should be used (should alwayss be the case, except for study purpose...
                    if (interceptor > interceptType::none)
                        status = rescueInterceptor(parameters, oldParameters, gradient, oldGradient,
                                                   lowerParameterBounds, upperParameterBounds, cost, subsequentFails,
                                                   iteration, f, batches[batchIdx], batchLogger.get(), reporter);

                    // If we still have a failure, stop optimization
                    if (status == functionEvaluationFailure)
                        return finish(cost, parameters, minibatchExitStatus::invalidNumber, reporter, batchLogger.get());

                } else {
                    // Cost function evaluation was succeful, so we can increase the step size
                    subsequentFails = std::max(subsequentFails - 1, 0);
                    learningRateUpdater->increaseLearningRate();

                    // Overwrite old parameters and old gradient, since they won't be needed any more
                    std::copy(gradient.begin(), gradient.end(), oldGradient.begin())
                    std::copy(parameters.begin(), parameters.end(), parameters.begin())
                }

                /* Update parameters after successful gradient evaluation */
                handleStep(parameters, oldParameters, gradient, lowerParameterBounds, upperParameterBounds,
                           cost, iteration, f, batches[batchIdx], batchLogger.get(), reporter);

            }

            // epoch finished, write the values in hdf5-file
            if (reporter)
                reporter->iterationFinished(parameters, cost, gradient);

            if (getVectorNorm(gradient) <= gradientNormThreshold) {
                // evaluate on full data set
                auto dataSpan = std::vector < BATCH_ELEMENT > (data.cbegin(), data.cend());
                evaluate(f, parameters, dataSpan, cost, gradient, epochLogger.get(), reporter);
                return finish(cost, parameters, minibatchExitStatus::gradientNormConvergence, reporter,
                              epochLogger.get());
            }
        }

        // evaluate on full data set
        auto dataSpan = std::vector < BATCH_ELEMENT > (data.cbegin(), data.cend());
        evaluate(f, parameters, dataSpan, cost, gradient, &logger, reporter);

        return finish(cost, parameters, minibatchExitStatus::maxEpochsExceeded, reporter, &logger);
    }

    FunctionEvaluationStatus evaluate(SummedGradientFunction<BATCH_ELEMENT> const& f,
                                      gsl::span<const double> parameters,
                                      std::vector<BATCH_ELEMENT> datasets,
                                      double &cost,
                                      gsl::span<double> gradient,
                                      Logger *logger,
                                      OptimizationReporter *reporter) const {
        if (reporter) {
            reporter->beforeCostFunctionCall(parameters);
            reporter->logger->setPrefix(logger->getPrefix());
        }

        double cpuTime = 0.0;
        auto status = f.evaluate(parameters, datasets, cost, gradient, logger, &cpuTime);

        // Normalize to batch size
        double batchSize = datasets.size();
        cost /= batchSize;
        for (auto &g : gradient)
            g /= batchSize;

        if (reporter) {
            reporter->cpuTimeIterationSec += cpuTime;
            reporter->cpuTimeTotalSec += cpuTime;
            reporter->afterCostFunctionCall(parameters, cost, gradient);
        }

        return status;
    }

    std::tuple<int, double, std::vector<double> > finish(double cost,
                                                         std::vector<double> const& parameters,
                                                         minibatchExitStatus status,
                                                         OptimizationReporter *reporter,
                                                         Logger *logger) {
        if (logger) {
            switch (status) {
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

        if (reporter)
            reporter->finished(cost, parameters, (int) status);

        return std::tuple<int, double, std::vector<double> >((int) status, cost, parameters);
    }

    /**
     * @brief Try to rescue optimization run at cost function failure
     *
     * @param parameters current parameter vector
     * @param oldParameters parameter vector before last step
     * @param gradient current cost function gradient
     * @param oldGradient cost function gradient before last step
     * @param cost new cost function value after interception
     * @param subsequentFails number of iterations during rescue interceptor
     * @param f Function to minize
     * @param data Full data set on which f will be evaluated
     * @param logger Logger instance for status messages
     * @param reporter OptimizationReporter instance for tracking progress
     * @return FunctionEvaluationStatus
     */
    FunctionEvaluationStatus rescueInterceptor(gsl::span<double> parameters,
                                               gsl::span<double> oldParameters,
                                               gsl::span<double> gradient,
                                               gsl::span<double> oldGradient,
                                               gsl::span<const double> lowerParameterBounds,
                                               gsl::span<const double> upperParameterBounds,
                                               double &cost,
                                               int &subsequentFails,
                                               int iteration,
                                               SummedGradientFunction<BATCH_ELEMENT> const& f,
                                               std::vector<BATCH_ELEMENT> datasets,
                                               Logger *logger,
                                               OptimizationReporter *reporter) {
        
        // initialize disgnostic variables 
        int maxSubsequentFails = 10;
        bool finalFail = false;
        bool initialFail = false;
        cost = NAN;
        FunctionEvaluationStatus status = functionEvaluationFailure;
        std::vector<double> parDifference(parameters.size(), 0.0);
        std::stringstream ss;
        
        if (reporter) {
            reporter->beforeCostFunctionCall(parameters);
            reporter->logger->setPrefix(logger->getPrefix());
        }
        
        // debug output
        if (logger) {
            std::vector<double> firstDifference = getVectorDifference(parameters[i], oldParameters[i]);
            std::stringstream first_ss;
            first_ss << " Interceptor is active! Former step size: " << getVectorNorm(firstDifference) << std::endl;
            logger->logmessage(LOGLVL_DEBUG, first_ss.str().c_str());
        }
        

        // Cost function evaluation failed: We need to intercept
        while (status == functionEvaluationFailure) {
            
            // If the objective function evaluation failed, we want to undo the step
            ++subsequentFails;
            parameterUpdater->undoLastStep();
            std::copy(oldGradient.begin(), oldGradient.end(), gradient.begin())
            std::copy(oldParameters.begin(), oldParameters.end(), parameters.begin())

            // note the previous parameter step for debugging
            parDifference = getVectorDifference(parameters[i], oldParameters[i]);
            
            // Check if there are NaNs in the parameter vector now (e.g., fail at first iteration)
            if (std::any_of(parameters.begin(), parameters.end(), [](double d) {return std::isnan(d);}))
                initialFail = true;
            if (subsequentFails >= maxSubsequentFails)
                finalFail = true;

            // If nothing helps and no cold restart wanted: cancel optimization
            if (initialFail or (finalFail and interceptor != interceptType::reduceStepAndRestart)) {
                logger->logmessage(LOGLVL_DEBUG, "Failure at intial point of optimization. Stopping.");
                return functionEvaluationFailure;
            }
            
            if (finalFail) {
                /* Reducing step size did not work. 
                 * Do a cold restart and take a very small step.
                 */
                ss << "Final failure. Rescue interceptor could not recover optimization run." << std::endl;
                subsequentFails = 0;
                parameterUpdater->clearCache();
                learningRateUpdater->setReductionFactor(1e-5);
            } else {
                /* We did not fail too often: we reduce the step size */
                 learningRateUpdater->reduceLearningRate();
                 ss << "Failure. Reducing learning rate and retry..." << std::endl;
            }
            
            // debug output
            ss << ": Interceptor, before new evalaluation: " << std::endl 
                << " New cost: " << cost 
                << ", new |g|2: " << getVectorNorm(gradient)
                << ", new LearningRate: " << learningRate 
                << ", real step length: " << getVectorNorm(parDifference)
                << ", Failures: " << subsequentFails 
                << std::endl;
            
            // Do the next step
            learningRate = learningRateUpdater->getCurrentLearningRate();
            parameterUpdater->updateParameters(learningRate, iteration, gradient, parameters, lowerParameterBounds,
                                               upperParameterBounds);

            // Re-evaluate the cost function and hope for the best
            status = evaluate(f, parameters, datasets, cost, gradient, logger, reporter);
            
            // get Difference for debugging
            parDifference = getVectorDifference(parameters[i], oldParameters[i]);
            
            // debug output
            ss << ": Interceptor, after new evalaluation: " << std::endl 
                << " New Cost: " << cost 
                << " new |g|2: " << getVectorNorm(gradient)
                << " new LearningRate: " << learningRate 
                << "real step length: " << getVectorNorm(parDifference) 
                << std::endl;
            if (logger)
                logger->logmessage(LOGLVL_DEBUG, ss.str().c_str());

        }

        return status;
    }

    /**
     * @brief Try to determine a good step length in descent direction
     *
     * @param parameters current parameter vector
     * @param oldParameters parameter vector before last step
     * @param gradient current cost function gradient
     * @param oldGradient cost function gradient before last step
     * @param cost new cost function value after interception
     * @param subsequentFails number of iterations during rescue interceptor
     * @param f Function to minize
     * @param data Full data set on which f will be evaluated
     * @param logger Logger instance for status messages
     * @param reporter OptimizationReporter instance for tracking progress
     * @return FunctionEvaluationStatus 
     */
    void handleStep(gsl::span<double> parameters,
                    gsl::span<double> oldParameters,
                    gsl::span<double> gradient,
                    gsl::span<const double> lowerParameterBounds,
                    gsl::span<const double> upperParameterBounds,
                    double cost,
                    int iteration,
                    SummedGradientFunction<BATCH_ELEMENT> const& f,
                    std::vector<BATCH_ELEMENT> datasets,
                    Logger *logger,
                    OptimizationReporter *reporter) {

        /* Retrieve step length and try a full step */
        double stepLength = learningRateUpdater->getCurrentLearningRate();
        parameterUpdater->updateParameters(stepLength, iteration, gradient, parameters, 
                                           lowerParameterBounds, upperParameterBounds);

        /* If no line search desired: that's it! */
        if (lineSearchSteps == 0)
            return;
        
        /* Define lambda function for step length evaluation  */
        std::function<double (double)> evalLineSearch = [&f, &datasets, iteration,
                   &parameters, &oldParameters, &gradient, 
                   &lowerParameterBounds, &upperParameterBounds, 
                   &logger, &reporter, this](double alpha) {
            
            /* Reset oldParameters and re-update with new step length */
            std::copy(oldParameters.begin(), oldParameters.end(), parameters.begin())
            parameterUpdater->undoLastStep();
            
            parameterUpdater->updateParameters(alpha, iteration, gradient, parameters, 
                                               lowerParameterBounds, upperParameterBounds);
            /* Write new cost funtion value and return */
            double newCost = NAN;
            evaluate(f, parameters, datasets, newCost, gsl::span<double>(), logger, reporter);
            return newCost;
        };
        
        /* From here, we check for a decreasing cost function,
         * otherwise, we already exited. 
         * 
         * The idea:
         * cost0 is the previous cost before the parameter update,
         * we want the update to be lower.
         * 
         * First compute a naive step as if there was now line-search -> get cost1
         * If cost1 > cost0: Try to improve -> get cost2
         * If also cost2 > cost0:
         *     If only short line-search is desired:  take min(cost1, cost2)
         *     If more longer line-search is desired: step into performLineSearch()
         * */
        
        double cost1 = evalLineSearch(stepLength);

        /* Return on improvement */
        if (cost1 <= cost) {
            logger->logmessage(LOGLVL_DEBUG, " Line-Search: Step was good right away...");
            return;
        }
        
        // Debugging output
        std::stringstream line_ss;
        std::stringstream parabola_ss;
        std::vector<double> parDifference = getVectorDifference(parameters[i], oldParameters[i]);

        /* No improvement: compute update direction */
        std::vector<double> direction(parameters.size(), NAN);
        for (unsigned int i = 0; i < gradient.size(); ++i)
            direction[i] = parameters[i] - oldParameters[i];
        double dirNorm = getVectorNorm(direction);
        for (unsigned int i = 0; i < gradient.size(); ++i)
            direction[i] /= dirNorm;
        
        /* Is the step direction a descent direction? */
        double dirGradient = getScalarProduct(direction, gradient);
        if (dirGradient > 0) {
            /* No descent direction, no hope for improvement:
             * Try to do something smart anyway */
            parabola_ss << " Line-Search: Not a descent direction! " << std::endl;
            
            /* Fit a parabola to decide whether a smaller or 
             * a bigger step seems more promising */
            double newStepLength = NAN;
            if (cost1 < cost - 2.0 * dirNorm * dirGradient) {
                parabola_ss << " Stepping longer now. " << std::endl;
                newStepLength = stepLength * 2.0;
            } else {
                parabola_ss << " Stepping shorter now. " << std::endl;
                newStepLength = stepLength / 2.0;
            }
            
            /* re-evaluate cost function */
            double cost2 = evalLineSearch(newStepLength);
            
            if (cost2 > cost1) {
                /* The parabola idea didn't work. Just admit the step, as it is */
                cost1 = evalLineSearch(stepLength);
                parabola_ss << " Step adaption did not work..." << std::endl;
            }
            parabola_ss << "   cost0: " << cost << " (steplength: " << getVectorNorm(parDifference)
                << "),   cost1: " << cost1 << " (steplength: " << stepLength
                << "),   cost2: " << cost2 << " (steplength: " << newStepLength
                << ") " << std::endl;
            if (logger)
                logger->logmessage(LOGLVL_DEBUG, parabola_ss.str().c_str());
            
            /* We tried all we could */
            return;
        }
        
        /* Original step was too big, but we're facing a descent direction 
         * Propose a new step based on a parabolic interpolation */
        double newStepLength = -0.5 * dirGradient * std::pow(stepLength, 2.0) /
                (cost1 - cost - dirGradient * stepLength);
        double cost2 = evalLineSearch(newStepLength);
        
        // Debugging output
        if (logger) {
            line_ss << " Line-Search: two steps were done, results are: " << std::endl
                << "   cost0: " << cost << " (steplength: " << getVectorNorm(parDifference)
                << "),   cost1: " << cost1 << " (steplength: " << stepLength
                << "),   cost2: " << cost2 << " (steplength: " << newStepLength
                << ")" << std::endl;
            logger->logmessage(LOGLVL_DEBUG, line_ss.str().c_str());
        }
        
        /* If we did improve, return, otherwise iterate */
        if (cost2 < cost) return;
        
        if (lineSearchSteps < 2) {
            /* No more iteration wanted, but 2nd try was better than 1st */
            if (cost2 <= cost1) return;
            
            /* 1st try was better than 2nd, use it */
            cost1 = evalLineSearch(stepLength);
        } else {
            /* No descent found and line search option is set: iterate! */
            std::stringstream line_ss;
            if (logger) {
                line_ss << " Line-Search: Need to go to third order approximation, looping... " << std::endl;
                logger->logmessage(LOGLVL_DEBUG, line_ss.str().c_str());
            }
            performLineSearch(stepLength,
                              newStepLength,
                              cost,
                              cost1,
                              cost2,
                              dirGradient,
                              evalLineSearch);
        }
    }

    /**
     * @brief Perform line search according to interpolation algo by [Dennis and Schnabel, 
     * Numerical Methods for Unconstrained Optimization and Non-linear Equations, 1993].
     *
     * @param alpha1 step length of first step
     * @param alpha2 step length of second step
     * @param cost objective function value before update
     * @param cost1 objective function value after first step
     * @param cost2 objective function value after second step
     * @param dirGradient alignment of gradient and step direction
     * @param costFunEvaluate objective function wo gradient
     */
    void performLineSearch(double alpha1,
                           double alpha2, 
                           double cost, 
                           double cost1, 
                           double cost2, 
                           double dirGradient,
                           std::function<double (double)> costFunEvaluate) {

        /* From here on, we will use cubic interpolation.
         * We need to compute the matrix-vector multiplication
         * 
         * a = 1/tmp_D * [ tmp_11, tmp_12 ] [tmp_v1]
         * b = 1/tmp_D * [ tmp_21, tmp_22 ] [tmp_v2]
         * 
         * in order to find the possible minimum at
         * 
         * alpha3 = -b + sqrt( ² - 3*a*dirGradient ) / (3*a)
         * 
         * Possibly, we have to iterrate this process. */

        /* Declare variables which will be needed outside the loop */
        double cost3 = NAN;
        
        for (int iStep = 2; iStep <= lineSearchSteps; ++iStep) {
            /* declare temporary variables */
            double tmp_D = std::pow(alpha2 * alpha1, 2.0) * (alpha2 - alpha1);
            double tmp_11 = std::pow(alpha1, 2.0);
            double tmp_12 = -1.0 * std::pow(alpha2, 2.0);
            double tmp_21 = -1.0 * std::pow(alpha1, 3.0);
            double tmp_22 = std::pow(alpha2, 3.0);
            double tmp_v1 = cost2 - cost - dirGradient * alpha2;
            double tmp_v2 = cost1 - cost - dirGradient * alpha1;
            double a = (tmp_11 * tmp_v1 + tmp_12 * tmp_v2) / tmp_D;
            double b = (tmp_21 * tmp_v1 + tmp_22 * tmp_v2) / tmp_D;
            
            /* Compute possible new step length */
            double alpha3 = (-b + std::sqrt(b*b - 3.0*a*dirGradient)) / (3.0*a);
            
            /* Evaluate line search function at alpha2 */
            double cost3 = costFunEvaluate(alpha3);
            
            /* If improvement, return */
            if (cost3 < cost) return;
            
            /* If no improvement, update values and re-iterate */
            if (iStep < lineSearchSteps) {
                cost1 = cost2;
                cost2 = cost3;
                alpha1 = alpha2;
                alpha2 = alpha3;
            }
        }
        
        /* No good step was found, but the max number 
         * of line search steps was reached */
        if (cost3 < std::min(cost1, cost2)) return;
        
        if (cost1 < cost2) {
            cost1 = costFunEvaluate(alpha1);
        } else {
            cost2 = costFunEvaluate(alpha2);
        }
    }
        
    
    std::unique_ptr<ParameterUpdater> parameterUpdater = std::make_unique<ParameterUpdaterVanilla>();

    // Set some default values
    interceptType interceptor = interceptType::reduceStepAndRestart;
    int lineSearchSteps = 0;
    int batchSize = 1;
    int maxEpochs = 1;
    double gradientNormThreshold = 0.0;
    double learningRate = 0.001;

    std::unique_ptr<LearningRateUpdater> learningRateUpdater = std::make_unique < LearningRateUpdater
            > (maxEpochs, learningRateInterp::linear);
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
std::unique_ptr<MinibatchOptimizer<BATCH_ELEMENT>> getMinibatchOptimizer(OptimizationOptions const& options) {
    auto optim = std::make_unique<MinibatchOptimizer<BATCH_ELEMENT>>();

    options.for_each<MinibatchOptimizer<BATCH_ELEMENT>*>(setMinibatchOption, optim.get());

    return optim;
}

std::tuple<int, double, std::vector<double> > runMinibatchOptimization(MinibatchOptimizationProblem<int> *problem);

/**
 * @brief Clip values to given element-wise bounds.
 * @param lowerBounds
 * @param upperBounds
 * @param x
 */
template<typename T>
void clipToBounds(gsl::span<const T> lowerBounds,
                  gsl::span<const T> upperBounds,
                  gsl::span<T> x) {
    if (lowerBounds.empty() && upperBounds.empty())
        return;

    RELEASE_ASSERT(lowerBounds.size() == upperBounds.size(), "");
    RELEASE_ASSERT(lowerBounds.size() == x.size(), "");

    for (int i = 0; static_cast<typename gsl::span<const T>::index_type>(i) < x.size(); ++i)
        x[i] = std::min(std::max(lowerBounds[i], x[i]), upperBounds[i]);
}

} // namespace parpe

#endif
