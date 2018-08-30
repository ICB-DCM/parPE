#ifndef PARPE_OPTIMIZATION_MINIBATCH_OPTIMIZATION_H
#define PARPE_OPTIMIZATION_MINIBATCH_OPTIMIZATION_H

#include <vector>

#include <cassert>
#include <cmath>
#include <random>
#include <algorithm>
#include <model.h>
#include <iostream>
#include <misc.h>
#include <optimizationOptions.h>
#include <optimizationProblem.h>

namespace parpe {

enum class minibatchExitStatus {
    gradientNormConvergence,
    maxEpochsExceeded,
    invalidNumber
};


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
 * Split a vector into batches of the given size. All but possibly the last one will be of Size batchSize.
 */
template<typename T>
std::vector<std::vector<T>> getBatches(gsl::span<const T> data, int batchSize) {
    int numBatches = ceil(static_cast<double>(data.size()) / batchSize);

    std::vector<std::vector<T>> batches(numBatches, std::vector<T>());
    for(int i = 0, batchIdx = -1; (unsigned) i < data.size(); ++i) {
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
double getVectorNorm(std::vector<double> const& v);


class ParameterUpdater {
public:
    virtual void updateParameters(gsl::span<const double> gradient, gsl::span<double> parameters) = 0;
    virtual ~ParameterUpdater() = default;
};

class ParameterUpdaterVanilla : public ParameterUpdater {
public:
    ParameterUpdaterVanilla() = default;
    ParameterUpdaterVanilla(double learningRate) : learningRate(learningRate) {}

    void updateParameters(gsl::span<const double> gradient, gsl::span<double> parameters) {
        int numParameters = gradient.size();
        for(int i = 0; i < numParameters; ++i)
            parameters[i] -= learningRate * gradient[i];
    }

    double learningRate = 1.0;
};

class ParameterUpdaterRmsProp : public ParameterUpdater {
public:
    ParameterUpdaterRmsProp() = default;
    ParameterUpdaterRmsProp(double learningRate) : learningRate(learningRate) {}

    void updateParameters(gsl::span<const double> gradient, gsl::span<double> parameters);

    int updates = 0;
    double learningRate = 1.0;
    double decayRate = 0.9;
    double eps = 1e-4;
    std::vector<double> gradientCache;
};



template<typename DATUM>
class MinibatchOptimizer {
public:
    std::tuple<int, double, std::vector<double> > optimize(
            SummedGradientFunction<DATUM> const& f,
            gsl::span<const DATUM> data,
            gsl::span<const double> initialParameters)
    {
        std::vector<double> parameters(initialParameters.begin(), initialParameters.end());
        std::vector<DATUM> shuffledData(data.begin(), data.end());

        int numParameters = f.numParameters();
        assert((unsigned) numParameters == parameters.size());
        std::vector<double> gradient(numParameters, 0.0);
        std::random_device rd;
        std::mt19937 rng(rd());
        minibatchExitStatus status = minibatchExitStatus::gradientNormConvergence;
        double cost = NAN;

        for(int epoch = 0; epoch < maxEpochs; ++epoch) {
            std::cout<<"Epoch "<<epoch<<std::endl;

            std::shuffle(shuffledData.begin(), shuffledData.end(), rng);
            auto batches = getBatches<DATUM>(shuffledData, batchSize);

            for(int batchIdx = 0; (unsigned) batchIdx < batches.size(); ++batchIdx) {
                auto status = f.evaluate(parameters, batches[batchIdx], cost, gradient, nullptr, nullptr);
                std::cout<<"\tBatch "<<batchIdx<<": p: "<<parameters
                        <<" Cost: "<<cost<<" Gradient:"<<gradient
                       << " |g|2: "<<getVectorNorm(gradient)
                       << " Batch: "<<batches[batchIdx]<<std::endl;
                // TODO: do something smarter
                if(status == functionEvaluationFailure) {
                    std::cout<<"Minibatch cost function evaluation failed."<<std::endl;
                    return std::tuple<int, double, std::vector<double> >((int)minibatchExitStatus::invalidNumber, cost, parameters);
                }

                parameterUpdater->updateParameters(gradient, parameters);
            }

            if(getVectorNorm(gradient) <= gradientNormThreshold) {
                std::cout<<"Convergence: gradientNormThreshold reached."<<std::endl;
                break;
            }

        }

        return std::tuple<int, double, std::vector<double> >((int)status, cost, parameters);

    }

    std::unique_ptr<ParameterUpdater> parameterUpdater = std::make_unique<ParameterUpdaterVanilla>();

    int batchSize = 1;
    int maxEpochs = 10;
    double gradientNormThreshold = 0.0;
};

void setMinibatchOption(const std::pair<const std::string, const std::string> &pair,
                        MinibatchOptimizer<int>* optimizer);


template<typename DATUM>
std::unique_ptr<MinibatchOptimizer<DATUM>> getMinibatchOptimizer(OptimizationOptions const& options) {
    auto optim = std::make_unique<MinibatchOptimizer<int>>();

    options.for_each<MinibatchOptimizer<int>*>(setMinibatchOption, optim.get());

    return optim;
}

std::tuple<int, double, std::vector<double> > runMinibatchOptimization(OptimizationProblem *problem);


} // namespace parpe

#endif
