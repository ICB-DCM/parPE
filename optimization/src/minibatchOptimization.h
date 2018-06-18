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

namespace parpe {

enum class minibatchExitStatus {
    gradientNormConvergence,
    maxEpochsExceeded,
};

/**
 * Split a vector into batches of the given size. All but possibly the last one will be of Size batchSize.
 */
template<typename DATUM>
std::vector<std::vector<DATUM>> getBatches(std::vector<DATUM> const& data, int batchSize) {
    int numBatches = ceil(static_cast<double>(data.size()) / batchSize);

    std::vector<std::vector<DATUM>> batches(numBatches, std::vector<DATUM>());
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
double getVectorNorm(std::vector<double> const& v) {
    double norm = 0.0;
    for(auto e : v)
        norm += std::pow(e, 2.0);
    norm = std::sqrt(norm);
    return norm;
}


class ParameterUpdater {
public:
    virtual void updateParameters(std::vector<double> const& gradient, std::vector<double>& parameters) = 0;
    virtual ~ParameterUpdater() = default;
};

class ParameterUpdaterVanilla : public ParameterUpdater {
public:
    ParameterUpdaterVanilla() = default;
    ParameterUpdaterVanilla(double learningRate) : learningRate(learningRate) {}

    void updateParameters(std::vector<double> const& gradient, std::vector<double>& parameters) {
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

    void updateParameters(std::vector<double> const& gradient, std::vector<double>& parameters) {
        if(++updates % 10 == 0) {
            std::cout<<"Adapting learning rate "<<learningRate;
            learningRate /= 2.0;
            std::cout<<" to "<<learningRate<<std::endl;
        }

        int numParameters = gradient.size();
        gradientCache.resize(numParameters);

        for(int i = 0; i < numParameters; ++i) {
            gradientCache[i] = decayRate * gradientCache[i]
                    + (1 - decayRate) * gradient[i] * gradient[i];

            parameters[i] += - learningRate * gradient[i] / (std::sqrt(gradientCache[i]) + eps);
        }

    }
    int updates = 0;
    double learningRate = 1.0;
    double decayRate = 0.9;
    double eps = 1e-4;
    std::vector<double> gradientCache;
};



template<typename DATUM>
class MinibatchOptimizer {
public:
    std::tuple<int, double, std::vector<double> > optimize(SummedGradientFunction<DATUM> const& f, std::vector<DATUM> data, int batchSize, std::vector<double> parameters) {
        int numParameters = f.numParameters();
        assert((unsigned) numParameters == parameters.size());
        std::vector<double> gradient(numParameters, 0.0);
        std::random_device rd;
        std::mt19937 rng(rd());
        minibatchExitStatus status = minibatchExitStatus::gradientNormConvergence;
        double cost = NAN;

        for(int epoch = 0; epoch < maxEpochs; ++epoch) {
            std::cout<<"Epoch "<<epoch<<std::endl;

            std::shuffle(data.begin(), data.end(), rng);
            auto batches = getBatches(data, batchSize);

            for(int batchIdx = 0; (unsigned) batchIdx < batches.size(); ++batchIdx) {
                f.evaluate(parameters, batches[batchIdx], cost, gradient);
                std::cout<<"\tBatch "<<batchIdx<<": p: "<<parameters
                        <<" Cost: "<<cost<<" Gradient:"<<gradient
                       << " |g|2: "<<getVectorNorm(gradient)
                       << " Batch: "<<batches[batchIdx]<<std::endl;
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

    int maxEpochs = 10;
    double gradientNormThreshold = 0.0;
};




}

#endif
