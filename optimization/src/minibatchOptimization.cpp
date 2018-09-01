#include "minibatchOptimization.h"

#include <parpeException.h>

namespace parpe {

double getVectorNorm(gsl::span<const double> v) {
    double norm = 0.0;
    for(auto e : v)
        norm += std::pow(e, 2.0);
    norm = std::sqrt(norm);
    return norm;
}

ParameterUpdaterRmsProp::ParameterUpdaterRmsProp(double learningRate) : learningRate(learningRate) {}

void ParameterUpdaterRmsProp::updateParameters(gsl::span<double const> gradient,
                                               gsl::span<double> parameters,
                                               gsl::span<const double> lowerBounds,
                                               gsl::span<const double> upperBounds)
{
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

    clipToBounds(lowerBounds, upperBounds, parameters);
}

void setMinibatchOption(const std::pair<const std::string, const std::string> &pair, MinibatchOptimizer<int> *optimizer) {
    const std::string &key = pair.first;
    const std::string &val = pair.second;

    if(key == "maxEpochs") {
        optimizer->maxEpochs = std::stoi(val);
    } else if(key == "batchSize") {
        optimizer->batchSize = std::stoi(val);
    } else if(key == "gradientNormThreshold") {
        optimizer->gradientNormThreshold = std::stod(val);
    } else if(key == "parameterUpdater") {
        if(val == "Vanilla") {
            // already default optimizer->parameterUpdater = std::make_unique<ParameterUpdaterVanilla>();
        } else if(val == "RmsProp" && !dynamic_cast<ParameterUpdaterRmsProp*>(optimizer->parameterUpdater.get())) {
            // this might have been set previously if there was an updater-specific option before
            optimizer->parameterUpdater = std::make_unique<ParameterUpdaterRmsProp>();
        } else {
            logmessage(LOGLVL_WARNING, "Ignoring unknown Minibatch parameterUpdater %s.", val.c_str());
        }
    } else if(key == "RmsProp-learningRate") {
        if(!dynamic_cast<ParameterUpdaterRmsProp*>(optimizer->parameterUpdater.get())) {
            optimizer->parameterUpdater = std::make_unique<ParameterUpdaterRmsProp>();
        }
        auto updater = dynamic_cast<ParameterUpdaterRmsProp*>(optimizer->parameterUpdater.get());
        RELEASE_ASSERT(updater, "");

        updater->learningRate = std::stod(val);
    } else if(key == "Vanilla-learningRate") {
        if(auto updater = dynamic_cast<ParameterUpdaterVanilla*>(optimizer->parameterUpdater.get())) {
            updater->learningRate = std::stod(val);
        }
    } else {
        logmessage(LOGLVL_WARNING, "Ignoring unknown optimization option %s.", key.c_str());
        return;
    }

    logmessage(LOGLVL_DEBUG, "Set optimization option %s to %s.", key.c_str(), val.c_str());
}

std::tuple<int, double, std::vector<double> > runMinibatchOptimization(MinibatchOptimizationProblem<int> *problem)
{
    auto minibatchOptimizer = getMinibatchOptimizer<int>(problem->getOptimizationOptions());

    auto costFun = problem->getGradientFunction();

    std::vector<double> initialParameters(costFun->numParameters());
    problem->fillInitialParameters(initialParameters);

    std::vector<double> lowerParameterBounds(costFun->numParameters());
    problem->fillParametersMin(lowerParameterBounds);

    std::vector<double> upperParameterBounds(costFun->numParameters());
    problem->fillParametersMax(upperParameterBounds);

    auto data = problem->getTrainingData();

    return minibatchOptimizer->optimize(*costFun, data, initialParameters,
                                        lowerParameterBounds, upperParameterBounds,
                                        problem->getReporter().get(),
                                        problem->logger.get());
}

ParameterUpdaterVanilla::ParameterUpdaterVanilla(double learningRate) : learningRate(learningRate) {}

void ParameterUpdaterVanilla::updateParameters(gsl::span<const double> gradient,
                                               gsl::span<double> parameters,
                                               gsl::span<const double> lowerBounds,
                                               gsl::span<const double> upperBounds)
{
    int numParameters = gradient.size();
    for(int i = 0; i < numParameters; ++i) {
        // logmessage(LOGLVL_DEBUG, "p_%d: %f - %f = %f", i, parameters[i], learningRate * gradient[i], parameters[i] - learningRate * gradient[i]);
        parameters[i] -= learningRate * gradient[i];
    }

    clipToBounds(lowerBounds, upperBounds, parameters);
}

}
