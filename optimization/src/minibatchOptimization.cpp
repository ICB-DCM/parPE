#include "minibatchOptimization.h"


namespace parpe {

double getVectorNorm(const std::vector<double> &v) {
    double norm = 0.0;
    for(auto e : v)
        norm += std::pow(e, 2.0);
    norm = std::sqrt(norm);
    return norm;
}

void ParameterUpdaterRmsProp::updateParameters(gsl::span<double const> gradient, gsl::span<double> parameters) {
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
    } else {
        logmessage(LOGLVL_WARNING, "Ignoring unknown optimization option %s.", key.c_str());
        return;
    }

    logmessage(LOGLVL_DEBUG, "Set optimization option %s to %s.", key.c_str(), val.c_str());
}

std::tuple<int, double, std::vector<double> > runMinibatchOptimization(OptimizationProblem *problem) {
    auto minibatchProblem = dynamic_cast<MinibatchOptimizationProblem<int>*>(problem);
    if(!minibatchProblem)
        throw ParPEException("Minibatch optimizer selected but given optimization problem cannot be solved by minibatch optimizer");

    auto minibatchOptimizer = getMinibatchOptimizer<int>(problem->getOptimizationOptions());
    std::vector<double> parameters(problem->costFun->numParameters());
    problem->fillInitialParameters(parameters);

    auto costFun = minibatchProblem->getGradientFunction();
    auto data = minibatchProblem->getTrainingData();

    return minibatchOptimizer->optimize(*costFun, data, parameters);
}

}
