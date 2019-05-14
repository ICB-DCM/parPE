#include <parpeoptimization/minibatchOptimization.h>

#include <parpecommon/parpeException.h>

namespace parpe {

double getScalarProduct(gsl::span<const double> v,
                        gsl::span<const double> w) {
    double scalarProduct = 0.0;
    for (unsigned int i = 0; i < v.size(); ++i)
        scalarProduct += v[i] * w[i];
    return scalarProduct;
}

double getVectorNorm(gsl::span<const double> v) {
    return std::sqrt(getScalarProduct(v, v));
}

std::vector<double> getVectorDifference(gsl::span<const double> v,
                                        gsl::span<const double> w) {
    Expects(v.size() == w.size());
    std::vector<double> difference(v.size(), 0.0);
    for (unsigned int i = 0; i < v.size(); ++i)
        difference[i] = v[i] - w[i];
    
    return difference;
}

void setMinibatchOption(const std::pair<const std::string, const std::string> &pair,
                        MinibatchOptimizer<int> *optimizer) {
    const std::string &key = pair.first;
    const std::string &val = pair.second;

    /* Get options from h5-file */
    if (key == "maxEpochs") {
        optimizer->maxEpochs = std::stoi(val);
    } else if (key == "batchSize") {
        optimizer->batchSize = std::stoi(val);
    } else if (key == "gradientNormThreshold") {
        optimizer->gradientNormThreshold = std::stod(val);
    } else if (key == "lineSearchSteps") {
        optimizer->lineSearchSteps = std::stoi(val);
    } else if (key == "rescueInterceptor") {
        if (val == "none" or val == "0") {
            optimizer->interceptor = parpe::interceptType::none;
        } else if (val == "reduceStep" or val == "1") {
            optimizer->interceptor = parpe::interceptType::reduceStep;
        } else if (val == "reduceStepAndRestart" or val == "2") {
            optimizer->interceptor = parpe::interceptType::reduceStepAndRestart;
        }
    } else if (key == "parameterUpdater") {
        if (val == "Vanilla") {
            // already default optimizer->parameterUpdater = std::make_unique<ParameterUpdaterVanilla>();
        } else if (val == "RmsProp" && !dynamic_cast<ParameterUpdaterRmsProp*>(optimizer->parameterUpdater.get())) {
            // this might have been set previously if there was an updater-specific option before
            optimizer->parameterUpdater = std::make_unique<ParameterUpdaterRmsProp>();
        } else if (val == "Adam" && !dynamic_cast<ParameterUpdaterAdam*>(optimizer->parameterUpdater.get())) {
            // this might have been set previously if there was an updater-specific option before
            optimizer->parameterUpdater = std::make_unique<ParameterUpdaterAdam>();
        } else if (val == "AdamClassic" && !dynamic_cast<ParameterUpdaterAdamClassic*>(optimizer->parameterUpdater.get())) {
            // this might have been set previously if there was an updater-specific option before
            optimizer->parameterUpdater = std::make_unique<ParameterUpdaterAdamClassic>();
        } else {
            logmessage(LOGLVL_WARNING, "Ignoring unknown Minibatch parameterUpdater %s.", val.c_str());
        }
    } else if (key == "learningRateInterpMode") {
        if (val == "linear") {
            optimizer->learningRateUpdater = std::make_unique < LearningRateUpdater
                    > (optimizer->maxEpochs, parpe::learningRateInterp::linear);
        } else if (val == "inverseLinear") {
            optimizer->learningRateUpdater = std::make_unique < LearningRateUpdater
                    > (optimizer->maxEpochs, parpe::learningRateInterp::inverseLinear);
        } else if (val == "logarithmic") {
            optimizer->learningRateUpdater = std::make_unique < LearningRateUpdater
                    > (optimizer->maxEpochs, parpe::learningRateInterp::logarithmic);
        }
    } else if (key == "startLearningRate") {
            optimizer->learningRateUpdater->setStartLearningRate(std::stod(val));
    } else if (key == "endLearningRate") {
            optimizer->learningRateUpdater->setEndLearningRate(std::stod(val));
    } else {
        logmessage(LOGLVL_WARNING, "Ignoring unknown optimization option %s.", key.c_str());
        return;
    }

    logmessage(LOGLVL_DEBUG, "Set optimization option %s to %s.", key.c_str(), val.c_str());
}

std::tuple<int, double, std::vector<double> > runMinibatchOptimization(MinibatchOptimizationProblem<int> *problem) {
    auto minibatchOptimizer = getMinibatchOptimizer<int>(problem->getOptimizationOptions());

    auto costFun = problem->getGradientFunction();

    std::vector<double> initialParameters(costFun->numParameters());
    problem->fillInitialParameters(initialParameters);

    std::vector<double> lowerParameterBounds(costFun->numParameters());
    problem->fillParametersMin(lowerParameterBounds);

    std::vector<double> upperParameterBounds(costFun->numParameters());
    problem->fillParametersMax(upperParameterBounds);

    auto data = problem->getTrainingData();

    return minibatchOptimizer->optimize(*costFun, data, initialParameters, lowerParameterBounds, upperParameterBounds,
                                        problem->getReporter().get(), problem->logger.get());
}


LearningRateUpdater::LearningRateUpdater(int maxEpochs, learningRateInterp learningRateInterpMode)
    :maxEpochs(maxEpochs),
      learningRateInterpMode(learningRateInterpMode)
{
}

void LearningRateUpdater::updateLearningRate(int currentEpoch) {

    // Depending on the interpolation mode the current learning rate computed must be...
    if (learningRateInterpMode == learningRateInterp::linear) {
        currentLearningRate = startLearningRate
                - (startLearningRate - endLearningRate) * ((double) currentEpoch) / ((double) maxEpochs - 1);
    } else if (learningRateInterpMode == learningRateInterp::inverseLinear) {
        currentLearningRate = 1 / startLearningRate
                - (1 / startLearningRate - 1 / endLearningRate) * ((double) currentEpoch) / ((double) maxEpochs - 1);
        currentLearningRate = 1 / currentLearningRate;
    } else if (learningRateInterpMode == learningRateInterp::logarithmic) {
        currentLearningRate = log(startLearningRate)
                - (log(startLearningRate) - log(endLearningRate)) * ((double) currentEpoch) / ((double) maxEpochs - 1);
        currentLearningRate = exp(currentLearningRate);
    }
}

double LearningRateUpdater::getCurrentLearningRate() {
    return currentLearningRate * reductionFactor;
}

void LearningRateUpdater::reduceLearningRate() {
    double learningRateReduction = 5.0;
    reductionFactor = reductionFactor / learningRateReduction;
}

void LearningRateUpdater::increaseLearningRate() {
    double learningRateIncrease = 1.3;
    reductionFactor = std::min(reductionFactor * learningRateIncrease, 1.0);
}

void LearningRateUpdater::setReductionFactor(double newReductionFactor) {
    reductionFactor = newReductionFactor;
}

void LearningRateUpdater::setMaxEpochs(int newMaxEpochs) {
    maxEpochs = newMaxEpochs;
}

void LearningRateUpdater::setStartLearningRate(double learningRate)
{
    startLearningRate = learningRate;
}

void LearningRateUpdater::setEndLearningRate(double learningRate)
{
    endLearningRate = learningRate;
}

void ParameterUpdaterRmsProp::initialize(unsigned int numParameters) {
    gradientNormCache.resize(numParameters);
    oldGradientNormCache.resize(numParameters);
    std::fill(gradientNormCache.begin(), gradientNormCache.end(), 0.0);
    std::fill(oldGradientNormCache.begin(), oldGradientNormCache.end(), 0.0);
}

void ParameterUpdaterRmsProp::updateParameters(double learningRate,
                                               int  /*iteration*/,
                                               gsl::span<double const> gradient,
                                               gsl::span<double> parameters,
                                               gsl::span<const double> lowerBounds,
                                               gsl::span<const double> upperBounds) {

    int numParameters = gradient.size();

    oldGradientNormCache = gradientNormCache;

    for (int i = 0; i < numParameters; ++i) {
        gradientNormCache[i] = decayRate * gradientNormCache[i] + (1 - decayRate) * gradient[i] * gradient[i];

        parameters[i] += -learningRate * gradient[i] / (std::sqrt(gradientNormCache[i]) + delta);
    }

    clipToBounds(lowerBounds, upperBounds, parameters);
}

void ParameterUpdaterRmsProp::undoLastStep() {
    // The cached gradient norm needs to be restored, since the new one is probably NaN
    gradientNormCache = oldGradientNormCache;
}

void ParameterUpdaterRmsProp::clearCache() {
    // Reset all cached information
    std::fill(gradientNormCache.begin(), gradientNormCache.end(), 0.0);
    std::fill(oldGradientNormCache.begin(), oldGradientNormCache.end(), 0.0);
}

void ParameterUpdaterMomentum::initialize(unsigned int numParameters) {
    momentum.resize(numParameters);
    oldMomentum.resize(numParameters);
    std::fill(momentum.begin(), momentum.end(), 0.0);
    std::fill(oldMomentum.begin(), oldMomentum.end(), 0.0);
}

void ParameterUpdaterMomentum::updateParameters(double learningRate,
                                               int  /*iteration*/,
                                               gsl::span<double const> gradient,
                                               gsl::span<double> parameters,
                                               gsl::span<const double> lowerBounds,
                                               gsl::span<const double> upperBounds) {

    int numParameters = gradient.size();

    oldMomentum = momentum;

    for (int i = 0; i < numParameters; ++i) {
        momentum[i] = decayRate * momentum[i] + (1 - decayRate) * gradient[i];
        
        momentumNormalize = std::max(getVectorNorm(momentum), 1.0)
        
        parameters[i] += -learningRate * momentum[i] / momentumNormalize;
    }

    clipToBounds(lowerBounds, upperBounds, parameters);
}

void ParameterUpdaterMomentum::undoLastStep() {
    // The cached gradient norm needs to be restored, since the new one is probably NaN
    momentum = oldMomentum;
}

void ParameterUpdaterMomentum::clearCache() {
    // Reset all cached information
    std::fill(momentum.begin(), momentum.end(), 0.0);
}


void ParameterUpdaterAdam::initialize(unsigned int numParameters) {
    gradientCache.resize(numParameters);
    gradientNormCache.resize(numParameters);
    oldGradientCache.resize(numParameters);
    oldGradientNormCache.resize(numParameters);
    std::fill(gradientNormCache.begin(), gradientNormCache.end(), 0.0);
    std::fill(oldGradientNormCache.begin(), oldGradientNormCache.end(), 0.0);
    std::fill(gradientCache.begin(), gradientCache.end(), 0.0);
    std::fill(oldGradientCache.begin(), oldGradientCache.end(), 0.0);
}

void ParameterUpdaterAdam::updateParameters(double learningRate,
                                            int iteration,
                                            gsl::span<double const> gradient,
                                            gsl::span<double> parameters,
                                            gsl::span<const double> lowerBounds,
                                            gsl::span<const double> upperBounds) {

    int numParameters = gradient.size();
    double tmpNumerator;
    double tmpDenominator;

    oldGradientNormCache = gradientNormCache;
    oldGradientCache = gradientCache;
    
    for (int i = 0; i < numParameters; ++i) {        
        // compute new steps from last gradient information
        gradientCache[i] = decayRateGradient * gradientCache[i] + (1 - decayRateGradient) * gradient[i];
        gradientNormCache[i] = decayRateGradientNorm * gradientNormCache[i]
                + (1 - decayRateGradientNorm) * gradient[i] * gradient[i];

        tmpNumerator = gradientCache[i] / (1 - std::pow(decayRateGradient, (double) iteration));
        tmpDenominator = std::sqrt(gradientNormCache[i] / (1 - std::pow(decayRateGradientNorm, (double) iteration)))
                + delta;

        parameters[i] += -learningRate * tmpNumerator / tmpDenominator;
    }

    clipToBounds(lowerBounds, upperBounds, parameters);
}

void ParameterUpdaterAdam::undoLastStep() {
    // The cached gradient norm needs to be restored, since the new one is probably NaN
    gradientNormCache = oldGradientNormCache;
    gradientCache = oldGradientCache;
}

void ParameterUpdaterAdam::clearCache() {
    // Reset all cached information
    std::fill(gradientNormCache.begin(), gradientNormCache.end(), 0.0);
    std::fill(oldGradientNormCache.begin(), oldGradientNormCache.end(), 0.0);
    std::fill(gradientCache.begin(), gradientCache.end(), 0.0);
    std::fill(oldGradientCache.begin(), oldGradientCache.end(), 0.0);
}

void ParameterUpdaterAdamClassic::initialize(unsigned int numParameters) {
    gradientCache.resize(numParameters);
    gradientNormCache.resize(numParameters);
    oldGradientCache.resize(numParameters);
    oldGradientNormCache.resize(numParameters);
    std::fill(gradientNormCache.begin(), gradientNormCache.end(), 0.0);
    std::fill(oldGradientNormCache.begin(), oldGradientNormCache.end(), 0.0);
    std::fill(gradientCache.begin(), gradientCache.end(), 0.0);
    std::fill(oldGradientCache.begin(), oldGradientCache.end(), 0.0);
}

void ParameterUpdaterAdamClassic::updateParameters(double learningRate,
                                            int iteration,
                                            gsl::span<double const> gradient,
                                            gsl::span<double> parameters,
                                            gsl::span<const double> lowerBounds,
                                            gsl::span<const double> upperBounds) {

    int numParameters = gradient.size();
    double tmpNumerator;
    double tmpDenominator;

    oldGradientNormCache = gradientNormCache;
    oldGradientCache = gradientCache;
    
    for (int i = 0; i < numParameters; ++i) {        
        // compute new steps from last gradient information
        gradientCache[i] = decayRateGradient * gradientCache[i] + (1 - decayRateGradient) * gradient[i];
        gradientNormCache[i] = decayRateGradientNorm * gradientNormCache[i]
                + (1 - decayRateGradientNorm) * gradient[i] * gradient[i];

        tmpNumerator = gradientCache[i] / (1 - std::pow(decayRateGradient, (double) iteration));
        tmpDenominator = std::sqrt(gradientNormCache[i] / (1 - std::pow(decayRateGradientNorm, (double) iteration)))
                + delta;

        parameters[i] += -learningRate * tmpNumerator / tmpDenominator;
    }

    clipToBounds(lowerBounds, upperBounds, parameters);
}

void ParameterUpdaterAdamClassic::undoLastStep() {
    // The cached gradient norm needs to be restored, since the new one is probably NaN
    gradientNormCache = oldGradientNormCache;
    gradientCache = oldGradientCache;
}

void ParameterUpdaterAdamClassic::clearCache() {
    // Reset all cached information
    std::fill(gradientNormCache.begin(), gradientNormCache.end(), 0.0);
    std::fill(oldGradientNormCache.begin(), oldGradientNormCache.end(), 0.0);
    std::fill(gradientCache.begin(), gradientCache.end(), 0.0);
    std::fill(oldGradientCache.begin(), oldGradientCache.end(), 0.0);
}


void ParameterUpdaterVanilla::updateParameters(double learningRate,
                                               int  /*iteration*/,
                                               gsl::span<const double> gradient,
                                               gsl::span<double> parameters,
                                               gsl::span<const double> lowerBounds,
                                               gsl::span<const double> upperBounds) {
    int numParameters = gradient.size();
    double delta = 1e-8;

    for (int i = 0; i < numParameters; ++i) {
        // logmessage(LOGLVL_DEBUG, "p_%d: %f - %f = %f", i, parameters[i], learningRate * gradient[i], parameters[i] - learningRate * gradient[i]);
        parameters[i] -= learningRate * gradient[i] / (getVectorNorm(gradient) + delta);
    }

    clipToBounds(lowerBounds, upperBounds, parameters);
}

void ParameterUpdaterVanilla::initialize(unsigned int numParameters) {}

void ParameterUpdaterVanilla::clearCache() {}

void ParameterUpdaterVanilla::undoLastStep() {}

}
