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

void setMinibatchOption(const std::pair<const std::string, const std::string> &pair, MinibatchOptimizer<int> *optimizer) {
    const std::string &key = pair.first;
    const std::string &val = pair.second;

    if(key == "maxEpochs") {
        optimizer->maxEpochs = std::stoi(val);
    } else if(key == "batchSize") {
        optimizer->batchSize = std::stoi(val);
    } else if(key == "gradientNormThreshold") {
        optimizer->gradientNormThreshold = std::stod(val);
    } else if(key == "rescueInterceptor") {
        optimizer->interceptor = std::stoi(val);
    } else if(key == "parameterUpdater") {
        if(val == "Vanilla") {
            // already default optimizer->parameterUpdater = std::make_unique<ParameterUpdaterVanilla>();
        } else if(val == "RmsProp" && !dynamic_cast<ParameterUpdaterRmsProp*>(optimizer->parameterUpdater.get())) {
            // this might have been set previously if there was an updater-specific option before
            optimizer->parameterUpdater = std::make_unique<ParameterUpdaterRmsProp>();
        } else if(val == "Adam" && !dynamic_cast<ParameterUpdaterAdam*>(optimizer->parameterUpdater.get())) {
            // this might have been set previously if there was an updater-specific option before
            optimizer->parameterUpdater = std::make_unique<ParameterUpdaterAdam>();
        } else {
            logmessage(LOGLVL_WARNING, "Ignoring unknown Minibatch parameterUpdater %s.", val.c_str());
        }
    } else if(key == "learningRateInterpMode") {
        if(val == "linear") {
		    optimizer->learningRateUpdater = std::make_unique<LearningRateUpdater>(optimizer->maxEpochs, parpe::learningRateInterp::linear );
        } else if(val == "inverseLinear") {
            optimizer->learningRateUpdater = std::make_unique<LearningRateUpdater>(optimizer->maxEpochs, parpe::learningRateInterp::inverseLinear );
        } else if(val == "logarithmic") {
			optimizer->learningRateUpdater = std::make_unique<LearningRateUpdater>(optimizer->maxEpochs, parpe::learningRateInterp::logarithmic );
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



/**
 * LearningRateUpdater
 * The LearningRateUpdater provides the possibility to reduce the learning rate per epoch
 * and makes it possible to adapt the learning rate accroding to succes or failure of
 * the ODE solver.
 */
void LearningRateUpdater::updateLearningRate(int currentEpoch) {

	// Depending on the interpolation mode the current learning rate computed must be...
	if(learningRateInterpMode == learningRateInterp::linear) {
    	currentLearningRate = startLearningRate - (startLearningRate - endLearningRate) * ((double) currentEpoch) / ((double) maxEpochs - 1);
	} else if(learningRateInterpMode == learningRateInterp::inverseLinear) {
		currentLearningRate = 1 / startLearningRate - (1 / startLearningRate - 1 / endLearningRate) * ((double) currentEpoch) / ((double) maxEpochs - 1);
		currentLearningRate = 1 / currentLearningRate;
	} else if(learningRateInterpMode == learningRateInterp::logarithmic) {
		currentLearningRate = log(startLearningRate) - (log(startLearningRate) - log(endLearningRate)) * ((double) currentEpoch) / ((double) maxEpochs - 1);
		currentLearningRate = exp(currentLearningRate);
	}
}

double LearningRateUpdater::getCurrentLearningRate() {
	return currentLearningRate * reductionFactor;
}

void LearningRateUpdater::reduceLearningRate() {
	reductionFactor = reductionFactor / 5;
}

void LearningRateUpdater::increaseLearningRate() {
	reductionFactor = std::min(reductionFactor * 1.3, 1.0);
}

void LearningRateUpdater::setReductionFactor(double newReductionFactor) {
	reductionFactor = newReductionFactor;
}

void LearningRateUpdater::setMaxEpochs(int newMaxEpochs) {
	maxEpochs = newMaxEpochs;
}


/* 
	Minibatch optimizer: RMSProp Updater
	The RMSProp updater currently takes two inputs:
		* start value of the learning rate
		* end value of the learning rate
*/
void ParameterUpdaterRmsProp::initialize(unsigned int numParameters)
{
    gradientNormCache.resize(numParameters);
    oldGradientNormCache.resize(numParameters);
}

void ParameterUpdaterRmsProp::updateParameters(double learningRate,
											   int iteration,
					                           gsl::span<double const> gradient,
                                               gsl::span<double> parameters,
                                               gsl::span<const double> lowerBounds,
                                               gsl::span<const double> upperBounds)
{

    int numParameters = gradient.size();
    
    oldGradientNormCache = gradientNormCache;
    
    for(int i = 0; i < numParameters; ++i) {
        gradientNormCache[i] = decayRate * gradientNormCache[i]
                + (1 - decayRate) * gradient[i] * gradient[i];

        parameters[i] += - learningRate * gradient[i] / (std::sqrt(gradientNormCache[i]) + delta);
    }

    clipToBounds(lowerBounds, upperBounds, parameters);
}

void ParameterUpdaterRmsProp::undoLastStep()
{
	// The cached gradient norm needs to be restored, since the new one is probably NaN
    gradientNormCache = oldGradientNormCache;
}

void ParameterUpdaterRmsProp::clearCache()
{
	// Reset all cached information
	int numParameters = gradientNormCache.size();
    for(int ip = 0; ip < numParameters; ++ip) {
        gradientNormCache[ip] = 0.0;
        oldGradientNormCache[ip] = 0.0;
    }
}



/* 
	Minibatch optimizer: Adam Updater
	The Adam updater currently takes two inputs:
		* start value of the learning rate
		* end value of the learning rate
*/
void ParameterUpdaterAdam::initialize(unsigned int numParameters)
{
    gradientCache.resize(numParameters);
    gradientNormCache.resize(numParameters);
    oldGradientCache.resize(numParameters);
    oldGradientNormCache.resize(numParameters);
}

void ParameterUpdaterAdam::updateParameters(double learningRate,
											int iteration,
					                        gsl::span<double const> gradient,
                                            gsl::span<double> parameters,
                                            gsl::span<const double> lowerBounds,
                                            gsl::span<const double> upperBounds)
{

    int numParameters = gradient.size();
    double tmpNumerator;
    double tmpDenominator;
    
    oldGradientNormCache = gradientNormCache;
    oldGradientCache = gradientCache;
    
    for(int i = 0; i < numParameters; ++i) {
        gradientCache[i] = decayRateGradient * gradientCache[i]
                + (1 - decayRateGradient) * gradient[i];
        gradientNormCache[i] = decayRateGradientNorm * gradientNormCache[i]
                + (1 - decayRateGradientNorm) * gradient[i] * gradient[i];

        tmpNumerator = gradientCache[i] / (1 - std::pow(decayRateGradient, (double)iteration));
        tmpDenominator = std::sqrt(gradientNormCache[i] / (1 - std::pow(decayRateGradientNorm, (double)iteration))) + delta;
        
        parameters[i] += - learningRate * tmpNumerator / tmpDenominator;
    }

    clipToBounds(lowerBounds, upperBounds, parameters);
}

void ParameterUpdaterAdam::undoLastStep()
{
	// The cached gradient norm needs to be restored, since the new one is probably NaN
    gradientNormCache = oldGradientNormCache;
    gradientCache = oldGradientCache;
}

void ParameterUpdaterAdam::clearCache()
{
	// Reset all cached information
	int numParameters = gradientNormCache.size();
    for(int ip = 0; ip < numParameters; ++ip) {
        gradientNormCache[ip] = 0.0;
        oldGradientNormCache[ip] = 0.0;
        gradientCache[ip] = 0.0;
        oldGradientCache[ip] = 0.0;
    }
}



/* 
	Minibatch optimizer: Vanilla SGD Updater
	The Vanilla SGD updater currently takes two inputs:
		* start value of the learning rate
		* end value of the learning rate
*/
void ParameterUpdaterVanilla::updateParameters(double learningRate,
											   int iteration,
                                               gsl::span<const double> gradient,
                                               gsl::span<double> parameters,
                                               gsl::span<const double> lowerBounds,
                                               gsl::span<const double> upperBounds)
{
    int numParameters = gradient.size();

    for(int i = 0; i < numParameters; ++i) {
        // logmessage(LOGLVL_DEBUG, "p_%d: %f - %f = %f", i, parameters[i], learningRate * gradient[i], parameters[i] - learningRate * gradient[i]);
        parameters[i] -= learningRate * gradient[i] / getVectorNorm(gradient);
    }

    clipToBounds(lowerBounds, upperBounds, parameters);
}

}
