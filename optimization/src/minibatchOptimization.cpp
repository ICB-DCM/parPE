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
    } else if(key == "parameterUpdater") {
        if(val == "Vanilla") {
            // already default optimizer->parameterUpdater = std::make_unique<ParameterUpdaterVanilla>();
        } else if(val == "RmsProp" && !dynamic_cast<ParameterUpdaterRmsProp*>(optimizer->parameterUpdater.get())) {
            // this might have been set previously if there was an updater-specific option before
            optimizer->parameterUpdater = std::make_unique<ParameterUpdaterRmsProp>();
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
    	currentLearningRate = startLearningRate - (startLearningRate - endLearningRate) * ((double) currentEpoch - 1) / ((double) maxEpochs - 1);
	} else if(learningRateInterpMode == learningRateInterp::inverseLinear) {
		currentLearningRate = 1 / startLearningRate - (1 / startLearningRate - 1 / endLearningRate) * ((double) currentEpoch - 1) / ((double) maxEpochs - 1);
		currentLearningRate = 1 / currentLearningRate;
	} else if(learningRateInterpMode == learningRateInterp::logarithmic) {
		currentLearningRate = log(startLearningRate) - (log(startLearningRate) - log(endLearningRate)) * ((double) currentEpoch - 1) / ((double) maxEpochs - 1);
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



/* 
	Minibatch optimizer: RMSProp Updater
	The RMSProp updater currently takes two inputs:
		* start value of the learning rate
		* end value of the learning rate
*/
void ParameterUpdaterRmsProp::updateParameters(double learningRate,
					                           gsl::span<double const> gradient,
                                               gsl::span<double> parameters,
                                               gsl::span<const double> lowerBounds,
                                               gsl::span<const double> upperBounds)
{

    int numParameters = gradient.size();
    gradientNormCache.resize(numParameters);

    for(int i = 0; i < numParameters; ++i) {
        gradientNormCache[i] = decayRate * gradientNormCache[i]
                + (1 - decayRate) * gradient[i] * gradient[i];

        parameters[i] += - learningRate * gradient[i] / (std::sqrt(gradientNormCache[i]) + delta);
    }

    clipToBounds(lowerBounds, upperBounds, parameters);
}



/* 
	Minibatch optimizer: Vanilla SGD Updater
	The Vanilla SGD updater currently takes two inputs:
		* start value of the learning rate
		* end value of the learning rate
*/
void ParameterUpdaterVanilla::updateParameters(double learningRate,
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
