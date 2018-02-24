#include "model.h"

#include <cmath>

namespace parpe {


template<typename X>
void Model<X>::evaluate(const double *parameters, const std::vector<X> &features, std::vector<double> &outputs) const {
    auto unusedGrad = std::vector<std::vector<double>>();
    evaluate(parameters, features, outputs, unusedGrad);

}

void LinearModel::evaluate(const double *parameters, const std::vector<std::vector<double> > &features, std::vector<double> &outputs, std::vector<std::vector<double> > &outputGradients) const {

    const int numObservations = features.size();
    const int numFeatures = features[0].size();
    const int numParams = numFeatures + 1;
    const int idxOffset = numParams - 1;

    for(int i = 0; i < numObservations; ++i) {
        outputs[i] = 0.0;
        for(int j = 0; j < numFeatures; ++j) {
            outputs[i] += features[i][j] * parameters[j];
        }
        outputs[i] += parameters[idxOffset];
    }

    if(outputGradients.size()) {
        // Simplify: [A, 1.0]
        // for each observation
        for(int i = 0; i < numObservations; ++i) {
            // for each parameter
            for(int j = 0; j < numParams - 1; ++j) {
                outputGradients[i][j] = 0.0;
                // for each feature
                outputGradients[i][j] += features[i][j];
            }
            outputGradients[i][idxOffset] = 1.0; // offset
        }
    }
}

FunctionEvaluationStatus LinearModelMSE::evaluate(const double * const parameters,
                                                std::vector<int> dataIndices, double &fval, double *gradient) const  {

    int numDatasets = dataIndices.size();

    // get data for indices
    std::vector<std::vector<double>> data(numDatasets);
    for(int i = 0; (unsigned) i < data.size(); ++i)
        data[i] = datasets[dataIndices[i]];

    // evaluate
    std::vector<double> outputs(numDatasets);
    std::vector<std::vector<double>>
            outputGradients(numDatasets,
                            std::vector<double>(numParameters(), NAN));
    lm.evaluate(parameters, data, outputs, outputGradients);

    // compute MSE
    fval = 0.0;
    for(int i = 0; i < numDatasets; ++i) {
        fval += std::pow(labels[dataIndices[i]] - outputs[i], 2.0);
    }
    fval /= numDatasets;

    // and MSE gradient
    if(gradient) {
        for(int p = 0; p < numParameters(); ++p) {
            gradient[p] = 0.0;
            for(int i = 0; i < numDatasets; ++i) {
                gradient[p] +=  -2.0 * outputGradients[i][p] * (labels[dataIndices[i]] - outputs[i]);
            }
            gradient[p] /= numDatasets;
        }
    }

    return functionEvaluationFailure;
}

}
