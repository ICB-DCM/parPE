#include "optimizerMinibatch.h"
#include "optimizationProblem.h"

#include <algorithm>
#include <cmath>

OptimizerMinibatch::OptimizerMinibatch() {}

int OptimizerMinibatch::optimize(OptimizationProblem *model, void *data,
                                 int numBatches) {
    // TODO random if NULL
    std::vector<double> parameters(model->initialParameters,
                                   model->numOptimizationParameters);
    std::vector<double> gradient(model->numOptimizationParameters);

    // TODO numDatasets to OptimizationProblem
    batches = getBatches(model->numConditions, numBatches);

    int maxEpochs = 1;

    for (int epoch = 0; epoch < maxEpochs; ++epoch) {
        for (int batch = 0; batch < numBatches; ++batch) {

            // evaluate on batch TODO: batch argument
            model->evaluateObjectiveFunction(parameters.data(), &cost,
                                             gradient.data(),
                                             batches[batch].data());

            // update parameters with gradient
            updateParameters(parameters, gradient);
        }

        // log optimizer iteration
        model->intermediateFunction();
    }
}

// pointer to vector?
std::vector<std::vector<int>> OptimizerMinibatch::getBatches(int numDatasets,
                                                             int numBatches) {
    std::vector<int> datasets(numDatasets);
    for (int i = 0; i < numDatasets; datasets[i] = i++)
        ;
    std::random_shuffle(datasets.begin(), datasets.end());

    std::vector<std::vector<int>> batches;

    int batchSize = ceil((double)numDatasets / numBatches);

    std::vector<int>::iterator batchBegin = datasets.begin();
    std::vector<int>::iterator batchEnd = batchBegin + batchSize;

    for (int i = 0; i < numBatches; ++i) {
        if (batchEnd > datasets.end())
            batchEnd = datasets.end();

        batches[i] = std::vector<int>(datasets.begin() + batchBegin, data.);

        batchBegin = batchEnd;
        batchEnd = batchBegin + batchSize;
    }

    return batches;
}
