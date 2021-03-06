#include <gtest/gtest.h>

#include <parpecommon/parpeConfig.h>
#include <parpeoptimization/minibatchOptimization.h>
#include <parpecommon/model.h>
#include <parpeoptimization/optimizationOptions.h>

#include "quadraticTestProblem.h"
#include "../parpecommon/testingMisc.h"
#include <parpecommon/misc.h>

#include <functional>
#include <random>
#include <algorithm>


TEST(MinibatchOptimization, CreatesBatches) {

    int numElements = 10;
    std::vector<int> input(numElements);
    std::iota(input.begin(), input.end(), 0);

    // single batch
    int batchSize = numElements;
    auto batchesAct = parpe::getBatches<int>(input, batchSize);
    EXPECT_EQ(1UL, batchesAct.size());
    EXPECT_TRUE(input == batchesAct[0]);

    // 2 batches, equal size
    batchSize = 5;
    batchesAct = parpe::getBatches<int>(input, batchSize);
    EXPECT_EQ(2UL, batchesAct.size());
    EXPECT_TRUE(std::vector<int>(input.begin(), input.begin() + batchSize) == batchesAct[0]);
    EXPECT_TRUE(std::vector<int>(input.begin() + batchSize, input.end()) == batchesAct[1]);

    // 2 batches, unequal
    batchSize = 6;
    batchesAct = parpe::getBatches<int>(input, batchSize);
    EXPECT_EQ(2UL, batchesAct.size());
    EXPECT_TRUE(std::vector<int>(input.begin(), input.begin() + batchSize) == batchesAct[0]);
    EXPECT_TRUE(std::vector<int>(input.begin() + batchSize, input.end()) == batchesAct[1]);
}

TEST(MinibatchOptimization, UpdatesParameters) {
    // Test whether the most simple parameter updater works reliably
    std::vector<double> gradient {3.0, 4.0};
    std::vector<double> parameters {2.0, 3.0};
    std::vector<double> lowerBounds {-5.0, -5.0};
    std::vector<double> upperBounds {5.0, 5.0};
    std::vector<double> parametersExp {1.7, 2.6};

    double learningRate = 0.5;
    int iteration = 1;
    double toleratedError = 0.0001;
    int numParameters = parameters.size();
    bool errored = false;

    parpe::ParameterUpdaterVanilla pu;
    pu.initialize(2);
    pu.updateParameters(learningRate, iteration, gradient, parameters, lowerBounds, upperBounds);

    std::transform(parametersExp.begin( ), parametersExp.end( ), parameters.begin( ), parametersExp.begin( ), std::minus<double>( ));

    for (int i = 0; i < numParameters; i++)
        if (!errored)
            if (parametersExp[i] > toleratedError or parametersExp[i] < -toleratedError)
                errored = true;

    EXPECT_TRUE(!errored);
}



class MinibatchOptimizationLinearModel : public ::testing::Test {

protected:
    void SetUp() override {
        generateRandomFeatures();

        dataIndices.resize(data.size());
        std::iota(dataIndices.begin(), dataIndices.end(), 0.0);

        // make predictions
        labels.resize(numDatasets);
        lm.evaluate(trueParameters, data, labels);
    }


    void generateRandomFeatures() {
        // generate data or feature vector
        std::uniform_real_distribution<double> unif(0, 10);
        std::mt19937 rng(rd());
        data.assign(numDatasets, std::vector<double>(trueParameters.size() - 1));
        for (int i = 0; i < numDatasets; ++i) {
            std::generate(data[i].begin(), data[i].end(), [&unif, &rng]() {return unif(rng);});
            //std::cout<<data[i];
        }

    }

    std::unique_ptr<parpe::LinearModelMSE> getLinearModelMSE() {
        // prepare model for optimization
        auto lm2 = std::make_unique < parpe::LinearModelMSE > (trueParameters.size());
        lm2->datasets = data;
        lm2->labels = labels;
        return lm2;
    }

    std::unique_ptr<parpe::OptimizationProblemImpl> getOptimizationProblem() {
        auto lm2 = getLinearModelMSE();

        auto sgf = std::make_unique < parpe::SummedGradientFunctionGradientFunctionAdapter<int>
                > (std::move(lm2), dataIndices);
        auto p = std::make_unique<parpe::OptimizationProblemImpl>();
        p->cost_fun_ = std::move(sgf);
        p->setParametersMin(std::vector<double>(trueParameters.size(), 0.0));
        p->setParametersMax(std::vector<double>(trueParameters.size(), 5.0));
        p->logger_ = std::make_unique<parpe::Logger>();
        return p;
    }


    void TearDown() override {
    }

    std::vector<double> trueParameters = { 3.0, 2.0 };
    int numDatasets = 10;
    int batchSize = 2;

    std::random_device rd;
    std::vector<std::vector<double>> data;

    parpe::LinearModel lm;
    std::vector<double> labels;

    std::vector<int> dataIndices;
};


TEST_F(MinibatchOptimizationLinearModel, CostWithTrueParametersIsZeroIndivdually) {
    // verify cost gradient with true parameters is 0
    auto lm2 = getLinearModelMSE();
    double mse = NAN;
    std::vector<double> gradient(trueParameters.size());
    for(int i = 0; i < numDatasets; ++i) {
        lm2->evaluate(trueParameters, i, mse, gradient, nullptr, nullptr);
        EXPECT_EQ(0.0, mse);
        EXPECT_TRUE(std::vector<double>(trueParameters.size(), 0.0) == gradient);
    }
}

TEST_F(MinibatchOptimizationLinearModel, CostWithTrueParametersIsZeroFull) {
    // verify cost gradient with true parameters is 0
    auto lm2 = getLinearModelMSE();
    double mse = NAN;
    std::vector<double> gradient(trueParameters.size());
    lm2->evaluate(trueParameters, dataIndices, mse, gradient, nullptr, nullptr);
    EXPECT_EQ(0.0, mse);
    EXPECT_TRUE(std::vector<double>(trueParameters.size(), 0.0) == gradient);
}

TEST_F(MinibatchOptimizationLinearModel, MinibatchSucceedFromOptimum) {
    // verify optimization succeeds with true parameters
    auto lm2 = getLinearModelMSE();
    parpe::MinibatchOptimizer<int> mb;
    mb.maxEpochs = 20;
    mb.batchSize = 2;
    std::vector<double> startingPoint = {3.0, 2.0};
    auto result = mb.optimize(*lm2, dataIndices, startingPoint, gsl::span<const double>(), gsl::span<const double>(), nullptr, nullptr);
    EXPECT_EQ((int)parpe::minibatchExitStatus::gradientNormConvergence, std::get<0>(result));
    EXPECT_EQ(0.0, std::get<1>(result));
    EXPECT_TRUE(trueParameters == std::get<2>(result));
}

TEST_F(MinibatchOptimizationLinearModel, LinearModelCheckCostGradient) {
    // use gradient checker
    auto p = getOptimizationProblem();

    for(int i = 0; i < 10; ++i)
        parpe::optimizationProblemGradientCheck(p.get(), 10, 1e-1);

    // TODO: check results automatically
}

#ifdef PARPE_ENABLE_IPOPT
#include <parpeoptimization/localOptimizationIpopt.h>
TEST_F(MinibatchOptimizationLinearModel, linearModelDoesBatchOptimizerSucceed) {
    // test batch optimizer
    auto p = getOptimizationProblem();

    parpe::OptimizerIpOpt o;
    auto oo = p->getOptimizationOptions();
    oo.maxOptimizerIterations = 20;
    p->setOptimizationOptions(oo);

    auto resultBatchOpt = o.optimize(p.get());
    EXPECT_NEAR(0.0, std::get<1>(resultBatchOpt), 1e-8);
    for(int i = 0; (unsigned) i < trueParameters.size(); ++i)
    EXPECT_NEAR(trueParameters[i], std::get<2>(resultBatchOpt)[i], 1e-6);
    // -> is identifiable and gradient okay
}
#endif

TEST_F(MinibatchOptimizationLinearModel, LinearModel) {
    // optimization/tests/unittests_optimization -sg minibatchOptimizationLinearModel -sn linearModel
    std::cout<<"True parameters "<<trueParameters<<std::endl;

    // try from non-optimal starting point
    std::vector<double> startingPoint = {2.0, 4};
    auto p = getOptimizationProblem();
    p->setInitialParameters(startingPoint);

    auto lm3 = std::make_unique<parpe::LinearModelMSE>(trueParameters.size());
    lm3->datasets = data;
    lm3->labels = labels;

    parpe::MinibatchOptimizer<int> mb;
    mb.maxEpochs = 100;
    //mb.parameterUpdater = std::make_unique<parpe::ParameterUpdaterVanilla>(0.02);
    mb.parameterUpdater = std::make_unique<parpe::ParameterUpdaterRmsProp>();
    mb.batchSize = batchSize;
    auto result = mb.optimize(*lm3, dataIndices, startingPoint,
            gsl::span<const double>(), gsl::span<const double>(),
            nullptr, nullptr);

    // TODO add some gaussian noise
    // std::normal_distribution<double> norm(0.0, 1.0);
    //for(auto &e : labels)
    //e += norm(rng);
    // TODO: add test with mockReporter (also for other optimizers)
}
