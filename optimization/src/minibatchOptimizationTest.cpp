#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>

#include <minibatchOptimization.h>
#include <model.h>
#include <optimizationOptions.h>
#include <quadraticTestProblem.h>
#include <testingMisc.h>
#include <misc.h>

#include <functional>
#include <random>
#include <algorithm>

// clang-format off
TEST_GROUP(minibatchOptimization){

    void setup() {
        mock().clear();
    }

    void teardown() {
        mock().checkExpectations();
        mock().clear();
    }
};
// clang-format on

TEST(minibatchOptimization, getBatches) {

    int numElements = 10;
    std::vector<int> input(numElements);
    std::iota(input.begin(), input.end(), 0);

    // single batch
    int batchSize = numElements;
    auto batchesAct = parpe::getBatches<int>(input, batchSize);
    CHECK_EQUAL(1, batchesAct.size());
    CHECK_TRUE(input == batchesAct[0]);

    // 2 batches, equal size
    batchSize = 5;
    batchesAct = parpe::getBatches<int>(input, batchSize);
    CHECK_EQUAL(2, batchesAct.size());
    CHECK_TRUE(std::vector<int>(input.begin(), input.begin() + batchSize) == batchesAct[0]);
    CHECK_TRUE(std::vector<int>(input.begin() + batchSize, input.end()) == batchesAct[1]);

    // 2 batches, unequal
    batchSize = 6;
    batchesAct = parpe::getBatches<int>(input, batchSize);
    CHECK_EQUAL(2, batchesAct.size());
    CHECK_TRUE(std::vector<int>(input.begin(), input.begin() + batchSize) == batchesAct[0]);
    CHECK_TRUE(std::vector<int>(input.begin() + batchSize, input.end()) == batchesAct[1]);
}


TEST(minibatchOptimization, updateParameters) {
    std::vector<double> gradient {1.0, 2.0};
    std::vector<double> parameters {2.0, 3.0};
    std::vector<double> parametersExp {1.5, 2.0};

    parpe::ParameterUpdaterVanilla pu(0.5);
    pu.updateParameters(gradient, parameters);
    CHECK_TRUE(parametersExp == parameters);
}

TEST(minibatchOptimization, updateParametersRMSProp) {
    std::vector<double> gradient {1.0, 2.0};
    std::vector<double> parameters {2.0, 3.0};
    //std::vector<double> parametersExp {1.5, 2.0};

    parpe::ParameterUpdaterRmsProp pu;
    pu.updateParameters(gradient, parameters);
    pu.updateParameters(gradient, parameters);
    //CHECK_TRUE(parametersExp == parameters);
}




// clang-format off
TEST_GROUP(minibatchOptimizationLinearModel){
    void setup() {
        mock().clear();

        generateRandomFeatures();

        dataIndices.resize(data.size());
        std::iota(dataIndices.begin(), dataIndices.end(), 0.0);

        // make predictions
        labels.resize(numDatasets);
        lm.evaluate(trueParameters, data, labels);
    }

    void generateRandomFeatures() {
        // generate data or feature vector
        std::uniform_real_distribution<double> unif(0,10);
        std::mt19937 rng(rd());
        data.assign(numDatasets,
                    std::vector<double>(trueParameters.size() - 1));
        for(int i = 0; i < numDatasets; ++i) {
            std::generate(data[i].begin(), data[i].end(), [&unif, &rng](){ return unif(rng); });
            //std::cout<<data[i];
        }

    }

    std::unique_ptr<parpe::LinearModelMSE> getLinearModelMSE() {
        // prepare model for optimization
        auto lm2 = std::make_unique<parpe::LinearModelMSE>(trueParameters.size());
        lm2->datasets = data;
        lm2->labels = labels;
        return lm2;
    }

    std::unique_ptr<parpe::OptimizationProblemImpl> getOptimizationProblem() {
        auto lm2 = getLinearModelMSE();

        auto sgf = std::make_unique<
                parpe::SummedGradientFunctionGradientFunctionAdapter<int>
                >(std::move(lm2), dataIndices);
        auto p = std::make_unique<parpe::OptimizationProblemImpl>();
        p->costFun = std::move(sgf);
        p->setParametersMin(std::vector<double>(trueParameters.size(), 0.0));
        p->setParametersMax(std::vector<double>(trueParameters.size(), 5.0));
        p->logger = std::make_unique<parpe::Logger>();
        return p;
    }

    void teardown() {
        mock().checkExpectations();
        mock().clear();
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
// clang-format on




TEST(minibatchOptimizationLinearModel, testCostWithTrueParametersIsZeroIndivdually) {
    // verify cost gradient with true parameters is 0
    auto lm2 = getLinearModelMSE();
    double mse = NAN;
    std::vector<double> gradient(trueParameters.size());
    for(int i = 0; i < numDatasets; ++i) {
        lm2->evaluate(trueParameters, i, mse, gradient, nullptr, nullptr);
        CHECK_EQUAL(0.0, mse);
        CHECK_TRUE(std::vector<double>(trueParameters.size(), 0.0) == gradient);
    }
}

TEST(minibatchOptimizationLinearModel, testCostWithTrueParametersIsZeroFull) {
    // verify cost gradient with true parameters is 0
    auto lm2 = getLinearModelMSE();
    double mse = NAN;
    std::vector<double> gradient(trueParameters.size());
    lm2->evaluate(trueParameters, dataIndices, mse, gradient, nullptr, nullptr);
    CHECK_EQUAL(0.0, mse);
    CHECK_TRUE(std::vector<double>(trueParameters.size(), 0.0) == gradient);
}

TEST(minibatchOptimizationLinearModel, testMinibatchSucceedFromOptimum) {
    // verify optimization succeeds with true parameters
    auto lm2 = getLinearModelMSE();
    parpe::MinibatchOptimizer<int> mb;
    mb.maxEpochs = 20;
    mb.batchSize = 2;
    std::vector<double> startingPoint = { 3.0, 2.0 };
    auto result = mb.optimize(*lm2, dataIndices, startingPoint);
    CHECK_EQUAL((int)parpe::minibatchExitStatus::gradientNormConvergence, std::get<0>(result));
    CHECK_EQUAL(0.0, std::get<1>(result));
    CHECK_TRUE(trueParameters == std::get<2>(result));
}

TEST(minibatchOptimizationLinearModel, linearModelCheckCostGradient) {
    // use gradient checker
    auto p = getOptimizationProblem();

    for(int i = 0; i < 10; ++i)
        parpe::optimizationProblemGradientCheck(p.get(), 10, 1e-1);

    // TODO: check results automatically
}

//#include <localOptimizationFsqp.h>
#include <localOptimizationIpopt.h>
TEST(minibatchOptimizationLinearModel, linearModelTestBatchOptimizerSucceeds) {
    // test batch optimizer
    auto p = getOptimizationProblem();

    //parpe::OptimizerFsqp o;
    parpe::OptimizerIpOpt o;
    auto oo = p->getOptimizationOptions();
    oo.maxOptimizerIterations = 20;
    p->setOptimizationOptions(oo);

    auto resultBatchOpt = o.optimize(p.get());
    DOUBLES_EQUAL(0.0, std::get<1>(resultBatchOpt), 1e-8);
    for(int i = 0; (unsigned) i < trueParameters.size(); ++i)
        DOUBLES_EQUAL(trueParameters[i], std::get<2>(resultBatchOpt)[i], 1e-6);
    // -> is identifiable and gradient okay
}


TEST(minibatchOptimizationLinearModel, linearModel) {
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
    auto result = mb.optimize(*lm3, dataIndices, startingPoint);


    // TODO add some gaussian noise
    // std::normal_distribution<double> norm(0.0, 1.0);
    //for(auto &e : labels)
    //e += norm(rng);
    // TODO: add test with mockReporter (also for other optimizers)
}
