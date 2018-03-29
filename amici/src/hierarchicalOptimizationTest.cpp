#include "hierarchicalOptimization.h"
#include "testingMisc.h"
#include "../../optimization/tests/quadraticTestProblem.h"
#include <amici/defines.h>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#define TESTFILE "../../../amici/tests/testhierarchical.h5"

// clang-format off
TEST_GROUP(hierarchicalOptimization){
    void setup() {
        mock().enable();
        mock().clear();
    }

    void teardown() {
        mock().checkExpectations();
        mock().clear();
    }
};
// clang-format on


TEST(hierarchicalOptimization, reader) {
    //     mappingToObservable = np.array([[ 0, 1, 0], [ 0, 2, 0], [ 1, 1, 1], [1, 2, 1], [1, 3, 1]])

    parpe::AnalyticalParameterHdf5Reader r(H5::H5File(TESTFILE, H5F_ACC_RDONLY),
                                           "/scalingParameterIndices",
                                           "/scalingParametersMapToObservables");

    auto exp1 = std::vector<int> {1, 2};
    CHECK_TRUE(exp1 == r.getConditionsForParameter(0));

    auto exp2 = std::vector<int> {1, 2, 3};
    CHECK_TRUE(exp2 == r.getConditionsForParameter(1));

    CHECK_THROWS(std::out_of_range, r.getObservablesForParameter(0, 0));

    auto exp3 = std::vector<int> {1};
    CHECK_TRUE(exp3 == r.getObservablesForParameter(1, 1));

    auto exp4 = std::vector<int> {0, 1};
    CHECK_TRUE(exp4 == r.getOptimizationParameterIndices());
}


class AnalyticalParameterProviderMock : public parpe::AnalyticalParameterProvider {
public:
    AnalyticalParameterProviderMock() = default;

    std::vector<int> getConditionsForParameter(int parameterIndex) const override {
        mock().actualCall("AnalyticalParameterProviderMock::getConditionsForParameter")
                .withIntParameter("parameterIndex", parameterIndex);

        return conditionsForParameter[parameterIndex];
    }

    std::vector<int> const& getObservablesForParameter(int parameterIndex, int conditionIdx) const override {
        mock().actualCall("AnalyticalParameterProviderMock::getObservablesForParameter")
                .withIntParameter("parameterIndex", parameterIndex)
                .withIntParameter("conditionIdx", conditionIdx);

        return mapping[parameterIndex].at(conditionIdx);
    }

    std::vector<int> getOptimizationParameterIndices() const override {
        mock().actualCall("AnalyticalParameterProviderMock::getOptimizationParameterIndices");

        return optimizationParameterIndices;
    }

    std::vector <std::vector<int>> conditionsForParameter;
    std::vector <int> optimizationParameterIndices;
    // x[scalingIdx][conditionIdx] -> std::vector of observableIndicies
    std::vector<std::map<int, std::vector<int>>> mapping;
};


class AmiciSummedGradientFunctionMock : public parpe::AmiciSummedGradientFunction<int> {
public:
    parpe::FunctionEvaluationStatus getModelOutputs(const double * const parameters,
                                                    std::vector<std::vector<double> > &modelOutput) const override {
        mock().actualCall("AmiciSummedGradientFunctionMock::getModelOutputs");

        modelOutput = this->modelOutput;
        lastParameters.assign(parameters, parameters + numParameters_);

        return parpe::functionEvaluationSuccess;
    }

    std::vector<std::vector<double>> getAllMeasurements() const override { return measurements; }

    parpe::FunctionEvaluationStatus evaluate(
            const double* const parameters,
            std::vector<int> datasets,
            double &fval,
            double* gradient) const override {

        mock().actualCall("AmiciSummedGradientFunctionMock::evaluate");

        return parpe::functionEvaluationSuccess;
    }

    amici::AMICI_parameter_scaling getParameterScaling(int parameterIndex) const override {
        return amici::AMICI_SCALING_LOG10;
    }

    parpe::FunctionEvaluationStatus evaluate(
            const double* const parameters,
            int dataset,
            double &fval,
            double* gradient) const override {

        mock().actualCall("AmiciSummedGradientFunctionMock::evaluate");

        return parpe::functionEvaluationSuccess;
    }

    int numParameters() const override {
        mock().actualCall("AmiciSummedGradientFunctionMock::numParameters");

        return numParameters_;
    }

    int numParameters_ = 4;
    int numConditions = 4;
    int numObservables = 3;
    int numTimepoints = 2;

    std::vector<std::vector<double> > modelOutput = {{1.0, 1.0, 1.0,
                                                      1.0, 1.0, 1.0},
                                                     {1.0, 1.0, 1.0,
                                                      1.0, 1.0, 1.0},
                                                     {1.0, 1.0, 1.0,
                                                      1.0, 1.0, 1.0},
                                                     {1.0, 1.0, 1.0,
                                                      1.0, 1.0, 1.0},};
    std::vector<std::vector<double>> measurements =  {{NAN, 1.0, 1.0,
                                                       1.0, 1.0, 1.0},
                                                      {2.0, 1.0, 1.0,
                                                       2.0, 1.0, NAN},
                                                      {2.0, 1.0, 1.0,
                                                       2.0, NAN, 1.0},
                                                      {1.0, 1.0, 1.0,
                                                       NAN, 1.0, 1.0},};
    mutable std::vector<double> lastParameters;
};


TEST(hierarchicalOptimization, hierarchicalOptimization) {
    mock().ignoreOtherCalls();

    auto fun = std::make_unique<AmiciSummedGradientFunctionMock>();
    auto fun2 = fun.get();
    auto r = std::make_unique<parpe::AnalyticalParameterHdf5Reader>(H5::H5File(TESTFILE, H5F_ACC_RDONLY),
                                                                    "/scalingParameterIndices",
                                                                    "/scalingParametersMapToObservables");
    auto scalingReader = r.get();
    auto o = std::make_unique<parpe::AnalyticalParameterHdf5Reader>(H5::H5File(TESTFILE, H5F_ACC_RDONLY),
                                                                    "/offsetParameterIndices",
                                                                    "/offsetParametersMapToObservables");
    parpe::HierachicalOptimizationWrapper w(std::move(fun), std::move(r), std::move(o),
                                            fun2->numConditions, fun2->numObservables, fun2->numTimepoints,
                                            parpe::ErrorModel::normal);

    CHECK_TRUE(w.numScalingFactors() == 2);

    std::vector<double> reducedParameters {3.0, 2.0};
    std::vector<double> fullParameters {3.0, 2.0, 1.5, 1.3}; // last 2 are scalings
    // scalings set to log10(1)
    std::vector<double> onesFullParameters {0.0, 0.0, 3.0, 2.0}; // last 2 are scalings

    std::vector<double> scalingDummy(w.numScalingFactors(), 0.0);
    std::vector<double> offsetDummy(w.numOffsetParameters(),0.0);
    CHECK_TRUE(onesFullParameters == parpe::spliceParameters(reducedParameters.data(), reducedParameters.size(),
                                                             w.getProportionalityFactorIndices(), w.getOffsetParameterIndices(),
                                                             scalingDummy, offsetDummy));

    // Ensure it is called with proper parameter vector:
    auto outputs = w.getUnscaledModelOutputs(fullParameters.data());
    CHECK_TRUE(onesFullParameters == fun2->lastParameters);

    auto s = w.computeAnalyticalScalings(0, outputs, fun2->measurements);
    CHECK_EQUAL(log10(2.0), s);

    applyOptimalScaling(0, 2.0, outputs, *scalingReader, fun2->numObservables, fun2->numTimepoints);
    // output has to be equal to measurement for all points scaled with this parameter
    CHECK_TRUE(outputs[1][0] == fun2->measurements[1][0]);
    CHECK_TRUE(outputs[1][3] == fun2->measurements[1][3]);
    CHECK_TRUE(outputs[2][0] == fun2->measurements[2][0]);
    CHECK_TRUE(outputs[2][3] == fun2->measurements[2][3]);

    // likelihood without offset must be 0 after scaling and if all other measurements/observables agree
    auto llh = w.computeNegLogLikelihood(fun2->measurements, outputs);
    double pi = atan(1)*4;
    double llhOffset = 0.5 * log(2 * pi) * 20;
    DOUBLES_EQUAL(llh - llhOffset, 0, 1e-10);

    //    w.computeAnalyticalScalings();
    //    w.evaluate();
    //    w.evaluateWithScalings();

    CHECK_EQUAL(2, w.numParameters());
}

TEST(hierarchicalOptimization, testNoAnalyticalParameters) {
    // Should only call fun::evaluate, nothing else

    // setup
    auto fun = std::make_unique<AmiciSummedGradientFunctionMock>();
    auto fun2 = fun.get();

    auto scalingProvider = std::make_unique<AnalyticalParameterProviderMock>();
    auto offsetProvider = std::make_unique<AnalyticalParameterProviderMock>();

    // for offsets and proportionality factors
    mock().expectNCalls(2, "AnalyticalParameterProviderMock::getOptimizationParameterIndices");

    parpe::HierachicalOptimizationWrapper w(std::move(fun), std::move(scalingProvider), std::move(offsetProvider),
                                            fun2->numConditions, fun2->numObservables, fun2->numTimepoints,
                                            parpe::ErrorModel::normal);


    mock().expectNCalls(1, "AmiciSummedGradientFunctionMock::evaluate");

    std::vector<double> parameters;
    double fval;
    w.evaluate(parameters.data(), fval, nullptr);
}


TEST(hierarchicalOptimization, testComputeAnalyticalScalings) {
    /* data
     * measurement  = data * 10
     * check scaling = 10
     * */
    int numObservables = 2;
    int numTimepoints = 2;
    int scalingIdx = 0;

    std::vector<std::vector<double> > modelOutputsUnscaled { {1.0, 2.0, 3.0, 4.0} };
    std::vector<std::vector<double> > measurements { {10.0, 20.0, 30.0, 40.0} };

    AnalyticalParameterProviderMock scalingProvider;
    scalingProvider.conditionsForParameter.push_back({0});
    scalingProvider.mapping.resize(1);
    scalingProvider.mapping[0][0] = {0};
    // unused: scalingProvider.optimizationParameterIndices;

    // TEST LIN
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getConditionsForParameter")
            .withIntParameter("parameterIndex", 0);
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getObservablesForParameter")
            .withIntParameter("parameterIndex", 0).withIntParameter("conditionIdx", 0);

    auto scaling = parpe::computeAnalyticalScalings(scalingIdx, amici::AMICI_SCALING_NONE,
                                                    modelOutputsUnscaled, measurements,
                                                    scalingProvider, numObservables, numTimepoints);
    CHECK_EQUAL(10.0, scaling);

    // TEST LOG10
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getConditionsForParameter")
            .withIntParameter("parameterIndex", 0);
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getObservablesForParameter")
            .withIntParameter("parameterIndex", 0).withIntParameter("conditionIdx", 0);

    scaling = parpe::computeAnalyticalScalings(scalingIdx, amici::AMICI_SCALING_LOG10,
                                               modelOutputsUnscaled, measurements,
                                               scalingProvider, numObservables, numTimepoints);
    CHECK_EQUAL(1.0, scaling);

    // TEST LOG10 NAN
    measurements[0][0] = NAN;
    modelOutputsUnscaled[0][0] = 2345; // not used
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getConditionsForParameter")
            .withIntParameter("parameterIndex", 0);
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getObservablesForParameter")
            .withIntParameter("parameterIndex", 0).withIntParameter("conditionIdx", 0);

    scaling = parpe::computeAnalyticalScalings(scalingIdx, amici::AMICI_SCALING_LOG10,
                                               modelOutputsUnscaled, measurements,
                                               scalingProvider, numObservables, numTimepoints);
    CHECK_EQUAL(1.0, scaling);

}


TEST(hierarchicalOptimization, testComputeAnalyticalOffsets) {
    /* data
     * measurement  = data + 10
     * check offset = 10
     * */
    int numObservables = 2;
    int numTimepoints = 2;
    int scalingIdx = 0;

    std::vector<std::vector<double> > modelOutputsUnscaled { {1.0, 2.0, 3.0, 4.0} };
    std::vector<std::vector<double> > measurements { {11.0, 12.0, 13.0, 14.0} };

    AnalyticalParameterProviderMock scalingProvider;
    scalingProvider.conditionsForParameter.push_back({0});
    scalingProvider.mapping.resize(1);
    scalingProvider.mapping[0][0] = {0};
    // unused: scalingProvider.optimizationParameterIndices;

    // TEST LIN
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getConditionsForParameter")
            .withIntParameter("parameterIndex", 0);
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getObservablesForParameter")
            .withIntParameter("parameterIndex", 0).withIntParameter("conditionIdx", 0);

    auto offset = parpe::computeAnalyticalOffsets(scalingIdx, amici::AMICI_SCALING_NONE,
                                                  modelOutputsUnscaled, measurements,
                                                  scalingProvider, numObservables, numTimepoints);
    CHECK_EQUAL(10.0, offset);

    // TEST LOG10
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getConditionsForParameter")
            .withIntParameter("parameterIndex", 0);
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getObservablesForParameter")
            .withIntParameter("parameterIndex", 0).withIntParameter("conditionIdx", 0);

    offset = parpe::computeAnalyticalOffsets(scalingIdx, amici::AMICI_SCALING_LOG10,
                                             modelOutputsUnscaled, measurements,
                                             scalingProvider, numObservables, numTimepoints);
    CHECK_EQUAL(1.0, offset);
}

TEST(hierarchicalOptimization, applyOptimalScaling) {
    int numObservables = 2;
    int numTimepoints = 2;
    int scalingIdx = 0;
    double scaling = 0.5;
    std::vector<std::vector<double> > modelOutputs { {2.0, 4.0, 6.0, 8.0} };
    std::vector<std::vector<double> > modelOutputsScaledExpected { {1.0, 4.0, 3.0, 8.0} };

    AnalyticalParameterProviderMock scalingProvider;
    scalingProvider.conditionsForParameter.push_back({0});
    scalingProvider.mapping.resize(1);
    // applies to all timepoints for observable 0 (== element 0 and 2)
    scalingProvider.mapping[0][0] = {0};
    // unused: scalingProvider.optimizationParameterIndices;

    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getConditionsForParameter")
            .withIntParameter("parameterIndex", 0);
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getObservablesForParameter")
            .withIntParameter("parameterIndex", 0).withIntParameter("conditionIdx", 0);

    parpe::applyOptimalScaling(scalingIdx, scaling, modelOutputs,
                               scalingProvider, numObservables, numTimepoints);

    CHECK_TRUE(modelOutputsScaledExpected == modelOutputs);
}


TEST(hierarchicalOptimization, applyOptimalOffset) {
    int numObservables = 2;
    int numTimepoints = 2;
    int offsetIdx = 0;
    double offset = 5;
    std::vector<std::vector<double> > modelOutputs { {-4.0, 4.0, -2.0, 8.0} };
    std::vector<std::vector<double> > modelOutputsScaledExpected { {1.0, 4.0, 3.0, 8.0} };

    AnalyticalParameterProviderMock offsetProvider;
    offsetProvider.conditionsForParameter.push_back({0});
    offsetProvider.mapping.resize(1);
    // applies to all timepoints for observable 0 (== element 0 and 2)
    offsetProvider.mapping[0][0] = {0};
    // unused: scalingProvider.optimizationParameterIndices;

    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getConditionsForParameter")
            .withIntParameter("parameterIndex", 0);
    mock().expectNCalls(1, "AnalyticalParameterProviderMock::getObservablesForParameter")
            .withIntParameter("parameterIndex", 0).withIntParameter("conditionIdx", 0);

    parpe::applyOptimalOffset(offsetIdx, offset, modelOutputs,
                              offsetProvider, numObservables, numTimepoints);

    CHECK_TRUE(modelOutputsScaledExpected == modelOutputs);
}


TEST(hierarchicalOptimization, testScaling) {
    CHECK_EQUAL(42.0, parpe::getUnscaledParameter(42.0, amici::AMICI_SCALING_NONE));
    CHECK_EQUAL(42.0, parpe::getScaledParameter(42.0, amici::AMICI_SCALING_NONE));

    CHECK_EQUAL(2.0, parpe::getScaledParameter(100.0, amici::AMICI_SCALING_LOG10));
    CHECK_EQUAL(100.0, parpe::getUnscaledParameter(2.0, amici::AMICI_SCALING_LOG10));

    CHECK_THROWS(parpe::ParPEException, parpe::getUnscaledParameter(42.0, amici::AMICI_SCALING_LN));
    CHECK_THROWS(parpe::ParPEException, parpe::getScaledParameter(42.0, amici::AMICI_SCALING_LN));

}

TEST(hierarchicalOptimization, spliceParameters) {
    std::vector<double> fullParametersExp {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};

    std::vector<double> reducedParameters {1.0, 5.0, 8.0};

    std::vector<int> proportionalityFactorIndices {2, 3, 7};
    std::vector<double> scalings {2.0, 3.0, 7.0};

    std::vector<int> offsetParameterIndices {0, 4, 6};
    std::vector<double> offsets {0.0, 4.0, 6.0};

    auto fullParametersAct = parpe::spliceParameters(reducedParameters.data(), reducedParameters.size(),
                                                     proportionalityFactorIndices, offsetParameterIndices,
                                                     scalings, offsets);

    CHECK_TRUE(fullParametersExp == fullParametersAct);
}

TEST(hierarchicalOptimization, spliceParametersNothingToDo) {
    std::vector<double> fullParametersExp {0.0, 1.0, 2.0};

    std::vector<double> reducedParameters {0.0, 1.0, 2.0};

    std::vector<int> proportionalityFactorIndices;
    std::vector<double> scalings;

    std::vector<int> offsetParameterIndices;
    std::vector<double> offsets;

    auto fullParametersAct = parpe::spliceParameters(reducedParameters.data(), reducedParameters.size(),
                                                     proportionalityFactorIndices, offsetParameterIndices,
                                                     scalings, offsets);

    CHECK_TRUE(fullParametersExp == fullParametersAct);
}


TEST(hierarchicalOptimization, filterParams) {
    // TODO
}

TEST(hierarchicalOptimization, problemWrapper) {
    // TODO test wrapper; need dataprovider?!
    //    mock().ignoreOtherCalls();
    //    parpe::QuadraticTestProblem problem;

    //    parpe::OptimizerIpOpt optimizer;
    //    auto result = optimizer.optimize(&problem);

    //    // check status, cost, parameter
    //    CHECK_EQUAL(0, std::get<0>(result));
    //    DOUBLES_EQUAL(42.0, std::get<1>(result), 1e-12);
    //    DOUBLES_EQUAL(-1.0, std::get<2>(result).at(0), 1e-12);
}
