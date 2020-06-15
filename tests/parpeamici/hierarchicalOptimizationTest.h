#include <parpeamici/hierarchicalOptimization.h>
#include <parpeamici/amiciMisc.h>
#include <parpecommon/parpeException.h>

#include "../parpeoptimization/quadraticTestProblem.h"
#include "../parpecommon/testingMisc.h"

#include <amici/defines.h>
#include <amici/misc.h>

#include <gtest/gtest.h>


// created by hierarchicalOptimizationTest.py form CMake
#define TESTFILE "testhierarchical.h5"

using ::testing::Mock;
using ::testing::_;
using ::testing::Ne;
using ::testing::Eq;
using ::testing::Return;
using ::testing::ReturnRef;
using ::testing::ReturnRefOfCopy;
using ::testing::SetArgReferee;
using ::testing::DoAll;

TEST(hierarchicalOptimization1, reader) {
    // mappingToObservable = np.array([[ 0, 1, 0], [ 0, 2, 0], [ 1, 1, 1], [1, 2, 1], [1, 3, 1]])

    parpe::AnalyticalParameterHdf5Reader r(H5::H5File(TESTFILE, H5F_ACC_RDONLY),
                                           "/scalingParameterIndices",
                                           "/scalingParametersMapToObservables");

    auto exp1 = std::vector<int> {1, 2};
    EXPECT_TRUE(exp1 == r.getConditionsForParameter(0));

    auto exp2 = std::vector<int> {1, 2, 3};
    EXPECT_TRUE(exp2 == r.getConditionsForParameter(1));

    EXPECT_THROW(r.getObservablesForParameter(0, 0), std::out_of_range);

    auto exp3 = std::vector<int> {1};
    EXPECT_TRUE(exp3 == r.getObservablesForParameter(1, 1));

    auto exp4 = std::vector<int> {0, 1};
    EXPECT_TRUE(exp4 == r.getOptimizationParameterIndices());
}


class AnalyticalParameterProviderMock
        : public parpe::AnalyticalParameterProvider {
public:
    AnalyticalParameterProviderMock() = default;

    MOCK_CONST_METHOD1(getConditionsForParameter,
                       std::vector<int>(int parameterIndex));
    MOCK_CONST_METHOD2(getObservablesForParameter,
                       std::vector<int> const&(int parameterIndex, int conditionIdx));
    MOCK_CONST_METHOD0(getOptimizationParameterIndices, std::vector<int>());
};


class AmiciSummedGradientFunctionMock : public parpe::AmiciSummedGradientFunction {
public:
    MOCK_CONST_METHOD4(getModelOutputs, parpe::FunctionEvaluationStatus(
                           gsl::span<double const> parameters,
                           std::vector<std::vector<double> > &modelOutput,
                           parpe::Logger *logger, double *cpuTime));

    MOCK_CONST_METHOD0(getAllMeasurements,  std::vector<std::vector<double>>());

    MOCK_CONST_METHOD0(getAllSigmas,  std::vector<std::vector<double>>());

    //    MOCK_CONST_METHOD6(evaluate, parpe::FunctionEvaluationStatus(
    //                           gsl::span<double const> parameters,
    //                           int dataset,
    //                           double &fval,
    //                           gsl::span<double> gradient,
    //                           parpe::Logger *logger, double *cpuTime));

    MOCK_CONST_METHOD6(evaluate, parpe::FunctionEvaluationStatus(
                           gsl::span<double const> parameters,
                           std::vector<int> datasets,
                           double &fval,
                           gsl::span<double> gradient,
                           parpe::Logger *logger, double *cpuTime));

    MOCK_CONST_METHOD1(getParameterScaling, amici::ParameterScaling(int parameterIndex));

    MOCK_CONST_METHOD0(numParameters,  int());
};


class hierarchicalOptimization : public ::testing::Test {

protected:
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
    std::vector<std::vector<double>> sigmas =  {{NAN, 1.0, 1.0,
                                                 1.0, 1.0, 1.0},
                                                {1.0, 1.0, 1.0,
                                                 1.0, 1.0, NAN},
                                                {1.0, 1.0, 1.0,
                                                 1.0, NAN, 1.0},
                                                {1.0, 1.0, 1.0,
                                                 NAN, 1.0, 1.0},};
};



TEST_F(hierarchicalOptimization, hierarchicalOptimization) {
    auto funUnqiue = std::make_unique<AmiciSummedGradientFunctionMock>();
    auto fun = funUnqiue.get();
    ON_CALL(*fun, numParameters()).WillByDefault(Return(numParameters_));
    ON_CALL(*fun, getParameterScaling(_))
            .WillByDefault(Return(amici::ParameterScaling::log10));

    auto scalingReaderUnique =
            std::make_unique<parpe::AnalyticalParameterHdf5Reader>(
                H5::H5File(TESTFILE, H5F_ACC_RDONLY),
                "/scalingParameterIndices",
                "/scalingParametersMapToObservables");
    auto scalingReader = scalingReaderUnique.get();

    auto offsetReaderUnique =
            std::make_unique<parpe::AnalyticalParameterHdf5Reader>(
                H5::H5File(TESTFILE, H5F_ACC_RDONLY),
                "/offsetParameterIndices",
                "/offsetParametersMapToObservables");

    auto sigmaReaderUnique =
            std::make_unique<parpe::AnalyticalParameterHdf5Reader>(
                H5::H5File(TESTFILE, H5F_ACC_RDONLY),
                "/sigmaParameterIndices",
                "/sigmaParametersMapToObservables");

    parpe::HierarchicalOptimizationWrapper hierarchicalOptimizationWrapper(
                std::move(funUnqiue),
                std::move(scalingReaderUnique), std::move(offsetReaderUnique), std::move(sigmaReaderUnique),
                numConditions, numObservables,
                parpe::ErrorModel::normal);

    EXPECT_EQ(2, hierarchicalOptimizationWrapper.numProportionalityFactors());

    std::vector<double> reducedParameters {3.0, 2.0};
    // last 2 are scalings
    std::vector<double> fullParameters {3.0, 2.0, 1.5, 1.3};
    // scalings set to log10(1)
    // last 2 are scalings
    std::vector<double> onesFullParameters {0.0, 0.0, 3.0, 2.0};

    std::vector<double> scalingDummy(
                hierarchicalOptimizationWrapper.numProportionalityFactors(), 0.0);
    std::vector<double> offsetDummy(
                hierarchicalOptimizationWrapper.numOffsetParameters(), 0.0);
    std::vector<double> sigmaDummy(
                hierarchicalOptimizationWrapper.numSigmaParameters(), 0.0);
    auto splicedParameter = parpe::spliceParameters(
                gsl::make_span(reducedParameters.data(),
                               reducedParameters.size()),
                hierarchicalOptimizationWrapper.getProportionalityFactorIndices(),
                hierarchicalOptimizationWrapper.getOffsetParameterIndices(),
                hierarchicalOptimizationWrapper.getSigmaParameterIndices(),
                scalingDummy, offsetDummy, sigmaDummy);
    EXPECT_EQ(onesFullParameters, splicedParameter);

    ON_CALL(*fun, getModelOutputs(_, _, _, _))
            .WillByDefault(DoAll(SetArgReferee<1>(modelOutput),
                                 Return(parpe::functionEvaluationSuccess)));

    // Ensure it is called with proper parameter vector:
    EXPECT_CALL(*fun, getModelOutputs(
                    gsl::span<const double>(onesFullParameters), _, _, _));

    auto outputs = hierarchicalOptimizationWrapper.getUnscaledModelOutputs(
                reducedParameters, nullptr, nullptr);
    Mock::VerifyAndClearExpectations(fun);

    auto s = parpe::getScaledParameter(
                parpe::computeAnalyticalScalings(
                    0, outputs, measurements,
                    *scalingReader, numObservables),
                amici::ParameterScaling::log10);
    EXPECT_EQ(log10(2.0), s);

    applyOptimalScaling(0, 2.0, outputs, *scalingReader, numObservables);
    // output has to be equal to measurement for all points scaled with this
    // parameter
    EXPECT_TRUE(outputs[1][0] == measurements[1][0]);
    EXPECT_TRUE(outputs[1][3] == measurements[1][3]);
    EXPECT_TRUE(outputs[2][0] == measurements[2][0]);
    EXPECT_TRUE(outputs[2][3] == measurements[2][3]);

    // likelihood without offset must be 0 after scaling and if all other
    // measurements/observables agree
    std::vector<std::vector<double>> sigmas(
                outputs.size(), std::vector<double>(outputs[0].size(), 1.0));
    const auto llh = parpe::computeNegLogLikelihood(measurements, outputs, sigmas);
    const double pi = atan(1) * 4;
    const double llhOffset = 0.5 * log(2 * pi) * 20;
    EXPECT_NEAR(llh - llhOffset, 0, 1e-10);

    //    w.computeAnalyticalScalings();
    //    w.evaluate();
    //    w.evaluateWithScalings();

    EXPECT_EQ(2, hierarchicalOptimizationWrapper.numParameters());
}

TEST_F(hierarchicalOptimization, testNoAnalyticalParameters) {
    // Should only call fun::evaluate, nothing else

    // setup
    auto fun = std::make_unique<AmiciSummedGradientFunctionMock>();
    auto funNonOwning = fun.get();

    ON_CALL(*fun, numParameters()).WillByDefault(Return(numParameters_));

    auto scalingProvider = std::make_unique<AnalyticalParameterProviderMock>();
    auto offsetProvider = std::make_unique<AnalyticalParameterProviderMock>();
    auto sigmaProvider = std::make_unique<AnalyticalParameterProviderMock>();

    EXPECT_CALL(*fun, numParameters()).Times(3);

    // for offsets and proportionality factors and sigmas (empty)
    EXPECT_CALL(*scalingProvider, getOptimizationParameterIndices());
    EXPECT_CALL(*offsetProvider, getOptimizationParameterIndices());
    EXPECT_CALL(*sigmaProvider, getOptimizationParameterIndices());

    parpe::HierarchicalOptimizationWrapper w(
                std::move(fun), std::move(scalingProvider),
                std::move(offsetProvider), std::move(sigmaProvider),
                numConditions, numObservables, parpe::ErrorModel::normal);

    EXPECT_CALL(*funNonOwning, evaluate(_, _, _, Eq(gsl::span<double>()), _, _));

    std::vector<double> parameters{3.0, 2.0, 1.5, 1.3};
    double fval;
    w.evaluate(parameters, fval, gsl::span<double>(), nullptr, nullptr);
}


TEST_F(hierarchicalOptimization, testComputeAnalyticalScalings) {
    /* data
     * measurement  = data * 10
     * check scaling = 10
     * */
    constexpr int numObservables = 2;
    // constexpr int numTimepoints = 2;
    constexpr int scalingIdx = 0;

    std::vector<std::vector<double> >
            modelOutputsUnscaled { {1.0, 2.0, 3.0, 4.0} };
    std::vector<std::vector<double> >
            measurements { {10.0, 20.0, 30.0, 40.0} };

    AnalyticalParameterProviderMock scalingProvider;

    // TEST LIN
    std::vector<int> res {0};
    ON_CALL(scalingProvider, getConditionsForParameter(0))
            .WillByDefault(Return(res));
    ON_CALL(scalingProvider, getObservablesForParameter(0, 0))
            .WillByDefault(ReturnRef(res));

    EXPECT_CALL(scalingProvider, getConditionsForParameter(0));
    EXPECT_CALL(scalingProvider, getObservablesForParameter(0, 0));

    const auto scaling = parpe::computeAnalyticalScalings(
                scalingIdx, modelOutputsUnscaled, measurements,
                scalingProvider, numObservables);
    const auto scaledScaling = parpe::getScaledParameter(
                scaling, amici::ParameterScaling::none);

    EXPECT_EQ(10.0, scaledScaling);

    // TEST LOG10
    EXPECT_CALL(scalingProvider, getConditionsForParameter(0));
    EXPECT_CALL(scalingProvider, getObservablesForParameter(0, 0));

    const auto scaling2 = parpe::computeAnalyticalScalings(
                scalingIdx, modelOutputsUnscaled, measurements,
                scalingProvider, numObservables);
    const auto scaledScaling2 = parpe::getScaledParameter(
                scaling2, amici::ParameterScaling::log10);

    EXPECT_EQ(1.0, scaledScaling2);

    // TEST LOG10 NAN
    measurements[0][0] = NAN;
    modelOutputsUnscaled[0][0] = 2345; // not used

    EXPECT_CALL(scalingProvider, getConditionsForParameter(0));
    EXPECT_CALL(scalingProvider, getObservablesForParameter(0, 0));

    const auto scaling3 = parpe::computeAnalyticalScalings(
                scalingIdx, modelOutputsUnscaled, measurements,
                scalingProvider, numObservables);
    const auto scaledScaling3 = parpe::getScaledParameter(
                scaling3, amici::ParameterScaling::log10);

    EXPECT_EQ(1.0, scaledScaling3);
}


TEST_F(hierarchicalOptimization, testComputeAnalyticalOffsets) {
    /* data
     * measurement  = data + 10
     * check offset = 10
     * */
    constexpr int numObservables = 2;
    // constexpr int numTimepoints = 2;
    constexpr int scalingIdx = 0;

    const std::vector<std::vector<double> >
            modelOutputsUnscaled { {1.0, 2.0, 3.0, 4.0} };
    const std::vector<std::vector<double> >
            measurements { {11.0, 12.0, 13.0, 14.0} };

    AnalyticalParameterProviderMock scalingProvider;

    // TEST LIN
    std::vector<int> res {0};
    ON_CALL(scalingProvider, getConditionsForParameter(0))
            .WillByDefault(Return(res));
    ON_CALL(scalingProvider, getObservablesForParameter(0, 0))
            .WillByDefault(ReturnRefOfCopy(res));

    EXPECT_CALL(scalingProvider, getConditionsForParameter(0));
    EXPECT_CALL(scalingProvider, getObservablesForParameter(0, 0));

    const auto offset = parpe::computeAnalyticalOffsets(
                scalingIdx, modelOutputsUnscaled, measurements,
                scalingProvider, numObservables);
    const auto scaledOffset = parpe::getScaledParameter(
                offset, amici::ParameterScaling::none);
    EXPECT_EQ(10.0, scaledOffset);

    // TEST LOG10
    EXPECT_CALL(scalingProvider, getConditionsForParameter(0));
    EXPECT_CALL(scalingProvider, getObservablesForParameter(0, 0));

    const auto offset2 = parpe::computeAnalyticalOffsets(
                scalingIdx, modelOutputsUnscaled, measurements,
                scalingProvider, numObservables);
    const auto scaledOffset2 = parpe::getScaledParameter(
                offset2, amici::ParameterScaling::log10);
    EXPECT_EQ(1.0, scaledOffset2);
}

TEST_F(hierarchicalOptimization, applyOptimalScaling) {
    constexpr int numObservables = 2;
    // constexpr int numTimepoints = 2;
    constexpr int scalingIdx = 0;
    constexpr double scaling = 0.5;
    const std::vector<std::vector<double> >
            modelOutputsScaledExpected { {1.0, 4.0, 3.0, 8.0} };
    std::vector<std::vector<double> > modelOutputs { {2.0, 4.0, 6.0, 8.0} };

    AnalyticalParameterProviderMock scalingProvider;

    std::vector<int> res {0};
    ON_CALL(scalingProvider, getConditionsForParameter(0))
            .WillByDefault(Return(res));
    // applies to all timepoints for observable 0 (== element 0 and 2)
    ON_CALL(scalingProvider, getObservablesForParameter(0, 0))
            .WillByDefault(ReturnRefOfCopy(res));

    EXPECT_CALL(scalingProvider, getConditionsForParameter(0));
    EXPECT_CALL(scalingProvider, getObservablesForParameter(0, 0));

    parpe::applyOptimalScaling(scalingIdx, scaling, modelOutputs,
                               scalingProvider, numObservables);

    EXPECT_EQ(modelOutputsScaledExpected, modelOutputs);
}


TEST_F(hierarchicalOptimization, applyOptimalOffset) {
    constexpr int numObservables = 2;
    // constexpr int numTimepoints = 2;
    constexpr int offsetIdx = 0;
    constexpr double offset = 5;
    const std::vector<std::vector<double> > modelOutputsScaledExpected { {1.0, 4.0, 3.0, 8.0} };
    std::vector<std::vector<double> > modelOutputs { {-4.0, 4.0, -2.0, 8.0} };


    AnalyticalParameterProviderMock offsetProvider;
    std::vector<int> res {0};
    ON_CALL(offsetProvider, getConditionsForParameter(0))
            .WillByDefault(Return(res));
    // applies to all timepoints for observable 0 (== element 0 and 2)
    ON_CALL(offsetProvider, getObservablesForParameter(0, 0))
            .WillByDefault(ReturnRefOfCopy(res));

    EXPECT_CALL(offsetProvider, getConditionsForParameter(0));
    EXPECT_CALL(offsetProvider, getObservablesForParameter(0, 0));

    parpe::applyOptimalOffset(offsetIdx, offset, modelOutputs,
                              offsetProvider, numObservables);

    EXPECT_EQ(modelOutputsScaledExpected, modelOutputs);
}


TEST_F(hierarchicalOptimization, testScaling) {
    EXPECT_EQ(42.0,
              amici::getUnscaledParameter(42.0, amici::ParameterScaling::none));
    EXPECT_EQ(42.0,
              parpe::getScaledParameter(42.0, amici::ParameterScaling::none));

    EXPECT_EQ(2.0,
              parpe::getScaledParameter(100.0, amici::ParameterScaling::log10));
    EXPECT_EQ(100.0,
              amici::getUnscaledParameter(2.0, amici::ParameterScaling::log10));

    EXPECT_DOUBLE_EQ(1.0,
                     parpe::getScaledParameter(std::exp(1.0),
                                               amici::ParameterScaling::ln));
    EXPECT_DOUBLE_EQ(std::exp(1),
                     amici::getUnscaledParameter(1.0,
                                                 amici::ParameterScaling::ln));
}

TEST(hierarchicalOptimization1, spliceParameters) {
    const std::vector<double>
            fullParametersExp {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};

    const std::vector<double> reducedParameters {1.0, 5.0, 8.0};

    const std::vector<int> proportionalityFactorIndices {2, 3, 7};
    const std::vector<double> scalings {2.0, 3.0, 7.0};

    const std::vector<int> offsetParameterIndices {0, 4};
    const std::vector<double> offsets {0.0, 4.0};

    const std::vector<int> sigmaParameterIndices {6};
    const std::vector<double> sigmas {6.0};

    auto fullParametersAct = parpe::spliceParameters(gsl::make_span(reducedParameters.data(), reducedParameters.size()),
                                                     proportionalityFactorIndices, offsetParameterIndices, sigmaParameterIndices,
                                                     scalings, offsets, sigmas);

    EXPECT_EQ(fullParametersExp, fullParametersAct);
}

TEST(hierarchicalOptimization1, spliceParametersNothingToDo) {
    const std::vector<double> fullParametersExp {0.0, 1.0, 2.0};

    const std::vector<double> reducedParameters {0.0, 1.0, 2.0};

    const std::vector<int> proportionalityFactorIndices;
    const std::vector<double> scalings;

    const std::vector<int> offsetParameterIndices;
    const std::vector<double> offsets;

    const std::vector<int> sigmaParameterIndices;
    const std::vector<double> sigmas;

    auto fullParametersAct = parpe::spliceParameters(
                reducedParameters, proportionalityFactorIndices,
                offsetParameterIndices, sigmaParameterIndices,
                scalings, offsets, sigmas);

    EXPECT_EQ(fullParametersExp, fullParametersAct);
}


TEST(hierarchicalOptimization1, fillFilteredParams) {
    const std::vector<double> resultExp {1.0, 2.0, 3.0, 4.0, 5.0};
    const std::vector<double>
            valuesToFilter {9.0, 1.0, 2.0, 9.0, 9.0, 3.0, 4.0, 5.0, 9.0};
    const std::vector<int> sortedIndicesToExclude {0, 3, 4, 8};

    auto resultAct = std::vector<double>(
                valuesToFilter.size() - sortedIndicesToExclude.size());
    parpe::fillFilteredParams(valuesToFilter, sortedIndicesToExclude,
                              resultAct);

    EXPECT_EQ(resultExp, resultAct);
}


TEST_F(hierarchicalOptimization, testWrappedFunIsCalledWithGradient) {
    // setup
    auto fun = std::make_unique<AmiciSummedGradientFunctionMock>();
    auto funNonOwning = fun.get();

    auto scalingProvider = std::make_unique<AnalyticalParameterProviderMock>();
    auto offsetProvider = std::make_unique<AnalyticalParameterProviderMock>();
    auto sigmaProvider = std::make_unique<AnalyticalParameterProviderMock>();
    auto scalingProviderNonOwning = scalingProvider.get();

    std::vector<int> res {0};
    std::vector<int> idxs {3};
    ON_CALL(*scalingProvider, getOptimizationParameterIndices())
            .WillByDefault(Return(idxs));
    ON_CALL(*scalingProvider, getConditionsForParameter(0))
            .WillByDefault(Return(res));
    // applies to all timepoints for observable 0 (== element 0 and 2)
    ON_CALL(*scalingProvider, getObservablesForParameter(0, 0))
            .WillByDefault(ReturnRefOfCopy(res));
    ON_CALL(*fun, numParameters()).WillByDefault(Return(numParameters_));
    ON_CALL(*fun, getModelOutputs(_, _, _, _))
            .WillByDefault(DoAll(SetArgReferee<1>(modelOutput),
                                 Return(parpe::functionEvaluationSuccess)));
    ON_CALL(*fun, getAllMeasurements()).WillByDefault(Return(measurements));
    ON_CALL(*fun, getAllSigmas()).WillByDefault(Return(sigmas));

    EXPECT_CALL(*fun, numParameters()).Times(2);
    EXPECT_CALL(*scalingProvider, getOptimizationParameterIndices());
    EXPECT_CALL(*offsetProvider, getOptimizationParameterIndices());
    EXPECT_CALL(*sigmaProvider, getOptimizationParameterIndices());

    parpe::HierarchicalOptimizationWrapper hierarchicalWrapper(
                std::move(fun), std::move(scalingProvider),
                std::move(offsetProvider), std::move(sigmaProvider),
                numConditions, numObservables,
                parpe::ErrorModel::normal);
    Mock::VerifyAndClearExpectations(funNonOwning);
    Mock::VerifyAndClearExpectations(scalingProviderNonOwning);
    //    Mock::VerifyAndClearExpectations(*offsetProvider);
    //    Mock::VerifyAndClearExpectations(*sigmaProvider);

    EXPECT_CALL(*funNonOwning, numParameters()).Times(1);

    const std::vector<double> parameters { 1.0, 2.0, 3.0, /*4.0*/ };
    EXPECT_EQ((unsigned) hierarchicalWrapper.numParameters(),
              parameters.size());

    std::vector<double> gradient(parameters.size());

    // ensure fun::evaluate is called with gradient
    EXPECT_CALL(*funNonOwning, numParameters());
    EXPECT_CALL(*funNonOwning, getModelOutputs(_, _, _, _));
    EXPECT_CALL(*scalingProviderNonOwning, getConditionsForParameter(0)).Times(2);
    EXPECT_CALL(*scalingProviderNonOwning, getObservablesForParameter(0, 0)).Times(2);
    EXPECT_CALL(*funNonOwning, evaluate(_, _, _, Ne(gsl::span<double>()), _, _));

    double fval;
    hierarchicalWrapper.evaluate(parameters, fval, gradient, nullptr, nullptr);

    Mock::VerifyAndClearExpectations(funNonOwning);
    Mock::VerifyAndClearExpectations(scalingProviderNonOwning);

    // test fun::evaluate is not called if no gradient (only get outputs)
    EXPECT_CALL(*funNonOwning, numParameters());
    EXPECT_CALL(*funNonOwning, getModelOutputs(_, _, _, _));
    EXPECT_CALL(*scalingProviderNonOwning, getConditionsForParameter(0))
            .Times(2);
    EXPECT_CALL(*scalingProviderNonOwning, getObservablesForParameter(0, 0))
            .Times(2);

    hierarchicalWrapper.evaluate(parameters, fval, gsl::span<double>(), nullptr, nullptr);
}

TEST(hierarchicalOptimization1, likelihoodOfMatchingData) {
    const std::vector<double> data {1.0, 2.0, 3.0};
    const std::vector<double> sigmas {1.0, 1.0, 1.0};

    const double pi = atan(1)*4;
    const double llhOffset = 0.5 * log(2 * pi);
    const double expected = llhOffset * static_cast<double>(data.size());

    auto actual = parpe::computeNegLogLikelihood(data, data, sigmas);
    EXPECT_EQ(expected, actual);
}


TEST_F(hierarchicalOptimization, problemWrapper) {
    //std::unique_ptr<parpe::OptimizationProblem> problem(new parpe::QuadraticTestProblem());
    //auto hCost = std::make_unique<parpe::HierarchicalOptimizationWrapper>();
    //auto wrappedFun = dynamic_cast<SummedGradientFunctionGradientFunctionAdapter<int>*>(wrappedProblem->costFun.get());

    //parpe::HierarchicalOptimizationProblemWrapper hw(std::move(problem), std::move(hCost));

    // TODO test wrapper; need dataprovider?!
    //    mock().ignoreOtherCalls();
    //    parpe::QuadraticTestProblem problem;

    //    parpe::OptimizerIpOpt optimizer;
    //    auto result = optimizer.optimize(&problem);

    //    // check status, cost, parameter
    //    EXPECT_EQ(0, std::get<0>(result));
    //    DOUBLES_EQUAL(42.0, std::get<1>(result), 1e-12);
    //    DOUBLES_EQUAL(-1.0, std::get<2>(result).at(0), 1e-12);
}
