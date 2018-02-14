#include "hierachicalOptimization.h"
#include "testingMisc.h"
#include "../../optimization/tests/quadraticTestProblem.h"

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#define TESTFILE "../../../amici/tests/testhierarchical.h5"

// clang-format off
TEST_GROUP(hierarchicalOptimization){
    void setup() {
    }

    void teardown() {
    }
};
// clang-format on

TEST(hierarchicalOptimization, reader) {
    //     mappingToObservable = np.array([[ 0, 1, 0], [ 0, 2, 0], [ 1, 1, 1], [1, 2, 1], [1, 3, 1]])

    parpe::ScalingFactorHdf5Reader r(H5::H5File(TESTFILE, H5F_ACC_RDONLY), "/");

    auto exp1 = std::vector<int> {1, 2};
    CHECK_TRUE(exp1 == r.getConditionsForScalingParameter(0));

    auto exp2 = std::vector<int> {1, 2, 3};
    CHECK_TRUE(exp2 == r.getConditionsForScalingParameter(1));

    CHECK_THROWS(std::out_of_range, r.getObservablesForScalingParameter(0, 0));

    auto exp3 = std::vector<int> {1};
    CHECK_TRUE(exp3 == r.getObservablesForScalingParameter(1, 1));

    auto exp4 = std::vector<int> {0, 1};
    CHECK_TRUE(exp4 == r.getProportionalityFactorIndices());
}


class AmiciSummedGradientFunctionDummy : public parpe::AmiciSummedGradientFunction<int> {
public:
    parpe::FunctionEvaluationStatus getModelOutputs(const double * const parameters,
                                                    std::vector<std::vector<double> > &modelOutput) const override {
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
        return parpe::functionEvaluationSuccess;
    }

    int numParameters() const override { return numParameters_; }

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

    auto fun = std::make_unique<AmiciSummedGradientFunctionDummy>();
    auto fun2 = fun.get();
    auto r = std::make_unique<parpe::ScalingFactorHdf5Reader>(H5::H5File(TESTFILE, H5F_ACC_RDONLY), "/");
    parpe::HierachicalOptimizationWrapper w(std::move(fun), std::move(r),
                                            fun->numConditions, fun->numObservables, fun->numTimepoints,
                                            parpe::ParameterTransformation::log10, parpe::ErrorModel::normal);

    CHECK_TRUE(w.numScalingFactors() == 2);

    std::vector<double> reducedParameters {3.0, 2.0};
    std::vector<double> fullParameters {3.0, 2.0, 1.5, 1.3}; // last 2 are scalings
    // scalings set to log10(1)
    std::vector<double> onesFullParameters {0.0, 0.0, 3.0, 2.0}; // last 2 are scalings

    std::vector<double> scalingDummy(w.numScalingFactors(), 0.0);
    CHECK_TRUE(onesFullParameters == w.spliceParameters(reducedParameters.data(), reducedParameters.size(), scalingDummy));

    // Ensure it is called with proper parameter vector:
    auto outputs = w.getUnscaledModelOutputs(fullParameters.data());
    CHECK_TRUE(onesFullParameters == fun2->lastParameters);

    auto s = w.computeAnalyticalScalings(0, outputs, fun2->measurements);
    CHECK_EQUAL(log10(2.0), s);

    w.applyOptimalScaling(0, 2.0, outputs);
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
