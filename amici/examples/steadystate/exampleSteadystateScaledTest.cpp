#include "testingMisc.h"
#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <misc.h>
#include <multiConditionProblem.h>
#include <localOptimizationIpopt.h>
#include <hierarchicalOptimization.h>

#include "wrapfunctions.h"

#include <amici/amici.h>

/* Tests using model model_steadystate_scaled */

// clang-format off
TEST_GROUP(steadystateProblemTests){
    /*
    const std::vector<double> t { 1.0e8 };
    const std::vector<double> k { 0.1, 0.4, 0.7, 1.0 };
    const std::vector<double> p { 1.0, 0.5, 0.4, 2.0, 0.1 };
    const std::vector<double> xSteadystateExp {0.456644592142075,
                                            0.437977375496898,
                                            0.033333333333333};
    */
    const int scalingParameterIdx = 5;
    const int offsetParameterIdx = 6;
    const int scaledObservableIdx = 3;
    const int offsettedObservableIdx = 4;
    const std::vector<double> t { 1.0e8 };
    // const std::vector<double> k { };
    //const std::vector<double> p { 1.0, 0.5, 0.4, 2.0, 0.1, 1.0, 1.0, 0.0 };
    const std::vector<double> x0 { 0.1, 0.4, 0.7 };
    const std::vector<double> xSteadystateExp {0.456644592142075,
                                               0.437977375496898,
                                               0.033333333333333};
    const std::vector<double> yExp {0.456644592142075,
                                    0.437977375496898,
                                    0.033333333333333,
                                    2.0 * 0.456644592142075,
                                    3.0 + 0.437977375496898,
                                    0.456644592142075};
    void setup(){

    }

    void teardown(){
    }
};
// clang-format on


TEST(steadystateProblemTests, testSteadystate) {
    /* Verify steadystate matches saved results and loglikelihood is correct */

    // verify steadystate
    auto model = getModel();
    model->setTimepoints(t);

    auto solver = model->getSolver();
    auto rdata = amici::runAmiciSimulation(*solver, nullptr, *model);

    // verify steadystate concentrations
    parpe::checkEqualArray(xSteadystateExp.data(),
                           rdata->x.data(),
                           xSteadystateExp.size(), 1e-5, 1e-5);

    // verify likelihood for matching measurement / simulation
    amici::ExpData edata {*model};
    edata.setObservedData(yExp);
    edata.setObservedDataStdDev(std::vector<double>(yExp.size(), 1.0));
    rdata = amici::runAmiciSimulation(*solver, &edata, *model);
    CHECK_EQUAL(rdata->status, AMICI_SUCCESS);
    DOUBLES_EQUAL(parpe::getLogLikelihoodOffset(edata.nt() * edata.nytrue), rdata->llh, 1e-5);
}

TEST(steadystateProblemTests, testSteadystateMultiCond) {
    auto model = getModel();
    auto modelNonOwning = model.get();
    auto p = model->getParameters();
    model->setTimepoints(t);
    model->setInitialStates(x0);
    //model->setParameters(p);

    parpe::MultiConditionDataProviderDefault dp(std::move(model), modelNonOwning->getSolver());

    dp.edata.push_back(amici::ExpData(*modelNonOwning));
    dp.edata[0].fixedParameters = modelNonOwning->getFixedParameters();
    dp.edata[0].setObservedData(yExp);
    dp.edata[0].setObservedDataStdDev(std::vector<double>(yExp.size(), 1.0));

    //parpe::AmiciSummedGradientFunction<int>(&dp, nullptr);
    parpe::MultiConditionProblem problem(&dp);
    double cost;
    problem.costFun->evaluate(p, cost, gsl::span<double>());
    DOUBLES_EQUAL(-parpe::getLogLikelihoodOffset(dp.edata[0].getObservedData().size()), cost, 1e-5);
}


TEST(steadystateProblemTests, testSteadystateHierarchical) {
    // introduce scaling parameters
    auto model = getModel();
    //model->setFixedParameters(k);
    model->setInitialStates(x0);
    //model->setParameters(p);
    model->setTimepoints(t);
    auto modelNonOwning = model.get();

    const double scalingExp = 2.0; // scaling parameter
    const double offsetExp = 2.0; // offset parameter
    const std::vector<double> pReduced { 1.0, 0.5, 0.4, 2.0,
                                         0.1, 1.0 };
    auto yScaledExp = yExp;
    yScaledExp[scaledObservableIdx] = scalingExp * yExp[0];
    yScaledExp[offsettedObservableIdx] = offsetExp + yExp[1];
    parpe::MultiConditionDataProviderDefault dp(std::move(model), modelNonOwning->getSolver());
    // x0?
    dp.edata.push_back(amici::ExpData(*modelNonOwning));
    dp.edata[0].fixedParameters = modelNonOwning->getFixedParameters();
    dp.edata[0].setObservedData(yScaledExp);
    dp.edata[0].setObservedDataStdDev(std::vector<double>(yExp.size(), 1.0));

    //parpe::MultiConditionProblem problem(&dp);

    auto scalings = std::make_unique<parpe::AnalyticalParameterProviderDefault>();
    scalings->conditionsForParameter.push_back({0});
    scalings->optimizationParameterIndices.push_back(scalingParameterIdx);
    // x[scalingIdx][conditionIdx] -> std::vector of observableIndicies
    scalings->mapping.resize(1);
    scalings->mapping[0][0] = {scaledObservableIdx};

    auto offsets = std::make_unique<parpe::AnalyticalParameterProviderDefault>();
    offsets->conditionsForParameter.push_back({0});
    offsets->optimizationParameterIndices.push_back(offsetParameterIdx);
    // x[scalingIdx][conditionIdx] -> std::vector of observableIndicies
    offsets->mapping.resize(1);
    offsets->mapping[0][0] = {offsettedObservableIdx};

    auto sigmas = std::make_unique<parpe::AnalyticalParameterProviderDefault>();

    auto gradFun = std::make_unique<parpe::AmiciSummedGradientFunction<int>>(&dp, nullptr, nullptr);
    parpe::HierarchicalOptimizationWrapper hier(std::move(gradFun),
                                               std::move(scalings),
                                               std::move(offsets),
                                               std::move(sigmas),
                                               dp.getNumberOfConditions(),
                                               modelNonOwning->nytrue,
                                               modelNonOwning->nt(),
                                               parpe::ErrorModel::normal);
    double cost;
    hier.evaluate(pReduced, cost, gsl::span<double>(), nullptr, nullptr);
    DOUBLES_EQUAL(-parpe::getLogLikelihoodOffset(dp.edata[0].getObservedData().size()), cost, 1e-5);

    const std::vector<double> pFull { 1.0, 0.5, 0.4, 2.0,
                                      0.1, scalingExp, offsetExp, 1.0 };
    hier.fun->evaluate(pFull, {0}, cost, gsl::span<double>(), nullptr, nullptr);
    DOUBLES_EQUAL(-parpe::getLogLikelihoodOffset(dp.edata[0].getObservedData().size()), cost, 1e-5);
}



TEST(steadystateProblemTests, testOptimizationHierarchical) {
    /* setup model & solver */
    // introduce scaling parameters
    auto model = getModel();
    //model->setFixedParameters(k);
    //model->setInitialStates(x0);
    //model->setParameters(p);
    model->setTimepoints(t);
    model->setParameterList({0, 1, 2, 3, 4, 5, 6, 7});
    auto modelNonOwning = model.get();

    auto solver = model->getSolver();
    solver->setSensitivityMethod(amici::SensitivityMethod::adjoint);

    /* generate scaled data */
    const double scalingExp = 2.0; // scaling parameter
    const double offsetExp = 2.0; // offset parameter
    const std::vector<double> pReduced { 1.0, 0.5, 0.4, 2.0,
                                         0.1, /*1.0,*/ /*1.0,*/ 1.0 };
    auto yScaledExp = yExp;
    yScaledExp[scaledObservableIdx] = scalingExp * yExp[0];
    yScaledExp[offsettedObservableIdx] = offsetExp + yExp[1];
    parpe::MultiConditionDataProviderDefault dp(std::move(model), std::move(solver));
    // x0?
    dp.edata.push_back(amici::ExpData(*modelNonOwning));
    dp.edata[0].fixedParameters = modelNonOwning->getFixedParameters();
    dp.edata[0].setObservedData(yScaledExp);
    dp.edata[0].setObservedDataStdDev(std::vector<double>(yScaledExp.size(), 1.0));
    //parpe::MultiConditionProblem problem(&dp);

    /* setup hierarchical optimization */
    // one scaling parameter
    auto scalings = std::make_unique<parpe::AnalyticalParameterProviderDefault>();
    scalings->conditionsForParameter.push_back({0});
    scalings->optimizationParameterIndices.push_back(scalingParameterIdx);
    // x[scalingIdx][conditionIdx] -> std::vector of observableIndicies
    scalings->mapping.resize(1);
    scalings->mapping[0][0] = {scaledObservableIdx};

    // analytical offset parameter
    auto offsets = std::make_unique<parpe::AnalyticalParameterProviderDefault>();
    offsets->conditionsForParameter.push_back({0});
    offsets->optimizationParameterIndices.push_back(offsetParameterIdx);
    // x[scalingIdx][conditionIdx] -> std::vector of observableIndicies
    offsets->mapping.resize(1);
    offsets->mapping[0][0] = {offsettedObservableIdx};

    auto sigmas = std::make_unique<parpe::AnalyticalParameterProviderDefault>();

    // create wrapper
    auto gradFun = std::make_unique<parpe::AmiciSummedGradientFunction<int>>(&dp, nullptr, nullptr);
    auto hier = std::make_unique<parpe::HierarchicalOptimizationWrapper>(std::move(gradFun),
                                                                        std::move(scalings),
                                                                        std::move(offsets),
                                                                        std::move(sigmas),
                                                                        dp.getNumberOfConditions(),
                                                                        modelNonOwning->nytrue,
                                                                        modelNonOwning->nt(),
                                                                        parpe::ErrorModel::normal);
    // evaluate and ensure scaling factor is computed so that y_mes = y_sim
    double cost;
    hier->evaluate(pReduced, cost, gsl::span<double>(), nullptr, nullptr);
    DOUBLES_EQUAL(-parpe::getLogLikelihoodOffset(dp.edata[0].getObservedData().size()), cost, 1e-5);

    const std::vector<double> pFull { 1.0, 0.5, 0.4, 2.0,
                                      0.1, scalingExp, offsetExp, 1.0 };
    hier->fun->evaluate(pFull, {0}, cost, gsl::span<double>(), nullptr, nullptr);
    DOUBLES_EQUAL(-parpe::getLogLikelihoodOffset(dp.edata[0].getObservedData().size()), cost, 1e-5);

    parpe::OptimizationProblemImpl problem(std::move(hier), std::make_unique<parpe::Logger>());
    //    std::vector<double> startingPoint = pReduced;
    //    for(auto& pp : startingPoint)
    //        pp += 1;
    //    problem.setInitialParameters(startingPoint);
    problem.setParametersMin({0, 0, 0, 0,
                              0, 0, 0});
    problem.setParametersMax({2, 2, 2, 2,
                              2, 2, 2});

    // TODO: need to switch on sensitivities
    parpe::OptimizerIpOpt optimizer;
    //auto result = optimizer.optimize(&problem);
    // check status, cost, parameter
    //CHECK_EQUAL(1, std::get<0>(result));
    //DOUBLES_EQUAL(-parpe::getLogLikelihoodOffset(dp.edata[0].my.size()), std::get<1>(result), 1e-5);
    //DOUBLES_EQUAL(-1.0, std::get<2>(result).at(0), 1e-8);
    //std::cout<<std::get<2>(result);

    // TODO: make identifiable
}
