#include "testingMisc.h"
#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include <misc.h>

#include "wrapfunctions.h"
#undef pi

#include <amici/amici.h>
#include "multiConditionProblem.h"
// clang-format off
TEST_GROUP(steadystateProblemTests){
    const std::vector<double> t { 1.0e8 };
    const std::vector<double> k { 0.1, 0.4, 0.7, 1.0 };
    const std::vector<double> p { 1.0, 0.5, 0.4, 2.0, 0.1 };
    const std::vector<double> xSteadystateExp {0.456644592142075,
                                            0.437977375496898,
                                            0.033333333333333};

    void setup(){

    }

    void teardown(){
    }
};
// clang-format on


TEST(steadystateProblemTests, testSteadystate) {
    // verify steadystate

    auto model = getModel();
    model->setFixedParameters(k);
    model->setParameters(p);
    model->setTimepoints(t);
    auto solver = model->getSolver();

    auto rdata = amici::runAmiciSimulation(*solver, nullptr, *model);

   parpe::checkEqualArray(xSteadystateExp.data(),
                          rdata->x.data(),
                          xSteadystateExp.size(), 1e-5, 1e-5);

   // verify likelihood for matching measurement / simulation
   amici::ExpData edata {*model};
   edata.my = xSteadystateExp;
   edata.sigmay.assign(edata.my.size(), 1.0);
   rdata = amici::runAmiciSimulation(*solver, &edata, *model);
   CHECK_EQUAL(rdata->status, AMICI_SUCCESS);
   DOUBLES_EQUAL(parpe::getLogLikelihoodOffset(edata.my.size()), rdata->llh, 1e-5);
}

TEST(steadystateProblemTests, testSteadystateMultiCond) {
    auto model = getModel();
    model->setTimepoints(t);

    parpe::MultiConditionDataProviderDefault dp(std::move(model));

    dp.p.push_back(p);
    dp.k.push_back(k);
    // need to like that until there is amici copy constructor
    dp.edata = std::vector<amici::ExpData>(1);
    dp.edata[0].my = xSteadystateExp;
    dp.edata[0].sigmay.assign(dp.edata[0].my.size(), 1.0);

    //parpe::AmiciSummedGradientFunction<int>(&dp, nullptr);
    parpe::MultiConditionProblem problem(&dp);
    double cost;
    problem.costFun->evaluate(p.data(), cost, nullptr);
    DOUBLES_EQUAL(-parpe::getLogLikelihoodOffset(dp.edata[0].my.size()), cost, 1e-5);
}


TEST(steadystateProblemTests, testSteadystateHierarchical) {
    //TODO
}
