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
    void setup(){

    }

    void teardown(){
    }
};
// clang-format on


TEST(steadystateProblemTests, testSteadystate) {
    // verify steadystate
    const std::vector<double> t { 1.0e8 };
    const std::vector<double> k { 0.1, 0.4, 0.7, 1.0 };
    const std::vector<double> p { 1.0, 0.5, 0.4, 2.0, 0.1 };

    auto model = getModel();
    model->setFixedParameters(k);
    model->setParameters(p);
    model->setTimepoints(t);
    auto solver = model->getSolver();

    auto rdata = amici::runAmiciSimulation(*solver, nullptr, *model);

    const std::vector<double> xSteadystateExp {0.456644592142075,
                                            0.437977375496898,
                                            0.033333333333333};
   parpe::checkEqualArray(xSteadystateExp.data(),
                          rdata->x.data(),
                          xSteadystateExp.size(), 1e-5, 1e-5);

   // verify likelihood for matching measurement / simulation
   amici::ExpData edata {*model};
   edata.my = xSteadystateExp;
   edata.sigmay.assign(edata.my.size(), 1.0);
   rdata = amici::runAmiciSimulation(*solver, &edata, *model);
   CHECK_EQUAL(rdata->status, AMICI_SUCCESS);
   DOUBLES_EQUAL(- parpe::getLikelihoodOffset(edata.my.size()), rdata->llh, 1e-5);
}

TEST(steadystateProblemTests, testSteadystateMultiCond) {
    // TODO
    //parpe::AmiciSummedGradientFunction<int>()
    //parpe::MultiConditionProblem p;
    // MultiCondition problem to interface -> hdf5 , non-hdf5 -> testing
}


TEST(steadystateProblemTests, testSteadystateHierarchical) {
    //TODO
}
