#include "simulationWorkerAmici.h"
#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "misc.h"
#include "testingMisc.h"
#include <amici/amici.h>

// required to be defined by file included below
#define NEW_OPTION_FILE "undefined"
#define HDFFILE "undefined"
#define HDFFILEWRITE "undefined"
#include "../tests/cpputest/testfunctions.h" // for Modell_Test

// clang-format off
TEST_GROUP(simulationWorkerAmici){
    void setup(){

    }

    void teardown(){
    }
};
// clang-format on


TEST(simulationWorkerAmici, testSerializeResultPackageMessage) {
//    amici::CVodeSolver s;
//    amici::Model m(1, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, amici::AMICI_O2MODE_NONE);

    parpe::JobResultAmiciSimulation results(std::make_unique<amici::ReturnData>(), 2.1);

    int msgSize = 0;
    auto buffer = std::unique_ptr<char[]>(
                amici::serializeToChar<parpe::JobResultAmiciSimulation>(results, &msgSize));

    parpe::JobResultAmiciSimulation resultsAct =
            amici::deserializeFromChar<parpe::JobResultAmiciSimulation>(buffer.get(), msgSize);

    CHECK_EQUAL(results.simulationTimeInSec, resultsAct.simulationTimeInSec);
}

TEST(simulationWorkerAmici, testSerializeWorkPackageMessage) {
    // serialize and deserialize workpackage with random content
    int nTheta = parpe::randInt(0, 5000);

    // setup userData
    amici::Model_Test m(0,0,0,0,0,0,0,0,0,0,0,0,0,0,amici::AMICI_O2MODE_NONE,std::vector<realtype>(nTheta),std::vector<realtype>(),std::vector<int>(),std::vector<realtype>(),std::vector<int>());
    amici::CVodeSolver expSolver;
    expSolver.setSensitivityMethod(amici::AMICI_SENSI_ASA);
    expSolver.setSensitivityOrder(amici::AMICI_SENSI_ORDER_FIRST);

    auto p = m.getParameters();
    parpe::fillArrayRandomDoubleSameInterval(-1e-8, 1e8, m.np(), p.data());
    m.setParameters(p);

    // setup extra payload
    int lenDataElements = 10;
    std::vector<double> expData(lenDataElements);
    for (int i = 0; i < lenDataElements; ++i)
        expData[i] = parpe::randDouble(-1e-8, 1e8);

    parpe::JobAmiciSimulation<std::vector<double>> expWp(&expSolver, &m, &expData);

    // serialize
    auto buffer = expWp.serialize();

    // deserialize
    std::vector<double> actData(lenDataElements);
    amici::Model_Test actModel(0,0,0,0,0,0,0,0,0,0,0,0,0,0,amici::AMICI_O2MODE_NONE,std::vector<realtype>(nTheta),std::vector<realtype>(),std::vector<int>(),std::vector<realtype>(),std::vector<int>());
    amici::CVodeSolver actSolver;

    parpe::JobAmiciSimulation<std::vector<double>> actWp(&actSolver, &actModel, &actData);
    actWp.deserialize(buffer.data(), buffer.size());

    // verify
    CHECK_EQUAL(actSolver.getSensitivityMethod(), expSolver.getSensitivityMethod());
    CHECK_EQUAL(actSolver.getSensitivityOrder(), expSolver.getSensitivityOrder());

    CHECK_TRUE(actModel.getParameters() == m.getParameters());
    for (int i = 0; i < lenDataElements; ++i) {
        DOUBLES_EQUAL(actData[i], expData[i], 0);
    }
}
