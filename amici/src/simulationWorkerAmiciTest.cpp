#include "simulationWorkerAmici.h"
#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "misc.h"
#include "testingMisc.h"
#include <amici_model.h>

// clang-format off
TEST_GROUP(simulationWorkerAmici){
    void setup(){

    }

    void teardown(){
    }
};
// clang-format on


TEST(simulationWorkerAmici, testSerializeResultPackageMessage) {
    amici::UserData u(1, 2, 3);
    amici::Model m(1, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, amici::AMICI_O2MODE_NONE);

    parpe::JobResultAmiciSimulation results(1, std::make_unique<amici::ReturnData>(&u, &m), 2.1);

    int msgSize = 0;
    auto buffer = std::unique_ptr<char[]>(
                amici::serializeToChar<parpe::JobResultAmiciSimulation>(&results, &msgSize));

    parpe::JobResultAmiciSimulation resultsAct =
            amici::deserializeFromChar<parpe::JobResultAmiciSimulation>(buffer.get(), msgSize);

    CHECK_EQUAL(results.simulationTimeInSec, resultsAct.simulationTimeInSec);
}

TEST(simulationWorkerAmici, testSerializeWorkPackageMessage) {
    // serialize and deserialize workpackage with random content
    int nTheta = randInt(0, 5000);

    // setup userData
    amici::UserData expUdata(nTheta, 4, 5);
    expUdata.sensi_meth = amici::AMICI_SENSI_ASA;
    expUdata.sensi = amici::AMICI_SENSI_ORDER_FIRST;
    parpe::fillArrayRandomDoubleSameInterval(-1e-8, 1e8, nTheta, expUdata.p);

    // setup extra payload
    int lenDataElements = 10;
    std::vector<double> expData(lenDataElements);
    for (int i = 0; i < lenDataElements; ++i)
        expData[i] = parpe::randDouble(-1e-8, 1e8);

    parpe::JobAmiciSimulation<std::vector<double>> expWp(&expUdata, &expData);

    // serialize
    auto buffer = expWp.serialize();

    // deserialize
    amici::UserData actUdata(nTheta, 4, 5);
    std::vector<double> actData(lenDataElements);
    parpe::JobAmiciSimulation<std::vector<double>> actWp(&actUdata, &actData);
    actWp.deserialize(buffer.data(), buffer.size());

    // verify
    CHECK_EQUAL(actUdata.sensi, expUdata.sensi);
    CHECK_EQUAL(actUdata.sensi_meth, expUdata.sensi_meth);

    for (int i = 0; i < nTheta; ++i) {
        DOUBLES_EQUAL(actUdata.p[i], expUdata.p[i], 0);
    }
    for (int i = 0; i < lenDataElements; ++i) {
        DOUBLES_EQUAL(actData[i], expData[i], 0);
    }
}
