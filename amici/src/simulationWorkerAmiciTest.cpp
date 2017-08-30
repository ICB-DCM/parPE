#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "misc.h"
#include "simulationWorkerAmici.h"
#include "testingMisc.h"

TEST_GROUP(simulationWorkerAmici){void setup(){

}

                                  void teardown(){

                                  }};

IGNORE_TEST(simulationWorkerAmici, testSerializeResultPackageMessage) {
//    // serialize and deserialize resultpackage with random content
//    int nTheta = randInt(0, 5000);

//    // generate random data
//    UserData udataExp;
//    udataExp.sensi_meth = AMICI_SENSI_ASA;
//    udataExp.sensi = AMICI_SENSI_ORDER_FIRST;

//    ReturnData rdataExp();
//    int statusExp = randInt(INT_MIN, INT_MAX);
//    double llhExp = randDouble(1e-8, 1e8);
//    rdataExp.llh = &llhExp;
//    rdataExp.sllh = new double[nTheta];
//    for (int i = 0; i < nTheta; ++i) {
//        rdataExp.sllh[i] = randDouble(1e-8, 1e8);
//    }

//    int resultPackageLength = JobResultAmiciSimulation::getLength(nTheta);
//    char *buffer = new char[resultPackageLength];
//    JobResultAmiciSimulation::serialize(&rdataExp, &udataExp, statusExp,
//                                        buffer);

//    // deserialize
//    JobResultAmiciSimulation rpAct;
//    rpAct.sllh = new double[nTheta];
//    rpAct.deserialize(buffer);
//    delete[] buffer;

//    CHECK_EQUAL(rpAct.status, statusExp);
//    DOUBLES_EQUAL(rpAct.llh, llhExp, 0);
//    for (int i = 0; i < nTheta; ++i) {
//        DOUBLES_EQUAL(rpAct.sllh[i], rdataExp.sllh[i], 0);
//    }

//    delete[] rdataExp.sllh;
//    delete[] rpAct.sllh;
}

TEST(simulationWorkerAmici, testSerializeWorkPackageMessage) {
    // serialize and deserialize workpackage with random content
    int nTheta = randInt(0, 5000);

    // generate data
    JobAmiciSimulation expWp;
    int lenDataElements = 10;
    expWp.lenData = lenDataElements * sizeof(double);
    expWp.data = alloca(expWp.lenData);
    for (int i = 0; i < lenDataElements; ++i)
        ((double *)expWp.data)[i] = randDouble(-1e-8, 1e8);
    expWp.sensitivityMethod = AMICI_SENSI_ASA;
    expWp.numSimulationParameters = nTheta;
    expWp.simulationParameters = new double[nTheta];
    for (int i = 0; i < nTheta; ++i) {
        expWp.simulationParameters[i] = randDouble(-1e-8, 1e8);
    }

    // serialize
    int workPackageLength =
        JobAmiciSimulation::getLength(nTheta, expWp.lenData);
    char *buffer = new char[workPackageLength];
    expWp.serialize(buffer);

    // deserialize
    JobAmiciSimulation actWp;
    actWp.data = alloca(expWp.lenData);
    actWp.simulationParameters = new double[nTheta];

    actWp.deserialize(buffer);
    delete[] buffer;

    // check
    CHECK_EQUAL(actWp.sensitivityMethod, expWp.sensitivityMethod);

    for (int i = 0; i < nTheta; ++i) {
        DOUBLES_EQUAL(actWp.simulationParameters[i],
                      expWp.simulationParameters[i], 0);
    }
    for (int i = 0; i < lenDataElements; ++i) {
        DOUBLES_EQUAL(((double *)actWp.data)[i], ((double *)expWp.data)[i], 0);
    }

    delete[] actWp.simulationParameters;
    delete[] expWp.simulationParameters;
}
