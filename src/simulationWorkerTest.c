#include "simulationWorker.h"
#include "testingMisc.h"
// TODO: test separetely
#include "../objectiveFunctionBenchmarkModel/dataprovider.h"
#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"


TEST_GROUP_C_SETUP(simulationWorker) {
}

TEST_GROUP_C_TEARDOWN(simulationWorker) {
}

TEST_C(simulationWorker, testSerializeResultPackageMessage) {
    // serialize and deserialize resultpackage with random content
    int nTheta = randInt(0, 5000);
    int resultPackageLength = getLengthResultPackageMessage(nTheta);

    char *buffer = alloca(resultPackageLength);

    resultPackageMessage rp;
    rp.status = randInt(INT_MIN, INT_MAX);
    rp.llh = randDouble(DBL_MIN, DBL_MAX);
    rp.sllh = alloca(nTheta * sizeof(*rp.sllh));
    for(int i = 0; i < nTheta; ++i) {
        rp.sllh[i] = randDouble(DBL_MIN, DBL_MAX);
    }

    serializeResultPackageMessage(rp, nTheta, buffer);

    int actStatus;
    double actLlh;
    double *actSllh = alloca(nTheta * sizeof(*actSllh));

    deserializeResultPackageMessage(buffer, nTheta, &actStatus, &actLlh, actSllh);

    CHECK_EQUAL_C_INT(actStatus, rp.status);
    CHECK_EQUAL_C_REAL(actLlh, rp.llh, 0);
    for(int i = 0; i < nTheta; ++i) {
        CHECK_EQUAL_C_REAL(actSllh[i], rp.sllh[i], 0);
    }
}

TEST_C(simulationWorker, testSerializeWorkPackageMessage) {
    // serialize and deserialize workpackage with random content
    int nTheta = randInt(0, 5000);
    int workPackageLength = getLengthWorkPackageMessage(nTheta);

    char *buffer = alloca(workPackageLength);
    workPackageMessage wp;
    wp.sensitivityMethod = randInt(INT_MIN, INT_MAX);
    datapath expPath;
    expPath.idxMultiStart = randInt(INT_MIN, INT_MAX);
    expPath.idxLocalOptimization = randInt(INT_MIN, INT_MAX);
    expPath.idxLocalOptimizationIteration = randInt(INT_MIN, INT_MAX);
    expPath.idxGenotype = randInt(INT_MIN, INT_MAX);
    expPath.idxExperiment = randInt(INT_MIN, INT_MAX);
    wp.data = &expPath;
    wp.lenData = sizeof(datapath);
    wp.theta = alloca(nTheta * sizeof(*wp.theta));
    for(int i = 0; i < nTheta; ++i) {
        wp.theta[i] = randDouble(-DBL_MIN, DBL_MAX);
    }

    serializeWorkPackageMessage(wp, nTheta, buffer);

    int actSensitivityMethod;
    datapath actPath;
    double *actTheta = alloca(nTheta * sizeof(*actTheta));

    deserializeWorkPackageMessage(buffer, nTheta, &actPath, actTheta, &actSensitivityMethod);

    CHECK_EQUAL_C_INT(actSensitivityMethod, wp.sensitivityMethod);
    CHECK_EQUAL_C_INT(actPath.idxMultiStart, expPath.idxMultiStart);
    CHECK_EQUAL_C_INT(actPath.idxLocalOptimization, expPath.idxLocalOptimization);
    CHECK_EQUAL_C_INT(actPath.idxLocalOptimizationIteration, expPath.idxLocalOptimizationIteration);
    CHECK_EQUAL_C_INT(actPath.idxGenotype, expPath.idxGenotype);
    CHECK_EQUAL_C_INT(actPath.idxExperiment, expPath.idxExperiment);

    for(int i = 0; i < nTheta; ++i) {
        CHECK_EQUAL_C_REAL(actTheta[i], wp.theta[i], 0);
    }
}
