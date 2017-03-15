#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

#include "queueworker.h"
// single out
#include "simulationworker.h"
#include "misc.h"

TEST_GROUP_C_SETUP(queueworker) {
}

TEST_GROUP_C_TEARDOWN(queueworker) {
}

TEST_C(queueworker, test_serializeResultPackageMessage) {
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

TEST_C(queueworker, test_serializeWorkPackageMessage) {
    // serialize and deserialize workpackage with random content
    int nTheta = randInt(0, 5000);
    int workPackageLength = getLengthWorkPackageMessage(nTheta);

    char *buffer = alloca(workPackageLength);
    workPackageMessage wp;
    wp.sensitivityMethod = randInt(INT_MIN, INT_MAX);
    wp.path.idxMultiStart = randInt(INT_MIN, INT_MAX);
    wp.path.idxLocalOptimization = randInt(INT_MIN, INT_MAX);
    wp.path.idxLocalOptimizationIteration = randInt(INT_MIN, INT_MAX);
    wp.path.idxGenotype = randInt(INT_MIN, INT_MAX);
    wp.path.idxExperiment = randInt(INT_MIN, INT_MAX);
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
    CHECK_EQUAL_C_INT(actPath.idxMultiStart, wp.path.idxMultiStart);
    CHECK_EQUAL_C_INT(actPath.idxLocalOptimization, wp.path.idxLocalOptimization);
    CHECK_EQUAL_C_INT(actPath.idxLocalOptimizationIteration, wp.path.idxLocalOptimizationIteration);
    CHECK_EQUAL_C_INT(actPath.idxGenotype, wp.path.idxGenotype);
    CHECK_EQUAL_C_INT(actPath.idxExperiment, wp.path.idxExperiment);

    for(int i = 0; i < nTheta; ++i) {
        CHECK_EQUAL_C_REAL(actTheta[i], wp.theta[i], 0);
    }
}
