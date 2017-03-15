#include "CppUTest/TestHarness_c.h"
#include "CppUTestExt/MockSupport_c.h"

// include .c because of static functions
#include "objectiveFunction.c"

// BEGIN required globals

#ifdef USE_MPE
// MPE event IDs for logging
const int mpe_event_begin_simulate, mpe_event_end_simulate;
const int mpe_event_begin_getrefs, mpe_event_end_getrefs;
const int mpe_event_begin_getdrugs, mpe_event_end_getdrugs;
const int mpe_event_begin_aggregate, mpe_event_end_aggregate;
#endif

// global mutex for HDF5 library calls
pthread_mutex_t mutexHDF;

// END required globals

/* TODO: test model/readSimulationUserData readSimulationExpData   getSimulationResults  */

TEST_GROUP_C_SETUP(objectivefunction) {
}

TEST_GROUP_C_TEARDOWN(objectivefunction) {
}

TEST_C(objectivefunction, test_reachedSteadyState) {
    int numStates = 4;
    double x[] = {1, 1, 1, 1};
    double xdot[] = {1e-6, 1e-6, 1e-6, 1e-6};
    int numTimepoints = 1;
    double tolerance = 1e-6;

    bool actual = reachedSteadyState(xdot, x, numTimepoints, numStates, tolerance);
    CHECK_EQUAL_C_BOOL(TRUE, actual);

    xdot[2] = 1e-5;
    actual = reachedSteadyState(xdot, x, numTimepoints, numStates, tolerance);
    CHECK_EQUAL_C_BOOL(FALSE, actual);
}
