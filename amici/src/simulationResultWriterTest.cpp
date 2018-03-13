#include <simulationResultWriter.h>
#include <hdf5Misc.h>
#include "testingMisc.h"
#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <amici/amici.h>
#include <amici/solver_cvodes.h>

#include "../tests/cpputest/testfunctions.h" // for Modell_Test

// clang-format off
TEST_GROUP(simulationResultWriter){
    void setup(){

    }

    void teardown(){
    }
};
// clang-format on


TEST(simulationResultWriter, testResultWriter) {
    // setup
    char tmpName[TMP_MAX];
    std::tmpnam(tmpName);

    parpe::SimulationResultWriter rw(tmpName, "/testResultWriter/");

    rw.saveLlh = true;
    rw.saveX = true;
    rw.saveYMes = true;
    rw.saveYSim = true;

    constexpr int numSimulations = 2;

    amici::Model_Test model( 3, 3, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0,
                            amici::AMICI_O2MODE_NONE,std::vector<realtype>(),
                            std::vector<realtype>(),std::vector<int>(),
                            std::vector<realtype>(),std::vector<int>());

    std::vector<double> timepoints {1.0, 2.0};
    model.setTimepoints(timepoints);

    amici::ExpData edata(model);
    std::vector<double> measurements {1.1, 2.1, 3.1, 4.1};
    CHECK_TRUE(measurements.size() == (unsigned) model.nytrue * model.nt());
    edata.setObservedData(measurements.data());

    amici::CVodeSolver solver;
    amici::ReturnData rdata(solver, &model);
    rdata.x.resize(model.nt() * model.nx);
    std::iota(rdata.x.begin(), rdata.x.end(), 0);

    rdata.llh = 1.2345;
    rdata.y.resize(measurements.size());
    std::iota(rdata.y.begin(), rdata.y.end(), 10);

    auto file = rw.reopenFile();

    // write
    rw.createDatasets(model, numSimulations);

    CHECK_TRUE(parpe::hdf5GroupExists(file.getId(), "/testResultWriter/"));
    CHECK_TRUE(parpe::hdf5DatasetExists(file.getId(), rw.llhPath.c_str()));
    CHECK_TRUE(parpe::hdf5DatasetExists(file.getId(), rw.xPath.c_str()));
    CHECK_TRUE(parpe::hdf5DatasetExists(file.getId(), rw.yMesPath.c_str()));
    CHECK_TRUE(parpe::hdf5DatasetExists(file.getId(), rw.ySimPath.c_str()));

    rw.saveSimulationResults(&edata, &rdata, 1);

    // verify
    std::vector<double> xAct(rdata.x.size());
    parpe::hdf5Read3DDoubleHyperslab(file.getId(), rw.xPath.c_str(),
                                     1, model.nt(), model.nxtrue,
                                     1, 0, 0, xAct.data());
    parpe::checkEqualArray(rdata.x.data(), xAct.data(), xAct.size(), 1e-16, 1e-16);

    std::vector<double> yMesAct(measurements.size());
    parpe::hdf5Read3DDoubleHyperslab(file.getId(), rw.yMesPath.c_str(),
                                     1, model.nt(), model.ny, 1, 0, 0, yMesAct.data());
    parpe::checkEqualArray(measurements.data(), yMesAct.data(), yMesAct.size(), 1e-16, 1e-16);

}

TEST(simulationResultWriter, testResultWriterNewExistingFile) {
    char tmpName[TMP_MAX];
    std::tmpnam(tmpName);

    // create file
    parpe::SimulationResultWriter rw1(tmpName, "/testResultWriter/");

    // append
    parpe::SimulationResultWriter rw2(tmpName, "/testResultWriter/");
}
