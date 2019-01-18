#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testingMisc.h"

#include <simulationResultWriter.h>
#include <hdf5Misc.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <amici/amici.h>
#include <amici/solver_cvodes.h>


// clang-format off
TEST_GROUP(simulationResultWriter){
    void setup(){

    }

    void teardown(){
    }
};
// clang-format on


TEST(simulationResultWriter, testResultWriter) {
    // setup ResultWriter
    char tmpName[TMP_MAX];
    std::tmpnam(tmpName);
    parpe::SimulationResultWriter rw(tmpName, "/testResultWriter/");

    rw.saveLlh = true;
    rw.saveX = true;
    rw.saveYMes = true;
    rw.saveYSim = true;

    // setup data
    constexpr int numSimulations = 2;
    constexpr int nx = 3;
    constexpr int nytrue = 2;
    const std::vector<double> timepoints {1.0, 2.0};

    amici::ExpData edata(nytrue, 0, 0, timepoints);
    std::vector<double> measurements {1.1, 2.1, 3.1, 4.1};
    CHECK_TRUE(measurements.size() == (unsigned) nytrue * timepoints.size());
    edata.setObservedData(measurements);

    amici::ReturnData rdata(timepoints, 0, 1, nx, nx, nx, nytrue, nytrue, 0, 0,
                            0, 0, 0, 0, timepoints.size(), 0,
                            std::vector<amici::ParameterScaling>(), amici::SecondOrderMode::none,
                            amici::SensitivityOrder::none, amici::SensitivityMethod::none);
    std::iota(rdata.x.begin(), rdata.x.end(), 0);
    rdata.llh = 1.2345;
    rdata.y.resize(measurements.size());
    std::iota(rdata.y.begin(), rdata.y.end(), 10);

    auto file = rw.reopenFile();

    // write
    rw.createDatasets(nytrue, nx, timepoints.size(), numSimulations);

    CHECK_TRUE(parpe::hdf5GroupExists(file.getId(), "/testResultWriter/"));
    CHECK_TRUE(parpe::hdf5DatasetExists(file.getId(), rw.llhPath.c_str()));
    CHECK_TRUE(parpe::hdf5DatasetExists(file.getId(), rw.xPath.c_str()));
    CHECK_TRUE(parpe::hdf5DatasetExists(file.getId(), rw.yMesPath.c_str()));
    CHECK_TRUE(parpe::hdf5DatasetExists(file.getId(), rw.ySimPath.c_str()));

    rw.saveSimulationResults(&edata, &rdata, 1);

    // verify
    auto xAct = parpe::hdf5Get3DDoubleHyperslab(
                file.getId(), rw.xPath.c_str(),
                1, timepoints.size(), nx,
                1, 0, 0);
    parpe::checkEqualArray(rdata.x.data(), xAct.data(), xAct.size(), 1e-16, 1e-16);

    auto yMesAct = parpe::hdf5Get3DDoubleHyperslab(
                file.getId(), rw.yMesPath.c_str(),
                1, timepoints.size(), nytrue, 1, 0, 0);
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
