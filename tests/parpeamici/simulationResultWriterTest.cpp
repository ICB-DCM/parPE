#include <gtest/gtest.h>

#include "../parpecommon/testingMisc.h"

#include <parpeamici/simulationResultWriter.h>
#include <parpecommon/hdf5Misc.h>

#include <vector>
#include <algorithm>
#include <numeric>

#include <amici/amici.h>
#include <amici/hdf5.h>
#include <amici/solver_cvodes.h>


TEST(SimulationResultWriter, ResultWriter) {
    // setup ResultWriter
    const char* tmpName = "parpeTest_testResultWriter.h5";
    auto _ = gsl::finally([tmpName] { remove(tmpName); });
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
    EXPECT_TRUE(measurements.size() == (unsigned) nytrue * timepoints.size());
    edata.setObservedData(measurements);

    amici::ReturnData rdata(
                timepoints,
                amici::ModelDimensions(nx, nx, nx, nx, 0, 0, 0, nytrue, nytrue,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       std::vector<int>(), 0, 0, 0, 0, 0, 0),
                0, 0, timepoints.size(), 0,
                std::vector<amici::ParameterScaling>(),
                amici::SecondOrderMode::none, amici::SensitivityOrder::none,
                amici::SensitivityMethod::none, amici::RDataReporting::full,
                true, true, 50);
    std::iota(rdata.x.begin(), rdata.x.end(), 0);
    rdata.llh = 1.2345;
    rdata.y.resize(measurements.size());
    std::iota(rdata.y.begin(), rdata.y.end(), 10);

    auto file = rw.reopenFile();

    // write
    rw.createDatasets(numSimulations);

    EXPECT_TRUE(parpe::hdf5GroupExists(file, "/testResultWriter/"));
    EXPECT_TRUE(file.nameExists(rw.llhPath));
    EXPECT_TRUE(file.nameExists(rw.xPath));
    EXPECT_TRUE(file.nameExists(rw.yMesPath));
    EXPECT_TRUE(file.nameExists(rw.ySimPath));
    EXPECT_TRUE(file.nameExists(rw.timePath));

    rw.saveSimulationResults(&edata, &rdata, 1);

    // verify
    hsize_t m, n;
    auto xAct = amici::hdf5::getDoubleDataset2D(file, rw.xPath + "/1", m, n);
    parpe::checkEqualArray(rdata.x.data(), xAct.data(), xAct.size(), 1e-16, 1e-16);

    auto yMesAct = amici::hdf5::getDoubleDataset2D(file, rw.yMesPath + "/1", m, n);
    parpe::checkEqualArray(measurements.data(), yMesAct.data(), yMesAct.size(), 1e-16, 1e-16);
}

TEST(SimulationResultWriter, ResultWriterNewExistingFile) {
    const char* tmpName = "parpeTest_testResultWriterNewExistingFile.h5";
    auto _ = gsl::finally([tmpName] { remove(tmpName); });

    // create file
    parpe::SimulationResultWriter rw1(tmpName, "/testResultWriter/");

    // append
    parpe::SimulationResultWriter rw2(tmpName, "/testResultWriter/");
}
