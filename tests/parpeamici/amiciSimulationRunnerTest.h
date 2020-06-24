#include <gtest/gtest.h>

#include <parpeamici/amiciSimulationRunner.h>
#include <parpecommon/misc.h>

#include "../parpecommon/testingMisc.h"

#include <amici/amici.h>


TEST(simulationWorkerAmici, testSerializeResultPackageMessage) {
    parpe::AmiciSimulationRunner::AmiciResultPackageSimple
            results = { 1.1, 2.345, std::vector<double>(1.0, 2.0),
                        std::vector<double>(3.0, 4),
                        std::vector<double>(3.3, 4),
                        std::vector<double>(1.0, 2), 10 };

    int msgSize = 0;
    auto buffer = std::unique_ptr<char[]>(
                amici::serializeToChar(results, &msgSize));

    parpe::AmiciSimulationRunner::AmiciResultPackageSimple resultsAct =
            amici::deserializeFromChar<
            parpe::AmiciSimulationRunner::AmiciResultPackageSimple>(
                buffer.get(), msgSize);

    EXPECT_EQ(resultsAct, results);
}
