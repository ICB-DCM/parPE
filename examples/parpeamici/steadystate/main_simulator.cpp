#include "steadyStateMultiConditionDataprovider.h"

#include <parpeamici/standaloneSimulator.h>
#include <parpecommon/misc.h>

#include <cstdio> // remove
#include <iostream>

std::unique_ptr<amici::Model> getModel();

int main(int argc, char **argv) {
    int status = EXIT_SUCCESS;

    parpe::initMpiIfNeeded(&argc, &argv);

    switch(argc)
    {
    case 6:
    {
        std::string dataFileName = argv[1];
        std::string dataFilePath = argv[2];
        std::string resultFileName = argv[3];
        std::string resultPath = argv[4];
        std::string simulationMode = argv[5];

        // TODO: testing-only remove result file
        remove(resultFileName.c_str());

        SteadyStateMultiConditionDataProvider dp(getModel(), dataFileName, dataFilePath + "/inputData");

        status = parpe::runSimulator(dp, simulationMode,
                                     dataFileName, dataFilePath,
                                     dataFileName, dataFilePath,
                                     resultFileName, resultPath);
        break;
    }
    case 8:
    {
        // simulate on test set: need optimizer result and test set data as inputs
        std::string conditionFileName = argv[1];
        std::string conditionFilePath = argv[2];
        std::string parameterFileName = argv[3];
        std::string parameterFilePath = argv[4];
        std::string resultFileName = argv[5];
        std::string resultPath = argv[6];
        std::string simulationMode = argv[7];

        // TODO: testing-only remove result file
        remove(resultFileName.c_str());

        SteadyStateMultiConditionDataProvider dp(getModel(), conditionFileName, conditionFilePath + "/inputData");

        status = parpe::runSimulator(dp, simulationMode,
                                     conditionFileName, conditionFilePath,
                                     parameterFileName, parameterFilePath,
                                     resultFileName, resultPath);
        break;
    }
    default:
        std::cerr<<"Error: wrong number of arguments.\n";
        std::cerr<<"Usage: ... CONDITION_FILE_NAME CONDITION_FILE_PATH [PARAMETER_FILE_NAME PARAMETER_FILE_PATH] OUTFILENAME OUTFILEPATH --at-optimum|--along-trajectory\n"; // |--parameter-matrix=PATH-UNSUPPORTED
        status = EXIT_FAILURE;
    }

    parpe::finalizeMpiIfNeeded();

    return status;
}
