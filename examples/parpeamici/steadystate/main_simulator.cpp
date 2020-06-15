#include <parpecommon/parpeConfig.h>

#include "steadyStateMultiConditionDataprovider.h"

#include <parpeamici/standaloneSimulator.h>
#include <parpecommon/misc.h>
#include <parpecommon/parpeConfig.h>

#include <cstdio> // remove
#include <iostream>
#include <stdexcept>

#ifdef PARPE_ENABLE_MPI
#include <mpi.h>
#endif

namespace amici::generic_model {
std::unique_ptr<amici::Model> getModel();
}


void printUsage() {
    std::cerr<<"Error: wrong number of arguments.\n";
    std::cerr<<"Usage: ... CONDITION_FILE_NAME CONDITION_FILE_PATH "
               "[PARAMETER_FILE_NAME PARAMETER_FILE_PATH] "
               "OUTFILENAME OUTFILEPATH "
               "--at-optimum|--along-trajectory "
               "--mpi|--nompi\n";
    // |--parameter-matrix=PATH-UNSUPPORTED
}

int main(int argc, char **argv) {
    int status = EXIT_SUCCESS;

    if(argc != 7 && argc != 9) {
        printUsage();
        return EXIT_FAILURE;
    }


    if(std::string(argv[argc -1]) == "--mpi") {
#ifdef PARPE_ENABLE_MPI
        MPI_Init(&argc, &argv);
#else
        throw std::runtime_error("parPE was built without MPI support.");
#endif
    } else if(std::string(argv[argc -1]) == "--nompi") {
        ;
    } else {
        printUsage();
        return EXIT_FAILURE;
    }

    if(argc == 7) {
        std::string dataFileName = argv[1];
        std::string dataFilePath = argv[2];
        std::string resultFileName = argv[3];
        std::string resultPath = argv[4];
        std::string simulationMode = argv[5];

        // TODO: testing-only remove result file
        if (parpe::getMpiRank() < 1)
            remove(resultFileName.c_str());

        SteadyStateMultiConditionDataProvider dp(
            amici::generic_model::getModel(), dataFileName,
            dataFilePath + "/inputData");

        status = parpe::runSimulator(dp, simulationMode,
                                     dataFileName, dataFilePath,
                                     dataFileName, dataFilePath,
                                     resultFileName, resultPath);
    } else if(argc == 9) {
        // simulate on test set: need optimizer result and test set data as inputs
        std::string conditionFileName = argv[1];
        std::string conditionFilePath = argv[2];
        std::string parameterFileName = argv[3];
        std::string parameterFilePath = argv[4];
        std::string resultFileName = argv[5];
        std::string resultPath = argv[6];
        std::string simulationMode = argv[7];

        // TODO: testing-only remove result file
        if (parpe::getMpiRank() < 1)
            remove(resultFileName.c_str());

        auto dpPath = conditionFilePath;
        {
            // check if this is a result file or a new input file
            // TODO: this should be handled cleaner
            auto file = parpe::hdf5OpenForReading(conditionFileName);
            if(parpe::hdf5GroupExists(
                        file.getId(),
                        (conditionFilePath + "/inputData").c_str()))
                dpPath = conditionFilePath + "/inputData";
        }

        SteadyStateMultiConditionDataProvider dp(
            amici::generic_model::getModel(), conditionFileName, dpPath);

        status = parpe::runSimulator(dp, simulationMode,
                                     conditionFileName, conditionFilePath,
                                     parameterFileName, parameterFilePath,
                                     resultFileName, resultPath);
    }

    parpe::finalizeMpiIfNeeded();

    return status;
}
