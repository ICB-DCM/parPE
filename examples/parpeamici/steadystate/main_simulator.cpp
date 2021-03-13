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
                 "PARAMETER_FILE_NAME PARAMETER_FILE_PATH "
                 "OUTFILENAME OUTFILEPATH "
                 "--at-optimum|--along-trajectory|--nominal "
                 "--mpi|--nompi --compute-inner|--nocompute-inner\n";
    // |--parameter-matrix=PATH-UNSUPPORTED
}

int main(int argc, char **argv) {
    int status = EXIT_SUCCESS;

    if(argc != 10) {
        printUsage();
        return EXIT_FAILURE;
    }

    bool computeInner;
    if(std::string(argv[argc -1]) == "--compute-inner") {
        computeInner = true;
    } else if(std::string(argv[argc -1]) == "--nocompute-inner") {
        computeInner = false;
    } else {
        printUsage();
        return EXIT_FAILURE;
    }

    if(std::string(argv[argc -2]) == "--mpi") {
#ifdef PARPE_ENABLE_MPI
        MPI_Init(&argc, &argv);
#else
        throw std::runtime_error("parPE was built without MPI support.");
#endif
    } else if(std::string(argv[argc -2]) == "--nompi") {
        ;
    } else {
        printUsage();
        return EXIT_FAILURE;
    }

    // simulate on test set: need optimizer result and test set data as inputs
    std::string conditionFileName = argv[1];
    std::string conditionFilePath = argv[2];
    std::string parameterFileName = argv[3];
    std::string parameterFilePath = argv[4];
    std::string resultFileName = argv[5];
    std::string resultPath = argv[6];
    std::string simulationMode = argv[7];

    SteadyStateMultiConditionDataProvider dp(
        amici::generic_model::getModel(),
        conditionFileName,
        conditionFilePath);

    status = parpe::runSimulator(dp, simulationMode,
                                 conditionFileName, conditionFilePath,
                                 parameterFileName, parameterFilePath,
                                 resultFileName, resultPath, computeInner);

    parpe::finalizeMpiIfNeeded();

    return status;
}
