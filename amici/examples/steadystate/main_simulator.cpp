#include "steadyStateMultiConditionDataprovider.h"
#include <standaloneSimulator.h>

#include <cstdio> // remove
#include <iostream>

#include <amici_model.h>

std::unique_ptr<amici::Model> getModel();

int main(int argc, char **argv) {
    if(argc != 6) {
        std::cerr<<"Error: wrong number of arguments.\n";
        std::cerr<<"Usage: ... INFILENAME INFILEPATH OUTFILENAME OUTFILEPATH --at-optimum|--along-trajectory|--parameter-matrix=PATH-UNSUPPORTED";
        return EXIT_FAILURE;
    }

    std::string dataFileName = argv[1];
    std::string dataFilePath = argv[2];
    std::string resultFileName = argv[3];
    std::string resultPath = argv[4];
    std::string simulationMode = argv[5];

    // TODO: testing-only remove result file
    remove(resultFileName.c_str());

    SteadyStateMultiConditionDataProvider dp(getModel(), dataFileName, dataFilePath + "/inputData");

    int status = parpe::runSimulator(dp, simulationMode, dataFileName, dataFilePath, resultFileName, resultPath);

    return status;
}
