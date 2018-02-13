#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>
#include <standaloneSimulator.h>
#include <optimizationOptions.h>

#include <cstdio>
#include <string>
#include <iostream>

#include <amici_model.h>

#include "steadyStateMultiConditionDataprovider.h"

enum class SimulatorOpType {finalParameters};

std::unique_ptr<amici::Model> getModel();

int run(parpe::StandaloneSimulator &sim, std::string simulationMode,
        std::string inFileName, std::string dataFilePath,
        std::string resultFileName, std::string resultPath,
        parpe::LoadBalancerMaster *loadBalancer) {

    if(simulationMode == "--at-optimum") {
        return parpe::runFinalParameters(sim, inFileName, resultFileName, resultPath, loadBalancer);
    } else if (simulationMode == "--along-trajectory") {
        return parpe::runAlongTrajectory(sim, inFileName, resultFileName, resultPath, loadBalancer);
    }

    return -1;
}


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

    int status = 0;
    SteadyStateMultiConditionDataProvider dp(getModel(), dataFileName, dataFilePath + "/inputData");
    parpe::StandaloneSimulator sim(&dp);

    int commSize = parpe::getMpiCommSize();

    if (commSize > 1) {
        if (parpe::getMpiRank() == 0) {
            parpe::LoadBalancerMaster loadBalancer;
            loadBalancer.run();
            status = run(sim, simulationMode, dataFileName, dataFilePath, resultFileName, resultPath, &loadBalancer);
            loadBalancer.terminate();
            loadBalancer.sendTerminationSignalToAllWorkers();
        } else {
            parpe::LoadBalancerWorker lbw;
            lbw.run([&sim](std::vector<char> &buffer, int jobId) {
                sim.messageHandler(buffer, jobId);
            });
        }
    } else {
        status = run(sim, simulationMode, dataFileName, dataFilePath, resultFileName, resultPath, nullptr);
    }

    return status;
}
