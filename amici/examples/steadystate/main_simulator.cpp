#include <cstdio>
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>
#include <standaloneSimulator.h>
#include <string>
#include <iostream>
#include <amici_model.h>
#include "steadyStateMultiConditionDataprovider.h"

std::unique_ptr<amici::Model> getModel();

int run(parpe::StandaloneSimulator &sim, std::string resultFileName, std::string resultPath,
        parpe::LoadBalancerMaster *loadBalancer) {

    // for each parameters in resultFileName
    auto model = getModel();
    std::vector<double> parameters(model->np(), 1);
    int errors = 0;
    errors += sim.run(resultFileName, resultPath, parameters, loadBalancer);
    return errors;
}

int main(int argc, char **argv) {
    if(argc != 5) {
        std::cerr<<"Error: wrong number of arguments.\n";
        std::cerr<<"Usage: ... INFILENAME INFILEPATH OUTFILENAME OUTFILEPATH";
        return EXIT_FAILURE;
    }

    std::string dataFileName = argv[1];
    std::string parameterPath = argv[2];
    std::string resultFileName = argv[3];
    std::string resultPath = argv[4];

    int status = 0;
    SteadyStateMultiConditionDataProvider dp(getModel(), dataFileName, "/inputData");
    parpe::StandaloneSimulator sim(&dp);

    int commSize = parpe::getMpiCommSize();

    if (commSize > 1) {
        if (parpe::getMpiRank() == 0) {
            parpe::LoadBalancerMaster loadBalancer;
            loadBalancer.run();
            status = run(sim, resultFileName, resultPath, &loadBalancer);
            loadBalancer.terminate();
            loadBalancer.sendTerminationSignalToAllWorkers();
        } else {
            parpe::LoadBalancerWorker lbw;
            lbw.run([&sim](std::vector<char> &buffer, int jobId) {
                sim.messageHandler(buffer, jobId);
            });
        }
    } else {
        status = run(sim, resultFileName, resultPath, nullptr);
    }

    return status;
}
