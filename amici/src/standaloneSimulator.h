#ifndef STANDALONESIMULATOR_H
#define STANDALONESIMULATOR_H

#include <amici.h>
#include <udata.h>
#include <edata.h>
#include <rdata.h>
#include <multiConditionProblem.h>

namespace parpe {

/**
 * @brief The StandaloneSimulator class is for running simulations for a given dataset
 * and given parameters in parallel or in sequential mode and saving the simulation results.
 */
class StandaloneSimulator
{
public:
    StandaloneSimulator(MultiConditionDataProvider *dp);

    int run(const std::string &resultFile, const std::string &resultPath, std::vector<double> parameters, LoadBalancerMaster *loadBalancer);

    void messageHandler(std::vector<char> &buffer, int jobId);

private:

    JobResultAmiciSimulation runSimulation(amici::UserData *udata, JobIdentifier path,
                                    int jobId);

    MultiConditionDataProvider *dataProvider = nullptr;
};


} // namespace parpe

#endif // STANDALONESIMULATOR_H
