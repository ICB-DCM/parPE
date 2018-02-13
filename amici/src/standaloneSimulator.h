#ifndef STANDALONESIMULATOR_H
#define STANDALONESIMULATOR_H

#include <multiConditionProblem.h>

#include <amici.h>
#include <edata.h>
#include <rdata.h>

namespace parpe {

/**
 * @brief The StandaloneSimulator class is for running simulations for a given dataset
 * and given parameters in parallel or in sequential mode and saving the simulation results.
 *
 * Command line interface should support:
 * ./simulate --at-optimum : use parameters from last iteration of all multi-start optimization runs
 *            --parameter-matrix : using arbitrary parameters from some matrix in hdf5 file
 *            --along-trajectory : use parameters along the optimization trajectory of all multi-start optimization runs
 */
class StandaloneSimulator
{
public:
    StandaloneSimulator(MultiConditionDataProvider *dp);

    int run(const std::string &resultFile, const std::string &resultPath, std::vector<double> parameters, LoadBalancerMaster *loadBalancer);

    void messageHandler(std::vector<char> &buffer, int jobId);

private:

    JobResultAmiciSimulation runSimulation(JobIdentifier path,
                                    int jobId);

    MultiConditionDataProvider *dataProvider = nullptr;
};

/**
 * @brief Read the final parameter set from parPE result file for the given local optimization index
 * @param startIndex
 * @param file
 * @return The final parameter vector
 */
std::vector<double> getFinalParameters(std::string startIndex, H5::H5File &file);

std::vector<std::vector<double>> getParameterTrajectory(std::string startIndex, H5::H5File &file);

int getNumStarts(H5::H5File file, std::string rootPath = "/");

int runFinalParameters(parpe::StandaloneSimulator &sim, std::string inFileName, std::string resultFileName, std::string resultPath,
        parpe::LoadBalancerMaster *loadBalancer);

int runAlongTrajectory(parpe::StandaloneSimulator &sim, std::string inFileName, std::string resultFileName, std::string resultPath,
        parpe::LoadBalancerMaster *loadBalancer);

} // namespace parpe

#endif // STANDALONESIMULATOR_H
