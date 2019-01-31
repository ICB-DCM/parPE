#ifndef STANDALONESIMULATOR_H
#define STANDALONESIMULATOR_H

#include <parpecommon/parpeConfig.h>
#include <parpeamici/multiConditionProblem.h>

#include <amici/amici.h>

namespace parpe {

/**
 * @brief The StandaloneSimulator class is for running simulations for a given dataset
 * and given parameters after optimization in parallel or in sequential mode and saving the simulation results.
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

    /**
     * @brief Run simulations for the given parameter and write results to an HDF5 file at the given location
     * @param resultFile Name of HDF5 output file. Will be created or appended.
     * @param resultPath HDF5 file root group name
     * @param optimizationParameters Parameters for simulation (results from hierarchical or standard optimization
     * @param loadBalancer LoadBalander instance for distributed memory parallel, or nullptr for shared memory parallel or sequential
     * @param inputFile File with simulation options and data used for optimization
     * @return Number of errors encountered
     */
    int run(const std::string &resultFile,
            const std::string &resultPath,
            std::vector<double> const& optimizationParameters,
            LoadBalancerMaster *loadBalancer,
            const H5::H5File &conditionFile, std::string conditionFilePath);

    void messageHandler(std::vector<char> &buffer, int jobId);

private:

    AmiciSimulationRunner::AmiciResultPackageSimple runSimulation(
            int conditionIdx,
            amici::Solver &solver,
            amici::Model &model);

    MultiConditionDataProvider *dataProvider = nullptr;

    /** Number of simulations to be sent to workers within one package (when running with MPI). */
    int maxSimulationsPerPackage = 8;
};


// enum class SimulatorOpType {finalParameters};

std::pair<int, double> getFunctionEvaluationWithMinimalCost(std::string const& datasetPath, H5::H5File const& file);


/**
 * @brief Read the final parameter set from parPE result file for the given local optimization index
 * @param startIndex
 * @param file
 * @return The final parameter vector
 */
std::vector<double> getFinalParameters(const std::string &startIndex, const H5::H5File &file);


std::vector<std::vector<double>> getParameterTrajectory(const std::string &startIndex, H5::H5File const& file);


int getNumStarts(const H5::H5File &file, const std::string &rootPath = "/");


int runFinalParameters(parpe::StandaloneSimulator &sim,
                       const std::string &conditionFileName,
                       const std::string &,
                       const std::string &parameterFileName,
                       const std::string &parameterFilePath,
                       const std::string &resultFileName,
                       const std::string &resultPath,
                       parpe::LoadBalancerMaster *loadBalancer);


int runAlongTrajectory(parpe::StandaloneSimulator &sim,
                       const std::string &conditionFileName,
                       const std::string &conditionFilePath,
                       const std::string &parameterFileName,
                       const std::string &parameterFilePath,
                       std::string const& resultFileName,
                       std::string const& resultPath,
                       parpe::LoadBalancerMaster *loadBalancer);


int runSimulator(MultiConditionDataProvider &dp,
                 std::string const& simulationMode,
                 const std::string &conditionFileName,
                 const std::string &conditionFilePath,
                 const std::string &parameterFileName,
                 const std::string &parameterFilePath,
                 std::string const& resultFileName,
                 std::string const& resultPath);

/**
 * @brief From the given parameter vector, extract outer optimization parameters, as defined in
 * the file HDF5 file parameterFile
 * @param fullParameters
 * @param parameterFile
 * @param parameterPath
 * @return
 */
std::vector<double> getOuterParameters(std::vector<double> const& fullParameters,
                                      H5::H5File const& parameterFile,
                                      std::string const& parameterPath);
} // namespace parpe

#endif // STANDALONESIMULATOR_H
