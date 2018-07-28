#ifndef PARPE_AMICI_SIMULATION_RESULT_WRITER_H
#define PARPE_AMICI_SIMULATION_RESULT_WRITER_H

#include <string>
#include <amici/amici.h>

#include <gsl/gsl-lite.hpp>
#include <H5Cpp.h>

namespace parpe {

/**
 * @brief The SimulationResultWriter class saves AMICI simulation results
 * for one or multiple conditions to an HDF5 file.
 */

class SimulationResultWriter {
public:

    SimulationResultWriter() = default;

    SimulationResultWriter(const H5::H5File &file, std::string rootPath);

    SimulationResultWriter(std::string const& hdf5FileName, std::string  rootPath);

    /**
     * @brief Create results datasets. Condition index is first dimension.
     * @param udata
     * @param edata
     * @param numSimulations
     */

    void createDatasets(const amici::Model &model,
                        int numberOfSimulations = 1);

    /**
     * @brief Save results for a single simulation to HDF5 file.
     * @param udata
     * @param edata
     * @param rdata
     * @param simulationIdx If >= 0: write results in the selected
     * position of the result data sets (-> createDatasets)
     */

    void saveSimulationResults(const amici::ExpData *edata,
                               const amici::ReturnData *rdata,
                               int simulationIdx);

    void saveMeasurements(gsl::span<const double> measurements, int nt, int nytrue, int simulationIdx) const;

    void saveModelOutputs(gsl::span<const double> outputs, int nt, int nytrue, int simulationIdx) const;

    void saveStates(gsl::span<const double> states, int nt, int nx, int simulationIdx) const;

    void saveLikelihood(double llh, int simulationIdx) const;

    H5::H5File reopenFile();

    bool saveX = false;
    bool saveLlh = false;
//    bool saveSllh = false;
    bool saveYSim = false;
    bool saveYMes = false;
//    bool saveP = false;
//    bool saveK = false;

    std::string yMesPath;
    std::string ySimPath;
    std::string xPath;
    std::string llhPath;

private:
    void updatePaths();

    std::string rootPath;
    H5::H5File file;
};

} // namespace parpe

#endif
