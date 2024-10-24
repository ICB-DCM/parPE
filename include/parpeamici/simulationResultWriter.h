#ifndef PARPE_AMICI_SIMULATION_RESULT_WRITER_H
#define PARPE_AMICI_SIMULATION_RESULT_WRITER_H

#include <amici/amici.h>
#include <string>

#include <H5Cpp.h>
#include <gsl/gsl-lite.hpp>

namespace parpe {

/**
 * @brief The SimulationResultWriter class saves AMICI simulation results
 * for one or multiple conditions to an HDF5 file.
 *
 * Will write to existing file or create new one. Will fail if the dataset
 * to be written already exist.
 *
 * Structure is
 *
 * ```
 * rootPath/
 * +- yMes/$conditionIdx [double: nt x ny] (nt may be condition-specific)
 * +- ySim/$conditionIdx [double: nt x ny] (nt may be condition-specific)
 * +- llh [double: nConditions]
 * +- ...
 * ```
 */

class SimulationResultWriter {

  public:
    SimulationResultWriter() = default;

    /**
     * @brief SimulationResultWriter
     * @param file HDF5 file object to write to
     * @param rootPath Path prefix inside HDF5 file
     */
    SimulationResultWriter(const H5::H5File& file, std::string rootPath);

    /**
     * @brief SimulationResultWriter
     * @param hdf5FileName HDF5 file to create or open for appending
     * @param rootPath Path prefix inside HDF5 file
     */
    SimulationResultWriter(
        std::string const& hdf5FileName,
        std::string rootPath);

    // Implement me
    SimulationResultWriter(SimulationResultWriter const&) = delete;

    ~SimulationResultWriter();

    /**
     * @brief Create results datasets.
     *
     * Must be called before first call to `save*`
     * @param numSimulations
     */
    void createDatasets(hsize_t numSimulations);

    void createDatasets(int numberOfSimulations = 1);

    /**
     * @brief Save results for a single simulation to HDF5 file.
     * @param edata
     * @param rdata
     * @param simulationIdx If >= 0: write results in the selected
     * position of the result data sets (-> createDatasets)
     */

    void saveSimulationResults(
        amici::ExpData const* edata,
        amici::ReturnData const* rdata,
        int simulationIdx);

    void
    saveTimepoints(gsl::span<double const> timepoints, int simulationIdx) const;

    void saveMeasurements(
        gsl::span<double const> measurements,
        int nt,
        int nytrue,
        int simulationIdx) const;

    void saveModelOutputs(
        gsl::span<double const> outputs,
        int nt,
        int nytrue,
        int simulationIdx) const;

    void saveStates(
        gsl::span<double const> states,
        int nt,
        int nx,
        int simulationIdx) const;

    void
    saveParameters(gsl::span<double const> parameters, int simulationIdx) const;

    void saveLikelihood(double llh, int simulationIdx) const;

    H5::H5File reopenFile();

    bool saveX = false;
    bool saveLlh = false;
    //    bool saveSllh = false;
    bool saveYSim = false;
    bool saveYMes = false;
    bool save_parameters_ = false;
    //    bool saveK = false;

    std::string yMesPath;
    std::string ySimPath;
    std::string xPath;
    std::string llhPath;
    std::string timePath;
    std::string parametersPath;

  private:
    void updatePaths();

    std::string rootPath;
    H5::H5File file;
};

} // namespace parpe

#endif
