#ifndef MULTICONDITIONDATAPROVIDER_H
#define MULTICONDITIONDATAPROVIDER_H

#include <hdf5Misc.h>
#include <optimizationOptions.h>

#include <amici/amici.h>
#include <gsl/gsl-lite.hpp>

#include <memory>
#include <string>
#include <vector>

#include <H5Cpp.h>


namespace parpe {

/**
 * @brief The MultiConditionDataProvider interface
 */
class MultiConditionDataProvider {
  public:

    virtual ~MultiConditionDataProvider() = default;

    /**
     * @brief Provides the number of conditions for which data is available and
     * simulations need to be run.
     * This is determined from the dimensions of the hdf5MeasurementPath
     * dataset.
     * @return Number of conditions
     */
    virtual int getNumberOfConditions() const = 0;

    virtual std::vector<int> getSimulationToOptimizationParameterMapping(int conditionIdx) const = 0;

    virtual void mapSimulationToOptimizationVariablesAddMultiply(
            int conditionIdx, gsl::span<double const> simulation, gsl::span<double> optimization, double coefficient = 1.0) const = 0;

    virtual void mapAndSetOptimizationToSimulationVariables(
            int conditionIdx, gsl::span<double const> optimization, gsl::span<double> simulation) const = 0;


    virtual amici::ParameterScaling getParameterScale(int optimizationParameterIndex) const = 0;


    virtual void updateSimulationParameters(int conditionIndex, gsl::span<double const> optimizationParams,
        amici::Model &model) const = 0;

    virtual std::unique_ptr<amici::ExpData> getExperimentalDataForCondition(int conditionIdx) const = 0;

    virtual std::vector<std::vector<double> > getAllMeasurements() const = 0;
    virtual std::vector<std::vector<double> > getAllSigmas() const = 0;

    /**
     * @brief Returns the number of optimization parameters of this problem
     * @return Number of parameters
     */
    virtual int getNumOptimizationParameters() const = 0;


    /**
     * @brief Returns a pointer to the underlying AMICI model
     * @return The model
     */
    virtual std::unique_ptr<amici::Model> getModel() const = 0;


    virtual std::unique_ptr<amici::Solver> getSolver() const = 0;

};


class MultiConditionDataProviderDefault : public MultiConditionDataProvider {
  public:
    MultiConditionDataProviderDefault(std::unique_ptr<amici::Model> model, std::unique_ptr<amici::Solver> solver);

    virtual ~MultiConditionDataProviderDefault() override = default;

    /**
     * @brief Provides the number of conditions for which data is available and
     * simulations need to be run.
     * This is determined from the dimensions of the hdf5MeasurementPath
     * dataset.
     * @return Number of conditions
     */
    virtual int getNumberOfConditions() const override;

    virtual std::vector<int> getSimulationToOptimizationParameterMapping(int conditionIdx) const override;

    virtual void mapSimulationToOptimizationVariablesAddMultiply(
            int conditionIdx, gsl::span<double const> simulation, gsl::span<double> optimization, double coefficient = 1.0) const override;

    virtual void mapAndSetOptimizationToSimulationVariables(
            int conditionIdx, gsl::span<double const> optimization, gsl::span<double> simulation) const override;


    virtual amici::ParameterScaling getParameterScale(int optimizationParameterIndex) const override;


    virtual void updateSimulationParameters(int conditionIndex, gsl::span<const double> optimizationParams,
        amici::Model &model) const override;

    virtual std::unique_ptr<amici::ExpData> getExperimentalDataForCondition(int conditionIdx) const override;

    virtual std::vector<std::vector<double> > getAllMeasurements() const override;
    virtual std::vector<std::vector<double> > getAllSigmas() const override;

    /**
     * @brief Returns the number of optimization parameters of this problem
     * @return Number of parameters
     */
    virtual int getNumOptimizationParameters() const override;


    /**
     * @brief Returns a pointer to the underlying AMICI model
     * @return The model
     */
    virtual std::unique_ptr<amici::Model> getModel() const override;


    virtual std::unique_ptr<amici::Solver> getSolver() const override;

    // TODO private
    std::vector<amici::ExpData> edata;

private:
    std::unique_ptr<amici::Model> model;
    std::unique_ptr<amici::Solver> solver;
};


/**
 * @brief The MultiConditionDataProvider class reads simulation data for
 * MultiConditionOptimizationProblem from a HDF5 file.
 *
 * This class assumes a certain layout of the underlying HDF5 file. Der dataset
 * names can be modified in hdf5*Path members.
 * Required dimensions:
 * * hdf5MeasurementPath, hdf5MeasurementSigmaPath: numObservables x
 * numConditions
 * * hdf5ConditionPath: numFixedParameters x numConditions
 * * hdf5AmiciOptionPath:
 * * hdf5ParameterPath:
 *
 * NOTE: The following dimensions are determined by the used AMICI model:
 * * numObservables := Model::ny
 * * numFixedParameters := Model::nk
 */

// TODO split; separate optimization from simulation
class MultiConditionDataProviderHDF5 : public MultiConditionDataProvider {
  public:
    MultiConditionDataProviderHDF5() = default;

    /**
     * @brief MultiConditionDataProvider
     * @param model A valid pointer to the amici::Model for which the data is to be provided.
     * @param hdf5Filename Path to the HDF5 file from which the data is to be read
     */
    MultiConditionDataProviderHDF5(std::unique_ptr<amici::Model> model, const std::string &hdf5Filename);

    /**
     * @brief See above.
     * @param model
     * @param hdf5Filename
     * @param rootPath The name of the HDF5 group under which the data is stored.
     */
    MultiConditionDataProviderHDF5(std::unique_ptr<amici::Model> model,
                                   std::string const& hdf5Filename,
                                   std::string const& rootPath);

    MultiConditionDataProviderHDF5(MultiConditionDataProviderHDF5 const&) = delete;

    virtual ~MultiConditionDataProviderHDF5() override = default;

    void openHdf5File(const std::string &hdf5Filename);

    /**
     * @brief Provides the number of conditions for which data is available and
     * simulations need to be run.
     * This is determined from the dimensions of the hdf5MeasurementPath
     * dataset.
     * @return Number of conditions
     */
    virtual int getNumberOfConditions() const override;

    virtual std::vector<int> getSimulationToOptimizationParameterMapping(int conditionIdx) const override;

    virtual void mapSimulationToOptimizationVariablesAddMultiply(
            int conditionIdx, gsl::span<double const> simulation, gsl::span<double> optimization, double coefficient = 1.0) const override;

    virtual void mapAndSetOptimizationToSimulationVariables(
            int conditionIdx, gsl::span<double const> optimization, gsl::span<double> simulation) const override;


    virtual amici::ParameterScaling getParameterScale(int optimizationParameterIndex) const override;

    /**
     * @brief Check if the data in the HDF5 file has consistent dimensions.
     * Aborts if not.
     */
    virtual void checkDataIntegrity() const;

    // void printInfo() const;

    virtual void readFixedSimulationParameters(int conditionIdx,
                                               double *buffer) const;

    virtual std::unique_ptr<amici::ExpData> getExperimentalDataForCondition(int conditionIdx) const override;

    std::vector<std::vector<double> > getAllMeasurements() const override;
    std::vector<std::vector<double> > getAllSigmas() const override;

    std::vector<double> getSigmaForConditionIndex(int conditionIdx) const;
    std::vector<double> getMeasurementForConditionIndex(int conditionIdx) const;

    /**
     * @brief getOptimizationParametersLowerBounds Get lower parameter bounds
     * NOTE: Currently the same bounds are assumed for kinetic parameters and
     * scaling parameters, ...
     * @param dataPath (not yet used)
     * @param buffer allocated memory to write parameter bounds
     */
    virtual void getOptimizationParametersLowerBounds(double *buffer) const;

    /**
     * @brief getOptimizationParametersUpperBounds Get upper parameter bounds
     * @param dataPath (not yet used)
     * @param buffer allocated memory to write parameter bounds
     */
    virtual void getOptimizationParametersUpperBounds(double *buffer) const;

    /**
     * @brief Returns the number of optimization parameters of this problem
     * @return Number of parameters
     */
    virtual int getNumOptimizationParameters() const override;


    /**
     * @brief Returns a pointer to a copy of the underlying AMICI model
     * as provided to the constructor
     * @return The model
     */
    virtual std::unique_ptr<amici::Model> getModel() const override;


    virtual std::unique_ptr<amici::Solver> getSolver() const override;

    /**
     * @brief Based on the array of optimization parameters, set the simulation
     * parameters in the given Model object to the ones for condition index.
     * @param conditionIndex
     * @param optimizationParams
     * @param udata
     */
    void updateSimulationParameters(int conditionIndex, gsl::span<const double> optimizationParams,
        amici::Model &model) const override;

    void copyInputData(const H5::H5File &target);

    /**
     * @brief Get the identifier of the used HDF5 file. Does not reopen. Do not close file.
     * @return The file ID
     */
    hid_t getHdf5FileId() const;

protected:
    void updateFixedSimulationParameters(int conditionIdx, amici::ExpData &edata) const;

    /**
     * @brief The model for which the data is to be read
     */
    std::unique_ptr<amici::Model> model;

    /**
     * @brief Absolute paths in the HDF5 file to the datasets
     * from which the respective data is to be read
     */
    std::string rootPath = "/";
    std::string hdf5MeasurementPath;
    std::string hdf5MeasurementSigmaPath;
    std::string hdf5ConditionPath;
    std::string hdf5ReferenceConditionPath;
    std::string hdf5AmiciOptionPath;
    std::string hdf5ParameterPath;
    std::string hdf5ParameterMinPath;
    std::string hdf5ParameterMaxPath;
    std::string hdf5ParameterScalingPath;
    std::string hdf5SimulationToOptimizationParameterMappingPath;

    /**
     * @brief HDF5 file handles for C++ and C API
     */
    H5::H5File file;

    std::unique_ptr<OptimizationOptions> optimizationOptions;
};

} // namespace parpe

#endif // MULTICONDITIONDATAPROVIDER_H
