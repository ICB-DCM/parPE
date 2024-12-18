#ifndef MULTICONDITIONDATAPROVIDER_H
#define MULTICONDITIONDATAPROVIDER_H

#include <parpecommon/hdf5Misc.h>
#include <parpeoptimization/optimizationOptions.h>

#include <amici/amici.h>

#include <gsl/gsl-lite.hpp>

#include <H5Cpp.h>

#include <memory>
#include <string>
#include <vector>

namespace parpe {

/**
 * @brief The MultiConditionDataProvider interface
 */
class MultiConditionDataProvider
{
  public:
    virtual ~MultiConditionDataProvider() = default;

    /**
     * @brief Provides the number of conditions for which data is available and
     * simulations need to be run.
     * @return Number of conditions
     */
    virtual int getNumberOfSimulationConditions() const = 0;

    /**
     * @brief Get mapping vector simulation_parameter_idx ->
     * optimization_parameter_idx for the given condition index.
     * @param conditionIdx
     * @return Mapping vector
     */
    virtual std::vector<int> getSimulationToOptimizationParameterMapping(
      int conditionIdx) const = 0;

    virtual void mapSimulationToOptimizationGradientAddMultiply(
      int conditionIdx,
      gsl::span<double const> simulation,
      gsl::span<double> optimization,
      gsl::span<const double> parameters,
      double coefficient = 1.0) const = 0;

    virtual void mapAndSetOptimizationToSimulationVariables(
      int conditionIdx,
      gsl::span<double const> optimization,
      gsl::span<double> simulation,
      gsl::span<amici::ParameterScaling> optimizationScale,
      gsl::span<amici::ParameterScaling> simulationScale) const = 0;

    /**
     * @brief Get the parameter scale for the given optimization parameter
     * @param simulationIdx
     * @return Parameter scale
     */
    virtual amici::ParameterScaling getParameterScaleOpt(
      int parameterIdx) const = 0;

    virtual std::vector<amici::ParameterScaling> getParameterScaleOpt()
      const = 0;

    /**
     * @brief Get the parameter scale vector for the given simulation
     * @param simulationIdx
     * @return
     */
    virtual std::vector<amici::ParameterScaling> getParameterScaleSim(
      int simulationIdx) const = 0;

    /**
     * @brief Get the parameter scale for the given parameter and simulation
     * @param simulationIdx
     * @return
     */
    virtual amici::ParameterScaling getParameterScaleSim(
      int simulationIdx,
      int modelParameterIdx) const = 0;

    virtual void updateSimulationParametersAndScale(
      int conditionIndex,
      gsl::span<double const> optimizationParams,
      amici::Model& model) const = 0;

    virtual std::unique_ptr<amici::ExpData> getExperimentalDataForCondition(
      int conditionIdx) const = 0;

    virtual std::vector<std::vector<double>> getAllMeasurements() const = 0;
    virtual std::vector<std::vector<double>> getAllSigmas() const = 0;

    /**
     * @brief Returns the number of optimization parameters of this problem
     * @return Number of parameters
     */
    virtual int getNumOptimizationParameters() const = 0;

    virtual std::vector<std::string> getProblemParameterIds() const = 0;

    /**
     * @brief Returns a pointer to the underlying AMICI model
     * @return The model
     */
    virtual std::unique_ptr<amici::Model> getModel() const = 0;

    virtual std::unique_ptr<amici::Solver> getSolver() const = 0;
};


/**
 * @brief In-memory data.
 *
 * !!Very limited implementation, currently only for testing!!
 */
class MultiConditionDataProviderDefault : public MultiConditionDataProvider
{
  public:
    MultiConditionDataProviderDefault(std::unique_ptr<amici::Model> model,
                                      std::unique_ptr<amici::Solver> solver);

    ~MultiConditionDataProviderDefault() override = default;

    /**
     * @brief Provides the number of conditions for which data is available and
     * simulations need to be run.
     * This is determined from the dimensions of the
     * MultiConditionDataProviderDefault::hdf5MeasurementPath
     * dataset.
     * @return Number of conditions
     */
    int getNumberOfSimulationConditions() const override;

    std::vector<int> getSimulationToOptimizationParameterMapping(
      int conditionIdx) const override;

    void mapSimulationToOptimizationGradientAddMultiply(
      int conditionIdx,
      gsl::span<double const> simulation,
      gsl::span<double> optimization,
      gsl::span<const double> parameters,
      double coefficient = 1.0) const override;

    void mapAndSetOptimizationToSimulationVariables(
      int conditionIdx,
      gsl::span<double const> optimization,
      gsl::span<double> simulation,
      gsl::span<amici::ParameterScaling> optimizationScale,
      gsl::span<amici::ParameterScaling> simulationScale) const override;

    std::vector<amici::ParameterScaling> getParameterScaleOpt()
      const override;

    amici::ParameterScaling getParameterScaleOpt(
      int optimizationParameterIndex) const override;

    amici::ParameterScaling getParameterScaleSim(
      int simulationIdx,
      int optimizationParameterIndex) const override;

    std::vector<amici::ParameterScaling> getParameterScaleSim(
      int) const override;

    void updateSimulationParametersAndScale(
      int conditionIndex,
      gsl::span<const double> optimizationParams,
      amici::Model& model) const override;

    std::unique_ptr<amici::ExpData> getExperimentalDataForCondition(
      int conditionIdx) const override;

    std::vector<std::vector<double>> getAllMeasurements()
      const override;
    std::vector<std::vector<double>> getAllSigmas() const override;

    /**
     * @brief Returns the number of optimization parameters of this problem
     * @return Number of parameters
     */
    int getNumOptimizationParameters() const override;

    /**
     * @brief Returns a pointer to the underlying AMICI model
     * @return The model
     */
    std::unique_ptr<amici::Model> getModel() const override;

    std::unique_ptr<amici::Solver> getSolver() const override;

    std::vector<std::string> getProblemParameterIds() const override;

    // TODO private
    std::vector<amici::ExpData> edata_;

  private:
    std::unique_ptr<amici::Model> model_;
    std::unique_ptr<amici::Solver> solver_;
};

/**
 * @brief The MultiConditionDataProvider class reads simulation data for
 * MultiConditionOptimizationProblem from a HDF5 file.
 *
 * This class assumes a certain layout of the underlying HDF5 file. Der dataset
 * names can be modified in hdf5*Path members.
 *
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
class MultiConditionDataProviderHDF5 : public MultiConditionDataProvider
{
  public:
    MultiConditionDataProviderHDF5() = default;

    /**
     * @brief MultiConditionDataProvider
     * @param model A valid pointer to the amici::Model for which the data is to
     * be provided.
     * @param hdf5Filename Path to the HDF5 file from which the data is to be
     * read
     */
    MultiConditionDataProviderHDF5(std::unique_ptr<amici::Model> model,
                                   const std::string& hdf5Filename);

    /**
     * @brief See above.
     * @param model
     * @param hdf5Filename
     * @param rootPath The name of the HDF5 group under which the data is
     * stored.
     */
    MultiConditionDataProviderHDF5(std::unique_ptr<amici::Model> model,
                                   std::string const& hdf5Filename,
                                   std::string const& rootPath);

    MultiConditionDataProviderHDF5(MultiConditionDataProviderHDF5 const&) =
      delete;

    ~MultiConditionDataProviderHDF5() override;

    /**
     * @brief Get the number of simulations required for objective function
     * evaluation. Currently, this amounts to the number
     * of conditions present in the data.
     * @return Number of conditions
     */
    int getNumberOfSimulationConditions() const override;

    /**
     * @brief Get index vector of length of model parameter with indices of
     * optimization parameters for the given condition.
     *
     * NOTE: This may contain -1 for parameter which are not mapped.
     * @param conditionIdx
     * @return
     */
    std::vector<int> getSimulationToOptimizationParameterMapping(
      int conditionIdx) const override;

    void mapSimulationToOptimizationGradientAddMultiply(
      int conditionIdx,
      gsl::span<double const> simulation,
      gsl::span<double> optimization,
      gsl::span<const double> parameters,
      double coefficient = 1.0) const override;

    void mapAndSetOptimizationToSimulationVariables(
      int conditionIdx,
      gsl::span<double const> optimization,
      gsl::span<double> simulation,
      gsl::span<amici::ParameterScaling> optimizationScale,
      gsl::span<amici::ParameterScaling> simulationScale) const override;

    std::vector<amici::ParameterScaling> getParameterScaleOpt()
      const override;

    amici::ParameterScaling getParameterScaleOpt(
      int parameterIdx) const override;

    std::vector<amici::ParameterScaling> getParameterScaleSim(
      int simulationIdx) const override;

    amici::ParameterScaling getParameterScaleSim(
      int simulationIdx,
      int modelParameterIdx) const override;

    /**
     * @brief Check if the data in the HDF5 file has consistent dimensions.
     * Aborts if not.
     */
    virtual void checkDataIntegrity() const;

    // void printInfo() const;

    virtual void readFixedSimulationParameters(int conditionIdx,
                                               gsl::span<double> buffer) const;

    std::unique_ptr<amici::ExpData> getExperimentalDataForCondition(
      int simulationIdx) const override;

    /**
     * @brief Get list of parameters w.r.t. which we need sensitivities
     * @param mapping Mapping from model simulation to objective parameters
     * @return AMICI's 'plist'
     */
    std::vector<int> getSensitivityParameterList(std::vector<int> const& mapping) const;

    std::vector<std::vector<double>> getAllMeasurements() const override;
    std::vector<std::vector<double>> getAllSigmas() const override;

    std::vector<double> getSigmaForSimulationIndex(int simulationIdx) const;
    std::vector<double> getMeasurementForSimulationIndex(
      int simulationIdx) const;

    /**
     * @brief Writes lower parameter bounds into the provided buffer
     * @param buffer allocated memory to write parameter bounds
     */
    virtual void getOptimizationParametersLowerBounds(
      gsl::span<double> buffer) const;

    /**
     * @brief Writes upper parameter bounds into the provided buffer
     * @param buffer allocated memory to write parameter bounds
     */
    virtual void getOptimizationParametersUpperBounds(
      gsl::span<double> buffer) const;

    /**
     * @brief Returns the number of optimization parameters of this problem
     * @return Number of parameters
     */
    int getNumOptimizationParameters() const override;

    /**
     * @brief Returns a pointer to a copy of the underlying AMICI model
     * as provided to the constructor
     * @return The model
     */
    std::unique_ptr<amici::Model> getModel() const override;

    std::unique_ptr<amici::Solver> getSolver() const override;

    /**
     * @brief Based on the array of optimization parameters, set the simulation
     * parameters in the given Model object to the ones for simulation index.
     * @param simulationIdx Index of the simulation condition for which to set
     * model parameters.
     * @param optimizationParams Problem parameters from which to extract
     * simulation parameters.
     * @param model Model on which to set parameter values and scale
     */
    void updateSimulationParametersAndScale(
      int simulationIdx,
      gsl::span<const double> optimizationParams,
      amici::Model& model) const override;

    void copyInputData(const H5::H5File& target);

    void getSimAndPreeqConditions(
      const int simulationIdx,
      int& preequilibrationConditionIdx,
      int& simulationConditionIdx
    ) const;

    std::vector<int> getReinitializationIndices(const int simulationIdx) const;

    /**
     * @brief Get a copy of the HDF5 file handle.
     * @return File handle
     */
    H5::H5File getHdf5File() const;

    void setModel(std::unique_ptr<amici::Model> model);

    std::vector<std::string> getProblemParameterIds() const override;

  protected:
    /**
     * @brief Update the constants in AMICI ExpData object. Reads a slab for the
     * given simulation from fixed parameters matrix.
     *
     * @param simulationIdx Index of the experimental condition for which the
     * parameters should be taken.
     * @param edata The object to be updated.
     */
    void updateFixedSimulationParameters(int simulationIdx,
                                         amici::ExpData& edata) const;


  private:
    /**
     * @brief The model for which the data is to be read
     */
    std::unique_ptr<amici::Model> model_;

    /**
     * @brief Absolute paths in the HDF5 file to the datasets
     * from which the respective data is to be read
     */
    std::string root_path_ = "/";
    std::string hdf5_measurement_path_;
    std::string hdf5_measurement_sigma_path_;
    std::string hdf5_condition_path_;
    std::string hdf5_reference_condition_path_;
    std::string hdf5_amici_options_path_;
    std::string hdf5_parameter_path_;
    std::string hdf5_parameter_min_path_;
    std::string hdf5_parameter_max_path_;
    std::string hdf5_parameter_scale_simulation_path_;
    std::string hdf5_parameter_scale_optimization_path_;
    std::string hdf5_simulation_to_optimization_parameter_mapping_path_;
    std::string hdf5_parameter_overrides_path;
    std::string hdf5_parameter_ids_path_;
    std::string hdf5_reinitialization_idxs_path_;

    /**
     * @brief HDF5 file handles for C++ and C API
     */
    H5::H5File file_;

    std::unique_ptr<OptimizationOptions> optimization_options_;
};

double
applyChainRule(double gradient,
               double parameter,
               amici::ParameterScaling oldScale,
               amici::ParameterScaling newScale);

} // namespace parpe

#endif // MULTICONDITIONDATAPROVIDER_H
