#ifndef MULTICONDITIONDATAPROVIDER_H
#define MULTICONDITIONDATAPROVIDER_H

#include <hdf5Misc.h>

#include <amici/amici.h>

#include <memory>
#include <string>
#include <vector>

#include <H5Cpp.h>


namespace parpe {

/** Struct to tell simulation workers which dataset they are operating on.
  */
struct JobIdentifier {
    /** current multistart batch (e.g. for crossvalidation) */
    int idxMultiStart = 0;

    /** current start index in multistart run */
    int idxLocalOptimization = 0;

    /** iteration of local solver or epoch for minibatch */
    int idxLocalOptimizationIteration = 0;

    // TODO int idxMiniBatch           /** current minibatch index */

    /** condition index (current data record) */
    // TODO Only this one is used for the moment
    int idxConditions = 0;

    void print() const;

    void sprint(char *buffer) const;
};


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
            int conditionIdx, const double *simulation, double *optimization, double coefficient = 1.0) const = 0;

    virtual void mapAndSetOptimizationToSimulationVariables(
            int conditionIdx, const double *optimization, double *simulation) const = 0;


    virtual amici::AMICI_parameter_scaling getParameterScale(int optimizationParameterIndex) const = 0;

    /**
     * @brief Update fixed model parameters for the specified condition.
     * @param conditionIdx
     * @param model The Model instance to be updated
     * @return Status, 0 on success, non-zero otherwise
     */
    virtual void updateFixedSimulationParameters(int conditionIdx,
                                                amici::Model &model) const = 0;

    virtual void updateSimulationParameters(int conditionIndex, const double *optimizationParams,
        amici::Model &model) const = 0;

    virtual std::unique_ptr<amici::ExpData> getExperimentalDataForCondition(int conditionIdx) const = 0;

    virtual std::vector<std::vector<double> > getAllMeasurements() const = 0;

    /**
     * @brief getOptimizationParametersLowerBounds Get lower parameter bounds
     * NOTE: Currently the same bounds are assumed for kinetic parameters and
     * scaling parameters, ...
     * @param dataPath (not yet used)
     * @param buffer allocated memory to write parameter bounds
     */
    virtual void getOptimizationParametersLowerBounds(double *buffer) const = 0;

    /**
     * @brief getOptimizationParametersUpperBounds Get upper parameter bounds
     * @param dataPath (not yet used)
     * @param buffer allocated memory to write parameter bounds
     */
    virtual void getOptimizationParametersUpperBounds(double *buffer) const = 0;

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

    /**
     * @brief Returns a new Model instance with options read from the HDF5
     * file.
     * Fixed parameters are set for the specified condition (variable parameters
     * are not).
     * @return A new Model instance.
     */
    virtual std::unique_ptr<amici::Model> getModelForCondition(int conditionIdx) const = 0;
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
 *
 * The vector of optimization variables is assumed to be [x_0, ...,
 * x_(numCommonParameter-1), conditionSpecificParameters].
 * conditionSpecificParameters := [cond0par0, cond0par1, ...,
 * cond0_par_(numConditionSpecificParametersPerSimulation-1),
 * cond_(numConditions-1)_(numConditionSpecificParametersPerSimulation-1) ]
 */

// TODO split; separate optimization from simulation
class MultiConditionDataProviderHDF5 : public MultiConditionDataProvider {
  public:
    MultiConditionDataProviderHDF5() = default;

    /**
     * @brief MultiConditionDataProvider
     * @param model A valid pointer to the amici::Model for which the data is to be provided.
     * The user is responsible for deleting the Model.
     * @param hdf5Filename Path to the HDF5 file from which the data is to be read
     */
    MultiConditionDataProviderHDF5(std::unique_ptr<amici::Model> model, std::string hdf5Filename);

    /**
     * @brief See above.
     * @param model
     * @param hdf5Filename
     * @param rootPath The name of the HDF5 group under which the data is stored.
     */
    MultiConditionDataProviderHDF5(std::unique_ptr<amici::Model> model, std::string hdf5Filename,
                               std::string rootPath);

    virtual ~MultiConditionDataProviderHDF5() = default;

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
            int conditionIdx, const double *simulation, double *optimization, double coefficient = 1.0) const override;

    virtual void mapAndSetOptimizationToSimulationVariables(
            int conditionIdx, const double *optimization, double *simulation) const override;


    virtual amici::AMICI_parameter_scaling getParameterScale(int optimizationParameterIndex) const override;

    /**
     * @brief Check if the data in the HDF5 file has consistent dimensions.
     * Aborts if not.
     */
    virtual void checkDataIntegrity() const;

    // void printInfo() const;


    /**
     * @brief Update fixed model parameters in of the passed UserData object for
     * the specified condition.
     * @param conditionIdx
     * @param udata The UserData instance to be updated
     */
    virtual void updateFixedSimulationParameters(int conditionIdx,
                                                amici::Model &model) const override;

    virtual std::unique_ptr<amici::ExpData> getExperimentalDataForCondition(int conditionIdx) const override;

    std::vector<std::vector<double> > getAllMeasurements() const override;

    /**
     * @brief getOptimizationParametersLowerBounds Get lower parameter bounds
     * NOTE: Currently the same bounds are assumed for kinetic parameters and
     * scaling parameters, ...
     * @param dataPath (not yet used)
     * @param buffer allocated memory to write parameter bounds
     */
    virtual void getOptimizationParametersLowerBounds(double *buffer) const override;

    /**
     * @brief getOptimizationParametersUpperBounds Get upper parameter bounds
     * @param dataPath (not yet used)
     * @param buffer allocated memory to write parameter bounds
     */
    virtual void getOptimizationParametersUpperBounds(double *buffer) const override;

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

    /**
     * @brief Returns a new Model instance with options read from the HDF5
     * file.
     * Fixed parameters are set for the specified condition (variable parameters
     * are not).
     * @return A new Model instance.
     */
    virtual std::unique_ptr<amici::Model> getModelForCondition(int conditionIdx) const override;

    /**
     * @brief Based on the array of optimization parameters, set the simulation
     * parameters in the given UserData object to the ones for condition index.
     * @param conditionIndex
     * @param optimizationParams
     * @param udata
     */
    void updateSimulationParameters(int conditionIndex, const double *optimizationParams,
        amici::Model &model) const override;

    void copyInputData(H5::H5File target);

    /**
     * @brief Get the identifier of the used HDF5 file. Does not reopen. Do not close file.
     * @return The file ID
     */
    hid_t getHdf5FileId() const;

protected:
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
    std::string hdf5AmiciOptionPath;
    std::string hdf5ParameterPath;
    std::string hdf5ParameterMinPath;
    std::string hdf5ParameterMaxPath;
    std::string hdf5ParameterScalingPath;

    /**
     * @brief HDF5 file handles for C++ and C API
     */
    H5::H5File file;
    hid_t fileId = 0;

};

} // namespace parpe



namespace boost {
namespace serialization {

template <class Archive>
void serialize(Archive &ar, parpe::JobIdentifier &d, const unsigned int version) {
    ar & d.idxMultiStart;
    ar & d.idxLocalOptimization;
    ar & d.idxLocalOptimizationIteration;
    ar & d.idxConditions;
}


} // namespace serialization
} // namespace boost

#endif // MULTICONDITIONDATAPROVIDER_H
