#ifndef MULTICONDITIONDATAPROVIDER_H
#define MULTICONDITIONDATAPROVIDER_H

#include <hdf5Misc.h>
#include <string>
#include <amici.h>

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

    void print();

    void sprint(char *buffer);
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

class MultiConditionDataProvider {
  public:
    MultiConditionDataProvider(amici::Model *model, std::string hdf5Filename);

    MultiConditionDataProvider(amici::Model *model, std::string hdf5Filename,
                               std::string rootPath);

    virtual ~MultiConditionDataProvider();

    /**
     * @brief Provides the number of conditions for which data is available and
     * simulations need to be run.
     * This is determined from the dimensions of the hdf5MeasurementPath
     * dataset.
     * @return Number of conditions
     */
    virtual int getNumberOfConditions() const;

    /**
     * @brief Number of model- oder optimization-parameters that are different
     * between the different conditions.
     * @return That number
     */
    virtual int getNumConditionSpecificParametersPerSimulation() const;

    /**
     * @brief Update fixed model parameters in of the passed UserData object for
     * the specified condition.
     * @param conditionIdx
     * @param udata The UserData instance to be updated
     * @return Status, 0 on success, non-zero otherwise
     */
    virtual int updateFixedSimulationParameters(int conditionIdx,
                                                amici::UserData *udata) const;

    virtual amici::ExpData *
    getExperimentalDataForCondition(int conditionIdx,
                                    const amici::UserData *udata) const;

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
    virtual int getNumOptimizationParameters() const;

    /**
     * @brief Get number of parameters common to all conditions.
     * @return Number of parameters
     */
    virtual int getNumCommonParameters() const;

    /**
     * @brief Returns a pointer to the underlying AMICI model
     * @return The model
     */
    virtual amici::Model *getModel() const;

    /**
     * @brief Returns a new UserData instance with options read from the HDF5
     * file.
     * Fixed and variable parameters are for no particular condition.
     * @return A new UserData instance.
     */
    virtual amici::UserData *getUserData() const;

    /**
     * @brief Returns a new UserData instance with options read from the HDF5
     * file.
     * Fixed parameters are set for the specified condition (variable parameters
     * are not).
     * @return A new UserData instance.
     */
    virtual amici::UserData *getUserDataForCondition(int conditionIdx) const;

    virtual int getIndexOfFirstConditionSpecificOptimizationParameter(
        int conditionIdx) const;

    void updateSimulationParameters(int conditionIndex, const double *optimizationParams,
        amici::UserData *udata) const;


    virtual void updateConditionSpecificSimulationParameters(
        int conditionIndex, const double *optimizationParams,
        amici::UserData *udata) const;

    hid_t getHdf5FileId() const;

  protected:
    MultiConditionDataProvider() = default;

    std::string hdf5MeasurementPath;
    std::string hdf5MeasurementSigmaPath;
    std::string hdf5ConditionPath;
    std::string hdf5AmiciOptionPath;
    std::string hdf5ParameterPath;

    hid_t fileId = 0;

    amici::Model *model = nullptr;
};

} // namespace parpe

#endif // MULTICONDITIONDATAPROVIDER_H
