#ifndef MULTICONDITIONDATAPROVIDER_H
#define MULTICONDITIONDATAPROVIDER_H

#include <hdf5Misc.h>
#include <string>
#include <udata.h>

class ExpData;

/** Struct to tell simulation workers which dataset they are operating on
  */
typedef struct JobIdentifier_tag {
    int idxMultiStart;                 /** current multistart batch (e.g. for crossvalidation) */
    int idxLocalOptimization;          /** current start index in multistart run */
    int idxLocalOptimizationIteration; /** iteration of local solver or epoch for minibatch */
    // TODO int idxMiniBatch           /** current minibatch index */
    int idxConditions;                 /** experiment index */ // TODO Only this one is used for the moment
} JobIdentifier;

void printJobIdentifier(JobIdentifier id);

void sprintJobIdentifier(char *buffer, JobIdentifier id);


/**
 * @brief The MultiConditionDataProvider class reads simulation data for MultiConditionOptimizationProblem from a HDF5 file.
 *
 * This class assumes a certain layout of the underlying HDF5 file. Der dataset names can be modified in hdf5*Path members.
 * Required dimensions:
 * * hdf5MeasurementPath, hdf5MeasurementSigmaPath: numObservables x numConditions
 * * hdf5ConditionPath: numFixedParameters x numConditions
 * * hdf5AmiciOptionPath:
 * * hdf5ParameterPath:
 *
 * NOTE: The following dimensions are determined by the used AMICI model:
 * * numObservables := UserData::ny
 * * numFixedParameters := UserData::nk
 *
 * The vector of optimization variables is assumed to be [x_0, ..., x_(numCommonParameter-1), conditionSpecificParameters].
 * conditionSpecificParameters := [cond0par0, cond0par1, ..., cond0_par_(numConditionSpecificParametersPerSimulation-1),
 * cond_(numConditions-1)_(numConditionSpecificParametersPerSimulation-1) ]
 */

class MultiConditionDataProvider
{
public:
    MultiConditionDataProvider(const char *hdf5Filename);

    MultiConditionDataProvider(const char *hdf5Filename, std::string rootPath);

    /**
     * @brief Provides the number of conditions for which data is available and simulations need to be run.
     * This is determined from the dimensions of the hdf5MeasurementPath dataset.
     * @return
     */
    virtual int getNumberOfConditions() const;

    /**
     * @brief Number of model- oder optimization-parameters that are different between the different conditions.
     * @return
     */

    virtual int getNumConditionSpecificParametersPerSimulation() const;

    virtual int updateFixedSimulationParameters(int conditionIdx, UserData *udata) const;

    virtual ExpData *getExperimentalDataForExperimentAndUpdateFixedParameters(int conditionIdx, UserData *udata) const;

    virtual ExpData *getExperimentalDataForCondition(int conditionIdx) const;

    /**
     * @brief getOptimizationParametersLowerBounds Get lower parameter bounds
     * NOTE: Currently the same bounds are assumed for kinetic parameters and scaling parameters, ...
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

    virtual int getNumOptimizationParameters() const;

    virtual int getNumCommonParameters() const;

    // TODO remove, since always need more info than pure model dimensions (e.g. nt)
    virtual UserData getModelDims() const;

    virtual UserData *getUserData() const;

    virtual UserData *getUserDataForCondition(int conditionIdx) const;

    virtual int getIndexOfFirstConditionSpecificOptimizationParameter(int conditionIdx) const;

    virtual void updateConditionSpecificSimulationParameters(int conditionIndex, const double *optimizationParams, UserData *udata) const;

    virtual ~MultiConditionDataProvider();

    std::string hdf5MeasurementPath;
    std::string hdf5MeasurementSigmaPath;
    std::string hdf5ConditionPath;
    std::string hdf5AmiciOptionPath;
    std::string hdf5ParameterPath;

    hid_t fileId = 0;

protected:
    MultiConditionDataProvider();
    UserData modelDims;
};

#endif // MULTICONDITIONDATAPROVIDER_H
