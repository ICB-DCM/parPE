#ifndef MULTICONDITIONDATAPROVIDER_H
#define MULTICONDITIONDATAPROVIDER_H

#include <hdf5Misc.h>
#include <udata.h>
#include <edata.h>

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
 * @brief The MultiConditionDataProvider class reads data for MultiConditionOptimizationProblem from a HDF5 file.
 */

class MultiConditionDataProvider
{
public:
    MultiConditionDataProvider(UserData udata, const char *hdf5Filename);

    virtual int getNumberOfConditions();

    virtual int getNumConditionSpecificParametersPerSimulation();

    virtual int updateFixedSimulationParameters(int conditionIdx, UserData *udata);

    virtual ExpData *getExperimentalDataForExperimentAndUpdateUserData(int conditionIdx, UserData *udata);

    /**
     * @brief getOptimizationParametersLowerBounds Get lower parameter bounds
     * NOTE: Currently the same bounds are assumed for kinetic parameters and scaling parameters, ...
     * @param dataPath (not yet used)
     * @param buffer allocated memory to write parameter bounds
     */
    virtual void getOptimizationParametersLowerBounds(double *buffer);

    /**
     * @brief getOptimizationParametersUpperBounds Get upper parameter bounds
     * @param dataPath (not yet used)
     * @param buffer allocated memory to write parameter bounds
     */
    virtual void getOptimizationParametersUpperBounds(double *buffer);

    virtual int getNumOptimizationParameters();

    virtual int getNumCommonParameters();

    virtual UserData getModelDims();

    virtual int getIndexOfFirstConditionSpecificOptimizationParameter(int conditionIdx);

    virtual ~MultiConditionDataProvider();


protected:
    UserData modelDims;
    hid_t file_id;

};

#endif // MULTICONDITIONDATAPROVIDER_H
