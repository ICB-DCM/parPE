#include "logger.h"
#include "src/ami_hdf5.h"
#include <string.h>
#include <assert.h>

hid_t file_id = 0;

loggerdata initResultHDFFile(const char *filename, const char *rootPath)
{
    loggerdata logstruct;
    //logstruct.file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    logstruct.file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); // overwrites TODO: backup
    logstruct.rootPath = malloc(sizeof(char) * (strlen(rootPath) + 1));
    strcpy(logstruct.rootPath, rootPath);


    // create data structures ... extendable. ..
    // create /multistarts/n/objectiveFunctionValue/

    printf("Init logging to '%s'\n", filename); fflush(stdout);

    return logstruct;
}


void closeResultHDFFile(loggerdata logstruct)
{
    free(logstruct.rootPath);
    herr_t status = H5Fclose(logstruct.file_id);

}



void flushLogger(loggerdata logstruct)
{
    herr_t status = H5Fflush(logstruct.file_id, H5F_SCOPE_LOCAL);
}

bool hdf5DatasetExists(hid_t file_id, const char *datasetName)
{
    return H5Lexists(file_id, datasetName, H5P_DEFAULT) > 0;
}

bool hdf5GroupExists(hid_t file_id, const char *groupName)
{
    // switch off error handler, check existance and reenable
    herr_t (*old_func)(void*);
    void *old_client_data;
    H5Eget_auto1(&old_func, &old_client_data);
    H5Eset_auto1(NULL, NULL);

    herr_t status = H5Gget_objinfo(file_id, groupName, 0, NULL);

    H5Eset_auto1(old_func, old_client_data);

    return status >= 0;
}

void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively)
{
    hid_t groupCreationPropertyList = H5P_DEFAULT;

    if(recursively) {
        groupCreationPropertyList = H5Pcreate(H5P_LINK_CREATE);
        herr_t status = H5Pset_create_intermediate_group (groupCreationPropertyList, 1);
    }

    hid_t group = H5Gcreate (file_id, groupPath, groupCreationPropertyList, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group);
}

void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char *datasetPath, int stride)
{
    int rank = 2;
    hsize_t initialDimensions[2]  = {stride, 0};
    hsize_t maximumDimensions[2]  = {stride, H5S_UNLIMITED};

    hid_t dataspace = H5Screate_simple(rank, initialDimensions, maximumDimensions);

    // need chunking for extendable dataset
    hsize_t chunkDimensions[2] = {stride, 1};
    hid_t datasetCreationProperty = H5Pcreate(H5P_DATASET_CREATE);
    herr_t status = H5Pset_chunk(datasetCreationProperty, rank, chunkDimensions);

    hid_t dataset = H5Dcreate2(file_id, datasetPath, H5T_NATIVE_DOUBLE, dataspace,
                               H5P_DEFAULT, datasetCreationProperty, H5P_DEFAULT);

    status = H5Dclose(dataset);
    status = H5Sclose(dataspace);
}

void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id, const char *datasetPath, const double *buffer)
{
    hid_t dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);

    // extend
    hid_t filespace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(filespace);
    assert(rank == 2);

    hsize_t currentDimensions[2];
    int status_n = H5Sget_simple_extent_dims(filespace, currentDimensions, NULL);

    hsize_t newDimensions[2]  = {currentDimensions[0], currentDimensions[1] + 1};
    herr_t status = H5Dset_extent(dataset, newDimensions);

    filespace = H5Dget_space(dataset);
    hsize_t offset[2] = {0, currentDimensions[1]};
    hsize_t slabsize[2] = {currentDimensions[0], 1};

    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, slabsize, NULL);

    hid_t memspace = H5Screate_simple(rank, slabsize, NULL);

    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, buffer);

    status = H5Dclose(dataset);
    status = H5Sclose(filespace);
    status = H5Sclose(memspace);
}

char *myStringCat(const char *first, const char *second)
{
    char *concatenation = malloc(sizeof(char) * (strlen(first) + strlen(second) + 1));
    strcpy(concatenation, first);
    strcat(concatenation, second);
    return concatenation;
}

void logLocalOptimizerObjectiveFunctionEvaluation(loggerdata logstruct, int numFunctionCalls, double *theta, double objectiveFunctionValue, const double *gradient, double timeElapsedInSeconds, int nTheta)
{
    // TODO need to know whoch multistart... need to provide logging struct? add to ipopt uzserdata
    char *fullGroupPath = myStringCat(logstruct.rootPath, "/trajectory/optimizerObjFunEvaluations/");

    hdf5CreateOrExtendAndWriteToDouble2DArray(logstruct.file_id, fullGroupPath, "negLogLikelihood", &objectiveFunctionValue, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(logstruct.file_id, fullGroupPath, "negLogLikelihoodGradient", gradient, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(logstruct.file_id, fullGroupPath, "p", theta, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(logstruct.file_id, fullGroupPath, "evalFTime", &timeElapsedInSeconds, 1);

    free(fullGroupPath);

    flushLogger(logstruct);
}

void logLocalOptimizerIteration(loggerdata logstruct, int numIterations, double *theta, double objectiveFunctionValue, const double *gradient, double timeElapsedInSeconds, int nTheta)
{
    // TODO need to know whoch multistart... need to provide logging struct? add to ipopt uzserdata
    char *fullGroupPath = myStringCat(logstruct.rootPath, "/trajectory/optimizerIterations/");

    hdf5CreateOrExtendAndWriteToDouble2DArray(logstruct.file_id, fullGroupPath, "negLogLikelihood", &objectiveFunctionValue, 1);
    hdf5CreateOrExtendAndWriteToDouble2DArray(logstruct.file_id, fullGroupPath, "negLogLikelihoodGradient", gradient, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(logstruct.file_id, fullGroupPath, "p", theta, nTheta);
    hdf5CreateOrExtendAndWriteToDouble2DArray(logstruct.file_id, fullGroupPath, "evalFTime", &timeElapsedInSeconds, 1);

    free(fullGroupPath);

    flushLogger(logstruct);
}

void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride)
{
    if(!hdf5GroupExists(file_id, parentPath)) {
        hdf5CreateGroup(file_id, parentPath, true);
    }

    char *fullDatasetPath = myStringCat(parentPath, datasetName);

    if(!hdf5DatasetExists(file_id, fullDatasetPath)) {
        hdf5CreateExtendableDouble2DArray(file_id, fullDatasetPath, stride);
    }

    hdf5Extend2ndDimensionAndWriteToDouble2DArray(file_id, fullDatasetPath, buffer);

    free(fullDatasetPath);
}
