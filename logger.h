#ifndef LOGGER_H
#define LOGGER_H

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdbool.h>

#define LOGGER_H_LOGFILE "mylog.log"

typedef struct {
    hid_t file_id;
    char *rootPath;

} loggerdata;


loggerdata initResultHDFFile(const char *filename, const char *rootPath);
void closeResultHDFFile(loggerdata logstruct);
void logLocalOptimizerObjectiveFunctionEvaluation(loggerdata logstruct, int numFunctionCalls, double *theta, double objectiveFunctionValue, const double *gradient, double timeElapsedInSeconds, int nTheta);
void logLocalOptimizerIteration(loggerdata logstruct, int numIterations, double *theta, double objectiveFunctionValue, const double *gradient, double timeElapsedInSeconds, int nTheta);
void flushLogger(loggerdata logstruct);

// TODO
bool hdf5DatasetExists(hid_t file_id, const char *datasetName);
bool hdf5GroupExists(hid_t file_id, const char *groupName);
void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively);
void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char* datasetPath, int stride);
void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id, const char *datasetPath, const double *buffer);
void hdf5CreateOrExtendAndWriteToDouble2DArray(hid_t file_id, const char *parentPath, const char *datasetName, const double *buffer, int stride);
char *myStringCat(const char *first, const char *second);
#endif
