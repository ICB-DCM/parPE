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
void writeObjectiveFunctionData(loggerdata logstruct, const int iter, const double obj_value, const double * const gradient, int size);
void writeOptimizationParameters(loggerdata logstruct, const int iter, const double * const theta, const int nTheta);
void writeObjectiveFunctionValue(loggerdata logstruct, const int iter, const double obj_value);
void writeObjectiveFunctionGradient(loggerdata logstruct, const int iter, const double * const gradient, int size);
void writeEvalFTime(loggerdata logstruct, const int iter, double timeElapsedInSeconds);
void flushLogger(loggerdata logstruct);

// TODO
bool hdf5DatasetExists(hid_t file_id, const char *datasetName);
bool hdf5GroupExists(hid_t file_id, const char *groupName);
void hdf5CreateGroup(hid_t file_id, const char *groupPath, bool recursively);
void hdf5CreateExtendableDouble2DArray(hid_t file_id, const char* datasetPath, int stride);
void hdf5Extend2ndDimensionAndWriteToDouble2DArray(hid_t file_id, const char *datasetPath, const double *buffer);
char *myStringCat(const char *first, const char *second);
#endif
