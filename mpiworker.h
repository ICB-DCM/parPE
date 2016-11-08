#ifndef MPI_WORKER_H
#define MPI_WORKER_H

#include "dataprovider.h"

#define MPI_TAG_EXIT_SIGNAL 0
#define MPI_WORKER_H_VERBOSE 0

typedef struct workPackageMessage_tag {
    datapath path;
    double *theta;
    int sensitivityMethod;
} workPackageMessage;

typedef struct resultPackageMessage_tag {
    int status;
    double llh;
    double *sllh;
} resultPackageMessage;

int getLengthWorkPackageMessage(int nTheta);

int getLengthResultPackageMessage(int nTheta);

void serializeWorkPackageMessage(workPackageMessage work, int nTheta, char *buffer);
void deserializeWorkPackageMessage(const char *msg, int nTheta, datapath *path, double *theta, int *sensitivityMethod);

void serializeResultPackageMessage(resultPackageMessage result, int nTheta, char *buffer);
void deserializeResultPackageMessage(char *buffer, int nTheta, int *status, double *llh, double *sllh);

void doWorkerWork();

resultPackageMessage handleWorkPackage(const char *buffer, UserData *udata, ReturnData **prdata, int tag);

void sendTerminationSignalToAllWorkers();

#endif
