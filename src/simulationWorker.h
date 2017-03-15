#ifndef SIMULATION_WORKER_H
#define SIMULATION_WORKER_H

// TODO:for datapath, remove
#include "../objectiveFunctionBenchmarkModel/dataprovider.h"

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


void doWorkerWork();

void handleWorkPackage(char *buffer, int jobId);

int getLengthWorkPackageMessage(int nTheta);

int getLengthResultPackageMessage(int nTheta);

void serializeWorkPackageMessage(workPackageMessage work, int nTheta, char *buffer);
void deserializeWorkPackageMessage(const char *msg, int nTheta, datapath *path, double *theta, int *sensitivityMethod);

void serializeResultPackageMessage(resultPackageMessage result, int nTheta, char *buffer);
void deserializeResultPackageMessage(char *buffer, int nTheta, int *status, double *llh, double *sllh);

#endif
