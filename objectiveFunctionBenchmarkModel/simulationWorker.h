#ifndef SIMULATION_WORKER_H
#define SIMULATION_WORKER_H

#include <include/udata.h>

typedef struct workPackageMessage_tag {
    /** Simulation data or dataset Id */
    int lenData;
    void* data;
    /** Simulation parameters */
    double *theta;
    int sensitivityMethod;
} workPackageMessage;

typedef struct resultPackageMessage_tag {
    /** simulation return status */
    int status;
    /** log likelihood */
    double llh;
    /** log likelihood gradient*/
    double *sllh;
} resultPackageMessage;


void doWorkerWork(UserData *_udata);

void handleWorkPackage(char *buffer, int jobId);

int getLengthWorkPackageMessage(int nTheta);

int getLengthResultPackageMessage(int nTheta);

void serializeWorkPackageMessage(workPackageMessage work, int nTheta, char *buffer);
void deserializeWorkPackageMessage(const char *msg, int nTheta, void *data, double *theta, int *sensitivityMethod);

void serializeResultPackageMessage(resultPackageMessage result, int nTheta, char *buffer);
void deserializeResultPackageMessage(char *buffer, int nTheta, int *status, double *llh, double *sllh);

#endif
