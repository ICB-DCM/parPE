#ifndef MPI_WORKER_H
#define MPI_WORKER_H

#include "include/amici.h"

typedef struct {
    int packageId;
    UserData udata;
    ExpData edata;
} workpackage;

void doWorkerWork();
workpackage receiveData();
void reportToMaster(ReturnData rdata);

void sendTerminationSignalToAllWorkers();


#endif
