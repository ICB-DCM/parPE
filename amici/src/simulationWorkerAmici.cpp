#include "simulationWorkerAmici.h"
#include <string.h>

void JobAmiciSimulation::toUserData(const char* buffer, UserData *udata, void *userData) {
    JobAmiciSimulation work;
    work.data = userData;
    work.simulationParameters = udata->p;
    work.deserialize(buffer);

    udata->sensi_meth = (AMI_sensi_meth) work.sensitivityMethod;
    udata->sensi = work.sensitivityMethod > 0 ? AMI_SENSI_ORDER_FIRST : AMI_SENSI_ORDER_NONE;
}

int JobAmiciSimulation::getLength(int numSimulationParameters, int sizeOfData)
{
    return sizeof(int) // user data size
            + sizeOfData // userdata
            + sizeof(int) // numTheta
            + sizeof(double) * numSimulationParameters // theta
            + sizeof(int); // sensi
}

void JobAmiciSimulation::serialize(char *buffer)
{
    size_t size = 0;

    size = sizeof(lenData);
    memcpy(buffer, &lenData, size);
    buffer += size;

    size = lenData;
    memcpy(buffer, data, size);
    buffer += size;

    size = sizeof(numSimulationParameters);
    memcpy(buffer, &numSimulationParameters, size);
    buffer += size;

    size = numSimulationParameters * sizeof(double);
    memcpy(buffer, simulationParameters, size);
    buffer += size;

    size = sizeof(int);
    memcpy(buffer, &sensitivityMethod, size);
    buffer += size;

}

void JobAmiciSimulation::deserialize(const char *msg)
{
    size_t size;

    size = sizeof(lenData);
    lenData = *(int *) msg;
    msg += size;

    size = lenData;
    memcpy(data, msg, size);
    msg += size;

    size = sizeof(numSimulationParameters);
    numSimulationParameters = *(int *) msg;
    msg += size;

    size = numSimulationParameters * sizeof(double);
    memcpy(simulationParameters, msg, size);
    msg += size;

    sensitivityMethod = *(int *) msg;
    size = sizeof(int);
    msg += size;
}

void JobResultAmiciSimulation::serialize(const ReturnData *rdata, const UserData *udata, int status, char *buffer)
{
    size_t size = 0;

    size = sizeof(int);
    memcpy(buffer, &udata->np, size);
    buffer += size;

    size = sizeof(int);
    memcpy(buffer, &status, size);
    buffer += size;

    size = sizeof(double);
    memcpy(buffer, &rdata->llh[0], size);
    buffer += size;

    size = sizeof(int);
    int sensiSize = udata->sensi_meth > 0 ? udata->np : 0;
    memcpy(buffer, &sensiSize, size);
    buffer += size;

    if(sensiSize) {
        size = sensiSize * sizeof(double);
        memcpy(buffer, rdata->sllh, size);
        buffer += size;
    }
}

void JobResultAmiciSimulation::deserialize(char *buffer)
{
    size_t size;

    numSimulationParameters = *(int *) buffer;
    size = sizeof(int);
    buffer += size;

    status = *(int *) buffer;
    size = sizeof(int);
    buffer += size;

    llh = *(double *) buffer;
    size = sizeof(double);
    buffer += size;

    int sensiSize = *(int *) buffer;
    size = sizeof(int);
    buffer += size;

    if(sensiSize) {
        size = sensiSize * sizeof(double);
        memcpy(sllh, buffer, size);
        buffer += size;
    }
}

int JobResultAmiciSimulation::getLength(int numSimulationParameters)
{
    return sizeof(int) // num Theta
            + sizeof(int) // status
            + sizeof(double) // llh
            + sizeof(int) // lenSensi
            + sizeof(double) * numSimulationParameters; // sensitivites
}
