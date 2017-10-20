#include "simulationWorkerAmici.h"
#include <cstring>

namespace parpe {

void JobAmiciSimulation::toUserData(const char *buffer, amici::UserData *udata,
                                    void *userData) {
    JobAmiciSimulation work;
    work.data = userData;
    work.simulationParameters = udata->p;
    work.deserialize(buffer);

    udata->sensi_meth = (amici::AMICI_sensi_meth)work.sensitivityMethod;
    udata->sensi = work.sensitivityMethod > 0 ? amici::AMICI_SENSI_ORDER_FIRST
                                              : amici::AMICI_SENSI_ORDER_NONE;
}

int JobAmiciSimulation::getLength(int numSimulationParameters, int sizeOfData) {
    return sizeof(int)                                // user data size
           + sizeOfData                               // userdata
           + sizeof(int)                              // numTheta
           + sizeof(double) * numSimulationParameters // theta
           + sizeof(int);                             // sensi
}

void JobAmiciSimulation::serialize(char *buffer) {
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

void JobAmiciSimulation::deserialize(const char *msg) {
    size_t size;

    size = sizeof(lenData);
    lenData = *(int *)msg;
    msg += size;

    size = lenData;
    memcpy(data, msg, size);
    msg += size;

    size = sizeof(numSimulationParameters);
    numSimulationParameters = *(int *)msg;
    msg += size;

    size = numSimulationParameters * sizeof(double);
    memcpy(simulationParameters, msg, size);
    msg += size;

    sensitivityMethod = *(int *)msg;
    size = sizeof(int);
    msg += size;
}

} // namespace parpe
