#ifndef SIMULATION_WORKER_H
#define SIMULATION_WORKER_H

#include <include/rdata.h>
#include <include/udata.h>

class JobAmiciSimulation {
public:
    /** Simulation data or dataset Id */
    int lenData;
    void* data;

    /** number of simulation parameters */
    int numSimulationParameters;

    /** Simulation parameters */
    double *simulationParameters;

    int sensitivityMethod;

    void serialize(char *buffer);

    void deserialize(const char *msg);

    static int getLength(int numSimulationParameters, int sizeOfData);

    static void toUserData(const char* buffer, UserData *udata, void *userData);
};

class JobResultAmiciSimulation {
public:
    /** number of simulation parameters */
    int numSimulationParameters;

    /** simulation return status */
    int status;

    /** log likelihood */
    double llh;

    /** log likelihood gradient*/
    double *sllh;

    static void serialize(const ReturnData *rdata, const UserData *udata, int status, char *buffer);
    void deserialize(char *buffer);


    static int getLength(int numSimulationParameters);

};
#endif
