#ifndef SIMULATION_WORKER_H
#define SIMULATION_WORKER_H

#include <include/rdata.h>
#include <include/udata.h>
#include <amici_serialization.h>

namespace parpe {

struct JobAmiciSimulation {
  public:
    /** Simulation data or dataset Id */
    int lenData;
    void *data;

    /** number of simulation parameters */
    int numSimulationParameters;

    /** Simulation parameters */
    double *simulationParameters;

    int sensitivityMethod;

    void serialize(char *buffer);

    void deserialize(const char *msg);

    static int getLength(int numSimulationParameters, int sizeOfData);

    static void toUserData(const char *buffer, UserData *udata, void *userData);
};

class JobResultAmiciSimulation {
public:
    JobResultAmiciSimulation() = default;
    JobResultAmiciSimulation(int status, ReturnData *rdata, double simulationTimeInSec)
        : status(status), rdata(rdata), simulationTimeInSec(simulationTimeInSec) {}

    ~JobResultAmiciSimulation() {
    }

    /** simulation return status */
    int status = 1;

    /** log likelihood */
    ReturnData *rdata = nullptr;

    double simulationTimeInSec = -1;
};

} // namespace parpe

namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, parpe::JobResultAmiciSimulation &d, const unsigned int version) {
    ar & d.status;
    if (Archive::is_loading::value) {
        d.rdata = new ReturnData();
    }
    ar & *d.rdata;
    ar & d.simulationTimeInSec;
}
} // namespace serialization
} // namespace boost

#endif

