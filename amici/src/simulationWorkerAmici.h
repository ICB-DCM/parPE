#ifndef SIMULATION_WORKER_H
#define SIMULATION_WORKER_H

#include <include/rdata.h>
#include <include/udata.h>
#include <amici_serialization.h>
#include <misc.h>

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

    static void toUserData(const char *buffer, amici::UserData *udata, void *userData);
};

class JobResultAmiciSimulation {
public:
    JobResultAmiciSimulation() = default;
    JobResultAmiciSimulation(int status, std::unique_ptr<amici::ReturnData> rdata, double simulationTimeInSec)
        : status(status), rdata(std::move(rdata)), simulationTimeInSec(simulationTimeInSec) {}

    ~JobResultAmiciSimulation() {
    }

    JobResultAmiciSimulation ( JobResultAmiciSimulation && other) {
        std::swap(status, other.status);
        std::swap(rdata, other.rdata);
        std::swap(simulationTimeInSec, other.simulationTimeInSec);
    }

    /** simulation return status */
    int status = 1;

    /** log likelihood */
    std::unique_ptr<amici::ReturnData> rdata;

    double simulationTimeInSec = -1;
};

} // namespace parpe

namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, parpe::JobResultAmiciSimulation &d, const unsigned int version) {
    ar & d.status;
    if (Archive::is_loading::value) {
        d.rdata = std::make_unique<amici::ReturnData>();
    }
    ar & *d.rdata.get();
    ar & d.simulationTimeInSec;
}
} // namespace serialization
} // namespace boost

#endif

