#ifndef SIMULATION_WORKER_H
#define SIMULATION_WORKER_H

#include <include/rdata.h>
#include <include/udata.h>
#include <amici_serialization.h>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <misc.h>
#include <udata.h>

namespace parpe {

/**
 * @brief The JobAmiciSimulation struct contains data to be sent to a worker to run a simulation
 */
template<typename USERDATA>
struct JobAmiciSimulation {
  public:
    JobAmiciSimulation(amici::UserData *udata, USERDATA *data)
        : udata(udata), data(data)
    {}
    amici::UserData *udata = nullptr;

    /** Simulation data or dataset Id */
    USERDATA *data = nullptr;

    std::vector<char> serialize() {
        std::string serialized;
        ::boost::iostreams::back_insert_device<std::string> inserter(serialized);
        ::boost::iostreams::stream<::boost::iostreams::back_insert_device<std::string>>
            s(inserter);
        ::boost::archive::binary_oarchive oar(s);
        oar << *data
            << udata->np
            << boost::serialization::make_array<double>(udata->p, udata->np)
            << udata->sensi_meth;

        s.flush();

        std::vector<char> buf(serialized.begin(), serialized.end());

        return buf;
    }

    void deserialize(const char* buffer, int size) {
        ::boost::iostreams::basic_array_source<char> device(buffer, size);
        ::boost::iostreams::stream<::boost::iostreams::basic_array_source<char>> s(
            device);
        ::boost::archive::binary_iarchive iar(s);
        iar >> *data;
        int numSimulationParameters;
        iar >> numSimulationParameters;
        assert(numSimulationParameters == udata->np);
        iar >> boost::serialization::make_array<double>(udata->p, numSimulationParameters);
        int sensitivityMethod;
        iar >> sensitivityMethod;
        udata->sensi_meth = (amici::AMICI_sensi_meth) sensitivityMethod;
        udata->sensi = sensitivityMethod > 0 ? amici::AMICI_SENSI_ORDER_FIRST
                                             : amici::AMICI_SENSI_ORDER_NONE;

    }
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

/**
 * For use with boost::serialization. Serialize raw C++ array.
 */
template <class Archive, typename T>
void archiveRawArray(Archive &ar, T **p, int &numElements) {
    if (Archive::is_loading::value) {
        ar &numElements;
        *p = numElements ? new T[numElements] : nullptr;
    } else {
        numElements = *p == nullptr ? 0 : numElements;
        ar &numElements;
    }
    ar &boost::serialization::make_array<T>(*p, numElements);
}

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

