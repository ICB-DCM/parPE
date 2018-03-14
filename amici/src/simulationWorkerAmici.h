#ifndef SIMULATION_WORKER_H
#define SIMULATION_WORKER_H

#include <misc.h>

#include <amici/rdata.h>
#include <amici/serialization.h>

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

namespace parpe {

/**
 * @brief The JobAmiciSimulation struct contains data to be sent to a worker to run a simulation
 */
template<typename USERDATA>
struct JobAmiciSimulation {
  public:
    JobAmiciSimulation(amici::Solver *solver, amici::Model *model, USERDATA *data)
        : solver(solver), model(model), data(data)
    {}
    amici::Solver *solver = nullptr;
    amici::Model *model = nullptr;

    /** Simulation data or dataset Id */
    USERDATA *data = nullptr;

    std::vector<char> serialize() {
        std::string serialized;
        ::boost::iostreams::back_insert_device<std::string> inserter(serialized);
        ::boost::iostreams::stream<::boost::iostreams::back_insert_device<std::string>>
            s(inserter);
        ::boost::archive::binary_oarchive oar(s);
        // older version boost::serialization can only serialize lvalues?!
        auto p = model->getParameters();
        auto sens = solver->getSensitivityMethod();
        oar << *data
            << p
            << sens;

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

        std::vector<double> parameters;
        iar >> parameters;
        model->setParameters(parameters);

        int sensitivityMethod;
        iar >> sensitivityMethod;
        solver->setSensitivityMethod(static_cast<amici::AMICI_sensi_meth>(sensitivityMethod));
        solver->setSensitivityOrder(sensitivityMethod > 0 ?
                                        amici::AMICI_SENSI_ORDER_FIRST
                                      : amici::AMICI_SENSI_ORDER_NONE);

    }
};

class JobResultAmiciSimulation {
public:
    JobResultAmiciSimulation() = default;
    JobResultAmiciSimulation(int status, std::unique_ptr<amici::ReturnData> rdata, double simulationTimeInSec)
        : status(status), rdata(std::move(rdata)), simulationTimeInSec(simulationTimeInSec) {}

    ~JobResultAmiciSimulation() = default;

    friend void swap(JobResultAmiciSimulation& first, JobResultAmiciSimulation& second) {
        using std::swap;
        swap(first.status, second.status);
        swap(first.simulationTimeInSec, second.simulationTimeInSec);
        swap(first.rdata, second.rdata);
    }

    JobResultAmiciSimulation ( JobResultAmiciSimulation && other) {
        swap(*this, other);
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

