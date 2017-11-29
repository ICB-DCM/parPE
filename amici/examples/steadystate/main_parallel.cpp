#include "hdf5Misc.h"
#include "steadystateProblemParallel.h"
#include <LoadBalancerMaster.h>
#include <LoadBalancerWorker.h>
#include <logging.h>
#include <mpi.h>
#include <cstdio>
#include <cstring>
#include <unistd.h>
#include <string>

/** @file
 *
 * This example demonstrates the use of the loadbalancer / queue for parallel
 * ODE simulation.
 */

void initMPI(int *argc, char ***argv);

int main(int argc, char **argv) {
    int status = 0;

    initMPI(&argc, &argv);

    if(argc != 2) {
        std::cerr<<"Error: wrong number of arguments. Exactly one argument for data file expected.";
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    std::string dataFileName = argv[1];

    parpe::initHDF5Mutex();

    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if (commSize == 1) {
        // run in serial mode
        SteadystateProblemParallel problem {nullptr, dataFileName};
        status = getLocalOptimum(&problem);

    } else {
        // run in parallel
        SteadystateProblemParallel problem {nullptr, dataFileName};

        int mpiRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        if (mpiRank == 0) {
            parpe::LoadBalancerMaster lbm;
            problem.loadBalancer = &lbm;
            lbm.run();

            status = parpe::getLocalOptimum(&problem);

            lbm.terminate();
            lbm.sendTerminationSignalToAllWorkers();
        } else {
            parpe::LoadBalancerWorker lbw;
            lbw.run([&problem](std::vector<char> &buffer, int jobId) { problem.messageHandler(buffer, jobId); });
        }
    }

    MPI_Finalize();

    return status;
}

void initMPI(int *argc, char ***argv) {
    int mpiErr = MPI_Init(argc, argv);
    if (mpiErr != MPI_SUCCESS) {
        parpe::logmessage(parpe::LOGLVL_CRITICAL, "Problem initializing MPI. Exiting.");
        exit(1);
    }

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    if (mpiRank == 0) {
        int commSize;
        MPI_Comm_size(MPI_COMM_WORLD, &commSize);

        parpe::logmessage(parpe::LOGLVL_INFO, "Running with %d MPI processes.", commSize);
    }
}
