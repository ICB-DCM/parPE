#include "optimizationApplication.h"
#include "logging.h"
#include "hdf5Misc.h"
#include "loadBalancerMaster.h"
#include <pthread.h>
#include <mpi.h>

OptimizationApplication::OptimizationApplication() : dataFileName(NULL), problem(NULL), resultWriter(NULL)
{

}

OptimizationApplication::OptimizationApplication(int argc, char **argv) : OptimizationApplication()
{
    // TODO: check if initialized already
    initMPI(&argc, &argv);

    initHDF5Mutex();

    if(argc == 2) {
        dataFileName = argv[1];
    }
}

void OptimizationApplication::initMPI(int *argc, char ***argv)
{
    int mpiErr = MPI_Init(argc, argv);
    if(mpiErr != MPI_SUCCESS) {
        logmessage(LOGLVL_CRITICAL, "Problem initializing MPI. Exiting.");
        exit(1);
    }

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    if(mpiRank == 0) {
        int commSize;
        MPI_Comm_size(MPI_COMM_WORLD, &commSize);

        logmessage(LOGLVL_INFO, "Running with %d MPI processes.", commSize);
    }
}


int OptimizationApplication::run()
{
    int status = 0;

    if(!dataFileName) {
        logmessage(LOGLVL_CRITICAL, "No input file provided. Must provide input file as first and only argument or set OptimizationApplication::inputFileName manually.");
        return 1;
    }
    initProblem(dataFileName, NULL); // TODO second argument

    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if(getMpiRank() == 0) {
        if(commSize > 1)
            loadBalancerStartMaster();

        status = runMaster();

        if(commSize > 1) {
            loadBalancerTerminate();
            sendTerminationSignalToAllWorkers();
        }

    } else {
        runWorker();
    }

    return status;
}

OptimizationApplication::~OptimizationApplication()
{
    destroyProblem();

    destroyHDF5Mutex();

    MPI_Finalize();
}

int OptimizationApplication::getMpiRank()
{
    int mpiRank = -1;

    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);

    if(mpiInitialized) {
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    }

    return mpiRank;
}

int OptimizationApplication::getMpiCommSize()
{
    int mpiCommSize = -1;

    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);

    if(mpiInitialized) {
        MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
    }

    return mpiCommSize;
}

