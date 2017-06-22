#include "optimizationApplication.h"
#include <mpi.h>
#include "logging.h"
#include "hdf5Misc.h"
#include <pthread.h>

OptimizationApplication::OptimizationApplication() : dataFileName(NULL), problem(NULL)
{

}

OptimizationApplication::OptimizationApplication(OptimizationProblem *problem, int argc, char **argv) : OptimizationApplication()
{
    this->problem = problem;

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
    if(!dataFileName) {
        logmessage(LOGLVL_CRITICAL, "No input file provided. Must provide input file as first and only argument or set OptimizationApplication::inputFileName manually.");
        return 1;
    }

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    char outfilefull[200];
    sprintf(outfilefull, resultFileName, mpiRank);
    initResultHDFFile(outfilefull, true); // option


    return 0;
}

OptimizationApplication::~OptimizationApplication()
{
    destroyHDF5Mutex();

    MPI_Finalize();
}

