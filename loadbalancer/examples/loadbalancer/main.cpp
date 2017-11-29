#include "LoadBalancerMaster.h"
#include "LoadBalancerWorker.h"
#include <cstdlib>
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

#define NUM_JOBS 1000

/*
 * Testing code for MPI load balancing.
 *
 */

/**
 * @brief master send a double to any of the workers, wait for completion,
 * verify result
 * @return number of errors
 */
int master() {
    parpe::LoadBalancerMaster lbm;
    lbm.run();

    int numJobs = NUM_JOBS;
    int numJobsFinished = 0;

    parpe::JobData jobdata[numJobs];

    // mutex to wait for simulations to finish
    pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

    for (int i = 0; i < numJobs; ++i) {
        parpe::JobData *job = &jobdata[i];
        job->jobDone = &numJobsFinished;
        job->jobDoneChangedCondition = &cond;
        job->jobDoneChangedMutex = &mutex;
        job->sendBuffer.resize(sizeof(double));
        *(double *)job->sendBuffer.data() = i;
        lbm.queueJob(job);
    }

    // wait for simulations to finish
    pthread_mutex_lock(&mutex);
    while (numJobsFinished < numJobs)
        pthread_cond_wait(&cond, &mutex);
    pthread_mutex_unlock(&mutex);
    pthread_mutex_destroy(&mutex);
    pthread_cond_destroy(&cond);

    // check results
    int errors = 0;
    for (int i = 0; i < numJobs; ++i) {
        double *buffer = (double *)(jobdata[i].recvBuffer);

        if (*buffer != 2 * i)
            printf("ERROR: %d was %f\n", i, *buffer);
        delete[] buffer;
    }

    lbm.terminate();
    lbm.sendTerminationSignalToAllWorkers();

    return errors;
}

class DuplicatingLoadBalancerWorker : public parpe::LoadBalancerWorker {
    /**
     * @brief messageHandler On the worker side, take the received value,
     * multiply by 2, return
     * @param buffer
     * @param size
     * @param jobId
     * @param userData
     */

    void messageHandler(std::vector<char> &buffer, int jobId) override {
        // read message
        double value = *reinterpret_cast<double *>(buffer.data());
        //    printf("Received %f\n", value);

        // sleep(1);

        // prepare result
        buffer.resize(sizeof(double));
        double *result = reinterpret_cast<double *>(buffer.data());
        *result = value * 2;
        //    printf("Sending %f\n", *result);
    }
};

void worker() {
    DuplicatingLoadBalancerWorker w;
    w.run();
}

int main(int argc, char **argv) {
    int status = 0;

    MPI_Init(&argc, &argv);

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if (mpiRank == 0) {
        status = master();
    } else {
        worker();
    }

    MPI_Finalize();

    return status;
}
