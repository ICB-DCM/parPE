#include "loadBalancerMaster.h"
#include "loadBalancerWorker.h"
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
    loadBalancerStartMaster();

    int numJobs = NUM_JOBS;
    int numJobsFinished = 0;

    JobData jobdata[numJobs];

    // mutex to wait for simulations to finish
    pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

    for (int i = 0; i < numJobs; ++i) {
        JobData *job = &jobdata[i];
        job->jobDone = &numJobsFinished;
        job->jobDoneChangedCondition = &cond;
        job->jobDoneChangedMutex = &mutex;
        job->lenSendBuffer = sizeof(double);
        job->sendBuffer = (char *)malloc(job->lenSendBuffer);
        *(double *)job->sendBuffer = i;
        loadBalancerQueueJob(job);
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
        free(buffer);
    }

    loadBalancerTerminate();
    sendTerminationSignalToAllWorkers();

    return errors;
}

/**
 * @brief messageHandler On the worker side, take the received value, multiply
 * by 2, return
 * @param buffer
 * @param size
 * @param jobId
 * @param userData
 */
void messageHandler(char **buffer, int *size, int jobId, void *userData) {
    // reuse allocated memory
    //    double *result = (double*) *buffer;
    //    *result *= 2;

    // read message
    double value = **((double **)buffer);
    free(*buffer);

    // sleep(1);

    // prepare result
    *size = sizeof(double);
    *buffer = (char *)malloc(*size);
    double *result = (double *)*buffer;
    *result = value * 2;
}

void worker() { loadBalancerWorkerRun(messageHandler, NULL); }

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
