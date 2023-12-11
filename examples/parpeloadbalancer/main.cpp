#include <parpeloadbalancer/loadBalancerMaster.h>
#include <parpeloadbalancer/loadBalancerWorker.h>

#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <mutex>
#include <condition_variable>

#include <mpi.h>

constexpr int NUM_JOBS = 1000;

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
    std::condition_variable cond;
    std::mutex mutex;

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
    std::unique_lock lock(mutex);
    cond.wait(lock, [&numJobsFinished, &numJobs]{
        return numJobsFinished == numJobs;});

    // check results
    int errors = 0;
    for (int i = 0; i < numJobs; ++i) {
        auto buffer = (double *)(jobdata[i].recvBuffer.data());

        if (*buffer != 2 * i)
            printf("ERROR: %d was %f\n", i, *buffer);
    }

    lbm.terminate();
    lbm.sendTerminationSignalToAllWorkers();

    return errors;
}

/**
 * @brief On the worker side, take the received value, multiply by 2, return
 * @param buffer
 * @param jobId
 */
void duplicatingMessageHandler(std::vector<char> &buffer, int  /*jobId*/) {
    // read message
    double value = *reinterpret_cast<double *>(buffer.data());
    //    printf("Received %f\n", value);

    // sleep(1);

    // prepare result
    buffer.resize(sizeof(double));
    auto result = reinterpret_cast<double *>(buffer.data());
    *result = value * 2;
    //    printf("Sending %f\n", *result);
}

void worker() {
    parpe::LoadBalancerWorker lbw;
    lbw.run(duplicatingMessageHandler);
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
