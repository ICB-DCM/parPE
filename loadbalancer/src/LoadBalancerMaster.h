#ifndef LOADBALANCERMASTER_H
#define LOADBALANCERMASTER_H

#include <mpi.h>
#include <pthread.h>
#include <queue>
#include <semaphore.h>

//#define MASTER_QUEUE_H_SHOW_COMMUNICATION

/** data to be sent to workers */
typedef struct JobData_tag {
    /** auto-assigned (unique number up to MAX_INT) */
    int jobId;

    /** size of data to send */
    int lenSendBuffer;
    /** data to send */
    char *sendBuffer;

    /** size of data to receive (set when job finished) */
    int lenRecvBuffer;
    /** data to receive (set when job finished) */
    char *recvBuffer;

    /** incremented by one, once the results have been received */
    int *jobDone;

    /** is signaled after jobDone has been incremented */
    pthread_cond_t *jobDoneChangedCondition;
    /** is locked to signal jobDoneChangedCondition condition  */
    pthread_mutex_t *jobDoneChangedMutex;
} JobData;

class LoadBalancerMaster {
  public:
    /**
     * @brief Start the load balancer using all available MPI processes.
     * Requires MPI to be initialized.
     */
    void run();

#ifndef QUEUE_MASTER_TEST
    static void assertMPIInitialized();
#endif

    /**
     * @brief Assign job ID and append to queue for sending to workers.
     * @param data
     */

    void queueJob(JobData *data);

    /**
     * @brief Stop the loadbalancer thread
     */
    void terminate();

    void sendTerminationSignalToAllWorkers();

    JobData initJobData(int lenSendBuffer, char *sendBuffer, int *jobDone,
                        pthread_cond_t *jobDoneChangedCondition,
                        pthread_mutex_t *jobDoneChangedMutex);

  protected:
    /**
     * @brief Thread entry point. This is run from run()
     * @param "this"
     * @return 0, always
     */
    static void *threadEntryPoint(void *vpLoadBalancerMaster);

    void *loadBalancerThreadRun();

    void freeEmptiedSendBuffers();

    int handleFinishedJobs();

    int getNextFreeWorkerIndex();

    /**
     * @brief getNextJob Pop oldest element from the queue and return.
     * @return The first queue element.
     */
    JobData *getNextJob();

    /**
     * @brief sendToWorker Send the given work package to the given worker and
     * track requests
     * @param workerIdx
     * @param queueElement
     */
    void sendToWorker(int workerIdx, JobData *data);

    /**
     * @brief receiveFinished Message received from worker, mark job as done.
     * @param workerID
     * @param jobID
     */
    int handleReply(MPI_Status *mpiStatus);

    bool isRunning = false;

    int numWorkers = 0;
    std::queue<JobData *> queue;
    int lastJobId = 0;
    // one for each worker, index is off by one from MPI rank
    // because no job is sent to master (rank 0)
    bool *workerIsBusy = nullptr;
    MPI_Request *sendRequests = nullptr; // TODO: option: free(sendbuffer)
    JobData **sentJobsData = nullptr;
    pthread_mutex_t mutexQueue = PTHREAD_MUTEX_INITIALIZER;
    sem_t semQueue = {};
    pthread_t queueThread = 0;
};

#endif // LOADBALANCERMASTER_H
