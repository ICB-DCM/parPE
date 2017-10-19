#ifndef LOADBALANCERMASTER_H
#define LOADBALANCERMASTER_H

#include <mpi.h>
#include <pthread.h>
#include <queue>
#include <semaphore.h>

//#define MASTER_QUEUE_H_SHOW_COMMUNICATION

namespace parPE {

/** data to be sent to workers */
struct JobData {
    JobData() = default;

    JobData(int lenSendBuffer, char *sendBuffer, int *jobDone,
            pthread_cond_t *jobDoneChangedCondition,
            pthread_mutex_t *jobDoneChangedMutex)
        : lenSendBuffer(lenSendBuffer), jobDone(jobDone),
          jobDoneChangedCondition(jobDoneChangedCondition),
          jobDoneChangedMutex(jobDoneChangedMutex) {
        this->sendBuffer = sendBuffer ? sendBuffer : new char[lenSendBuffer];
    }

    /** auto-assigned (unique number up to MAX_INT) */
    int jobId;

    /** size of data to send */
    int lenSendBuffer = 0;
    /** data to send */
    char *sendBuffer = nullptr;

    /** size of data to receive (set when job finished) */
    int lenRecvBuffer = 0;
    /** data to receive (set when job finished) */
    char *recvBuffer = nullptr;

    /** incremented by one, once the results have been received */
    int *jobDone = nullptr;

    /** is signaled after jobDone has been incremented */
    pthread_cond_t *jobDoneChangedCondition = nullptr;
    /** is locked to signal jobDoneChangedCondition condition  */
    pthread_mutex_t *jobDoneChangedMutex = nullptr;
};

class LoadBalancerMaster {
  public:
    /**
     * @brief Start the load balancer using all available MPI processes.
     * Requires MPI to be initialized.
     */
    void run();

    ~LoadBalancerMaster();

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

    bool isRunning() const;

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

    bool isRunning_ = false;

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

} // namespace parPE

#endif // LOADBALANCERMASTER_H
