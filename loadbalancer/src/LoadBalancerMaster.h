#ifndef LOADBALANCERMASTER_H
#define LOADBALANCERMASTER_H

#include <mpi.h>

#include <pthread.h>
#include <queue>
#include <semaphore.h>
#include <functional>

//#define MASTER_QUEUE_H_SHOW_COMMUNICATION

namespace parpe {

/** data to be sent to workers */
struct JobData {
    JobData() = default;

    JobData(int *jobDone,
            pthread_cond_t *jobDoneChangedCondition,
            pthread_mutex_t *jobDoneChangedMutex)
        : jobDone(jobDone),
          jobDoneChangedCondition(jobDoneChangedCondition),
          jobDoneChangedMutex(jobDoneChangedMutex) {
    }

    /** auto-assigned (unique number up to MAX_INT) */
    int jobId;

    /** data to send */
    std::vector<char> sendBuffer;

    /** data to receive (set when job finished) */
    std::vector<char> recvBuffer;

    /** incremented by one, once the results have been received */
    int *jobDone = nullptr;

    /** is signaled after jobDone has been incremented */
    pthread_cond_t *jobDoneChangedCondition = nullptr;
    /** is locked to signal jobDoneChangedCondition condition  */
    pthread_mutex_t *jobDoneChangedMutex = nullptr;

    /** callback when job is finished */
    std::function<void(JobData*)> callbackJobFinished = nullptr;
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
    static void assertMpiActive();
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

  private:
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
    std::vector<bool> workerIsBusy;
    std::vector<MPI_Request> sendRequests;
    std::vector<JobData *> sentJobsData;

    pthread_mutex_t mutexQueue = PTHREAD_MUTEX_INITIALIZER;

    sem_t semQueue = {};

    pthread_t queueThread = 0;
};

} // namespace parpe

#endif // LOADBALANCERMASTER_H
