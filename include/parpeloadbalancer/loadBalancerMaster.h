#ifndef LOADBALANCERMASTER_H
#define LOADBALANCERMASTER_H

#include <parpecommon/parpeConfig.h>

#include <pthread.h>
#include <queue>
#include <semaphore.h>
#include <functional>

#ifdef PARPE_ENABLE_MPI
#include <mpi.h>
#endif

namespace parpe {

/** Data to be sent to workers */
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
    int jobId = -1;

    /** data to send */
    std::vector<char> sendBuffer;

    /** data to receive (set when job finished) */
    std::vector<char> recvBuffer;

    /** incremented by one, once the results have been received (if set) */
    int *jobDone = nullptr;

    /** is signaled after jobDone has been incremented (if set) */
    pthread_cond_t *jobDoneChangedCondition = nullptr;
    /** is locked to signal jobDoneChangedCondition condition (if set) */
    pthread_mutex_t *jobDoneChangedMutex = nullptr;

    /** callback when job is finished (if set) */
    std::function<void(JobData*)> callbackJobFinished = nullptr;
};


#ifdef PARPE_ENABLE_MPI
/**
 * @brief The LoadBalancerMaster class sends jobs to workers, receives the
 * results and signals the client.
 */
class LoadBalancerMaster {
  public:
    LoadBalancerMaster() = default;

    LoadBalancerMaster(LoadBalancerMaster& other) = delete;

    LoadBalancerMaster& operator=(const LoadBalancerMaster& other) = delete;

    LoadBalancerMaster(LoadBalancerMaster &&other) noexcept = delete;

    LoadBalancerMaster const & operator=(LoadBalancerMaster &&fp) = delete;

    ~LoadBalancerMaster();

    /**
     * @brief Start the load balancer using all available MPI processes.
     * Requires MPI to be initialized.
     */
    void run();

#ifndef QUEUE_MASTER_TEST
    static void assertMpiActive();
#endif

    /**
     * @brief Assign job ID and append to queue for sending to workers.
     * @param data Data to be sent (user keeps ownership).
     */
    void queueJob(JobData *data);

    /**
     * @brief Stop the loadbalancer thread
     */
    void terminate();

    /**
     * @brief Send termination signal to all workers and wait for receive.
     */
    void sendTerminationSignalToAllWorkers();

    /**
     * @brief Returns whether we are ready to accept jobs (`run` was called, but
     * `terminate` was not).
     * @return true if running, false otherwise
     */
    bool isRunning() const;

    /**
     * @brief Get number of jobs in queue
     * @return Number of jobs in queue waiting to be sent
     */
    int getNumQueuedJobs() const;

  private:
    /**
     * @brief Thread entry point. This is run from run()
     * @param `this`
     * @return nullptr, always
     */
    static void *threadEntryPoint(void *vpLoadBalancerMaster);

    /**
     * @brief Main function of the load balancer thread.
     *
     * Called from threadEntryPoint. Will never return.
     */
    void loadBalancerThreadRun();

    /**
     * @brief Frees all send buffers after respective MPI messages have been
     * sent
     */
    void freeEmptiedSendBuffers();

    /**
     * @brief Check for finished jobs, receive their results and send next job
     * if jobs are waiting.
     * @return Index (not rank) of worker that has finished and not yet received
     * new work or NO_FREE_WORKER if no such worker.
     */
    int handleFinishedJobs();

    /**
     * @brief Get index (not rank) of next free worker.
     * @return That index or NO_FREE_WORKER if no such worker
     */
    int getNextFreeWorkerIndex();

    /**
     * @brief Pop oldest element from the queue and return.
     * @return The first queue element.
     */
    JobData *getNextJob();

    /**
     * @brief Send the given work package to the given worker and track
     * requests
     * @param workerIdx Index (not rank)
     * @param data Job data to send
     */
    void sendToWorker(int workerIdx, JobData *data);

    /**
     * @brief Handle the result message from a worker as indicated by mpiStatus.
     *
     * @param mpiStatus Receive the indicated message, mark job as done,
     * signal reception.
     */
    int handleReply(MPI_Status *mpiStatus);

    /**
     * @brief Check if jobs are waiting in queue and send to specified worker.
     *
     * Non-blocking.
     *
     * @param freeWorkerIndex Worker index to send job to. Assumed to be free.
     * @return Whether a job has been sent
     */
    bool sendQueuedJob(int freeWorkerIndex);

    /** MPI communicator we are working on */
    MPI_Comm mpiComm = MPI_COMM_WORLD;

    /** MPI data type for job and result packages */
    MPI_Datatype mpiJobDataType = MPI_BYTE;

    /** Indicates whether we are ready to handle jobs */
    bool isRunning_ = false;

    /** Number of workers we can send jobs to */
    int numWorkers = 0;

    /** Queue with jobs to be sent to workers */
    std::queue<JobData *> queue;

    /** Last assigned job ID used as MPI message tag */
    int lastJobId = 0;

    /** Keeps track of whether the respective worker is busy (received a job,
     * but has not sent the result) or is ready to accept a new job.
     *
     * Length is `numWorkers`. Index is off by one from MPI rank because no job
     * is sent to master (rank 0).
     */
    std::vector<bool> workerIsBusy;

    /** MPI requests for jobs sent asynchronously to workers. Used to track when
     * the respective send buffers can be freed. */
    std::vector<MPI_Request> sendRequests;

    /** Jobs that have been sent to workers. Required for handling replies and
     * signaling the client that processing has completed. */
    std::vector<JobData *> sentJobsData;

    /** Mutex to protect access to `queue`. */
    pthread_mutex_t mutexQueue = PTHREAD_MUTEX_INITIALIZER;

    /** Semaphore to limit queue length and avoid potentially huge memory
     * allocation for all send and receive buffers. Note that using this might
     * come with a decreasing performance due to frequent rescheduling
     */
    sem_t semQueue = {};

    /** Thread that runs the message dispatcher. */
    pthread_t queueThread = 0;

    /** Value to indicate that there is currently no known free worker. */
    constexpr static int NO_FREE_WORKER = -1;
};

#endif

} // namespace parpe

#endif // LOADBALANCERMASTER_H
