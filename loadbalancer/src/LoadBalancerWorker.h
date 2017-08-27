#ifndef LOADBALANCERWORKER_H
#define LOADBALANCERWORKER_H

#define MPI_TAG_EXIT_SIGNAL 0
#define QUEUE_WORKER_H_VERBOSE 0

class LoadBalancerWorker {
  public:
    LoadBalancerWorker() = default;

    /**
     * messageHandler is called by run when a message is received. The
     * message is contained in buffer.
     * @param buffer The message
     * @param size The size of the message in bytes
     * @param jobId is a message identifier, unique over the range of MAX_INT
     * messages
     */
    virtual void messageHandler(char **buffer, int *size, int jobId) = 0;

    void run();

  private:
    /**
     * @brief waitForAndHandleJobs
     * @return true: received termination signal
     */
    bool waitForAndHandleJobs();
};

#endif // LOADBALANCERWORKER_H
