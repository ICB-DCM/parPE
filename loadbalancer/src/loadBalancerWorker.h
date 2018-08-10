#ifndef LOADBALANCERWORKER_H
#define LOADBALANCERWORKER_H

#define MPI_TAG_EXIT_SIGNAL 0

#include <vector>
#include <functional>

namespace parpe {

class LoadBalancerWorker {
  public:
    LoadBalancerWorker() = default;

    /**
     * messageHandler is called by run when a message is received. The
     * message is contained in buffer.
     * @param buffer The message
     * @param jobId is a message identifier, unique over the range of MAX_INT
     * messages
     */
    using messageHandlerFunc = std::function<void (std::vector<char> &buffer, int jobId)>;

    void run(messageHandlerFunc messageHandler);

  private:
    /**
     * @brief waitForAndHandleJobs
     * @return true: received termination signal
     */
    bool waitForAndHandleJobs(messageHandlerFunc messageHandler);
};

} // namespace parpe

#endif // LOADBALANCERWORKER_H
