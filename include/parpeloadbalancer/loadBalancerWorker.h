#ifndef LOADBALANCERWORKER_H
#define LOADBALANCERWORKER_H

#include <parpecommon/parpeConfig.h>

#include <functional>
#include <vector>

constexpr int MPI_TAG_EXIT_SIGNAL = 0;

namespace parpe {

class LoadBalancerWorker {
  public:
    LoadBalancerWorker() = default;

    /**
     * Callback function for when a message is received.
     *
     * @param buffer The message
     * @param jobId A message identifier, unique over the range of MAX_INT
     * messages.
     */
    using messageHandlerFunc =
        std::function<void(std::vector<char>& buffer, int jobId)>;

    void run(messageHandlerFunc const& messageHandler);

  private:
    /**
     * @brief Probe for and dispatch the next incoming job
     * @return `true` if the termination signal was received, `false` otherwise.
     */
    bool waitForAndHandleJobs(messageHandlerFunc const& messageHandler);
};

} // namespace parpe

#endif // LOADBALANCERWORKER_H
