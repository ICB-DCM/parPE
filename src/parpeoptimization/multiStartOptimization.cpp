#include <boost/asio.hpp>

#include <parpeoptimization/multiStartOptimization.h>
#include <parpecommon/logging.h>
#include <parpecommon/parpeException.h>

#include <cstdlib>
#include <cstring>
#include <future>

namespace parpe {


MultiStartOptimization::MultiStartOptimization(
    MultiStartOptimizationProblem &problem,
    bool runParallel,
    int first_start_idx)
    : msProblem(problem),
      numberOfStarts(problem.getNumberOfStarts()),
      restartOnFailure(problem.restartOnFailure()),
      runParallel(runParallel),
      first_start_idx(first_start_idx)
{

}

void MultiStartOptimization::run() {
    if (runParallel)
        runMultiThreaded();
    else
        runSingleThreaded();
}

void MultiStartOptimization::runMultiThreaded() const
{
    // Determine thread pool size
    // (note that hardware_concurrency() may return 0)
    auto num_threads = std::max(std::thread::hardware_concurrency(), 1U);
    if(auto env = std::getenv("PARPE_NUM_PARALLEL_STARTS")) {
        num_threads = std::stoi(env);
    }
    num_threads = std::min(num_threads,
                           static_cast<unsigned int>(numberOfStarts));

    logmessage(loglevel::debug,
               "Running %d starts using %d threads",
               numberOfStarts, num_threads);

    boost::asio::thread_pool pool(num_threads);

    auto num_successful_starts = 0;
    auto num_finished_starts = 0;
    auto lastStartIdx = -1;

    // submit the minimum number of starts
    std::vector<std::future<std::pair<int, int>>> futures;
    futures.reserve(numberOfStarts);
    for (int start_idx = 0; start_idx < numberOfStarts; ++start_idx) {
        futures.push_back(
            boost::asio::post(
                pool,
                std::packaged_task<std::pair<int, int>()>([this, start_idx] {
                    return std::make_pair(start_idx, runStart(start_idx));
                })));
        ++lastStartIdx;
    }

    // Report finished runs and restart if necessary
    while ((restartOnFailure && num_successful_starts < numberOfStarts)
           || (!restartOnFailure && num_finished_starts < numberOfStarts)) {
        for (auto &future: futures) {
            // future value might have been retrieved before
            if(!future.valid()) {
                continue;
            }

            if(auto status = future.wait_for(std::chrono::milliseconds(1));
                status != std::future_status::ready) {
                continue;
            }

            ++num_finished_starts;
            auto [start_idx, retval] = future.get();

            if (retval == 0) {
                // success
                logmessage(loglevel::debug,
                           "Optimization #%d finished successfully",
                           start_idx);
                ++num_successful_starts;
            } else if (!restartOnFailure) {
                // failure, no new start
                logmessage(loglevel::debug,
                           "Optimization ms #%d finished "
                           "unsuccessfully. Not trying "
                           "new starting point.",
                           start_idx);
            } else {
                // failure, new start
                logmessage(loglevel::debug,
                           "Thread ms #%d finished unsuccessfully... "
                           "trying new starting point", start_idx);
                ++lastStartIdx;

                future = boost::asio::post(
                    pool,
                    std::packaged_task<std::pair<int, int>()>(
                        [this, start_idx=lastStartIdx] {
                            return std::make_pair(start_idx, runStart(start_idx));
                        }));
            }
        }
    }

    pool.join();

    logmessage(loglevel::debug, "Multi-start optimization finished.");
}

void MultiStartOptimization::runSingleThreaded()
{
    logmessage(loglevel::debug,
               "Starting runParallelMultiStartOptimization with %d starts sequentially",
               numberOfStarts);

    int ms = 0;
    int numSucceeded = 0;

    while(true) {
        if(restartOnFailure && numSucceeded == numberOfStarts)
            break;

        if(ms == numberOfStarts)
            break;

        auto result = runStart(ms);

        if(result) {
            logmessage(loglevel::debug,
                       "Start #%d finished successfully", ms);
            ++numSucceeded;
        } else {
            logmessage(loglevel::debug, "Start ms #%d finished "
                                        "unsuccessfully.",ms);
        }
        ++ms;
    }

    logmessage(loglevel::debug, "runParallelMultiStartOptimization finished");
}

void MultiStartOptimization::setRunParallel(bool runParallel)
{
    this->runParallel = runParallel;
}

int MultiStartOptimization::runStart(int start_idx) const
{
    logmessage(loglevel::debug,
               "Starting local optimization #%d", start_idx);

    auto problem = msProblem.getLocalProblem(first_start_idx + start_idx);
    return getLocalOptimum(problem.get());
}

} // namespace parpe
