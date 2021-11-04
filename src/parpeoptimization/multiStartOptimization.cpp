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

void MultiStartOptimization::runMultiThreaded()
{
    // Determine thread pool size
    auto num_threads = std::thread::hardware_concurrency();
    if(auto env = std::getenv("PARPE_NUM_PARALLEL_STARTS")) {
        num_threads = std::stod(env);
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
                    logmessage(loglevel::debug,
                               "Starting local optimization #%d", start_idx);

                    auto problem = this->msProblem.getLocalProblem(start_idx);
                    return std::make_pair(start_idx,
                                          getLocalOptimum(problem.get()));
                })));
        ++lastStartIdx;
    }

    // Report finished runs and restart if necessary
    while ((restartOnFailure && num_successful_starts < numberOfStarts)
           || (!restartOnFailure && num_finished_starts < numberOfStarts)) {
        for (auto &future: futures) {
            if(!future.valid()) {
                continue;
            }
            if(auto status = future.wait_for(std::chrono::milliseconds(1));
                status != std::future_status::ready) {
                continue;
            }
            ++num_finished_starts;

            auto [start_idx, retval] = future.get();
            if (retval == 0 || !restartOnFailure) {
                if (retval == 0) {
                    logmessage(loglevel::debug,
                               "Optimization #%d finished successfully",
                               start_idx);
                    ++num_successful_starts;
                } else {
                    logmessage(loglevel::debug,
                               "Optimization ms #%d finished "
                               "unsuccessfully. Not trying "
                               "new starting point.",
                               start_idx);
                }
            } else {
                logmessage(loglevel::warning,
                           "Thread ms #%d finished "
                           "unsuccessfully... trying new "
                           "starting point",
                           start_idx);
                ++lastStartIdx;

                future = boost::asio::post(
                    pool,
                    std::packaged_task<std::pair<int, int>()>(
                        [this, start_idx=lastStartIdx] {
                            logmessage(loglevel::debug,
                                       "Starting local optimization #%d",
                                       start_idx);

                            auto problem = msProblem.getLocalProblem(start_idx);
                            return std::make_pair(
                                start_idx, getLocalOptimum(problem.get()));
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

        auto problem = msProblem.getLocalProblem(first_start_idx + ms);
        auto result = getLocalOptimum(problem.get());
        if(result) {
            logmessage(loglevel::debug,
                       "Start #%d finished successfully", ms);
            ++numSucceeded;
        } else {
            logmessage(loglevel::debug, "Thread ms #%d finished "
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

std::vector<OptimizationProblem *>
MultiStartOptimization::createLocalOptimizationProblems() {
    std::vector<OptimizationProblem *> localProblems(numberOfStarts);

    for (int ms = 0; ms < numberOfStarts; ++ms) {
        localProblems[ms] = msProblem.getLocalProblem(first_start_idx + ms).release();
    }

    return localProblems;
}

} // namespace parpe
