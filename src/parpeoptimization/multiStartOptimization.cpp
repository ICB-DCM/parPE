#include <parpeoptimization/multiStartOptimization.h>
#include <parpecommon/logging.h>
#include <parpecommon/parpeException.h>

#include <cstdlib>
#include <cstring>
#include <pthread.h>
#include <unistd.h>
#include <cassert>

namespace parpe {


MultiStartOptimization::MultiStartOptimization(
        MultiStartOptimizationProblem &problem,
        bool runParallel)
    : msProblem(problem),
      numberOfStarts(problem.getNumberOfStarts()),
      restartOnFailure(problem.restartOnFailure()),
      runParallel(runParallel)
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
    logmessage(LOGLVL_DEBUG,
               "Starting runParallelMultiStartOptimization with %d starts",
               numberOfStarts);

    std::vector<pthread_t> localOptimizationThreads(numberOfStarts);

    std::vector<OptimizationProblem *> localProblems =
            createLocalOptimizationProblems();

    if(localProblems.size() != static_cast<std::vector<OptimizationProblem *>::size_type>(numberOfStarts)) {
        throw ParPEException("Number of problems does not match number of specific starts.");
    }

    pthread_attr_t threadAttr;
    pthread_attr_init(&threadAttr);
    pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);

    int lastStartIdx = -1;
    // launch threads for required number of starts
    for (int ms = 0; ms < numberOfStarts; ++ms) {
        ++lastStartIdx;

        logmessage(LOGLVL_DEBUG,
                   "Spawning thread for local optimization #%d (%d)",
                   lastStartIdx, ms);

        pthread_create(&localOptimizationThreads.at(ms), &threadAttr,
                       getLocalOptimumThreadWrapper, static_cast<void *>(localProblems[ms]));
    }

    int numCompleted = 0;

    while (numCompleted < numberOfStarts) {
        for (int ms = 0; ms < numberOfStarts; ++ms) {
            // problem still running?
            if (!localProblems[ms])
                continue;

            int *threadStatus = nullptr;
#ifndef __APPLE__
            // TODO(#84) pthread_tryjoin_np is not available on macOS. can replace easily by pthread_join, but this would only allow restarting failed threads rather late, so we disable the retry option for now.
            int joinStatus = pthread_tryjoin_np(localOptimizationThreads[ms],
                                                reinterpret_cast<void **>(&threadStatus));
#else
            int joinStatus = pthread_join(localOptimizationThreads[ms],
                                                reinterpret_cast<void **>(&threadStatus));
#endif
            if (joinStatus == 0) { // joined successful
                delete localProblems[ms];
                localProblems[ms] = nullptr;

                if (*threadStatus == 0 || !restartOnFailure) {
                    if (*threadStatus == 0) {
                        logmessage(LOGLVL_DEBUG,
                                   "Thread ms #%d finished successfully", ms);
                    } else {
                        logmessage(LOGLVL_DEBUG, "Thread ms #%d finished "
                                                 "unsuccessfully. Not trying "
                                                 "new starting point.",
                                   ms);
                    }
                    ++numCompleted;
                }
#ifndef __APPLE__
                else {
                    logmessage(LOGLVL_WARNING, "Thread ms #%d finished "
                                               "unsuccessfully... trying new "
                                               "starting point",
                               ms);
                    ++lastStartIdx;

                    localProblems[ms] = msProblem.getLocalProblem(lastStartIdx).release();
                    logmessage(
                                LOGLVL_DEBUG,
                                "Spawning thread for local optimization #%d (%d)",
                                lastStartIdx, ms);
                    pthread_create(&localOptimizationThreads[ms], &threadAttr,
                                   getLocalOptimumThreadWrapper,
                                   static_cast<void *>(localProblems[ms]));
                }
#endif
                delete threadStatus;
            }
        }

        sleep(1); // TODO: replace by condition via ThreadWrapper
    }

    logmessage(LOGLVL_DEBUG, "runParallelMultiStartOptimization finished");

    pthread_attr_destroy(&threadAttr);
}

void MultiStartOptimization::runSingleThreaded()
{
    logmessage(LOGLVL_DEBUG,
               "Starting runParallelMultiStartOptimization with %d starts sequentially",
               numberOfStarts);

    int ms = 0;
    int numSucceeded = 0;

    while(true) {
        if(restartOnFailure && numSucceeded == numberOfStarts)
            break;

        if(ms == numberOfStarts)
            break;

        auto problem = msProblem.getLocalProblem(ms);
        auto result = getLocalOptimum(problem.get());
        if(result) {
            logmessage(LOGLVL_DEBUG,
                       "Start #%d finished successfully", ms);
            ++numSucceeded;
        } else {
            logmessage(LOGLVL_DEBUG, "Thread ms #%d finished "
                                     "unsuccessfully.",ms);
        }
        ++ms;
    }

    logmessage(LOGLVL_DEBUG, "runParallelMultiStartOptimization finished");
}

void MultiStartOptimization::setRunParallel(bool runParallel)
{
    this->runParallel = runParallel;
}

std::vector<OptimizationProblem *>
MultiStartOptimization::createLocalOptimizationProblems() {
    std::vector<OptimizationProblem *> localProblems(numberOfStarts);

    for (int ms = 0; ms < numberOfStarts; ++ms) {
        localProblems[ms] = msProblem.getLocalProblem(ms).release();
    }

    return localProblems;
}

} // namespace parpe
