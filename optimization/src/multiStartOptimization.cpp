#include "multiStartOptimization.h"
#include "logging.h"
#include <cstdlib>
#include <cstring>
#include <pthread.h>
#include <unistd.h>

namespace parpe {

MultiStartOptimization::MultiStartOptimization(int numberOfStarts,
                                               bool restartOnFailure)
    : numberOfStarts(numberOfStarts), restartOnFailure(restartOnFailure) {}

int MultiStartOptimization::run() {
    return runSingleThreaded();
    return runMultiThreaded();
}

int MultiStartOptimization::runMultiThreaded()
{
    logmessage(LOGLVL_DEBUG,
               "Starting runParallelMultiStartOptimization with %d starts",
               numberOfStarts);

    std::vector<pthread_t> localOptimizationThreads(numberOfStarts);

    std::vector<OptimizationProblem *> localProblems =
        createLocalOptimizationProblems();

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
                       getLocalOptimumThreadWrapper, (void *)localProblems[ms]);
    }

    int numCompleted = 0;

    while (numCompleted < numberOfStarts) {
        for (int ms = 0; ms < numberOfStarts; ++ms) {
            // problem still running?
            if (!localProblems[ms])
                continue;

            int *threadStatus = 0;
            int joinStatus = pthread_tryjoin_np(localOptimizationThreads[ms],
                                                (void **)&threadStatus);

            if (joinStatus == 0) { // joined successful
                delete localProblems[ms];
                localProblems[ms] = NULL;

                if (*threadStatus == 0 || !restartOnFailure) {
                    if (*threadStatus == 0)
                        logmessage(LOGLVL_DEBUG,
                                   "Thread ms #%d finished successfully", ms);
                    else
                        logmessage(LOGLVL_DEBUG, "Thread ms #%d finished "
                                                 "unsuccessfully. Not trying "
                                                 "new starting point.",
                                   ms);
                    ++numCompleted;
                } else {
                    logmessage(LOGLVL_WARNING, "Thread ms #%d finished "
                                               "unsuccessfully... trying new "
                                               "starting point",
                               ms);
                    ++lastStartIdx;

                    localProblems[ms] = getLocalProblem(lastStartIdx).release();
                    logmessage(
                        LOGLVL_DEBUG,
                        "Spawning thread for local optimization #%d (%d)",
                        lastStartIdx, ms);
                    pthread_create(&localOptimizationThreads[ms], &threadAttr,
                                   getLocalOptimumThreadWrapper,
                                   (void *)localProblems[ms]);
                }
                delete threadStatus;
            }
        }

        sleep(0.1); // TODO: replace by condition via ThreadWrapper
    }

    logmessage(LOGLVL_DEBUG, "runParallelMultiStartOptimization finished");

    pthread_attr_destroy(&threadAttr);

    return 0;
}

int MultiStartOptimization::runSingleThreaded()
{
    logmessage(LOGLVL_DEBUG,
               "Starting runParallelMultiStartOptimization with %d starts sequentially",
               numberOfStarts);

    int ms = 0;
    int numSucceeded = 0;

    while(true) {
        if(restartOnFailure && numSucceeded == numberOfStarts)
            break;
        else if(ms == numberOfStarts)
            break;

        auto problem = getLocalProblem(ms);
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

    return 0;

}

std::unique_ptr<OptimizationProblem>
MultiStartOptimization::getLocalProblem(int multiStartIndex) {
    std::unique_ptr<OptimizationProblem> problem = getLocalProblemImpl(multiStartIndex);

    return problem;
}

std::vector<OptimizationProblem *>
MultiStartOptimization::createLocalOptimizationProblems() {
    std::vector<OptimizationProblem *> localProblems(numberOfStarts);

    for (int ms = 0; ms < numberOfStarts; ++ms) {
        localProblems[ms] = getLocalProblem(ms).release();
    }

    return localProblems;
}

} // namespace parpe
