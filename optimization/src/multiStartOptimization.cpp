#include "multiStartOptimization.h"

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "logging.h"

int runParallelMultiStartOptimization(
    OptimizationProblemGeneratorForMultiStart *problemGenerator,
    int numberOfStarts, bool restartOnFailure) {
    logmessage(LOGLVL_DEBUG,
               "Starting runParallelMultiStartOptimization with %d starts",
               numberOfStarts);

    pthread_t *localOptimizationThreads =
        (pthread_t *)alloca(numberOfStarts * sizeof(pthread_t));

    OptimizationProblem **localProblems =
        problemGenerator->createLocalOptimizationProblems(numberOfStarts);

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

        pthread_create(&localOptimizationThreads[ms], &threadAttr,
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

                    localProblems[ms] =
                        problemGenerator->getLocalProblem(lastStartIdx);
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

    delete[] localProblems;
    pthread_attr_destroy(&threadAttr);

    return 0;
}

OptimizationProblemGeneratorForMultiStart::
    OptimizationProblemGeneratorForMultiStart() {}

OptimizationProblem *OptimizationProblemGeneratorForMultiStart::getLocalProblem(
    int multiStartIndex) {
    OptimizationProblem *problem = getLocalProblemImpl(multiStartIndex);

    return problem;
}

OptimizationProblem **
OptimizationProblemGeneratorForMultiStart::createLocalOptimizationProblems(
    int numLocalOptimizations) {
    OptimizationProblem **localProblems =
        new OptimizationProblem *[numLocalOptimizations];

    for (int ms = 0; ms < numLocalOptimizations; ++ms) {
        localProblems[ms] = getLocalProblem(ms);
    }

    return localProblems;
}
