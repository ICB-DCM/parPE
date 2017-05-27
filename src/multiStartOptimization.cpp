#include "multiStartOptimization.h"

#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

#include "logging.h"


OptimizationProblem *getLocalProblem(optimizationProblemGeneratorForMultiStartFp problemGenerator,
                                     int multiStartIndex, void *userData) {
    OptimizationProblem *problem = problemGenerator(multiStartIndex, userData);

    if(!problem->initialParameters) {
        problem->initialParameters = new double [problem->numOptimizationParameters];
        getRandomStartingpoint(problem->parametersMin,
                               problem->parametersMax,
                               problem->numOptimizationParameters,
                               problem->initialParameters);
    }

    return problem;
}

OptimizationProblem **createLocalOptimizationProblems(optimizationProblemGeneratorForMultiStartFp problemGenerator,
                                                      int numLocalOptimizations, void *userData) {
    OptimizationProblem **localProblems = new OptimizationProblem*[numLocalOptimizations];

    for(int ms = 0; ms < numLocalOptimizations; ++ms) {
        localProblems[ms] = getLocalProblem(problemGenerator, ms, userData);
    }

    return localProblems;
}

int runParallelMultiStartOptimization(optimizationProblemGeneratorForMultiStartFp problemGenerator,
                                      int numberOfStarts, bool restartOnFailure, void *userData)
{
    logmessage(LOGLVL_DEBUG, "Starting runParallelMultiStartOptimization with %d starts", numberOfStarts);

    pthread_t *localOptimizationThreads = (pthread_t *) alloca(numberOfStarts * sizeof(pthread_t));

    OptimizationProblem **localProblems = createLocalOptimizationProblems(problemGenerator, numberOfStarts, userData);

    pthread_attr_t threadAttr;
    pthread_attr_init(&threadAttr);
    pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);

    int lastStartIdx = -1;

    // launch threads for required number of starts
    for(int ms = 0; ms < numberOfStarts; ++ms) {
        ++lastStartIdx;

        logmessage(LOGLVL_DEBUG, "Spawning thread for local optimization #%d (%d)", lastStartIdx, ms);

        pthread_create(&localOptimizationThreads[ms], &threadAttr, getLocalOptimumThreadWrapper, (void *)localProblems[ms]);
    }

    int numCompleted = 0;

    while(numCompleted < numberOfStarts) {
        for(int ms = 0; ms < numberOfStarts; ++ms) {
            if(!localOptimizationThreads[ms])
                continue;

            int *threadStatus = 0;
            int joinStatus = pthread_tryjoin_np(localOptimizationThreads[ms], (void**) &threadStatus);

            if(joinStatus == 0) { // joined successful
                if(*threadStatus == 0 || !restartOnFailure) {
                    if(*threadStatus == 0)
                        logmessage(LOGLVL_DEBUG, "Thread ms #%d finished successfully", ms);
                    else
                        logmessage(LOGLVL_DEBUG, "Thread ms #%d finished unsuccessfully. Not trying new starting point.", ms);
                    localOptimizationThreads[ms] = 0;
                    delete localProblems[ms];
                    ++numCompleted;
                } else {
                    delete localProblems[ms];
                    logmessage(LOGLVL_WARNING, "Thread ms #%d finished unsuccessfully... trying new starting point", ms);
                    ++lastStartIdx;

                    localProblems[ms] = getLocalProblem(problemGenerator, lastStartIdx, userData);
                    logmessage(LOGLVL_DEBUG, "Spawning thread for local optimization #%d (%d)", lastStartIdx, ms);
                    pthread_create(&localOptimizationThreads[ms], &threadAttr, getLocalOptimumThreadWrapper, (void *)localProblems[ms]);
                }
                delete threadStatus;
            }
        }

        sleep(0.1); // TODO: replace by condition via ThreadWrapper
    }
    delete[] localProblems;
    pthread_attr_destroy(&threadAttr);

    return 0;

}
