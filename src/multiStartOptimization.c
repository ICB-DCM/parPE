#include "multiStartOptimization.h"

#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

#include "misc.h"

MultiStartOptimization *multiStartOptimizationNew()
{
    MultiStartOptimization *ms = malloc(sizeof(*ms));
    memset(ms, 0, sizeof(*ms));

    return ms;

}


int runParallelMultiStartOptimization(MultiStartOptimization *multiStartOptimization)
{
    int numLocalOptimizations = multiStartOptimization->numberOfStarts;

    pthread_t *localOptimizationThreads = alloca(numLocalOptimizations * sizeof(pthread_t));

    OptimizationProblem localProblems[numLocalOptimizations]; // need to keep, since passed by ref to new thread

    pthread_attr_t threadAttr;
    pthread_attr_init(&threadAttr);
    pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);

    int lastStartIdx = -1;

    // launch threads for required number of starts
    for(int ms = 0; ms < numLocalOptimizations; ++ms) {
        localProblems[ms] = *multiStartOptimization->optimizationProblem;
        localProblems[ms].initialParameters = malloc(sizeof(double) * multiStartOptimization->optimizationProblem->numOptimizationParameters);
        multiStartOptimization->getInitialPoint(multiStartOptimization, ms, localProblems[ms].initialParameters);
        localProblems[ms].startIdx = ++lastStartIdx;

        logmessage(LOGLVL_DEBUG, "Spawning thread for local optimization #%d (%d)", lastStartIdx, ms);

        pthread_create(&localOptimizationThreads[ms], &threadAttr, getLocalOptimumThreadWrapper, (void *)&localProblems[ms]);
    }
    pthread_attr_destroy(&threadAttr);

    int numCompleted = 0;

    while(numCompleted < numLocalOptimizations) {
        for(int ms = 0; ms < numLocalOptimizations; ++ms) {
            if(!localOptimizationThreads[ms])
                continue;

            void *threadStatus = 0;
            int joinStatus = pthread_tryjoin_np(localOptimizationThreads[ms], &threadStatus);

            if(joinStatus == 0) { // joined successful
                if(*(int*)threadStatus == 0) {
                    logmessage(LOGLVL_DEBUG, "Thread ms #%d finished successfully", ms);
                    localOptimizationThreads[ms] = 0;
                    ++numCompleted;
                } else {
                    logmessage(LOGLVL_WARNING, "Thread ms #%d finished unsuccessfully... trying new starting point", ms);
                    localProblems[ms].startIdx = ++lastStartIdx;
                    logmessage(LOGLVL_DEBUG, "Spawning thread for local optimization #%d (%d)", lastStartIdx, ms);
                    pthread_create(&localOptimizationThreads[ms], &threadAttr, getLocalOptimumThreadWrapper, (void *)&localProblems[ms]);
                }
                free(threadStatus);
            }
        }

        sleep(1);
    }

    return 0;

}


