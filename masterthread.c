#include <alloca.h>
#include <unistd.h>

#include "masterthread.h"
#include "dataprovider.h"
#include "localoptimization.h"
#include "misc.h"

void *newMultiStartOptimization(void *multiStartIndexVP) {
    int multiStartIndex = *(int *) multiStartIndexVP;

    logmessage(LOGLVL_DEBUG, "Spawning thread for global optimization #%d", multiStartIndex);

    int numLocalOptimizations = getNumLocalOptimizationsForMultiStartRun(multiStartIndex);

    pthread_t *localOptimizationThreads = alloca(numLocalOptimizations * sizeof(pthread_t));

    pthread_attr_t threadAttr;
    pthread_attr_init(&threadAttr);
    pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);

    int ids[numLocalOptimizations];

    int lastStartIdx = -1;

    for(int ms = 0; ms < numLocalOptimizations; ++ms) {
        ids[ms] = multiStartIndex * 1000 + ++lastStartIdx;
        logmessage(LOGLVL_DEBUG, "Spawning thread for local optimization #%d.%d (%d)", multiStartIndex, lastStartIdx, ms);
        pthread_create(&localOptimizationThreads[ms], &threadAttr, newLocalOptimization, (void *)&ids[ms]);
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

                    ids[ms] = multiStartIndex * 1000 + ++lastStartIdx;
                    logmessage(LOGLVL_DEBUG, "Spawning thread for local optimization #%d.%d (%d)", multiStartIndex, lastStartIdx, ms);
                    pthread_create(&localOptimizationThreads[ms], &threadAttr, newLocalOptimization, (void *)&ids[ms]);
                }
                free(threadStatus);
            }
        }

        sleep(1);
    }

    logmessage(LOGLVL_DEBUG, "Leaving thread for global optimization #%d", multiStartIndex);

    return 0;
}


void *newLocalOptimization(void *idVP) {
    int id = *(int *) idVP;
    datapath datapath = {INT_MIN};

    datapath.idxMultiStart = id / 1000;
    datapath.idxLocalOptimization = id % 1000;

    // TODO pass options object, also add IpOpt options to config file
    logmessage(LOGLVL_DEBUG, "Starting newLocalOptimization #%d.%d", datapath.idxMultiStart, datapath.idxLocalOptimization);
    int *status = malloc(sizeof(int));
    *status = getLocalOptimum(datapath);
    logmessage(LOGLVL_DEBUG, "Finished newLocalOptimization #%d.%d", datapath.idxMultiStart, datapath.idxLocalOptimization);

    return status;
}
