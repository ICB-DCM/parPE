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
    for(int ms = 0; ms < numLocalOptimizations; ++ms) {
        ids[ms] = multiStartIndex * 1000 + ms;
        logmessage(LOGLVL_DEBUG, "Spawning thread for local optimization #%d.%d", multiStartIndex, ms);
        pthread_create(&localOptimizationThreads[ms], &threadAttr, newLocalOptimization, (void *)&ids[ms]);
    }
    pthread_attr_destroy(&threadAttr);

    for(int ms = 0; ms < numLocalOptimizations; ++ms) {
        pthread_join(localOptimizationThreads[ms], NULL);
        logmessage(LOGLVL_DEBUG, "Thread ms #%d finished", ms);
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
    getLocalOptimum(datapath);
    logmessage(LOGLVL_DEBUG, "Finished newLocalOptimization #%d.%d", datapath.idxMultiStart, datapath.idxLocalOptimization);

    // TODO need to acquire hdf5 mutex lock when running H5_term_library in at_exit
    return 0;
}
