#include <alloca.h>
#include <unistd.h>

#include "masterthread.h"
#include "dataprovider.h"
#include "localoptimization.h"

void *newMultiStartOptimization(void *multiStartIndexVP) {
    int multiStartIndex = *(int *) multiStartIndexVP;

    printf("Spawning thread for global optimization #%d\n", multiStartIndex);

    int numLocalOptimizations = getNumLocalOptimizationsForMultiStartRun(multiStartIndex);

    pthread_t *localOptimizationThreads = alloca(numLocalOptimizations * sizeof(pthread_t));

    pthread_attr_t threadAttr;
    pthread_attr_init(&threadAttr);
    pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);

    for(int ms = 0; ms < numLocalOptimizations; ++ms) {
        int id = multiStartIndex * 1000 + ms;
        printf("Spawning thread for local optimization #%d.%d\n", multiStartIndex, ms);
        pthread_create(&localOptimizationThreads[ms], &threadAttr, newLocalOptimization, (void *)&id);
    }
    pthread_attr_destroy(&threadAttr);

    for(int ms = 0; ms < numLocalOptimizations; ++ms) {
        pthread_join(localOptimizationThreads[ms], NULL);
        printf("Thread ms #%d finished\n", ms);
    }

    printf("Leaving thread for global optimization #%d\n", multiStartIndex);

    return 0;
}


void *newLocalOptimization(void *idVP) {
    int id = *(int *) idVP;
    datapath datapath;

    datapath.idxMultiStart = id / 1000;
    datapath.idxLocalOptimization = id % 1000;

    // TODO pass options object, also add IpOpt options to config file


    getLocalOptimum(datapath);
    printf("Finished newLocalOptimization #%d\n", id);

    return 0;
}
