#include "LoadBalancerMaster.h"

#include <cassert>
#include <climits>
#include <misc.h>

namespace parpe {

void LoadBalancerMaster::run() {
    if (isRunning_)
        return;

#ifndef QUEUE_MASTER_TEST
    assertMpiActive();
#endif

    int mpiCommSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiCommSize);
    assert(mpiCommSize > 1 &&
           "Need multiple MPI processes!"); // crashes otherwise

    numWorkers = mpiCommSize - 1;
    sentJobsData = new JobData *[numWorkers];
    workerIsBusy = new bool[numWorkers]();
    sendRequests = new MPI_Request[numWorkers];
    for (int i = 0; i < numWorkers; ++i) // have to initialize before can wait!
        sendRequests[i] = MPI_REQUEST_NULL;

    // Create semaphore to limit queue length
    // and avoid huge memory allocation for all send and receive buffers
    unsigned int queueMaxLength = mpiCommSize;
#ifdef SEM_VALUE_MAX
    queueMaxLength =
        SEM_VALUE_MAX < queueMaxLength ? SEM_VALUE_MAX : queueMaxLength;
#endif
    sem_init(&semQueue, 0, queueMaxLength);

    pthread_attr_t threadAttr;
    pthread_attr_init(&threadAttr);
    pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);
    pthread_create(&queueThread, &threadAttr, threadEntryPoint, this);
    pthread_attr_destroy(&threadAttr);

    isRunning_ = true;
}

LoadBalancerMaster::~LoadBalancerMaster()
{
    terminate();
}

#ifndef QUEUE_MASTER_TEST
void LoadBalancerMaster::assertMpiActive() {
    assert(getMpiActive());
}
#endif

void *LoadBalancerMaster::threadEntryPoint(void *vpLoadBalancerMaster) {
    auto master = static_cast<LoadBalancerMaster *>(vpLoadBalancerMaster);
    return master->loadBalancerThreadRun();
}

void *LoadBalancerMaster::loadBalancerThreadRun() {

    // dispatch queued work packages
    while (true) {

        // check if any job finished
        int lastFinishedWorkerIdx = handleFinishedJobs();

        freeEmptiedSendBuffers();

        // getNextFreeWorker
        int freeWorkerIndex = lastFinishedWorkerIdx;

        if (freeWorkerIndex < 0) {
            // no job finished recently, check free slots
            freeWorkerIndex = getNextFreeWorkerIndex();
        }

        if (freeWorkerIndex < 0) {
            // add cancellation point to avoid invalid reads in
            // loadBalancer.recvRequests
            pthread_testcancel();

            // all workers are busy, wait for next one to finish
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            freeWorkerIndex = handleReply(&status);
        }

        // found free worker, check for jobs to do
        JobData *currentQueueElement = getNextJob();

        if (currentQueueElement) {
            sendToWorker(freeWorkerIndex, currentQueueElement);
            sentJobsData[freeWorkerIndex] = currentQueueElement;
        }

        pthread_yield();
    };

    return nullptr;
}

void LoadBalancerMaster::freeEmptiedSendBuffers() {
    // free any emptied send buffers
    while (true) {
        int emptiedBufferIdx = -1;
        int anySendCompleted = 0;
        MPI_Testany(numWorkers, sendRequests, &emptiedBufferIdx,
                    &anySendCompleted, MPI_STATUS_IGNORE);

        if (anySendCompleted && emptiedBufferIdx != MPI_UNDEFINED) {
            delete[] sentJobsData[emptiedBufferIdx]->sendBuffer;
        } else {
            break;
        }
    }
}

int LoadBalancerMaster::handleFinishedJobs() {
    MPI_Status status;
    int finishedWorkerIdx = -1;

    // handle all finished jobs, if any
    while (true) {
        // add cancellation point to avoid invalid reads in
        // loadBalancer.recvRequests
        pthread_testcancel();

        int messageWaiting;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &messageWaiting,
                   &status);

        if (messageWaiting) {
            // some job is finished
            finishedWorkerIdx = handleReply(&status);
        } else {
            // there was nothing to be finished
            break;
        }
    }
    return finishedWorkerIdx;
}

int LoadBalancerMaster::getNextFreeWorkerIndex() {
    for (int i = 0; i < numWorkers; ++i) {
        if (workerIsBusy[i] == false) {
            return i;
        }
    }
    return -1;
}

JobData *LoadBalancerMaster::getNextJob() {

    pthread_mutex_lock(&mutexQueue);

    JobData *nextJob = nullptr;
    if (!queue.empty()) {
        nextJob = queue.front();
        queue.pop();
    }

    pthread_mutex_unlock(&mutexQueue);

    return nextJob;
}

void LoadBalancerMaster::sendToWorker(int workerIdx, JobData *data) {
    assert(workerIdx >= 0);
    assert(workerIdx < numWorkers);

    workerIsBusy[workerIdx] = true;

    int tag = data->jobId;
    int workerRank = workerIdx + 1;

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[36mSent job #%d to rank %d.\x1b[0m\n", tag, workerRank);
#endif

    MPI_Isend(data->sendBuffer, data->lenSendBuffer, MPI_BYTE, workerRank, tag,
              MPI_COMM_WORLD, &sendRequests[workerIdx]);

    sem_post(&semQueue);
}

void LoadBalancerMaster::queueJob(JobData *data) {
    assert(isRunning_);

    sem_wait(&semQueue);

    pthread_mutex_lock(&mutexQueue);

    if (lastJobId == INT_MAX) // Unlikely, but prevent overflow
        lastJobId = 0;

    data->jobId = ++lastJobId;

    queue.push(data);

    pthread_mutex_unlock(&mutexQueue);
}

void LoadBalancerMaster::terminate() {
    // avoid double termination
    pthread_mutex_lock(&mutexQueue);
    if (!isRunning_)
        return;
    isRunning_ = false;
    pthread_mutex_unlock(&mutexQueue);

    pthread_cancel(queueThread);
    // wait until canceled
    pthread_join(queueThread, nullptr);

    pthread_mutex_destroy(&mutexQueue);
    sem_destroy(&semQueue);

    if (sentJobsData)
        delete[] sentJobsData;
    if (sendRequests)
        delete[] sendRequests;
    if (workerIsBusy)
        delete[] workerIsBusy;

}

int LoadBalancerMaster::handleReply(MPI_Status *mpiStatus) {

    int workerIdx = mpiStatus->MPI_SOURCE - 1;
    JobData *data = sentJobsData[workerIdx];

    // allocate memory for result
    MPI_Get_count(mpiStatus, MPI_BYTE, &data->lenRecvBuffer);
    data->recvBuffer = new char[data->lenRecvBuffer];

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[32mReceiving result for job %d from %d (%dB)\x1b[0m\n",
           mpiStatus->MPI_TAG, mpiStatus->MPI_SOURCE, data->lenRecvBuffer);
#endif

    // receive
    MPI_Recv(data->recvBuffer, data->lenRecvBuffer, MPI_BYTE,
             mpiStatus->MPI_SOURCE, mpiStatus->MPI_TAG, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    workerIsBusy[workerIdx] = false;

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[32mReceived result for job %d from %d\x1b[0m\n",
           mpiStatus->MPI_TAG, mpiStatus->MPI_SOURCE);
#endif

    // signal job done
    if(data->jobDone)
        ++(*data->jobDone);

    // user-provided callback if specified
    if(data->callbackJobFinished)
        data->callbackJobFinished(data);

    pthread_mutex_lock(data->jobDoneChangedMutex);
    pthread_cond_signal(data->jobDoneChangedCondition);
    pthread_mutex_unlock(data->jobDoneChangedMutex);

    return workerIdx;
}

void LoadBalancerMaster::sendTerminationSignalToAllWorkers() {
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    MPI_Request reqs[commSize - 1];

    for (int i = 1; i < commSize; ++i) {
        reqs[i - 1] = MPI_REQUEST_NULL;
        MPI_Isend(MPI_BOTTOM, 0, MPI_INT, i, 0, MPI_COMM_WORLD, &reqs[i - 1]);
    }
    MPI_Waitall(commSize - 1, reqs, MPI_STATUS_IGNORE);
}

bool LoadBalancerMaster::isRunning() const
{
    return isRunning_;
}

} // namespace parpe
