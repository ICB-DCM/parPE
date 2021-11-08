#include <parpeloadbalancer/loadBalancerMaster.h>

#ifdef PARPE_ENABLE_MPI

#include <cassert>
#include <climits>
#include <chrono>

#include <parpecommon/misc.h>
#include <parpecommon/parpeException.h>

//#define MASTER_QUEUE_H_SHOW_COMMUNICATION 1

namespace parpe {

void LoadBalancerMaster::run() {
    if (isRunning_)
        return;

#ifndef QUEUE_MASTER_TEST
    assertMpiActive();
#endif

    int mpiCommSize;
    MPI_Comm_size(mpiComm, &mpiCommSize);

    if(mpiCommSize <= 2) {
        // crashes otherwise
        throw std::runtime_error("Need at least 2 MPI processes!");
    }

    numWorkers = mpiCommSize - 1;
    sentJobsData.resize(numWorkers, nullptr);
    workerIsBusy.resize(numWorkers, false);
    // have to initialize before can wait!
    sendRequests.resize(numWorkers, MPI_REQUEST_NULL);

    // Create semaphore to limit queue length
#ifdef SEM_VALUE_MAX
    unsigned int queueMaxLength = SEM_VALUE_MAX;
#else
    unsigned int queueMaxLength = UINT_MAX;
#endif
    sem_init(&semQueue, 0, queueMaxLength);

    queueThread = std::thread(&LoadBalancerMaster::loadBalancerThreadRun, this);

    isRunning_ = true;
}

LoadBalancerMaster::~LoadBalancerMaster()
{
    terminate();
    sem_destroy(&semQueue);
}

#ifndef QUEUE_MASTER_TEST
void LoadBalancerMaster::assertMpiActive() {
    Expects(getMpiActive());
}
#endif


void LoadBalancerMaster::loadBalancerThreadRun() {

    // dispatch queued work packages
    while (queue_thread_continue_) {
        int freeWorkerIndex = NO_FREE_WORKER;

        // empty send queue while there are free workers
        while(queue_thread_continue_ && (freeWorkerIndex = getNextFreeWorkerIndex()) >= 0
              && sendQueuedJob(freeWorkerIndex)) {}

        // check if any job finished
        handleFinishedJobs();

        freeEmptiedSendBuffers();
    };
}

void LoadBalancerMaster::freeEmptiedSendBuffers() {
    // free any emptied send buffers
    while (true) {
        int emptiedBufferIdx = MPI_UNDEFINED;
        int anySendCompleted = 0;
        MPI_Testany(sendRequests.size(), sendRequests.data(), &emptiedBufferIdx,
                    &anySendCompleted, MPI_STATUS_IGNORE);

        if (anySendCompleted && emptiedBufferIdx != MPI_UNDEFINED
                && sentJobsData[emptiedBufferIdx]) {
            /* By the time we check for send to be finished, we might have received the reply
             * already and the pointed-to object might have been already destroyed. This
             * is therefore set to nullptr when receiving the reply. */
            std::vector<char>().swap(sentJobsData[emptiedBufferIdx]->sendBuffer);
        } else {
            break;
        }
    }
}

int LoadBalancerMaster::handleFinishedJobs() {
    int finishedWorkerIdx = NO_FREE_WORKER;

    // handle all finished jobs, if any
    while (true) {
        // check for waiting incoming message
        MPI_Status status;
        int messageWaiting = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpiComm,
                   &messageWaiting, &status);

        if (messageWaiting) {
            // some job is finished, process that
            finishedWorkerIdx = handleReply(&status);

            // directly send new work if available
            if(sendQueuedJob(finishedWorkerIdx))
                finishedWorkerIdx = NO_FREE_WORKER; // not free anymore
        } else {
            // there was nothing to be finished
            break;
        }
    }
    return finishedWorkerIdx;
}

int LoadBalancerMaster::getNextFreeWorkerIndex() {
    for (unsigned int i = 0; i < workerIsBusy.size(); ++i) {
        if (!workerIsBusy[i])
            return i;
    }

    return NO_FREE_WORKER;
}

JobData *LoadBalancerMaster::getNextJob() {

    std::unique_lock lock(mutexQueue);

    JobData *nextJob = nullptr;
    if (!queue.empty()) {
        nextJob = queue.front();
        queue.pop();
    }

    return nextJob;
}

void LoadBalancerMaster::sendToWorker(int workerIdx, JobData *data) {
    Expects(workerIdx >= 0);
    Expects(workerIdx < numWorkers);

    workerIsBusy[workerIdx] = true;

    int tag = data->jobId;
    int workerRank = workerIdx + 1;

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[31mSending job #%d to rank %d (%luB).\x1b[0m\n", tag, workerRank, data->sendBuffer.size());
#endif

    MPI_Isend(data->sendBuffer.data(), data->sendBuffer.size(), mpiJobDataType,
              workerRank, tag,
              mpiComm, &sendRequests[workerIdx]);

    sem_post(&semQueue);
}

void LoadBalancerMaster::queueJob(JobData *data) {
    RELEASE_ASSERT(isRunning_, "Can't queue job while not running.");

    sem_wait(&semQueue);

    std::unique_lock lock(mutexQueue);

    if (lastJobId == INT_MAX) // Unlikely, but prevent overflow
        lastJobId = 0;

    data->jobId = ++lastJobId;

    queue.push(data);

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    int size = sizeof(*data) + data->sendBuffer.size() + data->recvBuffer.size();
    printf("\x1b[33mQueued job with size %dB. New queue length is %d.\x1b[0m\n", size, queue.size());
#endif

}

void LoadBalancerMaster::terminate() {
    if (!isRunning_) {
        // avoid double termination
        return;
    }
    isRunning_ = false;
    queue_thread_continue_ = false;
    queueThread.join();
}

int LoadBalancerMaster::handleReply(MPI_Status *mpiStatus) {

    int workerIdx = mpiStatus->MPI_SOURCE - 1;
    JobData *data = sentJobsData[workerIdx];
    sentJobsData[workerIdx] = nullptr;

    // allocate memory for result
    int lenRecvBuffer = 0;
    MPI_Get_count(mpiStatus, mpiJobDataType, &lenRecvBuffer);
    data->recvBuffer.resize(lenRecvBuffer);

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[32mReceiving result for job %d from %d (%luB)\x1b[0m\n",
           mpiStatus->MPI_TAG, mpiStatus->MPI_SOURCE, data->recvBuffer.size());
#endif

    // receive
    MPI_Recv(data->recvBuffer.data(), data->recvBuffer.size(), mpiJobDataType,
             mpiStatus->MPI_SOURCE, mpiStatus->MPI_TAG, mpiComm,
             MPI_STATUS_IGNORE);

    workerIsBusy[workerIdx] = false;

#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("\x1b[32mReceived result for job %d from %d\x1b[0m\n",
           mpiStatus->MPI_TAG, mpiStatus->MPI_SOURCE);
#endif

    // user-provided callback if specified
    if(data->callbackJobFinished)
        data->callbackJobFinished(data);


    // signal job done
    std::unique_lock<std::mutex> lock;
    if(data->jobDoneChangedMutex) {
        lock = std::unique_lock(*data->jobDoneChangedMutex);
    }
    if(data->jobDone)
        ++(*data->jobDone);
    if (data->jobDoneChangedCondition) {
        data->jobDoneChangedCondition->notify_all();
    }
    return workerIdx;
}

bool LoadBalancerMaster::sendQueuedJob(int freeWorkerIndex)
{
    if (freeWorkerIndex < 0)
        return false;

    JobData *currentQueueElement = getNextJob();

    if (currentQueueElement) {
        sendToWorker(freeWorkerIndex, currentQueueElement);
        sentJobsData[freeWorkerIndex] = currentQueueElement;
        return true;
    }
    return false;
}

void LoadBalancerMaster::sendTerminationSignalToAllWorkers() {
    int commSize;
    MPI_Comm_size(mpiComm, &commSize);

    MPI_Request reqs[commSize - 1];

    for (int i = 1; i < commSize; ++i) {
        reqs[i - 1] = MPI_REQUEST_NULL;
        MPI_Isend(MPI_BOTTOM, 0, MPI_INT, i, 0, mpiComm, &reqs[i - 1]);
    }
    MPI_Waitall(commSize - 1, reqs, MPI_STATUS_IGNORE);
}

bool LoadBalancerMaster::isRunning() const
{
#ifdef MASTER_QUEUE_H_SHOW_COMMUNICATION
    printf("LoadBalancerMaster::isRunning -> %d\n", isRunning_);
#endif
    return isRunning_;
}

int LoadBalancerMaster::getNumQueuedJobs() const
{
    std::unique_lock<std::mutex> lock(mutexQueue);
    return queue.size();
}

} // namespace parpe

#endif
