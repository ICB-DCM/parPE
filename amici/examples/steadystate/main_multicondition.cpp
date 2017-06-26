#include <stdio.h>
#include "steadystateProblemParallel.h"
#include <logging.h>
#include <loadBalancerWorker.h>
#include "hdf5Misc.h"
#include <unistd.h>
#include "SteadyStateMultiConditionProblem.h"
#include "wrapfunctions.h"
#include "simulationWorkerAmici.h"
#include "multiConditionProblemResultWriter.h"
#include "optimizationApplication.h"
/*
 * This example demonstrates the use of the loadbalancer / queue for parallel ODE simulation.
 */

void messageHandler(char **buffer, int *size, int jobId, void *userData);

class SteadystateApplication : public OptimizationApplication {
public:
    SteadystateApplication(int argc, char **argv) : OptimizationApplication(argc, argv) {}

    virtual void initProblem(const char *inFileArgument, const char *outFileArgument) {
        dataProvider = new SteadyStateMultiConditionDataProvider(inFileArgument);
        problem = new SteadyStateMultiConditionProblem(dataProvider);

        const char *outfilename = "testResultWriter_rank%03d.h5";
        char outfilefull[200];
        sprintf(outfilefull, outfilename, getMpiRank());
        JobIdentifier id = {0};
        resultWriter = new MultiConditionProblemResultWriter(problem, outfilefull, true, id);
        problem->resultWriter = resultWriter;
    }

    virtual void destroyProblem() {
        delete resultWriter;
        delete problem;
        delete dataProvider;
    }

    virtual void runWorker() {
        // TODO : move to base class; need wrapper; ParallelProblem interface?
        if(getMpiCommSize() > 1)
            loadBalancerWorkerRun(messageHandler, problem);
    }

    SteadyStateMultiConditionDataProvider *dataProvider;
};

class SteadystateLocalOptimizationApplication : public SteadystateApplication {
public:
    SteadystateLocalOptimizationApplication(int argc, char **argv) : SteadystateApplication(argc, argv) {}

    virtual int runMaster() {
            // Single optimization
            return getLocalOptimum(problem);
    }
};

class SteadystateMultiStartOptimizationApplication : public SteadystateApplication {
public:
    SteadystateMultiStartOptimizationApplication(int argc, char **argv) : SteadystateApplication(argc, argv) {}

    virtual int runMaster() {
        return getLocalOptimum(problem);

        int status = 0;
        // Multistart optimization
        OptimizationOptions options;
        options.maxOptimizerIterations = 1;
        options.numStarts = 1; // if numStarts > 1: need to use multiple MPI workers, otherwise simulation crashes due to CVODES threading issues

        MultiConditionProblemGeneratorForMultiStart generator;
        generator.options = &options;
        generator.resultWriter = reinterpret_cast<MultiConditionProblemResultWriter *>(problem->resultWriter);
        generator.dp = dataProvider;

        runParallelMultiStartOptimization(&generator, options.numStarts, options.retryOptimization);

        return status;
    }

};


int main(int argc, char **argv)
{
    int status = 0;

    SteadystateLocalOptimizationApplication app(argc, argv);
    status = app.run();

    return status;
}


void messageHandler(char** buffer, int *msgSize, int jobId, void *userData)
{
    //    int mpiRank;
    //    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    //    logmessage(LOGLVL_DEBUG, "Worker #%d: Job #%d received.", mpiRank, jobId);

    SteadyStateMultiConditionProblem *problem = (SteadyStateMultiConditionProblem *) userData;
    SteadyStateMultiConditionDataProvider *dataProvider = (SteadyStateMultiConditionDataProvider *)problem->getDataProvider();

    // unpack
    UserData udata = dataProvider->getModelDims(); // TODO get from buffer // TODO make sure this is full udata, not only model dins
    dataProvider->setupUserData(&udata);

    JobIdentifier path;
    JobAmiciSimulation::toUserData(*buffer, &udata, &path);
    free(*buffer);

    // work
    int status = 0;
    ReturnData *rdata = MultiConditionProblem::runAndLogSimulation(&udata, dataProvider, path, jobId, problem->resultWriter, &status);

    // pack & cleanup
    *msgSize = JobResultAmiciSimulation::getLength(udata.np);
    *buffer = (char*) malloc(*msgSize);
    JobResultAmiciSimulation::serialize(rdata, &udata, status, *buffer);

    delete rdata;

    /*

    SteadystateProblemParallel *problem = (SteadystateProblemParallel*) userData;
    UserData *udata = problem->udata;

    // unpack parameters
    int conditionIdx = (int) **buffer;
    int needGradient = (int) *(*buffer + sizeof(int));
    memcpy(udata->p, *buffer + 2 * sizeof(int), sizeof(double) * udata->np);
    free(*buffer);

    // read data for current conditions
    problem->readFixedParameters(conditionIdx);
    problem->readMeasurement(conditionIdx);
    problem->requireSensitivities(needGradient);

    // run simulation
    ReturnData *rdata = getSimulationResults(udata, problem->edata);

    // pack results
    *size = sizeof(double) * (udata->nplist + 1);
    *buffer = (char*) malloc(*size);
    double *doubleBuffer = (double*) *buffer;

    doubleBuffer[0] = rdata->llh[0];
    if(needGradient)
        for(int i = 0; i < udata->nplist; ++i)
            doubleBuffer[1 + i] = rdata->sllh[i];

    delete rdata;
    */
}
